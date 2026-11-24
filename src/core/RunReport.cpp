#include "bpfem/core/RunReport.hpp"

#include "bpfem/app/Application.hpp"
#include "bpfem/core/RunEnvironment.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iomanip>
#include <limits>
#include <locale>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#if __has_include("BpfemBuildProvenance.hpp")
#include "BpfemBuildProvenance.hpp"
#else
#define BPFEM_GIT_REVISION "unknown"
#define BPFEM_GIT_STATE "unknown"
#define BPFEM_SOURCE_SHA256 ""
#endif

#ifdef BPFEM_USE_MKL
#include <mkl.h>
#endif
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <bcrypt.h>
#pragma comment(lib, "bcrypt.lib")
#else
#include <openssl/evp.h>
#endif

namespace fem {

namespace {

std::string sha256(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot read for SHA-256: " + path.string());
    std::array<char, 64 * 1024> buffer{};
    std::array<unsigned char, 32> digest{};
    auto check = [&](bool ok) {
        if (!ok) throw std::runtime_error("Cannot compute SHA-256: " + path.string());
    };
#ifdef _WIN32
    struct Hash {
        BCRYPT_ALG_HANDLE algorithm = nullptr;
        BCRYPT_HASH_HANDLE handle = nullptr;
        ~Hash() {
            if (handle) BCryptDestroyHash(handle);
            if (algorithm) BCryptCloseAlgorithmProvider(algorithm, 0);
        }
    } hash;
    check(BCryptOpenAlgorithmProvider(&hash.algorithm, BCRYPT_SHA256_ALGORITHM, nullptr, 0) >= 0);
    check(BCryptCreateHash(hash.algorithm, &hash.handle, nullptr, 0, nullptr, 0, 0) >= 0);
#else
    std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> hash(EVP_MD_CTX_new(), EVP_MD_CTX_free);
    check(hash && EVP_DigestInit_ex(hash.get(), EVP_sha256(), nullptr) == 1);
#endif
    while (in.read(buffer.data(), buffer.size()) || in.gcount() > 0) {
#ifdef _WIN32
        check(BCryptHashData(hash.handle, reinterpret_cast<unsigned char*>(buffer.data()),
                             static_cast<ULONG>(in.gcount()), 0) >= 0);
#else
        check(EVP_DigestUpdate(hash.get(), buffer.data(), static_cast<std::size_t>(in.gcount())) == 1);
#endif
    }
    check(in.eof() && !in.bad());
#ifdef _WIN32
    check(BCryptFinishHash(hash.handle, digest.data(), static_cast<ULONG>(digest.size()), 0) >= 0);
#else
    unsigned int length = 0;
    check(EVP_DigestFinal_ex(hash.get(), digest.data(), &length) == 1 && length == digest.size());
#endif
    std::ostringstream result;
    result << std::hex << std::setfill('0');
    for (const auto byte : digest) result << std::setw(2) << static_cast<unsigned int>(byte);
    return result.str();
}

std::string precise(double value) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return out.str();
}

// Minimal JSON string escape for the values we actually emit (no embedded
// newlines / control chars expected, but quotes and backslashes do appear in
// Windows paths).
std::string jsonEscape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 2);
    for (char c : s) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20U) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    out += buf;
                } else {
                    out += c;
                }
        }
    }
    return out;
}

class JsonWriter {
public:
    explicit JsonWriter(std::ostream& os) : os_(os) {}

    void beginObject() {
        if (!suppressIndent_) sep();
        maybeIndent();
        os_ << "{";
        needComma_.push_back(false);
        indent_++;
        nl();
    }
    void endObject() {
        indent_--;
        nl();
        for (int i = 0; i < indent_; ++i) os_ << "  ";
        os_ << "}";
        needComma_.pop_back();
    }
    void beginArray() {
        if (!suppressIndent_) sep();
        maybeIndent();
        os_ << "[";
        needComma_.push_back(false);
        indent_++;
        nl();
    }
    void endArray() {
        indent_--;
        nl();
        for (int i = 0; i < indent_; ++i) os_ << "  ";
        os_ << "]";
        needComma_.pop_back();
    }

    void key(const std::string& k) {
        sep();
        for (int i = 0; i < indent_; ++i) os_ << "  ";
        os_ << "\"" << jsonEscape(k) << "\": ";
        suppressIndent_ = true;
    }
    void value(const std::string& s) { write_("\"" + jsonEscape(s) + "\""); }
    void value(int v)                { write_(std::to_string(v)); }
    void value(std::size_t v)        { write_(std::to_string(v)); }
    void value(double v)             {
        if (!std::isfinite(v)) {
            write_("null");
            return;
        }
        write_(precise(v));
    }
    void valueRaw(const std::string& v) { write_(v); }

private:
    void maybeIndent() {
        if (suppressIndent_) {
            suppressIndent_ = false;
            return;
        }
        for (int i = 0; i < indent_; ++i) os_ << "  ";
    }
    void sep() {
        if (!needComma_.empty() && needComma_.back()) {
            os_ << ",";
            nl();
        }
        if (!needComma_.empty()) {
            needComma_.back() = true;
        }
    }
    void write_(const std::string& s) {
        if (!suppressIndent_) {
            sep();
            for (int i = 0; i < indent_; ++i) os_ << "  ";
        }
        suppressIndent_ = false;
        os_ << s;
    }
    void nl() { os_ << "\n"; }

    std::ostream& os_;
    std::vector<bool> needComma_;
    int indent_ = 0;
    bool suppressIndent_ = false;
};

std::string optionValue(const RunReport& report, const std::string& key) {
    for (const auto& [k, v] : report.options) {
        if (k == key) {
            return v;
        }
    }
    return {};
}

}  // namespace

void captureRunInputs(RunReport& report, const Options& options, const char* executable) {
    report.aedtSha256 = sha256(options.aedt);
    report.meshSha256 = sha256(options.mesh);
    std::filesystem::path executablePath;
#ifdef _WIN32
    std::wstring path(32768, L'\0');
    const auto length = GetModuleFileNameW(nullptr, path.data(), static_cast<DWORD>(path.size()));
    if (length > 0 && length < path.size()) executablePath = path.substr(0, length);
#elif defined(__linux__)
    std::error_code error;
    executablePath = std::filesystem::read_symlink("/proc/self/exe", error);
#endif
    if (executablePath.empty()) executablePath = executable;
    if (std::filesystem::is_regular_file(executablePath)) {
        report.executableSha256 = sha256(executablePath);
    }
#ifdef BPFEM_USE_MKL
    report.mklMaxThreads = mkl_get_max_threads();
    report.mklPardisoMaxThreads = mkl_domain_get_max_threads(MKL_DOMAIN_PARDISO);
#endif
    if (const char* value = std::getenv("OMP_NUM_THREADS")) report.ompNumThreads = value;
}

void verifyRunInputs(const RunReport& report, const Options& options) {
    if (report.aedtSha256 != sha256(options.aedt) || report.meshSha256 != sha256(options.mesh)) {
        throw std::runtime_error("AEDT or mesh changed during project import");
    }
}

void fillReportFromEnvironment(RunReport& report, const RunEnvironment& env, const Options& options) {
    report.solverVersion = env.solverVersion;
    report.buildString = env.buildString;
    report.runtimeString = env.runtimeString;
    report.hostName = env.hostName;
    report.osDescription = env.osDescription;
    report.cpuDescription = env.cpuDescription;
    report.memoryDescription = env.memoryDescription;
    report.commandLine = env.commandLine;
    report.workingDirectory = env.workingDirectory;

    auto add = [&](const char* k, const std::string& v) {
        report.options.emplace_back(k, v);
    };
    add("aedt", options.aedt.u8string());
    add("mesh", options.mesh.u8string());
    add("outDir", options.outDir.u8string());
    add("maxSweepPoints", std::to_string(options.maxSweepPoints));
    add("maxIterations", std::to_string(options.maxIterations));
    add("basisOrder", std::to_string(options.basisOrder));
    add("fieldOutputOrder", std::to_string(options.fieldOutputOrder));
    add("tolerance", precise(options.tolerance));
    add("writeAllFields", options.writeAllFields ? "true" : "false");
    {
        std::string s = "direct";
        switch (options.sweepStrategy) {
            case SweepStrategy::Direct: s = "direct"; break;
            case SweepStrategy::Alps:   s = "alps"; break;
            case SweepStrategy::Awe:    s = "awe"; break;
            case SweepStrategy::Gawe:   s = "gawe"; break;
            case SweepStrategy::Mgawe:  s = "mgawe"; break;
            case SweepStrategy::Wcawe:  s = "wcawe"; break;
        }
        add("sweepStrategy", s);
    }
    add("alpsOrder", std::to_string(options.alpsOrder));
    add("alpsExpansionFrequencyHz", precise(options.alpsExpansionFrequencyHz));
    add("aweOrder", std::to_string(options.aweOrder));
    add("aweExpansionFrequencyHz", precise(options.aweExpansionFrequencyHz));
    add("gaweOrder", std::to_string(options.gaweOrder));
    add("gaweExpansionFrequencyHz", precise(options.gaweExpansionFrequencyHz));
    add("gaweDropTolerance", precise(options.gaweDropTolerance));
    add("mgawePointCount", std::to_string(options.mgawePointCount));
    add("mgaweOrder", std::to_string(options.mgaweOrder));
    add("mgaweDropTolerance", precise(options.mgaweDropTolerance));
    add("wcaweOrder", std::to_string(options.wcaweOrder));
    add("wcaweExpansionFrequencyHz", precise(options.wcaweExpansionFrequencyHz));
    add("wcaweDropTolerance", precise(options.wcaweDropTolerance));
    add("aweFieldPade", options.aweFieldPade ? "true" : "false");
    add("mgaweAdaptive", options.mgaweAdaptive ? "true" : "false");
    add("mgaweAdaptiveTrueResidual", options.mgaweAdaptiveTrueResidual ? "true" : "false");
    add("mgaweTrueResidualTolerance", precise(options.mgaweTrueResidualTolerance));
    add("mgawe2001", options.mgawe2001 ? "true" : "false");
    std::string expansions;
    for (const auto& [frequency, order] : options.mgaweExpansions) {
        if (!expansions.empty()) expansions += ',';
        expansions += precise(frequency) + ':' + std::to_string(order);
    }
    add("mgaweExpansions", expansions);
    add("wcaweCheckMoments", options.wcaweCheckMoments ? "true" : "false");
    add("wcaweResidualDiagnostics", options.wcaweResidualDiagnostics ? "true" : "false");
    add("portMethod",
        options.portMethod == PortMethod::Analytic    ? "analytic"
      : options.portMethod == PortMethod::Transfinite ? "tfe"
                                                       : "numerical");
    add("tfeModesPerPort", std::to_string(options.tfeModesPerPort));
    {
        std::string ls = "auto";
        switch (options.linearSolver) {
            case LinearSolverKind::Direct:   ls = "direct"; break;
            case LinearSolverKind::BiCGStab: ls = "bicgstab"; break;
            case LinearSolverKind::Gmres:    ls = "gmres"; break;
            case LinearSolverKind::Mumps:    ls = "mumps"; break;
            case LinearSolverKind::Auto:     ls = "auto"; break;
        }
        add("linearSolver", ls);
    }
    {
        std::string p = "none";
        switch (options.preconditioner) {
            case PreconditionerKind::Jacobi: p = "jacobi"; break;
            case PreconditionerKind::Ilu0:   p = "ilu0"; break;
            case PreconditionerKind::Schwarz: p = "schwarz"; break;
            case PreconditionerKind::None:   p = "none"; break;
        }
        add("preconditioner", p);
    }
    add("gmresRestart", std::to_string(options.gmresRestart));
    add("iterativeSpmv", options.iterativeMklSpmv ? "mkl" : "native");
    add("pardisoMode", options.pardisoOutOfCore ? "ooc" : "incore");
    add("mumpsMode", options.mumpsMode);
    add("mumpsBlrTolerance", precise(options.mumpsBlrTolerance));
    add("schwarzSubdomains", std::to_string(options.schwarzSubdomains));
    add("schwarzOverlap", std::to_string(options.schwarzOverlap));
    add("schwarzLocalSolver", options.schwarzLocalMumps ? "mumps" : "pardiso");
}

bool writeRunReport(const std::filesystem::path& path, const RunReport& report) {
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out.good()) {
        return false;
    }
    JsonWriter w(out);
    w.beginObject();
      w.key("schema");           w.value(report.schema);
      w.key("status");           w.value(report.status);
      w.key("solver_version");   w.value(report.solverVersion);
      w.key("build");            w.value(report.buildString);
      w.key("runtime");          w.value(report.runtimeString);
      w.key("host");             w.value(report.hostName);
      w.key("os");               w.value(report.osDescription);
      w.key("cpu");              w.value(report.cpuDescription);
      w.key("memory");           w.value(report.memoryDescription);
      w.key("cmdline");          w.value(report.commandLine);
      w.key("cwd");              w.value(report.workingDirectory);
      w.key("result_directory"); w.value(report.outputDirectory);
      w.key("inputs");
      w.beginObject();
        w.key("aedt_sha256"); w.value(report.aedtSha256);
        w.key("mesh_sha256"); w.value(report.meshSha256);
      w.endObject();
      w.key("provenance");
      w.beginObject();
        w.key("git_revision"); w.value(BPFEM_GIT_REVISION);
        w.key("git_worktree_state"); w.value(BPFEM_GIT_STATE);
        w.key("source_sha256"); w.value(BPFEM_SOURCE_SHA256);
        w.key("executable_sha256"); w.value(report.executableSha256);
      w.endObject();
      w.key("threads");
      w.beginObject();
        w.key("mkl_max_threads");
        if (report.mklMaxThreads > 0) w.value(report.mklMaxThreads);
        else w.value(std::numeric_limits<double>::quiet_NaN());
        w.key("mkl_pardiso_max_threads");
        if (report.mklPardisoMaxThreads > 0) w.value(report.mklPardisoMaxThreads);
        else w.value(std::numeric_limits<double>::quiet_NaN());
        w.key("omp_num_threads_env"); w.value(report.ompNumThreads);
      w.endObject();
      w.key("sweep_strategy");   w.value(report.sweepStrategyName);
      w.key("linear_solver_backend"); w.value(report.linearSolverBackend);
      w.key("offline_build_s"); w.value(report.offlineBuildSec);
      w.key("port_linearization_s"); w.value(report.portLinearizationSec);
      w.key("lanczos_operator_s"); w.value(report.lanczosOperatorSec);
      w.key("pole_decomposition_s"); w.value(report.poleDecompositionSec);
      w.key("orthogonalization_s"); w.value(report.orthogonalizationSec);
      w.key("rom_projection_s"); w.value(report.romProjectionSec);
      w.key("online_sweep_s"); w.value(report.onlineSweepSec);
      w.key("symbolic_analysis_count");
      w.value(static_cast<std::size_t>(report.symbolicAnalysisCount));
      w.key("symbolic_analysis_s"); w.value(report.symbolicAnalysisSec);
      w.key("numeric_factorization_count");
      w.value(static_cast<std::size_t>(report.numericFactorizationCount));
      w.key("numeric_factorization_s"); w.value(report.numericFactorizationSec);
      w.key("factorized_rhs_solve_count");
      w.value(static_cast<std::size_t>(report.factorizedRhsSolveCount));
      w.key("factorized_rhs_solve_s"); w.value(report.factorizedRhsSolveSec);
      w.key("factorized_solve_call_count");
      w.value(static_cast<std::size_t>(report.factorizedSolveCallCount));
      w.key("batch_rhs_max");
      w.value(static_cast<std::size_t>(report.batchRhsMax));
      w.key("estimated_solver_peak_kib");
      w.value(static_cast<std::size_t>(report.estimatedSolverPeakKiB));

      w.key("options");
      w.beginObject();
        for (const auto& [k, v] : report.options) {
            w.key(k);
            w.value(k == "outDir" && !report.outputDirectory.empty() ? report.outputDirectory : v);
        }
      w.endObject();

      w.key("phases");
      w.beginArray();
        for (const auto& p : report.phases) {
            w.beginObject();
              w.key("name");           w.value(p.name);
              w.key("elapsed_s");      w.value(p.elapsedSec);
              w.key("end_current_mb"); w.value(p.endCurrentMb);
              w.key("end_peak_mb");    w.value(p.endPeakMb);
            w.endObject();
        }
      w.endArray();

      w.key("sweep");
      w.beginArray();
        for (const auto& fp : report.sweep) {
            w.beginObject();
              w.key("freq_hz");     w.value(fp.frequencyHz);
              w.key("assemble_s");  w.value(fp.assembleSec);
              w.key("solve_s");     w.value(fp.solveSec);
              w.key("residual");    w.value(fp.residual);
              if (!fp.linearSolveStatus.empty()) {
                  w.key("absolute_residual"); w.value(fp.absoluteResidual);
                  w.key("linear_solve_status"); w.value(fp.linearSolveStatus);
              }
              w.key("iterations");  w.value(fp.iterations);
              w.key("s11_db");      w.value(fp.s11Db);
              w.key("s21_db");      w.value(fp.s21Db);
            w.endObject();
        }
      w.endArray();

      w.key("total_elapsed_s"); w.value(report.totalElapsedSec);
      w.key("peak_memory_mb");  w.value(report.peakMemoryMb);
    w.endObject();
    out << "\n";
    out.close();
    return out.good();
}

bool writeTimingReport(const std::filesystem::path& path, const RunReport& report) {
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out.good()) {
        return false;
    }

    double totalAssembleSec = 0.0;
    double totalSolveSec = 0.0;
    for (const auto& fp : report.sweep) {
        totalAssembleSec += fp.assembleSec;
        totalSolveSec += fp.solveSec;
    }

    JsonWriter w(out);
    w.beginObject();
      w.key("schema");              w.value(1);
      w.key("status");              w.value(report.status);
      w.key("algorithm");           w.value(report.sweepStrategyName.empty()
                                              ? optionValue(report, "sweepStrategy")
                                              : report.sweepStrategyName);
      w.key("basis_order");         w.value(optionValue(report, "basisOrder"));
      w.key("max_sweep_points");    w.value(optionValue(report, "maxSweepPoints"));
      w.key("linear_solver");       w.value(optionValue(report, "linearSolver"));
      w.key("pardiso_mode");        w.value(optionValue(report, "pardisoMode"));
      w.key("mumps_mode");          w.value(optionValue(report, "mumpsMode"));
      w.key("iterative_spmv");      w.value(optionValue(report, "iterativeSpmv"));
      w.key("schwarz_local_solver"); w.value(optionValue(report, "schwarzLocalSolver"));
      w.key("linear_solver_backend"); w.value(report.linearSolverBackend);
      w.key("preconditioner");      w.value(optionValue(report, "preconditioner"));
      w.key("offline_build_s");      w.value(report.offlineBuildSec);
      w.key("port_linearization_s"); w.value(report.portLinearizationSec);
      w.key("lanczos_operator_s");   w.value(report.lanczosOperatorSec);
      w.key("pole_decomposition_s"); w.value(report.poleDecompositionSec);
      w.key("orthogonalization_s");   w.value(report.orthogonalizationSec);
      w.key("rom_projection_s");      w.value(report.romProjectionSec);
      w.key("online_sweep_s");         w.value(report.onlineSweepSec);
      w.key("symbolic_analysis_count");
      w.value(static_cast<std::size_t>(report.symbolicAnalysisCount));
      w.key("symbolic_analysis_s");    w.value(report.symbolicAnalysisSec);
      w.key("numeric_factorization_count");
      w.value(static_cast<std::size_t>(report.numericFactorizationCount));
      w.key("numeric_factorization_s"); w.value(report.numericFactorizationSec);
      w.key("factorized_rhs_solve_count");
      w.value(static_cast<std::size_t>(report.factorizedRhsSolveCount));
      w.key("factorized_rhs_solve_s"); w.value(report.factorizedRhsSolveSec);
      w.key("factorized_solve_call_count");
      w.value(static_cast<std::size_t>(report.factorizedSolveCallCount));
      w.key("batch_rhs_max");
      w.value(static_cast<std::size_t>(report.batchRhsMax));
      w.key("total_elapsed_s");     w.value(report.totalElapsedSec);
      w.key("peak_memory_mb");      w.value(report.peakMemoryMb);
      w.key("frequency_points");    w.value(report.sweep.size());
      w.key("total_assemble_s");    w.value(totalAssembleSec);
      w.key("total_solve_s");       w.value(totalSolveSec);
      w.key("estimated_solver_peak_kib");
      w.value(static_cast<std::size_t>(report.estimatedSolverPeakKiB));
      w.key("phases");
      w.beginArray();
        for (const auto& p : report.phases) {
            w.beginObject();
              w.key("name");           w.value(p.name);
              w.key("elapsed_s");      w.value(p.elapsedSec);
              w.key("end_current_mb"); w.value(p.endCurrentMb);
              w.key("end_peak_mb");    w.value(p.endPeakMb);
            w.endObject();
        }
      w.endArray();
    w.endObject();
    out << "\n";
    out.close();
    return out.good();
}

}  // namespace fem
