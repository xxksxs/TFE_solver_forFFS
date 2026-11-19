#include "bpfem/core/RunReport.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <ios>
#include <iomanip>
#include <sstream>
#include <string>

namespace fem {

namespace {

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
        std::ostringstream oss;
        oss << std::setprecision(15) << v;
        write_(oss.str());
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

}  // namespace

std::string RunReport::optionValue(const std::string& key) const {
    for (const auto& [k, v] : options) {
        if (k == key) {
            return v;
        }
    }
    return {};
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
      w.key("sweep_strategy");   w.value(report.sweepStrategyName);
      w.key("linear_solver_backend"); w.value(report.linearSolverBackend);
      w.key("port_formulation"); w.value(report.portFormulation);
      w.key("linear_solver_family"); w.value(report.linearSolverFamily);
      w.key("active_preconditioner"); w.value(report.activePreconditioner);
      w.key("active_iterative_spmv"); w.value(report.activeIterativeSpmv);
      w.key("active_schwarz_local_solver"); w.value(report.activeSchwarzLocalSolver);
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
            w.value(v);
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
    out.flush();
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
                                              ? report.optionValue("sweepStrategy")
                                              : report.sweepStrategyName);
      w.key("basis_order");         w.value(report.optionValue("basisOrder"));
      w.key("max_sweep_points");    w.value(report.optionValue("maxSweepPoints"));
      w.key("linear_solver");       w.value(report.optionValue("linearSolver"));
      w.key("pardiso_mode");        w.value(report.optionValue("pardisoMode"));
      w.key("mumps_mode");          w.value(report.optionValue("mumpsMode"));
      w.key("iterative_spmv");      w.value(report.optionValue("iterativeSpmv"));
      w.key("schwarz_local_solver"); w.value(report.optionValue("schwarzLocalSolver"));
      w.key("linear_solver_backend"); w.value(report.linearSolverBackend);
      w.key("preconditioner");      w.value(report.optionValue("preconditioner"));
      w.key("port_formulation"); w.value(report.portFormulation);
      w.key("linear_solver_family"); w.value(report.linearSolverFamily);
      w.key("active_preconditioner"); w.value(report.activePreconditioner);
      w.key("active_iterative_spmv"); w.value(report.activeIterativeSpmv);
      w.key("active_schwarz_local_solver"); w.value(report.activeSchwarzLocalSolver);
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
    out.flush();
    return out.good();
}

}  // namespace fem
