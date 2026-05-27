#include "bpfem/core/RunReport.hpp"

#include "bpfem/app/Application.hpp"
#include "bpfem/core/RunEnvironment.hpp"

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
    add("aedt", options.aedt.string());
    add("mesh", options.mesh.string());
    add("outDir", options.outDir.string());
    add("maxSweepPoints", std::to_string(options.maxSweepPoints));
    add("maxIterations", std::to_string(options.maxIterations));
    add("basisOrder", std::to_string(options.basisOrder));
    add("fieldOutputOrder", std::to_string(options.fieldOutputOrder));
    add("tolerance", std::to_string(options.tolerance));
    add("writeAllFields", options.writeAllFields ? "true" : "false");
    add("sweepStrategy",
        options.sweepStrategy == SweepStrategy::Direct ? "direct" : "alps");
    add("alpsKrylovOrder", std::to_string(options.alpsKrylovOrder));
    add("alpsExpansionFrequencyHz", std::to_string(options.alpsExpansionFrequencyHz));
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
            case LinearSolverKind::Auto:     ls = "auto"; break;
        }
        add("linearSolver", ls);
    }
    {
        std::string p = "none";
        switch (options.preconditioner) {
            case PreconditionerKind::Jacobi: p = "jacobi"; break;
            case PreconditionerKind::Ilu0:   p = "ilu0"; break;
            case PreconditionerKind::None:   p = "none"; break;
        }
        add("preconditioner", p);
    }
    add("gmresRestart", std::to_string(options.gmresRestart));
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

}  // namespace fem
