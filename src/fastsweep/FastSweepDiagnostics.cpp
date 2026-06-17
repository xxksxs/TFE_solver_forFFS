#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

namespace fem::fastsweep {

namespace {

// 写出 JSON 字符串，并处理引号、反斜杠和换行转义。
void writeJsonString(std::ostream& out, const std::string& s) {
    out << '"';
    for (char ch : s) {
        if (ch == '"' || ch == '\\') {
            out << '\\' << ch;
        } else if (ch == '\n') {
            out << "\\n";
        } else {
            out << ch;
        }
    }
    out << '"';
}

}  // namespace

// 计算 lossless 二端口中 |S11|^2+|S21|^2 偏离 1 的最大值。
double maxPassivityError(const std::vector<SParameterPoint>& points) {
    double maxErr = 0.0;
    for (const auto& p : points) {
        const double s11 = std::abs(p.s11);
        const double s21 = std::abs(p.s21);
        maxErr = std::max(maxErr, std::abs(s11 * s11 + s21 * s21 - 1.0));
    }
    return maxErr;
}

// 将 fast-sweep 统一诊断信息写成 diagnostics.json。
bool writeDiagnosticsJson(const std::filesystem::path& path,
                          const FastSweepDiagnostics& d) {
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out.good()) {
        return false;
    }
    out << std::setprecision(17);
    out << "{\n";
    out << "  \"schema\": 1,\n";
    out << "  \"algorithm\": ";
    writeJsonString(out, d.algorithm);
    out << ",\n";
    out << "  \"expansion_frequencies_hz\": [";
    for (std::size_t i = 0; i < d.expansionFrequenciesHz.size(); ++i) {
        if (i != 0) {
            out << ", ";
        }
        out << d.expansionFrequenciesHz[i];
    }
    out << "],\n";
    out << "  \"requested_order\": " << d.requestedOrder << ",\n";
    out << "  \"rom_dimension\": " << d.romDimension << ",\n";
    out << "  \"retained_columns\": " << d.retainedColumns << ",\n";
    out << "  \"deflated_columns\": " << d.deflatedColumns << ",\n";
    out << "  \"basis_orthogonality_error\": " << d.basisOrthogonalityError << ",\n";
    out << "  \"pade_input_pivot_ratio\": " << d.padeInputPivotRatio << ",\n";
    out << "  \"pade_output_pivot_ratio\": " << d.padeOutputPivotRatio << ",\n";
    out << "  \"wcawe_moment_reconstruction_error\": "
        << d.wcaweMomentReconstructionError << ",\n";
    out << "  \"reduced_solve_succeeded\": "
        << (d.reducedSolveSucceeded ? "true" : "false") << ",\n";
    out << "  \"max_passivity_error\": " << d.maxPassivityError << "\n";
    out << "}\n";
    out.flush();
    return out.good();
}

}  // namespace fem::fastsweep
