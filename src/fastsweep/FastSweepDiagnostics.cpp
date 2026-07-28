#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

namespace fem::fastsweep {

namespace {

// 写出 JSON 字符串，并处理引号、反斜杠和换行转义。
void writeJsonString(std::ostream& out, const std::string& value) {
    out << '"';
    for (char character : value) {
        if (character == '"' || character == '\\') {
            out << '\\' << character;
        } else if (character == '\n') {
            out << "\\n";
        } else {
            out << character;
        }
    }
    out << '"';
}

}  // namespace

// 计算 lossless 二端口中 |S11|^2+|S21|^2 偏离 1 的最大值。
double maxPassivityError(const std::vector<SParameterPoint>& points) {
    double maxError = 0.0;
    for (const auto& point : points) {
        const double s11 = std::abs(point.s11);
        const double s21 = std::abs(point.s21);
        maxError = std::max(
            maxError, std::abs(s11 * s11 + s21 * s21 - 1.0));
    }
    return maxError;
}

// 将 fast-sweep 统一诊断信息写成 diagnostics.json。
bool writeDiagnosticsJson(const std::filesystem::path& path,
                          const FastSweepDiagnostics& diagnostics) {
    std::ofstream output(path, std::ios::out | std::ios::trunc);
    if (!output.good()) {
        return false;
    }
    const auto& d = diagnostics;
    output << std::setprecision(17);
    output << "{\n";
    output << "  \"schema\": 1,\n";
    output << "  \"algorithm\": ";
    writeJsonString(output, d.algorithm);
    output << ",\n";
    output << "  \"expansion_frequencies_hz\": [";
    for (std::size_t index = 0; index < d.expansionFrequenciesHz.size(); ++index) {
        if (index != 0) {
            output << ", ";
        }
        output << d.expansionFrequenciesHz[index];
    }
    output << "],\n";
    output << "  \"requested_order\": " << d.requestedOrder << ",\n";
    output << "  \"rom_dimension\": " << d.romDimension << ",\n";
    output << "  \"retained_columns\": " << d.retainedColumns << ",\n";
    output << "  \"deflated_columns\": " << d.deflatedColumns << ",\n";
    output << "  \"basis_orthogonality_error\": "
           << d.basisOrthogonalityError << ",\n";
    output << "  \"pade_input_pivot_ratio\": "
           << d.padeInputPivotRatio << ",\n";
    output << "  \"pade_output_pivot_ratio\": "
           << d.padeOutputPivotRatio << ",\n";
    output << "  \"wcawe_moment_reconstruction_error\": "
           << d.wcaweMomentReconstructionError << ",\n";
    output << "  \"wcawe_recurrence_residual_max\": "
           << d.wcaweRecurrenceResidualMax << ",\n";
    output << "  \"wcawe_basis_relation_residual_max\": "
           << d.wcaweBasisRelationResidualMax << ",\n";
    output << "  \"wcawe_orthogonality_error\": "
           << d.wcaweOrthogonalityError << ",\n";
    output << "  \"wcawe_min_u_diagonal\": "
           << d.wcaweMinUpperDiagonalAbs << ",\n";
    output << "  \"wcawe_u_diagonal_ratio\": "
           << d.wcaweUpperDiagonalRatio << ",\n";
    output << "  \"wcawe_triangular_solve_count\": "
           << d.wcaweTriangularSolveCount << ",\n";
    output << "  \"wcawe_breakdown_order\": "
           << d.wcaweBreakdownOrder << ",\n";
    output << "  \"wcawe_termination_reason\": ";
    writeJsonString(output, d.wcaweTerminationReason);
    output << ",\n";
    output << "  \"wcawe_moment_matching_error\": "
           << d.wcaweMomentMatchingError << ",\n";
    output << "  \"lanczos_biorthogonality_error\": "
           << d.lanczosBiorthogonalityError << ",\n";
    output << "  \"lanczos_tridiagonal_leakage\": "
           << d.lanczosTridiagonalLeakage << ",\n";
    output << "  \"lanczos_final_coupling\": "
           << d.lanczosFinalCoupling << ",\n";
    output << "  \"moment_matching_error\": "
           << d.momentMatchingError << ",\n";
    output << "  \"lanczos_recurrence_residual\": "
           << d.lanczosRecurrenceResidual << ",\n";
    output << "  \"lookahead_count\": " << d.lookAheadCount << ",\n";
    output << "  \"selective_reorthogonalization_count\": "
           << d.selectiveReorthogonalizationCount << ",\n";
    output << "  \"pole_residue_reconstruction_error\": "
           << d.poleResidueReconstructionError << ",\n";
    output << "  \"spurious_pole_count\": "
           << d.spuriousPoleCount << ",\n";
    output << "  \"factorized_rhs_solve_count\": "
           << d.factorizedRhsSolveCount << ",\n";
    output << "  \"factorized_solve_call_count\": "
           << d.factorizedSolveCallCount << ",\n";
    output << "  \"batch_rhs_max\": " << d.batchRhsMax << ",\n";
    output << "  \"port_linearization_s\": "
           << d.portLinearizationSec << ",\n";
    output << "  \"lanczos_operator_s\": "
           << d.lanczosOperatorSec << ",\n";
    output << "  \"orthogonalization_s\": "
           << d.orthogonalizationSec << ",\n";
    output << "  \"pole_decomposition_s\": "
           << d.poleDecompositionSec << ",\n";
    output << "  \"lanczos_breakdown_detected\": "
           << (d.lanczosBreakdownDetected ? "true" : "false") << ",\n";
    output << "  \"reduced_solve_succeeded\": "
           << (d.reducedSolveSucceeded ? "true" : "false") << ",\n";
    output << "  \"max_passivity_error\": "
           << d.maxPassivityError << "\n";
    output << "}\n";
    output.flush();
    return output.good();
}

}  // namespace fem::fastsweep
