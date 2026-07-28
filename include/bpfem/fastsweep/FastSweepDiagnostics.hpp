#pragma once

#include "bpfem/core/Types.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace fem::fastsweep {

struct FastSweepDiagnostics {
    std::string algorithm;
    std::vector<double> expansionFrequenciesHz;
    int requestedOrder = 0;
    int romDimension = 0;
    int retainedColumns = 0;
    int deflatedColumns = 0;
    double basisOrthogonalityError = 0.0;
    double padeInputPivotRatio = 0.0;
    double padeOutputPivotRatio = 0.0;
    // 兼容旧结果读取器；WCAWE v2 中该字段等于 basis relation residual。
    double wcaweMomentReconstructionError = 0.0;
    double wcaweRecurrenceResidualMax = 0.0;
    double wcaweBasisRelationResidualMax = 0.0;
    double wcaweOrthogonalityError = 0.0;
    double wcaweMinUpperDiagonalAbs = 0.0;
    double wcaweUpperDiagonalRatio = 0.0;
    std::uint64_t wcaweTriangularSolveCount = 0;
    int wcaweBreakdownOrder = 0;
    std::string wcaweTerminationReason;
    double wcaweMomentMatchingError = 0.0;
    double lanczosBiorthogonalityError = 0.0;
    double lanczosTridiagonalLeakage = 0.0;
    double lanczosFinalCoupling = 0.0;
    double momentMatchingError = 0.0;
    double lanczosRecurrenceResidual = 0.0;
    int lookAheadCount = 0;
    int selectiveReorthogonalizationCount = 0;
    double poleResidueReconstructionError = 0.0;
    int spuriousPoleCount = 0;
    std::uint64_t factorizedRhsSolveCount = 0;
    std::uint64_t factorizedSolveCallCount = 0;
    std::uint64_t batchRhsMax = 0;
    double portLinearizationSec = 0.0;
    double lanczosOperatorSec = 0.0;
    double orthogonalizationSec = 0.0;
    double poleDecompositionSec = 0.0;
    bool lanczosBreakdownDetected = false;
    bool reducedSolveSucceeded = true;
    double maxPassivityError = 0.0;
};

// 计算一组 S 参数点的最大无源性偏差。
double maxPassivityError(const std::vector<SParameterPoint>& points);

// 将 fast-sweep 诊断结构写成 JSON 文件。
bool writeDiagnosticsJson(const std::filesystem::path& path,
                          const FastSweepDiagnostics& diagnostics);

}  // namespace fem::fastsweep
