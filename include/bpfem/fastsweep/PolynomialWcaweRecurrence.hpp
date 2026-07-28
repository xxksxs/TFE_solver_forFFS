#pragma once

#include "bpfem/fastsweep/LinearAlgebra.hpp"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace fem::fastsweep {

using WcaweLinearOperator =
    std::function<std::vector<Complex>(const std::vector<Complex>&)>;
using WcaweRhsCoefficientFunction =
    std::function<std::vector<Complex>(std::size_t order)>;
using WcaweSolveFunction =
    std::function<std::vector<Complex>(const std::vector<Complex>& rhs)>;

struct PolynomialWcaweModel {
    std::size_t fullDimension = 0;
    std::size_t matrixDegree = 0;
    std::size_t rhsDegree = 0;
    WcaweLinearOperator applyA0;
    std::vector<WcaweLinearOperator> applyHigherOrder;
    WcaweRhsCoefficientFunction rhsCoefficientAt;
};

struct WcaweBuildOptions {
    int requestedOrder = 12;
    double breakdownTolerance = 1.0e-12;
    double recurrenceTolerance = 1.0e-10;
    double orthogonalityTolerance = 1.0e-10;
    int reorthogonalizationPasses = 2;
};

enum class WcaweTerminationReason {
    ReachedRequestedOrder,
    HappyBreakdown,
    NearSingularUpperTriangularFactor,
    SolverFailure,
    NonFiniteValue
};

struct WcaweStepDiagnostics {
    int order = 0;
    double candidateNormBeforeMgs = 0.0;
    double diagonalAbs = 0.0;
    double diagonalRatio = 0.0;
    double recurrenceResidual = 0.0;
    double basisRelationResidual = 0.0;
    double orthogonalityError = 0.0;
    int reorthogonalizationPasses = 0;
    std::uint64_t triangularSolveCount = 0;
};

struct WcaweBuildResult {
    std::vector<std::vector<Complex>> basis;
    std::vector<Complex> upperTriangularU;
    int requestedOrder = 0;
    int achievedOrder = 0;
    WcaweTerminationReason terminationReason =
        WcaweTerminationReason::ReachedRequestedOrder;
    std::string terminationMessage;
    std::vector<WcaweStepDiagnostics> steps;
    double maxRecurrenceResidual = 0.0;
    double maxBasisRelationResidual = 0.0;
    double maxOrthogonalityError = 0.0;
    double minUpperDiagonalAbs = 0.0;
    double upperDiagonalRatio = 0.0;
    double momentMatchingError = 0.0;
    double orthogonalizationSec = 0.0;
    std::uint64_t triangularSolveCount = 0;
};

class PolynomialWcaweRecurrence {
public:
    // 按 Slone 2003 式 (7)-(9)逐阶生成正交 WCAWE 基及上三角系数 U。
    static WcaweBuildResult build(
        const PolynomialWcaweModel& model,
        const WcaweSolveFunction& solveAtExpansion,
        const WcaweBuildOptions& options = {});

    // 组装第 n 阶论文校正 RHS；n 使用一基索引，basis/U 必须已包含前 n-1 阶。
    static std::vector<Complex> assembleCorrectedRhs(
        const PolynomialWcaweModel& model,
        const std::vector<std::vector<Complex>>& basis,
        const std::vector<Complex>& upperTriangularU,
        int stride,
        int n,
        double singularTolerance,
        std::uint64_t* triangularSolveCount = nullptr);
};

// 返回稳定的英文终止原因，供日志与 JSON 使用。
const char* wcaweTerminationReasonName(WcaweTerminationReason reason);

}  // namespace fem::fastsweep
