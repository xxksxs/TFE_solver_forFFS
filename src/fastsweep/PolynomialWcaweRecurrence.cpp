#include "bpfem/fastsweep/PolynomialWcaweRecurrence.hpp"

#include "bpfem/fastsweep/UpperTriangularBlockAction.hpp"
#include "bpfem/fastsweep/WellConditionedBasisBuilder.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace fem::fastsweep {

namespace {

// 判断全阶复向量是否全部为有限值。
bool isFiniteVector(const std::vector<Complex>& values) {
    for (const Complex value : values) {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
            return false;
        }
    }
    return true;
}

// 校验算子返回的向量长度和有限性，避免错误回调污染后续递推。
void validateVector(const std::vector<Complex>& values,
                    std::size_t expectedSize,
                    const char* context) {
    if (values.size() != expectedSize) {
        throw std::runtime_error(std::string(context) + " returned the wrong vector size");
    }
    if (!isFiniteVector(values)) {
        throw std::runtime_error(std::string(context) + " returned NaN or Inf");
    }
}

// 将 scale*source 累加到 target，供多个论文校正来源复用同一个全阶 RHS。
void addScaled(std::vector<Complex>& target,
               const std::vector<Complex>& source,
               Complex scale) {
    for (std::size_t i = 0; i < target.size(); ++i) {
        target[i] += scale * source[i];
    }
}

// 返回固定步长 U 的指定对角元绝对值。
double diagonalAbs(const std::vector<Complex>& upper, int stride, int index) {
    return std::abs(
        upper[static_cast<std::size_t>(index) * static_cast<std::size_t>(stride)
              + static_cast<std::size_t>(index)]);
}

// 计算 A0*candidate=rhs 的相对递推残差。
double recurrenceResidual(const PolynomialWcaweModel& model,
                          const std::vector<Complex>& candidate,
                          const std::vector<Complex>& rhs) {
    std::vector<Complex> residual = model.applyA0(candidate);
    validateVector(residual, model.fullDimension, "WCAWE A0 operator");
    for (std::size_t i = 0; i < residual.size(); ++i) {
        residual[i] -= rhs[i];
    }
    return norm2(residual) / std::max(norm2(rhs), 1.0e-300);
}

// 校验通用多项式模型和构建选项的基本契约。
void validateInputs(const PolynomialWcaweModel& model,
                    const WcaweSolveFunction& solveAtExpansion,
                    const WcaweBuildOptions& options) {
    if (model.fullDimension == 0 || model.matrixDegree < 1
        || !model.applyA0 || !model.rhsCoefficientAt || !solveAtExpansion) {
        throw std::invalid_argument(
            "PolynomialWcaweRecurrence: incomplete polynomial model");
    }
    if (model.applyHigherOrder.size() != model.matrixDegree) {
        throw std::invalid_argument(
            "PolynomialWcaweRecurrence: matrix degree does not match operator count");
    }
    for (const auto& op : model.applyHigherOrder) {
        if (!op) {
            throw std::invalid_argument(
                "PolynomialWcaweRecurrence: matrix coefficient operator is empty");
        }
    }
    if (options.requestedOrder < 1
        || !std::isfinite(options.breakdownTolerance)
        || options.breakdownTolerance <= 0.0
        || !std::isfinite(options.recurrenceTolerance)
        || options.recurrenceTolerance <= 0.0
        || !std::isfinite(options.orthogonalityTolerance)
        || options.orthogonalityTolerance <= 0.0
        || options.reorthogonalizationPasses < 1
        || options.reorthogonalizationPasses > 2) {
        throw std::invalid_argument(
            "PolynomialWcaweRecurrence: invalid build options");
    }
}

}  // namespace

// 把 WCAWE 终止枚举转换成稳定字符串，避免输出依赖编译器枚举格式。
const char* wcaweTerminationReasonName(WcaweTerminationReason reason) {
    switch (reason) {
        case WcaweTerminationReason::ReachedRequestedOrder:
            return "reached_requested_order";
        case WcaweTerminationReason::HappyBreakdown:
            return "happy_breakdown";
        case WcaweTerminationReason::NearSingularUpperTriangularFactor:
            return "near_singular_upper_triangular_factor";
        case WcaweTerminationReason::SolverFailure:
            return "solver_failure";
        case WcaweTerminationReason::NonFiniteValue:
            return "non_finite_value";
    }
    return "unknown";
}

// 按论文式 (7) 组装一个校正 RHS，并通过式 (9) 的连续子块回代读取历史 U。
std::vector<Complex> PolynomialWcaweRecurrence::assembleCorrectedRhs(
    const PolynomialWcaweModel& model,
    const std::vector<std::vector<Complex>>& basis,
    const std::vector<Complex>& upperTriangularU,
    int stride,
    int n,
    double singularTolerance,
    std::uint64_t* triangularSolveCount) {
    if (n < 1 || n > stride || static_cast<int>(basis.size()) != n - 1) {
        throw std::invalid_argument(
            "PolynomialWcaweRecurrence: inconsistent order in corrected RHS");
    }
    if (n == 1) {
        std::vector<Complex> rhs = model.rhsCoefficientAt(0);
        validateVector(rhs, model.fullDimension, "WCAWE RHS coefficient");
        return rhs;
    }

    std::vector<Complex> rhs(model.fullDimension, Complex(0.0, 0.0));
    const int maxRhsOrder =
        std::min(static_cast<int>(model.rhsDegree), n - 1);
    for (int m = 1; m <= maxRhsOrder; ++m) {
        const auto action = UpperTriangularBlockAction::applyToLastUnitVector(
            upperTriangularU, stride, n - 1, n, m, 1,
            singularTolerance, triangularSolveCount);
        std::vector<Complex> coefficient =
            model.rhsCoefficientAt(static_cast<std::size_t>(m));
        validateVector(coefficient, model.fullDimension, "WCAWE RHS coefficient");
        addScaled(rhs, coefficient, action.front());
    }

    std::vector<Complex> firstOrder =
        model.applyHigherOrder.front()(basis[static_cast<std::size_t>(n - 2)]);
    validateVector(firstOrder, model.fullDimension, "WCAWE A1 operator");
    addScaled(rhs, firstOrder, Complex(-1.0, 0.0));

    const int maxMatrixOrder =
        std::min(static_cast<int>(model.matrixDegree), n - 1);
    for (int m = 2; m <= maxMatrixOrder; ++m) {
        const auto action = UpperTriangularBlockAction::applyToLastUnitVector(
            upperTriangularU, stride, n - 1, n, m, 2,
            singularTolerance, triangularSolveCount);
        std::vector<Complex> combination(model.fullDimension, Complex(0.0, 0.0));
        for (int column = 0; column < n - m; ++column) {
            addScaled(combination, basis[static_cast<std::size_t>(column)],
                      action[static_cast<std::size_t>(column)]);
        }
        std::vector<Complex> higherOrder =
            model.applyHigherOrder[static_cast<std::size_t>(m - 1)](combination);
        validateVector(higherOrder, model.fullDimension, "WCAWE higher-order operator");
        addScaled(rhs, higherOrder, Complex(-1.0, 0.0));
    }
    return rhs;
}

// 逐阶执行论文校正、A0 回代和 MGS，生成下一阶会继续使用的 V 与 U。
WcaweBuildResult PolynomialWcaweRecurrence::build(
    const PolynomialWcaweModel& model,
    const WcaweSolveFunction& solveAtExpansion,
    const WcaweBuildOptions& options) {
    validateInputs(model, solveAtExpansion, options);

    WcaweBuildResult result;
    result.requestedOrder = options.requestedOrder;
    result.steps.reserve(static_cast<std::size_t>(options.requestedOrder));
    result.minUpperDiagonalAbs = std::numeric_limits<double>::infinity();

    WellConditionedBasisBuilder builder;
    builder.reset(options.requestedOrder, options.breakdownTolerance,
                  options.reorthogonalizationPasses);

    for (int n = 1; n <= options.requestedOrder; ++n) {
        const std::uint64_t solvesBefore = result.triangularSolveCount;
        std::vector<Complex> rhs;
        try {
            rhs = assembleCorrectedRhs(
                model, builder.basis(), builder.upperTriangularU(),
                builder.stride(), n, options.breakdownTolerance,
                &result.triangularSolveCount);
        } catch (const UpperTriangularBreakdown& error) {
            if (n == 1) {
                throw;
            }
            result.terminationReason =
                WcaweTerminationReason::NearSingularUpperTriangularFactor;
            result.terminationMessage = error.what();
            break;
        }

        std::vector<Complex> candidate;
        try {
            candidate = solveAtExpansion(rhs);
        } catch (const std::exception& error) {
            if (n == 1) {
                throw;
            }
            result.terminationReason = WcaweTerminationReason::SolverFailure;
            result.terminationMessage = error.what();
            break;
        }
        if (candidate.size() != model.fullDimension || !isFiniteVector(candidate)) {
            if (n == 1) {
                throw std::runtime_error(
                    "PolynomialWcaweRecurrence: invalid first solve result");
            }
            result.terminationReason = WcaweTerminationReason::NonFiniteValue;
            result.terminationMessage =
                "solver returned a wrong-sized or non-finite candidate";
            break;
        }

        const double residual = recurrenceResidual(model, candidate, rhs);
        const auto orthogonalizationStarted = std::chrono::steady_clock::now();
        const WcaweAppendResult append = builder.appendCandidate(candidate);
        result.orthogonalizationSec += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - orthogonalizationStarted).count();
        if (!append.accepted) {
            if (n == 1) {
                throw std::runtime_error(
                    "PolynomialWcaweRecurrence: first WCAWE candidate broke down");
            }
            result.terminationReason =
                append.candidateNormBeforeMgs <= options.breakdownTolerance
                    ? WcaweTerminationReason::HappyBreakdown
                    : WcaweTerminationReason::NearSingularUpperTriangularFactor;
            result.terminationMessage =
                "MGS produced a near-zero U diagonal at order " + std::to_string(n);
            break;
        }

        const double diagonal = append.diagonalAbs;
        result.minUpperDiagonalAbs = std::min(result.minUpperDiagonalAbs, diagonal);
        double maxDiagonal = 0.0;
        for (int index = 0; index < n; ++index) {
            maxDiagonal = std::max(
                maxDiagonal,
                diagonalAbs(builder.upperTriangularU(), builder.stride(), index));
        }
        result.upperDiagonalRatio =
            maxDiagonal / std::max(result.minUpperDiagonalAbs, 1.0e-300);

        WcaweStepDiagnostics step;
        step.order = n;
        step.candidateNormBeforeMgs = append.candidateNormBeforeMgs;
        step.diagonalAbs = diagonal;
        step.diagonalRatio = result.upperDiagonalRatio;
        step.recurrenceResidual = residual;
        step.basisRelationResidual = append.basisRelationResidual;
        step.orthogonalityError = append.orthogonalityError;
        step.reorthogonalizationPasses = append.reorthogonalizationPasses;
        step.triangularSolveCount =
            result.triangularSolveCount - solvesBefore;
        result.steps.push_back(step);
        result.maxRecurrenceResidual =
            std::max(result.maxRecurrenceResidual, residual);
        result.maxBasisRelationResidual =
            std::max(result.maxBasisRelationResidual, append.basisRelationResidual);
        result.maxOrthogonalityError =
            std::max(result.maxOrthogonalityError, append.orthogonalityError);
    }

    result.achievedOrder = builder.dimension();
    result.upperTriangularU = builder.upperTriangularU();
    result.basis = builder.takeBasis();
    if (result.achievedOrder == options.requestedOrder) {
        result.terminationReason = WcaweTerminationReason::ReachedRequestedOrder;
        result.terminationMessage.clear();
    }
    if (!std::isfinite(result.minUpperDiagonalAbs)) {
        result.minUpperDiagonalAbs = 0.0;
    }
    result.momentMatchingError =
        std::max(result.maxRecurrenceResidual, result.maxBasisRelationResidual);
    return result;
}

}  // namespace fem::fastsweep
