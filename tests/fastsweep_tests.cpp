#include "bpfem/fastsweep/GalerkinReducedModel.hpp"
#include "bpfem/fastsweep/LanczosPadeModel.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PadeApproximant.hpp"
#include "bpfem/fastsweep/PolynomialMomentRecurrence.hpp"
#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"
#include "bpfem/fastsweep/PolynomialWcaweRecurrence.hpp"
#include "bpfem/fastsweep/UpperTriangularBlockAction.hpp"
#include "bpfem/fastsweep/WellConditionedBasisBuilder.hpp"
#include "bpfem/fem/EdgeTopology.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/FactorizedSolveSession.hpp"
#include "bpfem/linalg/IFactorizedSparseSolver.hpp"
#include "bpfem/linalg/MklPardisoSolver.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Complex = std::complex<double>;

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void requireNear(double actual, double expected, double tolerance, const char* message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message);
    }
}

void requireNear(Complex actual, Complex expected, double tolerance, const char* message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message);
    }
}

// 记录分解和 RHS 调用次数，用于验证 fast-sweep 的求解生命周期。
class CountingFactorizedSolver final
    : public fem::linalg::ISparseSolver,
      public fem::linalg::IFactorizedSparseSolver {
public:
    fem::SolveResult solve(const fem::SparseMatrix&,
                           const std::vector<Complex>& rhs,
                           const fem::linalg::SolverConfig&) override {
        ++legacySolveCount;
        return {rhs, 1, 0.0};
    }

    void factorize(const fem::SparseMatrix&,
                   const fem::linalg::SolverConfig&) override {
        factorized = true;
        ++factorizeCount;
    }

    fem::SolveResult solveFactorized(
        const std::vector<Complex>& rhs,
        const fem::linalg::SolverConfig&) override {
        if (!factorized) {
            throw std::runtime_error("solveFactorized called before factorize");
        }
        ++factorizedSolveCount;
        return {rhs, 1, 0.0};
    }

    // 模拟一次调用处理多个 RHS，并分别统计向量数与调用数。
    std::vector<fem::SolveResult> solveFactorizedBatch(
        const std::vector<std::vector<Complex>>& rhsList,
        const fem::linalg::SolverConfig&) override {
        if (!factorized) {
            throw std::runtime_error("solveFactorizedBatch called before factorize");
        }
        ++batchSolveCallCount;
        batchRhsMax = std::max(batchRhsMax, static_cast<int>(rhsList.size()));
        factorizedSolveCount += static_cast<int>(rhsList.size());
        std::vector<fem::SolveResult> results;
        results.reserve(rhsList.size());
        for (const auto& rhs : rhsList) {
            results.push_back({rhs, 1, 0.0});
        }
        return results;
    }

    void clearFactorization() override { factorized = false; }

    fem::linalg::FactorizationStatistics factorizationStatistics() const override {
        fem::linalg::FactorizationStatistics stats;
        stats.numericFactorizationCount = factorizeCount;
        stats.factorizedRhsSolveCount = factorizedSolveCount;
        stats.factorizedSolveCallCount = batchSolveCallCount;
        stats.batchRhsMax = batchRhsMax;
        return stats;
    }

    const char* name() const override { return "counting-factorized"; }

    int factorizeCount = 0;
    int factorizedSolveCount = 0;
    int legacySolveCount = 0;
    int batchSolveCallCount = 0;
    int batchRhsMax = 0;
    bool factorized = false;
};

// 模拟不支持数值分解复用的迭代后端。
class CountingLegacySolver final : public fem::linalg::ISparseSolver {
public:
    fem::SolveResult solve(const fem::SparseMatrix&,
                           const std::vector<Complex>& rhs,
                           const fem::linalg::SolverConfig&) override {
        ++solveCount;
        return {rhs, 1, 0.0};
    }

    const char* name() const override { return "counting-legacy"; }

    int solveCount = 0;
};

// 验证单点、多点和迭代回退路径的分解/回代次数。
void testFactorizedSolveSessionCounts() {
    fem::SparseMatrix matrix(1);
    matrix.add(0, 0, Complex(2.0, 0.0));
    const std::vector<Complex> rhs = {Complex(1.0, 0.0)};

    CountingFactorizedSolver wcaweSolver;
    {
        fem::linalg::FactorizedSolveSession session(wcaweSolver, matrix);
        for (int order = 0; order < 12; ++order) {
            (void)session.solve(rhs);
        }
    }
    require(wcaweSolver.factorizeCount == 1, "WCAWE must factorize once");
    require(wcaweSolver.factorizedSolveCount == 12, "WCAWE must backsolve 12 RHS vectors");
    require(wcaweSolver.legacySolveCount == 0, "WCAWE unexpectedly used legacy solve");

    CountingFactorizedSolver batchSolver;
    {
        fem::linalg::FactorizedSolveSession session(batchSolver, matrix);
        const auto results = session.solveBatch({rhs, rhs});
        require(results.size() == 2, "Batch session result count mismatch");
    }
    require(batchSolver.factorizeCount == 1, "Batch session must factorize once");
    require(batchSolver.factorizedSolveCount == 2, "Batch RHS vector count mismatch");
    require(batchSolver.batchSolveCallCount == 1, "Batch solve call count mismatch");
    require(batchSolver.batchRhsMax == 2, "Batch maximum RHS count mismatch");

    CountingFactorizedSolver mgaweSolver;
    for (int point = 0; point < 3; ++point) {
        fem::linalg::FactorizedSolveSession session(mgaweSolver, matrix);
        for (int order = 0; order < 4; ++order) {
            (void)session.solve(rhs);
        }
    }
    require(mgaweSolver.factorizeCount == 3, "MGAWE must factorize once per expansion point");
    require(mgaweSolver.factorizedSolveCount == 12, "MGAWE must backsolve 12 RHS vectors");

    CountingLegacySolver iterativeSolver;
    fem::linalg::FactorizedSolveSession fallback(iterativeSolver, matrix);
    require(!fallback.reusesFactorization(), "Legacy solver incorrectly reports factor reuse");
    for (int order = 0; order < 4; ++order) {
        (void)fallback.solve(rhs);
    }
    require(iterativeSolver.solveCount == 4, "Legacy solver fallback count mismatch");
}

#ifdef BPFEM_USE_MKL
// 验证真实 PARDISO 的 phase=22/33 解、残差和模式复用状态。
void testPardisoFactorizedLifecycle() {
    fem::SparseMatrix matrix(2);
    matrix.add(0, 0, Complex(4.0, 0.0));
    matrix.add(0, 1, Complex(1.0, 0.0));
    matrix.add(1, 0, Complex(1.0, 0.0));
    matrix.add(1, 1, Complex(3.0, 0.0));

    fem::MklPardisoSolver solver;
    solver.factorize(matrix);
    const auto first = solver.solveFactorized({Complex(1.0, 0.0), Complex(2.0, 0.0)});
    const auto second = solver.solveFactorized({Complex(2.0, 0.0), Complex(1.0, 0.0)});
    const auto batch = solver.solveFactorizedBatch({
        {Complex(1.0, 0.0), Complex(2.0, 0.0)},
        {Complex(2.0, 0.0), Complex(1.0, 0.0)}
    });
    require(batch.size() == 2, "PARDISO batch result count mismatch");
    requireNear(batch[0].field[0], first.field[0], 1.0e-12,
                "PARDISO batch first solution mismatch");
    requireNear(batch[1].field[1], second.field[1], 1.0e-12,
                "PARDISO batch second solution mismatch");
    const auto firstAx = matrix.multiply(first.field);
    const auto secondAx = matrix.multiply(second.field);
    const double firstResidual = std::hypot(
        std::abs(firstAx[0] - Complex(1.0, 0.0)),
        std::abs(firstAx[1] - Complex(2.0, 0.0))) / std::sqrt(5.0);
    const double secondResidual = std::hypot(
        std::abs(secondAx[0] - Complex(2.0, 0.0)),
        std::abs(secondAx[1] - Complex(1.0, 0.0))) / std::sqrt(5.0);
    require(firstResidual < 1.0e-12, "PARDISO first factorized residual too large");
    require(secondResidual < 1.0e-12, "PARDISO second factorized residual too large");

    auto stats = solver.factorizationStatistics();
    require(stats.symbolicAnalysisCount == 1, "PARDISO symbolic analysis count mismatch");
    require(stats.numericFactorizationCount == 1, "PARDISO numeric factorization count mismatch");
    require(stats.factorizedRhsSolveCount == 4, "PARDISO RHS solve count mismatch");
    require(stats.factorizedSolveCallCount == 3, "PARDISO solve call count mismatch");
    require(stats.batchRhsMax == 2, "PARDISO batch maximum mismatch");

    fem::SparseMatrix samePattern(2);
    samePattern.add(0, 0, Complex(5.0, 0.0));
    samePattern.add(0, 1, Complex(1.0, 0.0));
    samePattern.add(1, 0, Complex(1.0, 0.0));
    samePattern.add(1, 1, Complex(4.0, 0.0));
    solver.factorize(samePattern);
    stats = solver.factorizationStatistics();
    require(stats.symbolicAnalysisCount == 1, "Same pattern repeated symbolic analysis");
    require(stats.numericFactorizationCount == 2, "Same pattern refactorization count mismatch");

    fem::SparseMatrix changedPattern(2);
    changedPattern.add(0, 0, Complex(2.0, 0.0));
    changedPattern.add(1, 1, Complex(3.0, 0.0));
    solver.factorize(changedPattern);
    stats = solver.factorizationStatistics();
    require(stats.symbolicAnalysisCount == 2, "Changed pattern did not trigger analysis");

    solver.clearFactorization();
    bool threw = false;
    try {
        (void)solver.solveFactorized({Complex(1.0, 0.0), Complex(1.0, 0.0)});
    } catch (const std::runtime_error&) {
        threw = true;
    }
    require(threw, "PARDISO solveFactorized must reject missing factors");
}
#endif

void testPolynomialMomentRecurrence() {
    std::vector<std::vector<Complex>> rhs = {
        {Complex(1.0, 0.0)},
        {Complex(0.0, 0.0)},
        {Complex(0.0, 0.0)}
    };
    std::vector<fem::fastsweep::PolynomialMomentRecurrence::LinearOperator> ops;
    ops.push_back([](const std::vector<Complex>& x) {
        return std::vector<Complex>{x[0]};
    });
    auto solve = [](const std::vector<Complex>& b) {
        return std::vector<Complex>{b[0] / Complex(2.0, 0.0)};
    };

    const auto moments =
        fem::fastsweep::PolynomialMomentRecurrence::generatePolynomialMoments(rhs, ops, solve);
    require(moments.size() == 3, "PolynomialMomentRecurrence returned wrong count");
    requireNear(moments[0][0], Complex(0.5, 0.0), 1.0e-14, "moment 0 mismatch");
    requireNear(moments[1][0], Complex(-0.25, 0.0), 1.0e-14, "moment 1 mismatch");
    requireNear(moments[2][0], Complex(0.125, 0.0), 1.0e-14, "moment 2 mismatch");
}

// 用中心有限差分验证 P1 线性化的 A1*x 与 b1 导数符号和尺度。
void testLosslessP1FiniteDifference() {
    fem::FEMAssembler::AffineSystem affine;
    affine.M = fem::SparseMatrix(2);
    affine.M.add(0, 0, Complex(2.0, 0.0));
    affine.M.add(1, 1, Complex(3.0, 0.0));
    fem::fastsweep::PolynomialPortMomentBuilder::LosslessLinearization linearization;
    linearization.matrixAtExpansion = fem::SparseMatrix(2);
    linearization.matrixAtExpansion.add(0, 0, Complex(4.0, 0.0));
    linearization.matrixAtExpansion.add(0, 1, Complex(0.5, 0.0));
    linearization.matrixAtExpansion.add(1, 0, Complex(0.5, 0.0));
    linearization.matrixAtExpansion.add(1, 1, Complex(5.0, 0.0));
    linearization.lambdaScale = 1.7;
    linearization.portAdmittanceFirstCoefficients = {0.25};
    linearization.rhs0 = {Complex(1.0, 0.2), Complex(-0.3, 0.1)};
    linearization.rhs1 = {Complex(0.05, -0.02), Complex(-0.04, 0.03)};
    const std::vector<std::vector<Complex>> portVectors = {
        {Complex(1.0, 0.0), Complex(0.5, 0.0)}
    };
    const std::vector<Complex> x = {Complex(0.7, -0.1), Complex(-0.2, 0.3)};
    const auto derivative =
        fem::fastsweep::PolynomialPortMomentBuilder::applyLosslessFirstOrderMatrix(
            affine, linearization, portVectors, x);
    const auto base = linearization.matrixAtExpansion.multiply(x);
    constexpr double step = 1.0e-6;
    std::vector<Complex> plus(2);
    std::vector<Complex> minus(2);
    for (std::size_t i = 0; i < 2; ++i) {
        plus[i] = base[i] + step * derivative[i];
        minus[i] = base[i] - step * derivative[i];
        const Complex finiteDifference = (plus[i] - minus[i]) / (2.0 * step);
        const double relativeError = std::abs(finiteDifference - derivative[i])
            / std::max(std::abs(derivative[i]), 1.0e-12);
        require(relativeError <= 1.0e-7, "P1 A1*x finite-difference mismatch");

        const Complex rhsPlus = linearization.rhs0[i] + step * linearization.rhs1[i];
        const Complex rhsMinus = linearization.rhs0[i] - step * linearization.rhs1[i];
        const Complex rhsDerivative = (rhsPlus - rhsMinus) / (2.0 * step);
        require(std::abs(rhsDerivative - linearization.rhs1[i])
                    / std::max(std::abs(linearization.rhs1[i]), 1.0e-12) <= 1.0e-7,
                "P1 b1 finite-difference mismatch");
    }
}

void testPadeApproximant() {
    const std::vector<Complex> moments = {
        Complex(1.0, 0.0),
        Complex(2.0, 0.0),
        Complex(4.0, 0.0)
    };
    const auto pade = fem::fastsweep::PadeApproximant::build(moments, 0, 1);
    require(pade.numeratorOrder() == 0, "Pade numerator order mismatch");
    require(pade.denominatorOrder() == 1, "Pade denominator order mismatch");
    requireNear(pade.denominator()[1], Complex(-2.0, 0.0), 1.0e-14,
                "Pade denominator coefficient mismatch");
    requireNear(pade.evaluate(Complex(0.1, 0.0)), Complex(1.25, 0.0), 1.0e-14,
                "Pade evaluation mismatch");
}

// 验证二阶双边 Lanczos 在完整二维空间中精确恢复传递函数和状态。
void testLanczosPadeModel() {
    using Vector = std::vector<Complex>;
    auto applyRight = [](const Vector& x) {
        return Vector{
            Complex(0.2, 0.0) * x[0] + Complex(0.1, 0.0) * x[1],
            Complex(-0.4, 0.0) * x[0] + Complex(0.3, 0.0) * x[1]
        };
    };
    auto applyTranspose = [](const Vector& x) {
        return Vector{
            Complex(0.2, 0.0) * x[0] + Complex(-0.4, 0.0) * x[1],
            Complex(0.1, 0.0) * x[0] + Complex(0.3, 0.0) * x[1]
        };
    };
    const Vector rightStart = {Complex(1.0, 0.0), Complex(2.0, 0.0)};
    const Vector output = {Complex(0.7, 0.0), Complex(-0.4, 0.0)};

    fem::fastsweep::LanczosPadeModel model;
    model.build(rightStart, output, applyRight, applyTranspose, 2, 1.0e-13);
    const Complex parameter(0.35, 0.0);
    Vector exact = rightStart;
    std::vector<Complex> system = {
        Complex(1.0 - 0.35 * 0.2, 0.0), Complex(-0.35 * 0.1, 0.0),
        Complex(0.35 * 0.4, 0.0), Complex(1.0 - 0.35 * 0.3, 0.0)
    };
    require(fem::fastsweep::denseSolve(system, exact, 2),
            "Lanczos reference solve failed");
    const Complex expected = output[0] * exact[0] + output[1] * exact[1];
    requireNear(model.evaluate(parameter), expected, 1.0e-11,
                "Lanczos-Pade transfer mismatch");
    const auto reconstructed = model.reconstruct(parameter);
    requireNear(reconstructed[0], exact[0], 1.0e-11,
                "Lanczos-Pade state component 0 mismatch");
    requireNear(reconstructed[1], exact[1], 1.0e-11,
                "Lanczos-Pade state component 1 mismatch");
    require(model.dimension() == 2, "Lanczos-Pade dimension mismatch");
    require(model.biorthogonalityError() < 1.0e-11,
            "Lanczos biorthogonality error too large");
}

// 验证标准三项递推的双正交性、投影三对角性和前 2q 阶矩匹配。
void testStandardBilateralLanczosInvariants() {
    using Vector = std::vector<Complex>;
    const std::vector<Complex> matrix = {
        {0.20, 0.01}, {0.10, -0.02}, {-0.05, 0.01}, {0.03, 0.00},
        {-0.40, 0.03}, {0.30, 0.00}, {0.08, -0.01}, {-0.02, 0.01},
        {0.07, 0.00}, {-0.09, 0.02}, {0.25, -0.02}, {0.11, 0.00},
        {-0.01, 0.01}, {0.06, 0.00}, {-0.12, 0.03}, {0.18, -0.01}
    };
    auto applyRight = [&](const Vector& x) {
        Vector out(4, Complex(0.0, 0.0));
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                out[static_cast<std::size_t>(i)] +=
                    matrix[static_cast<std::size_t>(i * 4 + j)]
                    * x[static_cast<std::size_t>(j)];
            }
        }
        return out;
    };
    auto applyTranspose = [&](const Vector& x) {
        Vector out(4, Complex(0.0, 0.0));
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                out[static_cast<std::size_t>(i)] +=
                    matrix[static_cast<std::size_t>(j * 4 + i)]
                    * x[static_cast<std::size_t>(j)];
            }
        }
        return out;
    };
    const Vector rightStart = {{1.0, 0.1}, {2.0, -0.2}, {-1.0, 0.3}, {0.5, 0.0}};
    const Vector output = {{0.7, 0.0}, {-0.4, 0.1}, {0.2, -0.1}, {0.3, 0.0}};

    fem::fastsweep::LanczosPadeModel model;
    constexpr int q = 3;
    model.build(rightStart, output, applyRight, applyTranspose, q,
                1.0e-13, true, true);
    require(model.dimension() == q, "Standard Lanczos dimension mismatch");
    require(model.biorthogonalityError() < 1.0e-11,
            "Standard Lanczos W^T V error too large");
    require(model.recurrenceResidual() < 1.0e-11,
            "Standard Lanczos recurrence residual too large");
    require(model.validateMomentMatching(applyRight, 2 * q) < 1.0e-11,
            "Standard Lanczos 2q moment matching failed");

    const auto& rightBasis = model.rightBasis();
    const auto& leftBasis = model.leftBasis();
    const auto& reduced = model.reducedOperator();
    for (int i = 0; i < q; ++i) {
        for (int j = 0; j < q; ++j) {
            const Complex projected = fem::fastsweep::bilinear(
                leftBasis[static_cast<std::size_t>(i)],
                applyRight(rightBasis[static_cast<std::size_t>(j)]));
            requireNear(projected, reduced[static_cast<std::size_t>(i * q + j)],
                        1.0e-11, "Standard Lanczos T=W^TGV mismatch");
        }
    }
#ifdef BPFEM_USE_MKL
    require(model.poleResidueReady(), "Pole-residue decomposition was not enabled");
    require(model.poleResidueReconstructionError() < 1.0e-12,
            "Pole-residue reconstruction error too large");
#endif
}

// 验证 happy breakdown、失败 look-ahead 与成功最小 look-ahead 的状态记录。
void testLanczosBreakdownPaths() {
    using Vector = std::vector<Complex>;
    const Vector start = {Complex(1.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)};
    const Vector output = start;

    auto identity = [](const Vector& x) { return x; };
    fem::fastsweep::LanczosPadeModel happy;
    happy.build(start, output, identity, identity, 3, 1.0e-13);
    require(happy.breakdownDetected() && happy.happyBreakdown(),
            "Happy breakdown was not identified");
    require(happy.dimension() == 1, "Happy breakdown retained an invalid column");

    const std::vector<Complex> failedMatrix = {
        0.0, 0.0, 1.0,
        1.0, 0.0, 0.0,
        0.0, 0.0, 0.0
    };
    auto multiply = [](const std::vector<Complex>& matrix,
                       const Vector& x,
                       bool transpose) {
        Vector out(3, Complex(0.0, 0.0));
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                const int index = transpose ? j * 3 + i : i * 3 + j;
                out[static_cast<std::size_t>(i)] +=
                    matrix[static_cast<std::size_t>(index)] * x[static_cast<std::size_t>(j)];
            }
        }
        return out;
    };
    fem::fastsweep::LanczosPadeModel failed;
    failed.build(start, output,
                 [&](const Vector& x) { return multiply(failedMatrix, x, false); },
                 [&](const Vector& x) { return multiply(failedMatrix, x, true); },
                 3, 1.0e-13);
    require(failed.breakdownDetected() && !failed.happyBreakdown(),
            "Ordinary coupling breakdown was not identified");
    require(failed.lookAheadCount() == 1, "Failed look-ahead count mismatch");

    const std::vector<Complex> cycleMatrix = {
        0.0, 0.0, 1.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0
    };
    fem::fastsweep::LanczosPadeModel recovered;
    recovered.build(start, output,
                    [&](const Vector& x) { return multiply(cycleMatrix, x, false); },
                    [&](const Vector& x) { return multiply(cycleMatrix, x, true); },
                    2, 1.0e-13);
    require(recovered.lookAheadCount() == 1, "Successful look-ahead was not recorded");
    require(recovered.dimension() == 2, "Look-ahead did not recover a second basis pair");
}

void testWellConditionedBasisBuilder() {
    fem::fastsweep::WellConditionedBasisBuilder builder;
    builder.reset(2, 1.0e-14, 2);
    const auto first = builder.appendCandidate(
        {Complex(2.0, 0.0), Complex(0.0, 0.0)});
    const auto second = builder.appendCandidate(
        {Complex(2.0, 0.0), Complex(1.0e-8, 0.0)});
    require(first.accepted, "WCAWE builder rejected first candidate");
    require(second.accepted, "WCAWE builder rejected independent candidate");
    require(builder.dimension() == 2, "WCAWE builder dimension mismatch");
    require(builder.upperTriangularU().size() == 4,
            "WCAWE U does not use the requested fixed stride");
    requireNear(builder.upperTriangularU()[2], Complex(0.0, 0.0), 0.0,
                "WCAWE U has a nonzero lower-triangular entry");
    require(second.orthogonalityError < 1.0e-10,
            "WCAWE orthogonality error too large");
    require(second.basisRelationResidual < 1.0e-12,
            "WCAWE candidate=V*U relation error too large");
}

// 计算小型行主序稠密矩阵乘积，作为 P_U 显式参考实现。
std::vector<Complex> multiplyDenseMatrices(const std::vector<Complex>& left,
                                           const std::vector<Complex>& right,
                                           int size) {
    std::vector<Complex> product(
        static_cast<std::size_t>(size * size), Complex(0.0, 0.0));
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            for (int inner = 0; inner < size; ++inner) {
                product[static_cast<std::size_t>(row * size + column)] +=
                    left[static_cast<std::size_t>(row * size + inner)]
                    * right[static_cast<std::size_t>(inner * size + column)];
            }
        }
    }
    return product;
}

// 用逐列 LU 求出小型矩阵逆，仅用于校验生产代码没有显式求逆且乘积顺序正确。
std::vector<Complex> invertDenseReference(const std::vector<Complex>& matrix,
                                          int size) {
    std::vector<Complex> inverse(
        static_cast<std::size_t>(size * size), Complex(0.0, 0.0));
    for (int column = 0; column < size; ++column) {
        std::vector<Complex> work = matrix;
        std::vector<Complex> rhs(static_cast<std::size_t>(size), Complex(0.0, 0.0));
        rhs[static_cast<std::size_t>(column)] = Complex(1.0, 0.0);
        require(fem::fastsweep::denseSolve(work, rhs, size),
                "Explicit U subblock reference solve failed");
        for (int row = 0; row < size; ++row) {
            inverse[static_cast<std::size_t>(row * size + column)] =
                rhs[static_cast<std::size_t>(row)];
        }
    }
    return inverse;
}

// 验证式 (9) 的非单位复杂 U 子块乘积，重点覆盖逆序回代。
void testUpperTriangularBlockAction() {
    constexpr int stride = 5;
    std::vector<Complex> upper(25, Complex(0.0, 0.0));
    for (int row = 0; row < stride; ++row) {
        for (int column = row; column < stride; ++column) {
            upper[static_cast<std::size_t>(row * stride + column)] =
                row == column
                    ? Complex(1.5 + 0.4 * row, 0.1 * (row + 1))
                    : Complex(0.07 * (row + 1) * (column + 2),
                              -0.03 * (column - row));
        }
    }

    constexpr int n = 5;
    constexpr int m = 3;
    constexpr int w = 1;
    constexpr int blockSize = n - m;
    std::vector<Complex> product = {
        Complex(1.0, 0.0), Complex(0.0, 0.0),
        Complex(0.0, 0.0), Complex(1.0, 0.0)
    };
    for (int t = w; t <= m; ++t) {
        std::vector<Complex> block(4, Complex(0.0, 0.0));
        for (int row = 0; row < blockSize; ++row) {
            for (int column = 0; column < blockSize; ++column) {
                block[static_cast<std::size_t>(row * blockSize + column)] =
                    upper[static_cast<std::size_t>((t - 1 + row) * stride
                                                   + t - 1 + column)];
            }
        }
        product = multiplyDenseMatrices(
            product, invertDenseReference(block, blockSize), blockSize);
    }

    std::uint64_t solveCount = 0;
    const auto actual =
        fem::fastsweep::UpperTriangularBlockAction::applyToLastUnitVector(
            upper, stride, n - 1, n, m, w, 1.0e-14, &solveCount);
    requireNear(actual[0], product[1], 1.0e-13,
                "P_U first component mismatch");
    requireNear(actual[1], product[3], 1.0e-13,
                "P_U second component mismatch");
    require(solveCount == 3, "P_U triangular solve count mismatch");
}

// 计算小型行主序复矩阵与向量的乘积，供通用多项式 WCAWE 单测使用。
std::vector<Complex> applyDenseMatrix(const std::vector<Complex>& matrix,
                                      const std::vector<Complex>& vector,
                                      int size) {
    std::vector<Complex> result(static_cast<std::size_t>(size), Complex(0.0, 0.0));
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            result[static_cast<std::size_t>(row)] +=
                matrix[static_cast<std::size_t>(row * size + column)]
                * vector[static_cast<std::size_t>(column)];
        }
    }
    return result;
}

// 计算向量相对正交基张成空间的投影残差。
double projectionResidual(
    const std::vector<Complex>& vector,
    const std::vector<std::vector<Complex>>& basis) {
    std::vector<Complex> residual(vector);
    for (const auto& basisVector : basis) {
        const Complex coefficient = fem::fastsweep::hdot(basisVector, vector);
        for (std::size_t row = 0; row < residual.size(); ++row) {
            residual[row] -= coefficient * basisVector[row];
        }
    }
    return fem::fastsweep::norm2(residual)
        / std::max(fem::fastsweep::norm2(vector), 1.0e-300);
}

// 验证论文式 (7)、式 (8)、U=I 退化以及真正 WCAWE 的矩子空间匹配。
void testPolynomialWcaweRecurrence() {
    constexpr int size = 6;
    constexpr int order = 4;
    std::vector<Complex> a0(36, Complex(0.0, 0.0));
    std::vector<Complex> a1(36, Complex(0.0, 0.0));
    std::vector<Complex> a2(36, Complex(0.0, 0.0));
    for (int row = 0; row < size; ++row) {
        a0[static_cast<std::size_t>(row * size + row)] =
            Complex(2.0 + 0.3 * row, 0.05 * row);
        a0[static_cast<std::size_t>(row * size + (row + 1) % size)] +=
            Complex(0.08, -0.02);
        a1[static_cast<std::size_t>(((row + 1) % size) * size + row)] =
            Complex(0.35 + 0.03 * row, 0.02);
        a1[static_cast<std::size_t>(row * size + row)] +=
            Complex(-0.07 * (row + 1), 0.01);
        a2[static_cast<std::size_t>(((row + 2) % size) * size + row)] =
            Complex(0.04 + 0.01 * row, -0.015);
    }
    const std::vector<Complex> b0 = {
        {2.0, 0.2}, {0.7, -0.1}, {-0.4, 0.3},
        {1.1, 0.0}, {0.2, -0.5}, {-0.8, 0.1}
    };
    const std::vector<Complex> b1 = {
        {0.1, -0.2}, {0.5, 0.1}, {0.3, -0.4},
        {-0.2, 0.2}, {0.7, 0.0}, {0.4, -0.1}
    };
    const std::vector<Complex> b2 = {
        {-0.05, 0.1}, {0.2, -0.1}, {0.15, 0.05},
        {0.3, -0.2}, {-0.25, 0.1}, {0.1, 0.2}
    };
    auto rhsAt = [&](std::size_t coefficientOrder) {
        if (coefficientOrder == 0) return b0;
        if (coefficientOrder == 1) return b1;
        if (coefficientOrder == 2) return b2;
        return std::vector<Complex>(size, Complex(0.0, 0.0));
    };
    auto solve = [&](const std::vector<Complex>& rhs) {
        std::vector<Complex> matrix = a0;
        std::vector<Complex> solution = rhs;
        require(fem::fastsweep::denseSolve(matrix, solution, size),
                "Small WCAWE A0 solve failed");
        return solution;
    };
    std::vector<fem::fastsweep::WcaweLinearOperator> operators;
    operators.push_back([&](const std::vector<Complex>& vector) {
        return applyDenseMatrix(a1, vector, size);
    });
    operators.push_back([&](const std::vector<Complex>& vector) {
        return applyDenseMatrix(a2, vector, size);
    });

    fem::fastsweep::PolynomialWcaweModel model;
    model.fullDimension = size;
    model.matrixDegree = 2;
    model.rhsDegree = 2;
    model.applyA0 = [&](const std::vector<Complex>& vector) {
        return applyDenseMatrix(a0, vector, size);
    };
    model.applyHigherOrder = operators;
    model.rhsCoefficientAt = rhsAt;
    fem::fastsweep::WcaweBuildOptions options;
    options.requestedOrder = order;
    options.breakdownTolerance = 1.0e-13;
    const auto result =
        fem::fastsweep::PolynomialWcaweRecurrence::build(model, solve, options);
    require(result.achievedOrder == order,
            "Polynomial WCAWE did not reach the requested order");
    require(result.maxRecurrenceResidual < 1.0e-11,
            "WCAWE equation (7) residual too large");
    require(result.maxBasisRelationResidual < 1.0e-12,
            "WCAWE equation (8) residual too large");
    require(result.maxOrthogonalityError < 1.0e-11,
            "WCAWE basis orthogonality error too large");

    std::vector<fem::fastsweep::PolynomialMomentRecurrence::LinearOperator>
        momentOperators = operators;
    const auto moments =
        fem::fastsweep::PolynomialMomentRecurrence::generatePolynomialMoments(
            order, rhsAt, momentOperators, solve);
    for (const auto& moment : moments) {
        require(projectionResidual(moment, result.basis) < 1.0e-10,
                "WCAWE basis failed traditional moment subspace matching");
    }

    std::vector<Complex> identityU(
        static_cast<std::size_t>(order * order), Complex(0.0, 0.0));
    for (int index = 0; index < order; ++index) {
        identityU[static_cast<std::size_t>(index * order + index)] =
            Complex(1.0, 0.0);
    }
    for (int n = 1; n <= order; ++n) {
        const std::vector<std::vector<Complex>> history(
            moments.begin(), moments.begin() + n - 1);
        const auto rhs =
            fem::fastsweep::PolynomialWcaweRecurrence::assembleCorrectedRhs(
                model, history, identityU, order, n, 1.0e-14);
        const auto identityCandidate = solve(rhs);
        std::vector<Complex> difference(identityCandidate);
        for (std::size_t row = 0; row < difference.size(); ++row) {
            difference[row] -= moments[static_cast<std::size_t>(n - 1)][row];
        }
        require(fem::fastsweep::norm2(difference)
                    / std::max(fem::fastsweep::norm2(identityCandidate), 1.0e-300)
                < 1.0e-12,
                "WCAWE U=I did not recover traditional AWE");
    }
}

// 用论文第 III-A 节对应的反例证明“归一化历史向量直接递推”不能替代 U 校正。
void testWcaweRejectsNaiveNormalizedRecurrence() {
    constexpr int size = 3;
    const std::vector<Complex> b0 = {
        Complex(2.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)
    };
    const std::vector<Complex> b1 = {
        Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(1.0, 0.0)
    };
    auto applyIdentity = [](const std::vector<Complex>& vector) { return vector; };
    auto applyA1 = [](const std::vector<Complex>& vector) {
        return std::vector<Complex>{
            Complex(0.0, 0.0), vector[0], Complex(0.0, 0.0)
        };
    };
    auto rhsAt = [&](std::size_t order) {
        if (order == 0) return b0;
        if (order == 1) return b1;
        return std::vector<Complex>(size, Complex(0.0, 0.0));
    };

    fem::fastsweep::PolynomialWcaweModel model;
    model.fullDimension = size;
    model.matrixDegree = 1;
    model.rhsDegree = 1;
    model.applyA0 = applyIdentity;
    model.applyHigherOrder = {applyA1};
    model.rhsCoefficientAt = rhsAt;
    fem::fastsweep::WcaweBuildOptions options;
    options.requestedOrder = 2;
    options.breakdownTolerance = 1.0e-14;
    const auto correct =
        fem::fastsweep::PolynomialWcaweRecurrence::build(
            model, applyIdentity, options);

    const std::vector<Complex> trueSecondMoment = {
        Complex(0.0, 0.0), Complex(-2.0, 0.0), Complex(1.0, 0.0)
    };
    require(projectionResidual(trueSecondMoment, correct.basis) < 1.0e-12,
            "Correct WCAWE failed the naive-recurrence counterexample");

    fem::fastsweep::WellConditionedBasisBuilder naive;
    naive.reset(2, 1.0e-14, 2);
    require(naive.appendCandidate(b0).accepted,
            "Naive counterexample rejected its first vector");
    std::vector<Complex> naiveRhs = b1;
    const auto a1v1 = applyA1(naive.basis().front());
    for (int row = 0; row < size; ++row) {
        naiveRhs[static_cast<std::size_t>(row)] -=
            a1v1[static_cast<std::size_t>(row)];
    }
    require(naive.appendCandidate(naiveRhs).accepted,
            "Naive counterexample rejected its second vector");
    require(projectionResidual(trueSecondMoment, naive.basis()) > 0.1,
            "Counterexample cannot distinguish naive normalized recurrence");
}

// 验证 U 对角 breakdown 后立即停止，不允许跳过失败列继续访问下一阶 P_U。
void testWcaweBreakdownStopsRecurrence() {
    int solveCount = 0;
    fem::fastsweep::PolynomialWcaweModel model;
    model.fullDimension = 3;
    model.matrixDegree = 1;
    model.rhsDegree = 0;
    model.applyA0 = [](const std::vector<Complex>& vector) { return vector; };
    model.applyHigherOrder = {
        [](const std::vector<Complex>& vector) { return vector; }
    };
    model.rhsCoefficientAt = [](std::size_t order) {
        return order == 0
            ? std::vector<Complex>{Complex(1.0, 0.0), Complex(0.0, 0.0),
                                   Complex(0.0, 0.0)}
            : std::vector<Complex>(3, Complex(0.0, 0.0));
    };
    fem::fastsweep::WcaweBuildOptions options;
    options.requestedOrder = 4;
    options.breakdownTolerance = 1.0e-13;
    const auto result = fem::fastsweep::PolynomialWcaweRecurrence::build(
        model,
        [&](const std::vector<Complex>& rhs) {
            ++solveCount;
            return rhs;
        },
        options);
    require(result.achievedOrder == 1,
            "WCAWE breakdown retained a dependent candidate");
    require(result.terminationReason
                == fem::fastsweep::WcaweTerminationReason::
                    NearSingularUpperTriangularFactor,
            "WCAWE breakdown termination reason mismatch");
    require(solveCount == 2,
            "WCAWE continued solving after the first breakdown");
}

fem::PortMode makePortMode(int faceId, int dof) {
    fem::PortMode mode;
    mode.faceId = faceId;
    mode.cutoffWavenumberSquared = 0.0;
    mode.couplingWeights.push_back({dof, 1.0});
    fem::PortQuadraturePoint qp;
    qp.modeFieldValue = {1.0, 0.0, 0.0};
    qp.normal = {0.0, 0.0, 1.0};
    qp.weight = 1.0;
    mode.quadrature.push_back(qp);
    return mode;
}

void testGalerkinReducedModel() {
    fem::ProjectDefinition project;
    project.ports.push_back({1, 10, 1, true, 1.0, 0.0});
    project.ports.push_back({2, 20, 1, false, 0.0, 0.0});

    fem::FEMAssembler::AffineSystem affine;
    affine.K = fem::SparseMatrix(2);
    affine.K.add(0, 0, Complex(1.0, 0.0));
    affine.K.add(1, 1, Complex(1.0, 0.0));
    affine.M = fem::SparseMatrix(2);
    affine.portCoupling = {{{0, 1.0}}, {{1, 1.0}}};
    affine.portCutoffSquared = {0.0, 0.0};
    affine.portFaceIds = {10, 20};
    affine.projectPortIndex = {0, 1};
    affine.isExcitationMode = {true, false};
    affine.lossless = true;

    fem::Mesh mesh;
    fem::EdgeTopology topology(mesh, 0);
    fem::PortModeSolver portModeSolver(mesh, topology);
    portModeSolver.setPrecomputed(10, makePortMode(10, 0));
    portModeSolver.setPrecomputed(20, makePortMode(20, 1));

    std::vector<std::vector<Complex>> basis = {
        {Complex(1.0, 0.0), Complex(0.0, 0.0)},
        {Complex(0.0, 0.0), Complex(1.0, 0.0)}
    };
    std::vector<std::vector<Complex>> portVectors = {
        {Complex(1.0, 0.0), Complex(0.0, 0.0)},
        {Complex(0.0, 0.0), Complex(1.0, 0.0)}
    };

    fem::fastsweep::GalerkinReducedModel model;
    model.build(project, affine, portModeSolver, std::move(basis), portVectors);
    require(model.dimension() == 2, "Galerkin ROM dimension mismatch");
    require(model.basisOrthogonalityError() < 1.0e-14,
            "Galerkin basis orthogonality error too large");
    const auto sp = model.evaluate(1.0e9);
    require(std::isfinite(sp.s11.real()) && std::isfinite(sp.s11.imag()),
            "Galerkin S11 is not finite");
    require(std::isfinite(sp.s21.real()) && std::isfinite(sp.s21.imag()),
            "Galerkin S21 is not finite");
    const auto field = model.reconstructField(1.0e9);
    require(field.size() == 2, "Galerkin reconstructed field size mismatch");
}

}  // namespace

int main() {
    testFactorizedSolveSessionCounts();
#ifdef BPFEM_USE_MKL
    testPardisoFactorizedLifecycle();
#endif
    testPolynomialMomentRecurrence();
    testLosslessP1FiniteDifference();
    testPadeApproximant();
    testLanczosPadeModel();
    testStandardBilateralLanczosInvariants();
    testLanczosBreakdownPaths();
    testWellConditionedBasisBuilder();
    testUpperTriangularBlockAction();
    testPolynomialWcaweRecurrence();
    testWcaweRejectsNaiveNormalizedRecurrence();
    testWcaweBreakdownStopsRecurrence();
    testGalerkinReducedModel();
    return 0;
}
