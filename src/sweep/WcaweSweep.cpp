#include "bpfem/sweep/WcaweSweep.hpp"

#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"
#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"
#include "bpfem/linalg/FactorizedSolveSession.hpp"
#include "bpfem/linalg/IFactorizedSparseSolver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem::sweep {

// 构造 WCAWE 扫频器，并把阶数和 breakdown 阈值规整到可用范围。
WcaweSweep::WcaweSweep(const ProjectDefinition& project,
                       const FEMAssembler& assembler,
                       const PortModeSolver& portModeSolver,
                       WcaweOptions options)
    : project_(project),
      assembler_(assembler),
      portModeSolver_(portModeSolver),
      options_(options) {
    options_.order = std::max(1, options_.order);
    if (!std::isfinite(options_.dropTolerance) || options_.dropTolerance <= 0.0) {
        options_.dropTolerance = 1.0e-12;
    }
}

// 执行论文一致的 WCAWE 扫频，并写出递推、U、正交性和 ROM 诊断。
SweepResult WcaweSweep::run(const std::vector<double>& frequencies,
                            const SweepContext& ctx) {
    const double expansionHz = options_.expansionFrequencyHz > 0.0
        ? options_.expansionFrequencyHz
        : (frequencies.empty()
            ? 0.5 * (project_.sweep.startHz + project_.sweep.endHz)
            : 0.5 * (frequencies.front() + frequencies.back()));

    ctx.log.info("WCAWE expansion frequency: "
                 + std::to_string(expansionHz / 1.0e9) + " GHz");
    ctx.log.info("WCAWE target basis order: " + std::to_string(options_.order));
    std::ostringstream dropTolerance;
    dropTolerance << std::scientific << std::setprecision(3)
                  << options_.dropTolerance;
    ctx.log.info("WCAWE U diagonal breakdown tolerance: " + dropTolerance.str());

    linalg::SolverConfig solverConfig;
    solverConfig.maxIterations = ctx.linearMaxIterations;
    solverConfig.tolerance = ctx.linearTolerance;
    const auto offlineStarted = std::chrono::steady_clock::now();
    const int dimension = buildOffline(expansionHz, ctx.solver, solverConfig);
    const double offlineBuildSec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - offlineStarted).count();
    ctx.log.info("WCAWE ROM dimension: " + std::to_string(dimension));
    ctx.log.info("WCAWE termination: "
                 + std::string(fastsweep::wcaweTerminationReasonName(
                     wcaweBuild_.terminationReason)));
    if (!wcaweBuild_.terminationMessage.empty()) {
        ctx.log.warn("WCAWE stopped early: " + wcaweBuild_.terminationMessage);
    }

    if (!ctx.outputDirectory.empty()) {
        const auto conditionPath = ctx.outputDirectory / "basis_condition.csv";
        if (writeBasisConditionCsv(conditionPath)) {
            ctx.log.info("Wrote " + conditionPath.string());
        } else {
            ctx.log.warn("Failed to write " + conditionPath.string());
        }
    }

    SweepResult output;
    output.offlineBuildSec = offlineBuildSec;
    output.orthogonalizationSec = orthogonalizationSec_;
    output.romProjectionSec = romProjectionSec_;
    output.points.reserve(frequencies.size());
    const auto onlineStarted = std::chrono::steady_clock::now();
    for (std::size_t index = 0; index < frequencies.size(); ++index) {
        const double frequencyHz = frequencies[index];
        if (index == 0 || (index + 1) % 50 == 0
            || index + 1 == frequencies.size()) {
            ctx.log.info("WCAWE evaluating point "
                         + std::to_string(index + 1) + "/"
                         + std::to_string(frequencies.size()) + " @ "
                         + std::to_string(frequencyHz / 1.0e9) + " GHz");
        }
        output.points.push_back(evaluate(frequencyHz));
        output.lastFrequencyHz = frequencyHz;
    }
    output.onlineSweepSec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - onlineStarted).count();
    if (!frequencies.empty()) {
        output.lastEdgeDofs = reconstructField(frequencies.back());
    }

    if (!ctx.outputDirectory.empty()) {
        fastsweep::FastSweepDiagnostics diagnostics;
        diagnostics.algorithm = name();
        diagnostics.expansionFrequenciesHz = {expansionHz};
        diagnostics.requestedOrder = options_.order;
        diagnostics.romDimension = romDim_;
        diagnostics.retainedColumns = romDim_;
        diagnostics.deflatedColumns = deflatedColumns_;
        diagnostics.basisOrthogonalityError = model_.basisOrthogonalityError();
        diagnostics.wcaweMomentReconstructionError =
            wcaweBuild_.maxBasisRelationResidual;
        diagnostics.wcaweRecurrenceResidualMax =
            wcaweBuild_.maxRecurrenceResidual;
        diagnostics.wcaweBasisRelationResidualMax =
            wcaweBuild_.maxBasisRelationResidual;
        diagnostics.wcaweOrthogonalityError =
            wcaweBuild_.maxOrthogonalityError;
        diagnostics.wcaweMinUpperDiagonalAbs =
            wcaweBuild_.minUpperDiagonalAbs;
        diagnostics.wcaweUpperDiagonalRatio =
            wcaweBuild_.upperDiagonalRatio;
        diagnostics.wcaweTriangularSolveCount =
            wcaweBuild_.triangularSolveCount;
        diagnostics.wcaweBreakdownOrder =
            wcaweBuild_.achievedOrder < wcaweBuild_.requestedOrder
                ? wcaweBuild_.achievedOrder + 1
                : 0;
        diagnostics.wcaweTerminationReason =
            fastsweep::wcaweTerminationReasonName(
                wcaweBuild_.terminationReason);
        diagnostics.wcaweMomentMatchingError =
            wcaweBuild_.momentMatchingError;
        if (const auto* factorized =
                dynamic_cast<const linalg::IFactorizedSparseSolver*>(&ctx.solver)) {
            const auto statistics = factorized->factorizationStatistics();
            diagnostics.factorizedRhsSolveCount =
                statistics.factorizedRhsSolveCount;
            diagnostics.factorizedSolveCallCount =
                statistics.factorizedSolveCallCount;
            diagnostics.batchRhsMax = statistics.batchRhsMax;
        }
        diagnostics.reducedSolveSucceeded = model_.lastSolveSucceeded();
        diagnostics.maxPassivityError =
            fastsweep::maxPassivityError(output.points);
        const auto diagnosticsPath =
            ctx.outputDirectory / "diagnostics.json";
        if (fastsweep::writeDiagnosticsJson(diagnosticsPath, diagnostics)) {
            ctx.log.info("Wrote " + diagnosticsPath.string());
        } else {
            ctx.log.warn("Failed to write " + diagnosticsPath.string());
        }
    }
    (void)ctx.onFieldSolved;
    return output;
}

// 将工程 P1 线性化适配成通用多项式模型，并按式 (7)-(9)构造 WCAWE 基。
int WcaweSweep::buildOffline(double expansionFrequencyHz,
                             linalg::ISparseSolver& solver,
                             const linalg::SolverConfig& solverConfig) {
    if (project_.ports.size() < 2) {
        throw std::runtime_error(
            "WcaweSweep: at least two project ports are required");
    }

    ready_ = false;
    expansionFrequencyHz_ = expansionFrequencyHz;
    numProjectPorts_ = static_cast<int>(project_.ports.size());
    affine_ = assembler_.buildAffineSystem();
    fullDim_ = static_cast<int>(affine_.K.size());
    numVirtualPorts_ = static_cast<int>(affine_.portCoupling.size());
    if (numVirtualPorts_ == 0) {
        throw std::runtime_error(
            "WcaweSweep: AffineSystem produced zero port entries");
    }
    if (!affine_.lossless) {
        throw std::runtime_error(
            "WcaweSweep: lossy materials are not supported. Use --sweep direct.");
    }

    const auto portVectors =
        fastsweep::PolynomialPortMomentBuilder::buildPortVectors(
            affine_, affine_.K.size());
    const auto linearization =
        fastsweep::PolynomialPortMomentBuilder::buildLosslessLinearization(
            project_, assembler_, portModeSolver_, affine_,
            expansionFrequencyHz, portVectors);

    fastsweep::PolynomialWcaweModel polynomial;
    polynomial.fullDimension = linearization.rhs0.size();
    polynomial.matrixDegree = 1;
    polynomial.rhsDegree = 1;
    polynomial.applyA0 =
        [&linearization](const std::vector<Complex>& vector) {
            return linearization.matrixAtExpansion.multiply(vector);
        };
    polynomial.applyHigherOrder.push_back(
        [this, &linearization, &portVectors](
            const std::vector<Complex>& vector) {
            return fastsweep::PolynomialPortMomentBuilder::
                applyLosslessFirstOrderMatrix(
                    affine_, linearization, portVectors, vector);
        });
    polynomial.rhsCoefficientAt =
        [&linearization](std::size_t order) {
            if (order == 0) {
                return linearization.rhs0;
            }
            if (order == 1) {
                return linearization.rhs1;
            }
            return std::vector<Complex>(
                linearization.rhs0.size(), Complex(0.0, 0.0));
        };

    linalg::FactorizedSolveSession solveSession(
        solver, linearization.matrixAtExpansion, solverConfig);
    auto solveAtExpansion =
        [&solveSession](const std::vector<Complex>& rhs) {
            SolveResult result = solveSession.solve(rhs);
            if (result.field.size() != rhs.size()
                || !std::isfinite(result.residual)) {
                throw std::runtime_error(
                    "WcaweSweep: expansion solve returned an invalid result");
            }
            return std::move(result.field);
        };

    fastsweep::WcaweBuildOptions buildOptions;
    buildOptions.requestedOrder = options_.order;
    buildOptions.breakdownTolerance = options_.dropTolerance;
    buildOptions.recurrenceTolerance = 1.0e-10;
    buildOptions.orthogonalityTolerance = 1.0e-10;
    buildOptions.reorthogonalizationPasses = 2;
    wcaweBuild_ = fastsweep::PolynomialWcaweRecurrence::build(
        polynomial, solveAtExpansion, buildOptions);
    orthogonalizationSec_ = wcaweBuild_.orthogonalizationSec;

    romDim_ = wcaweBuild_.achievedOrder;
    deflatedColumns_ = options_.order - romDim_;
    if (romDim_ == 0) {
        throw std::runtime_error(
            "WcaweSweep: paper WCAWE recurrence produced no usable basis");
    }

    const auto projectionStarted = std::chrono::steady_clock::now();
    model_.build(project_, affine_, portModeSolver_,
                 std::move(wcaweBuild_.basis), portVectors);
    romProjectionSec_ = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - projectionStarted).count();
    ready_ = true;
    return romDim_;
}

// 在降阶空间内求解指定频率的小型 ROM 线性系统。
std::vector<WcaweSweep::Complex> WcaweSweep::solveReduced(
    double frequencyHz) const {
    return model_.solveReduced(frequencyHz);
}

// 由 ROM 解提取指定频率的二端口 S 参数。
SParameterPoint WcaweSweep::evaluate(double frequencyHz) const {
    return model_.evaluate(frequencyHz);
}

// 将指定频率的 ROM 解提升回全阶边自由度，用于最后一个频点的场输出。
std::vector<WcaweSweep::Complex> WcaweSweep::reconstructField(
    double frequencyHz) const {
    return model_.reconstructField(frequencyHz);
}

// 写出每一阶的条件代理、U 对角元以及论文式 (7)、式 (8)残差。
bool WcaweSweep::writeBasisConditionCsv(
    const std::filesystem::path& path) const {
    std::ofstream output(path, std::ios::out | std::ios::trunc);
    if (!output.good()) {
        return false;
    }
    output
        << "order,awe_moment_condition_proxy,wcawe_basis_condition_proxy,"
           "u_diagonal_abs,u_diagonal_ratio,recurrence_residual,"
           "basis_relation_residual,orthogonality_error,"
           "triangular_solve_count\n";
    output << std::setprecision(17);
    for (const auto& step : wcaweBuild_.steps) {
        output << step.order << ','
               << step.diagonalRatio << ','
               << 1.0 + step.orthogonalityError << ','
               << step.diagonalAbs << ','
               << step.diagonalRatio << ','
               << step.recurrenceResidual << ','
               << step.basisRelationResidual << ','
               << step.orthogonalityError << ','
               << step.triangularSolveCount << '\n';
    }
    output.flush();
    return output.good();
}

}  // namespace fem::sweep
