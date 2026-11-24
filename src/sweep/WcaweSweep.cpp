#include "bpfem/sweep/WcaweSweep.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"
#include "bpfem/fastsweep/GalerkinReducedModel.hpp"
#include "bpfem/fastsweep/PolynomialMomentRecurrence.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PortModeUtilities.hpp"
#include "bpfem/linalg/FactorizedSolveSession.hpp"
#include "bpfem/linalg/IFactorizedSparseSolver.hpp"
#include "bpfem/post/ResultExtractor.hpp"

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
        if (!writeBasisConditionCsv(conditionPath)) {
            throw std::runtime_error("Failed to write required WCAWE output " + conditionPath.string());
        }
        ctx.log.info("Wrote " + conditionPath.string());
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
        diagnostics.model = "wavenumber_taylor_polynomial";
        diagnostics.expansionFrequenciesHz = {expansionHz};
        diagnostics.requestedOrder = options_.order;
        diagnostics.romDimension = romDim_;
        diagnostics.retainedColumns = romDim_;
        diagnostics.deflatedColumns = deflatedColumns_;
        diagnostics.basisOrthogonalityError = wcaweBuild_.maxOrthogonalityError;
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
        diagnostics.reducedSolveSucceeded = true;
        diagnostics.fullMomentProjectionError = fullMomentProjectionError_;
        diagnostics.maxPassivityError = fastsweep::maxPassivityError(output.points);
        if (diagnosticAffine_ && diagnosticTaylor_) {
            diagnostics.residualSamples.reserve(frequencies.size());
            for (const double frequency : frequencies) {
                diagnostics.residualSamples.push_back(fastsweep::losslessResiduals(project_, portModeSolver_,
                    *diagnosticAffine_, frequency, reconstructField(frequency), &*diagnosticTaylor_));
            }
        }
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

// 将工程高阶端口多项式适配成通用多项式模型，并按式 (7)-(9)构造 WCAWE 基。
int WcaweSweep::buildOffline(double expansionFrequencyHz,
                             linalg::ISparseSolver& solver,
                             const linalg::SolverConfig& solverConfig) {
    ready_ = false;
    if (project_.ports.size() < 2) {
        throw std::runtime_error(
            "WcaweSweep: at least two project ports are required");
    }

    fullMomentProjectionError_ = -1.0;
    diagnosticAffine_.reset();
    diagnosticTaylor_.reset();
    auto affine = assembler_.buildAffineSystem();
    if (affine.portCoupling.empty()) {
        throw std::runtime_error(
            "WcaweSweep: AffineSystem produced zero port entries");
    }
    if (!affine.lossless) {
        throw std::runtime_error(
            "WcaweSweep: lossy materials are not supported. Use --sweep direct.");
    }

    const auto portVectors =
        fastsweep::PolynomialPortMomentBuilder::buildPortVectors(
            affine, affine.K.size());
    using Builder = fastsweep::PolynomialPortMomentBuilder;
    // 即使只生成一或两个基向量，在线方程仍保留 k² 的二次项。
    auto series = Builder::buildTaylorSeries(project_, portModeSolver_, affine,
        expansionFrequencyHz, std::max(3, options_.order), Builder::Variable::WavenumberShift);
    std::vector<Complex> rhs0;
    const auto a0 = assembler_.assemble(expansionFrequencyHz, rhs0);
    if (std::all_of(rhs0.begin(), rhs0.end(), [](Complex value) { return value == Complex{}; })) {
        throw std::runtime_error("WCAWE at " + std::to_string(expansionFrequencyHz)
            + " Hz: zero excitation RHS at the expansion frequency");
    }

    fastsweep::PolynomialWcaweModel polynomial;
    polynomial.fullDimension = rhs0.size();
    polynomial.matrixDegree = series.lambda.size() - 1;
    polynomial.rhsDegree = polynomial.matrixDegree;
    polynomial.applyA0 = [&a0](const std::vector<Complex>& vector) { return a0.multiply(vector); };
    for (std::size_t r = 1; r <= polynomial.matrixDegree; ++r) {
        polynomial.applyHigherOrder.push_back([&, r](const std::vector<Complex>& vector) {
            return Builder::applyMatrixCoefficient(affine, series, portVectors, r, vector);
        });
    }
    polynomial.rhsCoefficientAt = [&](std::size_t order) {
        return order == 0 ? rhs0 : Builder::rhsCoefficient(series, portVectors, rhs0.size(), order);
    };

    fastsweep::WcaweBuildOptions buildOptions;
    buildOptions.requestedOrder = options_.order;
    buildOptions.breakdownTolerance = options_.dropTolerance;
    buildOptions.recurrenceTolerance = 1.0e-10;
    buildOptions.orthogonalityTolerance = 1.0e-10;
    buildOptions.reorthogonalizationPasses = 2;

    // A backsolve accepted at a looser tolerance could fail the corrected recurrence.
    auto expansionSolveConfig = solverConfig;
    expansionSolveConfig.tolerance = std::min(solverConfig.tolerance, buildOptions.recurrenceTolerance);
    linalg::FactorizedSolveSession solveSession(solver, a0, expansionSolveConfig);
    bool firstRhs = true;
    auto solveAtExpansion =
        [&](const std::vector<Complex>& rhs) {
            const bool isFirstRhs = firstRhs;
            firstRhs = false;
            try {
                SolveResult result = solveSession.solve(rhs);
                linalg::requireConverged(result, solver.name());
                return std::move(result.field);
            } catch (const std::exception& error) {
                // The recurrence propagates the first solve failure directly.
                if (isFirstRhs) {
                    throw std::runtime_error("WCAWE at " + std::to_string(expansionFrequencyHz)
                        + " Hz, RHS 1: " + error.what());
                }
                throw;
            }
        };

    wcaweBuild_ = fastsweep::PolynomialWcaweRecurrence::build(
        polynomial, solveAtExpansion, buildOptions);
    orthogonalizationSec_ = wcaweBuild_.orthogonalizationSec;
    using Termination = fastsweep::WcaweTerminationReason;
    const auto termination = wcaweBuild_.terminationReason;
    if (termination == Termination::SolverFailure || termination == Termination::ToleranceExceeded
        || termination == Termination::NonFiniteValue) {
        throw std::runtime_error("WCAWE at " + std::to_string(expansionFrequencyHz)
            + " Hz, RHS " + std::to_string(wcaweBuild_.achievedOrder + 1)
            + ": " + wcaweBuild_.terminationMessage);
    }

    romDim_ = wcaweBuild_.achievedOrder;
    deflatedColumns_ = options_.order - romDim_;
    if (romDim_ == 0) {
        throw std::runtime_error(
            "WcaweSweep: paper WCAWE recurrence produced no usable basis");
    }

    if (options_.checkMoments) {
        const auto moments = fastsweep::PolynomialMomentRecurrence::generatePolynomialMoments(
            static_cast<std::size_t>(romDim_), polynomial.rhsCoefficientAt,
            polynomial.applyHigherOrder, solveAtExpansion);
        fullMomentProjectionError_ = 0.0;
        for (const auto& moment : moments) {
            auto remainder = moment;
            for (const auto& v : wcaweBuild_.basis) {
                const auto coefficient = fastsweep::hdot(v, moment);
                for (std::size_t i = 0; i < moment.size(); ++i) remainder[i] -= coefficient * v[i];
            }
            fullMomentProjectionError_ = std::max(fullMomentProjectionError_,
                fastsweep::norm2(remainder) / std::max(fastsweep::norm2(moment), 1.0e-300));
        }
        if (fullMomentProjectionError_ > 1.0e-8) {
            throw std::runtime_error("WCAWE failed the independent full-order moment projection check");
        }
    }
    const auto projectionStarted = std::chrono::steady_clock::now();
    projectPolynomial(affine, series, a0, rhs0);
    expansionFrequencyHz_ = series.expansionHz;
    expansionLambda_ = series.lambda0;
    if (options_.residualDiagnostics) {
        portCutoffsSquared_ = affine.portCutoffSquared;
        diagnosticAffine_ = std::move(affine);
        diagnosticTaylor_ = std::move(series);
    } else {
        portCutoffsSquared_ = std::move(affine.portCutoffSquared);
    }
    romProjectionSec_ = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - projectionStarted).count();
    ready_ = true;
    return romDim_;
}

// Algorithm 1 投影使用 V^T；生成正交基时的 MGS 内积使用 V^H。
void WcaweSweep::projectPolynomial(
    const FEMAssembler::AffineSystem& affine,
    const fastsweep::PolynomialPortMomentBuilder::PortTaylorSeries& series,
    const SparseMatrix& matrixAtExpansion,
    const std::vector<Complex>& rhsAtExpansion) {
    const auto& basis = wcaweBuild_.basis;
    const auto q = basis.size();
    const auto coefficientCount = series.lambda.size();
    reducedMatrices_.assign(coefficientCount, std::vector<Complex>(q * q));
    reducedRhs_.assign(coefficientCount, std::vector<Complex>(q));
    std::vector<Complex> reducedMass(q * q);
    for (std::size_t j = 0; j < q; ++j) {
        const auto a0v = matrixAtExpansion.multiply(basis[j]);
        const auto mv = affine.M.multiply(basis[j]);
        reducedRhs_[0][j] = fastsweep::bilinear(basis[j], rhsAtExpansion);
        for (std::size_t i = 0; i <= j; ++i) {
            const auto a0Entry = fastsweep::bilinear(basis[i], a0v);
            const auto massEntry = fastsweep::bilinear(basis[i], mv);
            reducedMatrices_[0][i * q + j] = reducedMatrices_[0][j * q + i] = a0Entry;
            reducedMass[i * q + j] = reducedMass[j * q + i] = massEntry;
        }
    }
    for (std::size_t r = 1; r < coefficientCount; ++r) {
        for (std::size_t entry = 0; entry < reducedMass.size(); ++entry) {
            reducedMatrices_[r][entry] = -series.lambda[r] * reducedMass[entry];
        }
    }

    const auto dominant = fastsweep::dominantVirtualPortsByProject(
        portModeSolver_, affine, static_cast<int>(project_.ports.size()));
    for (auto& observation : observationVectors_) {
        observation.assign(q, Complex{});
    }
    // 每个端口只投影一次，随后为各阶系数加入降阶后的秩一贡献。
    for (std::size_t port = 0; port < affine.portCoupling.size(); ++port) {
        std::vector<Complex> projectedPort(q);
        for (std::size_t i = 0; i < q; ++i) {
            for (const auto& [dof, weight] : affine.portCoupling[port]) {
                projectedPort[i] += weight * basis[i][static_cast<std::size_t>(dof)];
            }
        }
        for (std::size_t r = 1; r < coefficientCount; ++r) {
            const Complex admittance(0.0, series.admittance[port][r]);
            const Complex excitation = series.excitation[port][r];
            for (std::size_t i = 0; i < q; ++i) {
                reducedRhs_[r][i] += excitation * projectedPort[i];
                const Complex rowFactor = admittance * projectedPort[i];
                for (std::size_t j = 0; j < q; ++j) {
                    reducedMatrices_[r][i * q + j] += rowFactor * projectedPort[j];
                }
            }
        }
        for (std::size_t observed = 0; observed < observationVectors_.size(); ++observed) {
            if (dominant[observed] == static_cast<int>(port)) {
                observationVectors_[observed] = projectedPort;
            }
        }
    }
}

// 在降阶空间内求解指定频率的小型 ROM 线性系统。
std::vector<WcaweSweep::Complex> WcaweSweep::solveReduced(
    double frequencyHz) const {
    if (!ready_) {
        throw std::runtime_error("WcaweSweep::solveReduced called before build");
    }
    const double k = 2.0 * pi * frequencyHz / c0;
    const bool sameBranch = std::all_of(portCutoffsSquared_.begin(), portCutoffsSquared_.end(),
        [&](double cutoff) { return (k * k > cutoff) == (expansionLambda_ > cutoff); });
    if (!sameBranch) {
        throw std::runtime_error("WCAWE has no polynomial expansion on the requested port cutoff branch at "
                                 + std::to_string(frequencyHz) + " Hz");
    }

    const double coordinate = frequencyHz / expansionFrequencyHz_ - 1.0;
    auto matrix = reducedMatrices_.back();
    auto rhs = reducedRhs_.back();
    for (std::size_t r = reducedMatrices_.size() - 1; r > 0; --r) {
        for (std::size_t entry = 0; entry < matrix.size(); ++entry) {
            matrix[entry] = matrix[entry] * coordinate + reducedMatrices_[r - 1][entry];
        }
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            rhs[i] = rhs[i] * coordinate + reducedRhs_[r - 1][i];
        }
    }
    if (!fastsweep::denseSolve(matrix, rhs, romDim_)) {
        throw std::runtime_error("WcaweSweep::solveReduced: reduced system singular at "
                                 + std::to_string(frequencyHz) + " Hz");
    }
    return rhs;
}

// 由 ROM 解提取指定频率的二端口 S 参数。
SParameterPoint WcaweSweep::evaluate(double frequencyHz) const {
    if (!ready_) {
        return SParameterPoint{frequencyHz};
    }
    const auto coordinates = solveReduced(frequencyHz);
    return ResultExtractor(project_, portModeSolver_).extractFromPortProjections(
        frequencyHz, fastsweep::bilinear(observationVectors_[0], coordinates),
        fastsweep::bilinear(observationVectors_[1], coordinates));
}

// 将指定频率的 ROM 解提升回全阶边自由度，用于最后一个频点的场输出。
std::vector<WcaweSweep::Complex> WcaweSweep::reconstructField(
    double frequencyHz) const {
    const auto coordinates = solveReduced(frequencyHz);
    const auto& basis = wcaweBuild_.basis;
    std::vector<Complex> field(basis.front().size());
    for (std::size_t j = 0; j < basis.size(); ++j) {
        for (std::size_t i = 0; i < field.size(); ++i) {
            field[i] += basis[j][i] * coordinates[j];
        }
    }
    return field;
}

// 写出各阶 U 对角元、正交性误差及论文式 (7)、式 (8) 的残差。
bool WcaweSweep::writeBasisConditionCsv(
    const std::filesystem::path& path) const {
    std::ofstream output(path, std::ios::out | std::ios::trunc);
    if (!output.good()) {
        return false;
    }
    output
        << "order,u_diagonal_abs,u_diagonal_ratio,recurrence_residual,"
           "basis_relation_residual,orthogonality_error,"
           "triangular_solve_count\n";
    output << std::setprecision(17);
    for (const auto& step : wcaweBuild_.steps) {
        output << step.order << ','
               << step.diagonalAbs << ','
               << step.diagonalRatio << ','
               << step.recurrenceResidual << ','
               << step.basisRelationResidual << ','
               << step.orthogonalityError << ','
               << step.triangularSolveCount << '\n';
    }
    output.close();
    return output.good();
}

}  // namespace fem::sweep
