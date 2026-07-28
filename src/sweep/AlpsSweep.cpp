#include "bpfem/sweep/AlpsSweep.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"
#include "bpfem/fastsweep/LinearAlgebra.hpp"
#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"
#include "bpfem/fastsweep/PortModeUtilities.hpp"
#include "bpfem/linalg/FactorizedSolveSession.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fem::sweep {
namespace {

using Complex = std::complex<double>;
using Vector = std::vector<Complex>;

// 从增广状态中取出完整场分量。
Vector fieldPart(const Vector& state, std::size_t fullDimension) {
    if (state.size() != fullDimension + 1) {
        throw std::runtime_error("AlpsSweep: invalid augmented state dimension");
    }
    return Vector(state.begin(), state.begin() + static_cast<std::ptrdiff_t>(fullDimension));
}

// 把完整场向量扩充一个标量状态，用于表示 b(t)=b0+t*b1。
Vector augmentedField(const Vector& field, Complex scalar) {
    Vector state = field;
    state.push_back(scalar);
    return state;
}

}  // namespace

// 保存依赖引用；大规模矩阵直到 run/buildOffline 时才装配。
AlpsSweep::AlpsSweep(const ProjectDefinition& project,
                     const FEMAssembler& assembler,
                     const PortModeSolver& portModeSolver,
                     AlpsOptions options)
    : project_(project), assembler_(assembler), portModeSolver_(portModeSolver), options_(options) {}

// 构造文档约定的初始展开点集合，并用最近局部 Padé 模型完成扫频。
SweepResult AlpsSweep::run(const std::vector<double>& frequencies, const SweepContext& ctx) {
    const auto offlineStarted = std::chrono::steady_clock::now();
    prepareAffineModel();
    std::vector<double> expansionFrequencies = {
        options_.expansionFrequencyHz > 0.0
            ? options_.expansionFrequencyHz
            : (frequencies.empty()
                ? 0.5 * (project_.sweep.startHz + project_.sweep.endHz)
                : 0.5 * (frequencies.front() + frequencies.back()))
    };
    expansionFrequencyHz_ = expansionFrequencies.front();
    const int localOrder = std::max(1, options_.order);
    ctx.log.info("ALPS model: standard single-point two-sided Lanczos-Pade with P1 port linearization");
    ctx.log.info("ALPS scalar Padé order: " + std::to_string(localOrder));
    ctx.log.info("ALPS expansion frequency: "
                 + std::to_string(expansionFrequencyHz_ / 1.0e9) + " GHz");

    linalg::SolverConfig solverCfg;
    solverCfg.maxIterations = ctx.linearMaxIterations;
    solverCfg.tolerance = ctx.linearTolerance;
    localModels_.clear();
    localModels_.reserve(expansionFrequencies.size());
    retainedColumns_ = 0;
    deflatedColumns_ = 0;
    romDim_ = 0;
    portLinearizationSec_ = 0.0;
    lanczosOperatorSec_ = 0.0;
    orthogonalizationSec_ = 0.0;
    poleDecompositionSec_ = 0.0;
    romProjectionSec_ = 0.0;
    for (double expansionHz : expansionFrequencies) {
        LocalPadeModel local = buildLocalModel(
            expansionHz, localOrder, ctx.solver, solverCfg);
        retainedColumns_ += local.input.dimension() + local.output.dimension();
        deflatedColumns_ += 2 * localOrder
            - local.input.dimension() - local.output.dimension();
        romDim_ += local.input.dimension();
        localModels_.push_back(std::move(local));
    }
    ready_ = !localModels_.empty();
    const double offlineBuildSec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - offlineStarted).count();
    ctx.log.info("ALPS total input Padé dimension: " + std::to_string(romDim_));
    ctx.log.info("ALPS retained Lanczos columns: " + std::to_string(retainedColumns_));
    ctx.log.info("ALPS deflated Lanczos columns: " + std::to_string(deflatedColumns_));

    SweepResult out;
    out.offlineBuildSec = offlineBuildSec;
    out.portLinearizationSec = portLinearizationSec_;
    out.lanczosOperatorSec = lanczosOperatorSec_;
    out.poleDecompositionSec = poleDecompositionSec_;
    out.orthogonalizationSec = orthogonalizationSec_;
    out.romProjectionSec = romProjectionSec_;
    out.points.reserve(frequencies.size());
    const auto onlineStarted = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        const double frequencyHz = frequencies[i];
        if (i == 0 || (i + 1) % 50 == 0 || i + 1 == frequencies.size()) {
            ctx.log.info("ALPS evaluating point " + std::to_string(i + 1) + "/"
                         + std::to_string(frequencies.size()) + " @ "
                         + std::to_string(frequencyHz / 1.0e9) + " GHz");
        }
        out.points.push_back(evaluate(frequencyHz));
        out.lastFrequencyHz = frequencyHz;
    }
    out.onlineSweepSec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - onlineStarted).count();
    if (!frequencies.empty()) {
        out.lastEdgeDofs = reconstructField(frequencies.back());
    }

    if (!ctx.outputDirectory.empty()) {
        fastsweep::FastSweepDiagnostics diag;
        diag.algorithm = name();
        diag.expansionFrequenciesHz = expansionFrequencies;
        diag.requestedOrder = options_.order;
        diag.romDimension = romDim_;
        diag.retainedColumns = retainedColumns_;
        diag.deflatedColumns = deflatedColumns_;
        diag.maxPassivityError = numProjectPorts_ == 2
            ? fastsweep::maxPassivityError(out.points) : 0.0;
        diag.portLinearizationSec = portLinearizationSec_;
        diag.lanczosOperatorSec = lanczosOperatorSec_;
        diag.orthogonalizationSec = orthogonalizationSec_;
        diag.poleDecompositionSec = poleDecompositionSec_;
        if (auto* factorized = dynamic_cast<linalg::IFactorizedSparseSolver*>(&ctx.solver)) {
            const auto stats = factorized->factorizationStatistics();
            diag.factorizedRhsSolveCount = stats.factorizedRhsSolveCount;
            diag.factorizedSolveCallCount = stats.factorizedSolveCallCount;
            diag.batchRhsMax = stats.batchRhsMax;
        }
        for (const auto& local : localModels_) {
            diag.reducedSolveSucceeded = diag.reducedSolveSucceeded
                && local.input.lastSolveSucceeded() && local.output.lastSolveSucceeded();
            diag.lanczosBiorthogonalityError = std::max(
                diag.lanczosBiorthogonalityError,
                std::max(local.input.biorthogonalityError(),
                         local.output.biorthogonalityError()));
            diag.lanczosTridiagonalLeakage = std::max(
                diag.lanczosTridiagonalLeakage,
                std::max(local.input.tridiagonalLeakage(),
                         local.output.tridiagonalLeakage()));
            const double coupling = std::min(
                local.input.finalCouplingMagnitude(), local.output.finalCouplingMagnitude());
            if (diag.lanczosFinalCoupling == 0.0 || coupling < diag.lanczosFinalCoupling) {
                diag.lanczosFinalCoupling = coupling;
            }
            diag.lanczosBreakdownDetected = diag.lanczosBreakdownDetected
                || local.input.breakdownDetected() || local.output.breakdownDetected();
            diag.momentMatchingError = std::max(
                diag.momentMatchingError,
                std::max(local.input.momentMatchingError(),
                         local.output.momentMatchingError()));
            diag.lanczosRecurrenceResidual = std::max(
                diag.lanczosRecurrenceResidual,
                std::max(local.input.recurrenceResidual(),
                         local.output.recurrenceResidual()));
            diag.lookAheadCount += local.input.lookAheadCount()
                + local.output.lookAheadCount();
            diag.selectiveReorthogonalizationCount +=
                local.input.selectiveReorthogonalizationCount()
                + local.output.selectiveReorthogonalizationCount();
            diag.poleResidueReconstructionError = std::max(
                diag.poleResidueReconstructionError,
                std::max(local.input.poleResidueReconstructionError(),
                         local.output.poleResidueReconstructionError()));
            if (!frequencies.empty()) {
                const auto parameterAt = [&](double frequencyHz) {
                    const double k0 = 2.0 * pi * frequencyHz / c0;
                    return (k0 * k0 - local.lambda0) / local.lambdaScale;
                };
                const double firstParameter = parameterAt(frequencies.front());
                const double lastParameter = parameterAt(frequencies.back());
                const double parameterMin = std::min(firstParameter, lastParameter);
                const double parameterMax = std::max(firstParameter, lastParameter);
                diag.spuriousPoleCount += local.input.countSuspectedSpuriousPoles(
                    parameterMin, parameterMax);
                diag.spuriousPoleCount += local.output.countSuspectedSpuriousPoles(
                    parameterMin, parameterMax);
            }
        }
        const auto path = ctx.outputDirectory / "diagnostics.json";
        if (fastsweep::writeDiagnosticsJson(path, diag)) {
            ctx.log.info("Wrote " + path.string());
        } else {
            ctx.log.warn("Failed to write " + path.string());
        }
    }
    (void)ctx.onFieldSolved;
    return out;
}

// 装配公共仿射数据，并找出用于 S11/S21 的主导虚拟端口。
void AlpsSweep::prepareAffineModel() {
    ready_ = false;
    numProjectPorts_ = static_cast<int>(project_.ports.size());
    if (numProjectPorts_ < 2) {
        throw std::runtime_error("AlpsSweep: at least two ports are required for S11/S21");
    }
    affine_ = assembler_.buildAffineSystem();
    fullDim_ = static_cast<int>(affine_.K.size());
    numVirtualPorts_ = static_cast<int>(affine_.portCoupling.size());
    if (numVirtualPorts_ == 0) {
        throw std::runtime_error("AlpsSweep: AffineSystem produced zero port entries");
    }
    if (!affine_.lossless) {
        throw std::runtime_error(
            "AlpsSweep: lossy materials are not supported; use --sweep direct");
    }
    portVectors_ = fastsweep::PolynomialPortMomentBuilder::buildPortVectors(
        affine_, affine_.K.size());
    const auto dominant = fastsweep::dominantVirtualPortsByProject(
        portModeSolver_, affine_, numProjectPorts_);
    inputVirtualPort_ = dominant.empty() ? -1 : dominant[0];
    outputVirtualPort_ = dominant.size() < 2 ? -1 : dominant[1];
    if (inputVirtualPort_ < 0 || outputVirtualPort_ < 0) {
        throw std::runtime_error("AlpsSweep: failed to resolve dominant input/output modes");
    }
}

// 在一个展开点构造增广 P1 算子和 S11/S21 两个双边 Lanczos 模型。
AlpsSweep::LocalPadeModel AlpsSweep::buildLocalModel(
    double expansionFrequencyHz,
    int localOrder,
    linalg::ISparseSolver& solver,
    const linalg::SolverConfig& solverConfig) {
    const auto linearizationStarted = std::chrono::steady_clock::now();
    auto linearization = fastsweep::PolynomialPortMomentBuilder::buildLosslessLinearization(
        project_, assembler_, portModeSolver_, affine_, expansionFrequencyHz, portVectors_);
    portLinearizationSec_ += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - linearizationStarted).count();

    linalg::FactorizedSolveSession solveSession(
        solver, linearization.matrixAtExpansion, solverConfig, false);
    const auto expansionSolve = solveSession.solve(linearization.rhs0);
    if (expansionSolve.field.size() != static_cast<std::size_t>(fullDim_)) {
        throw std::runtime_error("AlpsSweep: expansion solve returned the wrong field size");
    }
    const Vector rightStart = augmentedField(expansionSolve.field, Complex(1.0, 0.0));

    // 右移动 G=-L0^-1*L1；最后一个 Lanczos 步只需要该算子。
    auto applyRight = [&](const Vector& state) {
        const auto started = std::chrono::steady_clock::now();
        const Vector field = fieldPart(state, static_cast<std::size_t>(fullDim_));
        Vector rhs = fastsweep::PolynomialPortMomentBuilder::applyLosslessFirstOrderMatrix(
            affine_, linearization, portVectors_, field);
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            rhs[i] = -rhs[i] + linearization.rhs1[i] * state.back();
        }
        auto result = augmentedField(solveSession.solve(rhs).field, Complex(0.0, 0.0));
        lanczosOperatorSec_ += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        return result;
    };

    // 同时计算 G*v 与 G^T*w；PARDISO 后端把两个 RHS 合并为一次 phase=33。
    auto applyPair = [&](const Vector& rightState, const Vector& leftState) {
        const auto started = std::chrono::steady_clock::now();
        const Vector rightField = fieldPart(
            rightState, static_cast<std::size_t>(fullDim_));
        Vector rightRhs = fastsweep::PolynomialPortMomentBuilder::applyLosslessFirstOrderMatrix(
            affine_, linearization, portVectors_, rightField);
        for (std::size_t i = 0; i < rightRhs.size(); ++i) {
            rightRhs[i] = -rightRhs[i] + linearization.rhs1[i] * rightState.back();
        }
        Vector leftRhs = fieldPart(leftState, static_cast<std::size_t>(fullDim_));
        std::vector<Vector> rhsBatch;
        rhsBatch.reserve(2);
        rhsBatch.push_back(std::move(rightRhs));
        rhsBatch.push_back(std::move(leftRhs));
        auto solved = solveSession.solveBatch(rhsBatch);
        Vector leftImage = fastsweep::PolynomialPortMomentBuilder::applyLosslessFirstOrderMatrix(
            affine_, linearization, portVectors_, solved[1].field);
        for (auto& value : leftImage) {
            value = -value;
        }
        auto result = std::make_pair(
            augmentedField(solved[0].field, Complex(0.0, 0.0)),
            augmentedField(leftImage,
                           fastsweep::bilinear(linearization.rhs1, solved[1].field)));
        lanczosOperatorSec_ += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        return result;
    };

    const Vector inputOutput = augmentedField(
        portVectors_[static_cast<std::size_t>(inputVirtualPort_)], Complex(0.0, 0.0));
    const Vector transmissionOutput = augmentedField(
        portVectors_[static_cast<std::size_t>(outputVirtualPort_)], Complex(0.0, 0.0));
    LocalPadeModel local;
    local.expansionFrequencyHz = expansionFrequencyHz;
    local.lambda0 = linearization.lambda0;
    local.lambdaScale = linearization.lambdaScale;

    const double operatorBefore = lanczosOperatorSec_;
    const auto lanczosStarted = std::chrono::steady_clock::now();
    local.input.buildPaired(rightStart, inputOutput, applyRight, applyPair,
                            localOrder, options_.dropTolerance, true);
    local.output.buildPaired(rightStart, transmissionOutput, applyRight, applyPair,
                             localOrder, options_.dropTolerance, false);
    const double totalLanczosSec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - lanczosStarted).count();
    const double poleSec = local.input.poleDecompositionSec()
        + local.output.poleDecompositionSec();
    poleDecompositionSec_ += poleSec;
    romProjectionSec_ += poleSec;
    orthogonalizationSec_ += std::max(
        0.0, totalLanczosSec - (lanczosOperatorSec_ - operatorBefore) - poleSec);
    return local;
}

// 构造一个显式单展开点模型，保持测试和离线调用接口稳定。
int AlpsSweep::buildOffline(double expansionFrequencyHz,
                            linalg::ISparseSolver& solver,
                            const linalg::SolverConfig& solverConfig) {
    prepareAffineModel();
    expansionFrequencyHz_ = expansionFrequencyHz;
    const int localOrder = std::max(1, options_.order);
    retainedColumns_ = 0;
    deflatedColumns_ = 0;
    romDim_ = 0;
    portLinearizationSec_ = 0.0;
    lanczosOperatorSec_ = 0.0;
    orthogonalizationSec_ = 0.0;
    poleDecompositionSec_ = 0.0;
    romProjectionSec_ = 0.0;
    localModels_.clear();
    LocalPadeModel local = buildLocalModel(
        expansionFrequencyHz, localOrder, solver, solverConfig);
    retainedColumns_ = local.input.dimension() + local.output.dimension();
    deflatedColumns_ = 2 * localOrder - retainedColumns_;
    romDim_ = local.input.dimension();
    localModels_.push_back(std::move(local));
    ready_ = true;
    return romDim_;
}

// 用 q 与 q-2 阶 Padé 序列差异选择误差最小的局部模型。
const AlpsSweep::LocalPadeModel& AlpsSweep::nearestModel(double frequencyHz) const {
    if (!ready_ || localModels_.empty()) {
        throw std::runtime_error("AlpsSweep: no local Padé model is ready");
    }
    if (localModels_.size() == 1) {
        return localModels_.front();
    }
    const double k0 = 2.0 * pi * frequencyHz / c0;
    const double lambda = k0 * k0;
    const LocalPadeModel* best = &localModels_.front();
    double bestScore = std::numeric_limits<double>::infinity();
    double bestDistance = std::numeric_limits<double>::infinity();
    for (const auto& local : localModels_) {
        const Complex parameter((lambda - local.lambda0) / local.lambdaScale, 0.0);
        const double score = std::max(
            local.input.relativeTruncationError(parameter),
            local.output.relativeTruncationError(parameter));
        const double distance = std::abs(local.expansionFrequencyHz - frequencyHz);
        if (score < bestScore || (score == bestScore && distance < bestDistance)) {
            best = &local;
            bestScore = score;
            bestDistance = distance;
        }
    }
    return *best;
}

// 评估最近局部 Padé 场投影，并按端口功率归一化定义提取 S11/S21。
SParameterPoint AlpsSweep::evaluate(double frequencyHz) const {
    SParameterPoint point;
    point.frequencyHz = frequencyHz;
    if (!ready_) {
        return point;
    }
    const auto& local = nearestModel(frequencyHz);
    const double k0 = 2.0 * pi * frequencyHz / c0;
    const double lambda = k0 * k0;
    const Complex parameter((lambda - local.lambda0) / local.lambdaScale, 0.0);
    const Complex rawInput = local.input.evaluate(parameter);
    const Complex rawOutput = local.output.evaluate(parameter);

    const PortMode& inputMode = fastsweep::virtualPortMode(
        portModeSolver_, affine_, inputVirtualPort_);
    const PortMode& outputMode = fastsweep::virtualPortMode(
        portModeSolver_, affine_, outputVirtualPort_);
    const double inputNormalization = powerNormalizationFactor(inputMode, frequencyHz);
    const double outputNormalization = powerNormalizationFactor(outputMode, frequencyHz);
    if (inputNormalization <= 0.0 || outputNormalization <= 0.0) {
        return point;
    }
    const Complex incoming = rawInput / inputNormalization;
    const Complex outgoing = rawOutput / outputNormalization;
    const auto& excitation = project_.ports[0];
    const Complex incident = excitation.excited
        ? std::polar(std::sqrt(std::max(excitation.magnitudeW, 0.0)),
                     excitation.phaseDeg * pi / 180.0)
        : Complex(1.0, 0.0);
    if (incident == Complex(0.0, 0.0)) {
        return point;
    }
    point.s11 = (incoming - incident) / incident;
    point.s21 = outgoing / incident;
    return point;
}

// 重构最近输入模型的完整场；增广状态最后一个标量不写入 VTU。
std::vector<AlpsSweep::Complex> AlpsSweep::reconstructField(double frequencyHz) const {
    const auto& local = nearestModel(frequencyHz);
    const double k0 = 2.0 * pi * frequencyHz / c0;
    const double lambda = k0 * k0;
    const Complex parameter((lambda - local.lambda0) / local.lambdaScale, 0.0);
    return fieldPart(local.input.reconstruct(parameter), static_cast<std::size_t>(fullDim_));
}

}  // namespace fem::sweep
