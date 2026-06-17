#include "bpfem/sweep/WcaweSweep.hpp"

#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"
#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem::sweep {

// 构造 WCAWE 扫频器，并把阶数和丢弃阈值规整到可用范围。
WcaweSweep::WcaweSweep(const ProjectDefinition& project,
                       const FEMAssembler& assembler,
                       const PortModeSolver& portModeSolver,
                       WcaweOptions options)
    : project_(project),
      assembler_(assembler),
      portModeSolver_(portModeSolver),
      options_(options),
      basisBuilder_(options.dropTolerance) {
    options_.order = std::max(1, options_.order);
    if (options_.dropTolerance <= 0.0) {
        options_.dropTolerance = 1.0e-12;
      }
}

// 执行完整 WCAWE 扫频：离线生成良条件基，在线逐频点评估 S 参数，并写出诊断文件。
SweepResult WcaweSweep::run(const std::vector<double>& frequencies, const SweepContext& ctx) {
    const double expansionHz = options_.expansionFrequencyHz > 0.0
        ? options_.expansionFrequencyHz
        : (frequencies.empty()
            ? 0.5 * (project_.sweep.startHz + project_.sweep.endHz)
            : 0.5 * (frequencies.front() + frequencies.back()));

    ctx.log.info("WCAWE expansion frequency: " + std::to_string(expansionHz / 1.0e9) + " GHz");
    ctx.log.info("WCAWE target basis order: " + std::to_string(options_.order));
    std::ostringstream dropTol;
    dropTol << std::scientific << std::setprecision(3) << options_.dropTolerance;
    ctx.log.info("WCAWE MGS/R diagonal drop tolerance: " + dropTol.str());

    linalg::SolverConfig solverCfg;
    solverCfg.maxIterations = ctx.linearMaxIterations;
    solverCfg.tolerance = ctx.linearTolerance;
    const int dim = buildOffline(expansionHz, ctx.solver, solverCfg);
    ctx.log.info("WCAWE ROM dimension: " + std::to_string(dim));
    ctx.log.info("WCAWE deflated candidate columns: " + std::to_string(deflatedColumns_));

    if (!ctx.outputDirectory.empty()) {
        const auto conditionPath = ctx.outputDirectory / "basis_condition.csv";
        if (writeBasisConditionCsv(conditionPath)) {
            ctx.log.info("Wrote " + conditionPath.string());
        } else {
            ctx.log.warn("Failed to write " + conditionPath.string());
        }
    }

    SweepResult out;
    out.points.reserve(frequencies.size());
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        const double f = frequencies[i];
        if (i == 0 || (i + 1) % 50 == 0 || i + 1 == frequencies.size()) {
            ctx.log.info("WCAWE evaluating point " + std::to_string(i + 1) + "/"
                         + std::to_string(frequencies.size())
                         + " @ " + std::to_string(f / 1.0e9) + " GHz");
        }
        out.points.push_back(evaluate(f));
        out.lastFrequencyHz = f;
    }
    if (!frequencies.empty()) {
        out.lastEdgeDofs = reconstructField(frequencies.back());
    }

    if (!ctx.outputDirectory.empty()) {
        fastsweep::FastSweepDiagnostics diag;
        diag.algorithm = name();
        diag.expansionFrequenciesHz = {expansionHz};
        diag.requestedOrder = options_.order;
        diag.romDimension = romDim_;
        diag.retainedColumns = romDim_;
        diag.deflatedColumns = deflatedColumns_;
        diag.basisOrthogonalityError = model_.basisOrthogonalityError();
        diag.wcaweMomentReconstructionError = basisBuilder_.maxMomentReconstructionError();
        diag.reducedSolveSucceeded = model_.lastSolveSucceeded();
        diag.maxPassivityError = fastsweep::maxPassivityError(out.points);
        const auto diagnosticsPath = ctx.outputDirectory / "diagnostics.json";
        if (fastsweep::writeDiagnosticsJson(diagnosticsPath, diag)) {
            ctx.log.info("Wrote " + diagnosticsPath.string());
        } else {
            ctx.log.warn("Failed to write " + diagnosticsPath.string());
        }
    }
    (void)ctx.onFieldSolved;
    return out;
}

// 在展开频点处生成 AWE 矩向量，经过 WCAWE 正交化后构造统一 Galerkin ROM。
int WcaweSweep::buildOffline(double expansionFrequencyHz,
                             linalg::ISparseSolver& solver,
                             const linalg::SolverConfig& solverConfig) {
    if (project_.ports.size() < 2) {
        throw std::runtime_error("WcaweSweep: at least two project ports are required");
    }

    expansionFrequencyHz_ = expansionFrequencyHz;
    numProjectPorts_ = static_cast<int>(project_.ports.size());
    affine_ = assembler_.buildAffineSystem();
    fullDim_ = static_cast<int>(affine_.K.size());
    numVirtualPorts_ = static_cast<int>(affine_.portCoupling.size());
    if (numVirtualPorts_ == 0) {
        throw std::runtime_error("WcaweSweep: AffineSystem produced zero port entries");
    }
    if (!affine_.lossless) {
        throw std::runtime_error("WcaweSweep: lossy materials (sigma != 0) are not supported. Use --sweep direct.");
    }

    const auto portVectors =
        fastsweep::PolynomialPortMomentBuilder::buildPortVectors(affine_, affine_.K.size());

    basisBuilder_.clear();
    const auto aweMoments = buildAweMoments(expansionFrequencyHz, solver, solverConfig, portVectors);
    for (const auto& moment : aweMoments) {
        basisBuilder_.append(moment);
    }

    romDim_ = basisBuilder_.dimension();
    deflatedColumns_ = basisBuilder_.deflatedColumns();
    if (romDim_ == 0) {
        throw std::runtime_error("WcaweSweep: well-conditioned basis collapsed to zero dimension");
    }

    model_.build(project_, affine_, portModeSolver_, basisBuilder_.basis(), portVectors);
    ready_ = true;
    return romDim_;
}

// 调用公共矩递推模块，生成传统 AWE 矩向量，供 WCAWE 基构造器重新组合。
std::vector<std::vector<WcaweSweep::Complex>> WcaweSweep::buildAweMoments(
    double expansionFrequencyHz,
    linalg::ISparseSolver& solver,
    const linalg::SolverConfig& solverConfig,
    const std::vector<std::vector<Complex>>& portVectors) const {
    return fastsweep::PolynomialPortMomentBuilder::generateLosslessMoments(
        project_, assembler_, portModeSolver_, affine_, expansionFrequencyHz,
        std::max(1, options_.order), solver, solverConfig, portVectors);
}

// 在降阶空间内求解指定频率的小型 ROM 线性系统。
std::vector<WcaweSweep::Complex> WcaweSweep::solveReduced(double frequencyHz) const {
    return model_.solveReduced(frequencyHz);
}

// 由 ROM 解提取指定频率的二端口 S 参数。
SParameterPoint WcaweSweep::evaluate(double frequencyHz) const {
    return model_.evaluate(frequencyHz);
}

// 将指定频率的 ROM 解提升回全阶边自由度，用于最后一个频点的场输出。
std::vector<WcaweSweep::Complex> WcaweSweep::reconstructField(double frequencyHz) const {
    return model_.reconstructField(frequencyHz);
}

// 写出 WCAWE 条件数诊断表，记录传统矩基与正交基的稳定性差异。
bool WcaweSweep::writeBasisConditionCsv(const std::filesystem::path& path) const {
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out.good()) {
        return false;
    }
    out << "order,awe_condition_proxy,wcawe_condition_proxy,r_diagonal_abs,"
           "orthogonality_error,moment_reconstruction_error\n";
    out << std::setprecision(17);
    for (const auto& r : basisBuilder_.conditionRecords()) {
        out << r.order << ','
            << r.aweConditionProxy << ','
            << r.wcaweConditionProxy << ','
            << r.rDiagonalAbs << ','
            << r.orthogonalityError << ','
            << r.momentReconstructionError << '\n';
    }
    out.flush();
    return out.good();
}

}  // namespace fem::sweep
