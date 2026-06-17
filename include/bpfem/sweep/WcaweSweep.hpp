#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fastsweep/GalerkinReducedModel.hpp"
#include "bpfem/fastsweep/WellConditionedBasisBuilder.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"
#include "bpfem/sweep/ISweepStrategy.hpp"

#include <complex>
#include <filesystem>
#include <vector>

namespace fem::sweep {

struct WcaweOptions {
    // Target basis size. The builder may stop earlier if MGS detects a
    // near-dependent moment direction.
    int order = 12;

    // Expansion frequency in Hz. 0 means "use band center".
    double expansionFrequencyHz = 0.0;

    // Deflation threshold for the MGS residual / R diagonal.
    double dropTolerance = 1.0e-12;
};

class WcaweSweep : public ISweepStrategy {
public:
    using Complex = std::complex<double>;

    // 绑定工程、装配器、端口模式求解器和 WCAWE 配置。
    WcaweSweep(const ProjectDefinition& project,
               const FEMAssembler& assembler,
               const PortModeSolver& portModeSolver,
               WcaweOptions options = {});

    // 执行完整 WCAWE 扫频流程，并返回所有频点的 S 参数。
    SweepResult run(const std::vector<double>& frequencies,
                    const SweepContext& ctx) override;

    // 返回命令行和诊断文件使用的算法名。
    const char* name() const override { return "wcawe"; }

    // 离线阶段：生成良条件基并构建 Galerkin ROM。
    int buildOffline(double expansionFrequencyHz,
                     linalg::ISparseSolver& solver,
                     const linalg::SolverConfig& solverConfig);

    // 在线阶段：评估指定频率的 S 参数。
    SParameterPoint evaluate(double frequencyHz) const;

    // 在线阶段：重构指定频率的全阶场。
    std::vector<Complex> reconstructField(double frequencyHz) const;

    // 返回当前 ROM 维度。
    int dimension() const { return romDim_; }

    // 返回离线正交化中丢弃的候选列数。
    int deflatedColumns() const { return deflatedColumns_; }

    // 返回 WCAWE 离线阶段是否已经完成。
    bool ready() const { return ready_; }

    // 写出 WCAWE 条件数诊断 CSV。
    bool writeBasisConditionCsv(const std::filesystem::path& path) const;

private:
    // 生成传统 AWE 矩向量，随后由 WellConditionedBasisBuilder 变换为良条件基。
    std::vector<std::vector<Complex>> buildAweMoments(
        double expansionFrequencyHz,
        linalg::ISparseSolver& solver,
        const linalg::SolverConfig& solverConfig,
        const std::vector<std::vector<Complex>>& portVectors) const;

    // 求解指定频率的 reduced 坐标。
    std::vector<Complex> solveReduced(double frequencyHz) const;

    const ProjectDefinition& project_;
    const FEMAssembler& assembler_;
    const PortModeSolver& portModeSolver_;
    WcaweOptions options_;

    bool ready_ = false;
    int romDim_ = 0;
    int fullDim_ = 0;
    int numProjectPorts_ = 0;
    int numVirtualPorts_ = 0;
    int deflatedColumns_ = 0;
    double expansionFrequencyHz_ = 0.0;

    FEMAssembler::AffineSystem affine_;
    fastsweep::WellConditionedBasisBuilder basisBuilder_;
    fastsweep::GalerkinReducedModel model_;
};

}  // namespace fem::sweep
