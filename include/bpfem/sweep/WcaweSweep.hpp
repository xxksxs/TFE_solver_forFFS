#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fastsweep/PolynomialPortMomentBuilder.hpp"
#include "bpfem/fastsweep/FastSweepDiagnostics.hpp"
#include "bpfem/fastsweep/PolynomialWcaweRecurrence.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"
#include "bpfem/sweep/ISweepStrategy.hpp"

#include <array>
#include <complex>
#include <filesystem>
#include <vector>

namespace fem::sweep {

struct WcaweOptions {
    // 论文 WCAWE 单展开点递推的目标基阶数。
    bool checkMoments = false;
    bool residualDiagnostics = false;
    int order = 12;

    // 展开频率，单位 Hz；取 0 时使用扫频带中心。
    double expansionFrequencyHz = 0.0;

    // U 对角元和 MGS 后残差的相对 breakdown 阈值。
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

    // 执行论文一致的 WCAWE 离线构建、在线扫频和诊断输出。
    SweepResult run(const std::vector<double>& frequencies,
                    const SweepContext& ctx) override;

    // 返回命令行和诊断文件使用的算法名。
    const char* name() const override { return "wcawe"; }

    // 使用单展开点多项式逐阶生成 WCAWE 基并构建约化方程。
    int buildOffline(double expansionFrequencyHz,
                     linalg::ISparseSolver& solver,
                     const linalg::SolverConfig& solverConfig);

    // 在线阶段评估指定频率的 S 参数。
    SParameterPoint evaluate(double frequencyHz) const;

    // 在线阶段重构指定频率的全阶场。
    std::vector<Complex> reconstructField(double frequencyHz) const;

    // 返回当前 ROM 维度。
    int dimension() const { return romDim_; }

    // 返回因 breakdown 未生成的目标列数量。
    int deflatedColumns() const { return deflatedColumns_; }

    // 返回 WCAWE 离线阶段是否已经完成。
    bool ready() const { return ready_; }

    // 写出每阶 U 对角、递推残差和正交性诊断 CSV。
    bool writeBasisConditionCsv(const std::filesystem::path& path) const;

private:
    void projectPolynomial(
        const FEMAssembler::AffineSystem& affine,
        const fastsweep::PolynomialPortMomentBuilder::PortTaylorSeries& series,
        const SparseMatrix& matrixAtExpansion,
        const std::vector<Complex>& rhsAtExpansion);

    // 求解指定频率的 reduced 坐标。
    std::vector<Complex> solveReduced(double frequencyHz) const;

    const ProjectDefinition& project_;
    const FEMAssembler& assembler_;
    const PortModeSolver& portModeSolver_;
    WcaweOptions options_;

    double fullMomentProjectionError_ = -1.0;
    bool ready_ = false;
    int romDim_ = 0;
    int deflatedColumns_ = 0;
    double expansionFrequencyHz_ = 0.0;
    double expansionLambda_ = 0.0;
    double orthogonalizationSec_ = 0.0;
    double romProjectionSec_ = 0.0;

    std::vector<double> portCutoffsSquared_;
    std::optional<FEMAssembler::AffineSystem> diagnosticAffine_;
    std::optional<fastsweep::PolynomialPortMomentBuilder::PortTaylorSeries> diagnosticTaylor_;
    fastsweep::WcaweBuildResult wcaweBuild_;
    std::vector<std::vector<Complex>> reducedMatrices_;
    std::vector<std::vector<Complex>> reducedRhs_;
    std::array<std::vector<Complex>, 2> observationVectors_;
};

}  // namespace fem::sweep
