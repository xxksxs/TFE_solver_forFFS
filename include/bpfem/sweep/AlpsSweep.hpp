#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fastsweep/LanczosPadeModel.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/sweep/ISweepStrategy.hpp"

#include <complex>
#include <vector>

namespace fem::sweep {

// ALPS 使用端口一阶线性化和双边 Lanczos 隐式构造 Padé 传递函数。
struct AlpsOptions {
    int order = 12;                   // 每个标量 [q-1/q] Padé 模型的直接阶数
    int maxRestarts = 1;             // 1 为稳定单点模式；3 预留给残差驱动多点模式
    double dropTolerance = 1.0e-12;  // Lanczos 击穿判据
    double expansionFrequencyHz = 0.0;  // 0 表示使用频带中心
};

class AlpsSweep : public ISweepStrategy {
public:
    using Complex = std::complex<double>;

    // 保存装配器和端口求解器引用，离线阶段才执行大规模稀疏求解。
    AlpsSweep(const ProjectDefinition& project,
              const FEMAssembler& assembler,
              const PortModeSolver& portModeSolver,
              AlpsOptions options = {});

    // 在指定点或频带中心构造单点 Lanczos-Padé 模型，再评估全部频点。
    SweepResult run(const std::vector<double>& frequencies,
                    const SweepContext& ctx) override;

    // 返回命令行和诊断文件使用的算法标识。
    const char* name() const override { return "alps"; }

    // 构造一个指定展开点的模型，供单元测试和显式单点模式使用。
    int buildOffline(double expansionFrequencyHz,
                     linalg::ISparseSolver& solver,
                     const linalg::SolverConfig& solverConfig);

    // 由最近的局部 Padé 模型计算一个频点的 S11/S21。
    SParameterPoint evaluate(double frequencyHz) const;

    // 由最近的局部输入模型恢复完整棱自由度场。
    std::vector<Complex> reconstructField(double frequencyHz) const;

    // 返回所有局部输入 Padé 模型的维度总和。
    int dimension() const { return romDim_; }

    // 返回单点模式的展开频率；自动模式返回频带中心。
    double expansionFrequency() const { return expansionFrequencyHz_; }

    // 指示离线模型是否已构造完成。
    bool ready() const { return ready_; }

    // 返回成功保留的左右 Lanczos 列数总和。
    int retainedColumns() const { return retainedColumns_; }

    // 返回因击穿提前停止而未生成的列数。
    int deflatedColumns() const { return deflatedColumns_; }

private:
    struct LocalPadeModel {
        double expansionFrequencyHz = 0.0;
        double lambda0 = 0.0;
        double lambdaScale = 1.0;
        fastsweep::LanczosPadeModel input;
        fastsweep::LanczosPadeModel output;
    };

    // 装配频率无关系统并解析主导输入、输出虚拟端口。
    void prepareAffineModel();

    // 在一个展开点构造 P1 增广系统和两个标量 Padé 模型。
    LocalPadeModel buildLocalModel(double expansionFrequencyHz,
                                   int localOrder,
                                   linalg::ISparseSolver& solver,
                                   const linalg::SolverConfig& solverConfig);

    // 返回唯一局部模型；保留多模型选择逻辑供后续残差驱动扩展使用。
    const LocalPadeModel& nearestModel(double frequencyHz) const;

    const ProjectDefinition& project_;
    const FEMAssembler& assembler_;
    const PortModeSolver& portModeSolver_;
    AlpsOptions options_;

    bool ready_ = false;
    int romDim_ = 0;
    double expansionFrequencyHz_ = 0.0;
    int fullDim_ = 0;
    int numProjectPorts_ = 0;
    int numVirtualPorts_ = 0;
    int inputVirtualPort_ = -1;
    int outputVirtualPort_ = -1;
    int retainedColumns_ = 0;
    int deflatedColumns_ = 0;
    double portLinearizationSec_ = 0.0;
    double lanczosOperatorSec_ = 0.0;
    double orthogonalizationSec_ = 0.0;
    double poleDecompositionSec_ = 0.0;
    double romProjectionSec_ = 0.0;

    FEMAssembler::AffineSystem affine_;
    std::vector<std::vector<Complex>> portVectors_;
    std::vector<LocalPadeModel> localModels_;
};

}  // namespace fem::sweep
