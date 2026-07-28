#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"

#include <complex>
#include <vector>

namespace fem::fastsweep {

class PolynomialPortMomentBuilder {
public:
    using Complex = std::complex<double>;

    // ALPS/WCAWE 共用的无损 P1 系统：A(t)=A0+t*A1，b(t)=b0+t*b1。
    struct LosslessLinearization {
        SparseMatrix matrixAtExpansion;
        std::vector<Complex> rhs0;
        std::vector<Complex> rhs1;
        std::vector<double> portAdmittanceFirstCoefficients;
        double lambda0 = 0.0;
        double lambdaScale = 1.0;
    };

    // 将 AffineSystem 中的端口耦合向量提升为全阶复向量列表。
    static std::vector<std::vector<Complex>> buildPortVectors(
        const FEMAssembler::AffineSystem& affine,
        std::size_t fullDimension);

    // 生成无损端口边界问题的 AWE 矩向量，供 AWE/GAWE/MGAWE/WCAWE 共用。
    // 围绕展开点构造文档 P1 端口仿射近似，参数 t=(k0^2-lambda0)/lambdaScale。
    static LosslessLinearization buildLosslessLinearization(
        const ProjectDefinition& project,
        const FEMAssembler& assembler,
        const PortModeSolver& portModeSolver,
        const FEMAssembler::AffineSystem& affine,
        double expansionFrequencyHz,
        const std::vector<std::vector<Complex>>& portVectors);

    // 计算一阶矩阵系数 A1 乘向量，供 ALPS 左右移动算子共用。
    static std::vector<Complex> applyLosslessFirstOrderMatrix(
        const FEMAssembler::AffineSystem& affine,
        const LosslessLinearization& linearization,
        const std::vector<std::vector<Complex>>& portVectors,
        const std::vector<Complex>& x);

    static std::vector<std::vector<Complex>> generateLosslessMoments(
        const ProjectDefinition& project,
        const FEMAssembler& assembler,
        const PortModeSolver& portModeSolver,
        const FEMAssembler::AffineSystem& affine,
        double expansionFrequencyHz,
        int momentCount,
        linalg::ISparseSolver& solver,
        const linalg::SolverConfig& solverConfig,
        const std::vector<std::vector<Complex>>& portVectors);
};

}  // namespace fem::fastsweep
