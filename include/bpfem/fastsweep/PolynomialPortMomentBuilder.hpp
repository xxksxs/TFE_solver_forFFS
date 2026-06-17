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

    // 将 AffineSystem 中的端口耦合向量提升为全阶复向量列表。
    static std::vector<std::vector<Complex>> buildPortVectors(
        const FEMAssembler::AffineSystem& affine,
        std::size_t fullDimension);

    // 生成无损端口边界问题的 AWE 矩向量，供 AWE/GAWE/MGAWE/WCAWE 共用。
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
