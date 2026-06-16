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

    static std::vector<std::vector<Complex>> buildPortVectors(
        const FEMAssembler::AffineSystem& affine,
        std::size_t fullDimension);

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
