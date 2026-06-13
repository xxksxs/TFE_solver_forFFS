#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/sweep/ISweepStrategy.hpp"
#include "bpfem/sweep/MgaweSweep.hpp"

#include <complex>
#include <vector>

namespace fem::sweep {

struct GaweOptions {
    // Number of local AWE moment directions used to build the Galerkin basis.
    int order = 12;

    // Expansion frequency in Hz. 0 means "use band center".
    double expansionFrequencyHz = 0.0;

    // Modified Gram-Schmidt deflation threshold.
    double dropTolerance = 1.0e-10;
};

class GaweSweep : public ISweepStrategy {
public:
    using Complex = std::complex<double>;

    GaweSweep(const ProjectDefinition& project,
              const FEMAssembler& assembler,
              const PortModeSolver& portModeSolver,
              GaweOptions options = {});

    SweepResult run(const std::vector<double>& frequencies,
                    const SweepContext& ctx) override;

    const char* name() const override { return "gawe"; }

    int buildOffline(double expansionFrequencyHz,
                     linalg::ISparseSolver& solver,
                     const linalg::SolverConfig& solverConfig);
    SParameterPoint evaluate(double frequencyHz) const;
    std::vector<Complex> reconstructField(double frequencyHz) const;

    int dimension() const { return delegate_.dimension(); }
    int retainedColumns() const { return delegate_.retainedColumns(); }
    int deflatedColumns() const { return delegate_.deflatedColumns(); }
    bool ready() const { return delegate_.ready(); }

private:
    GaweOptions options_;
    MgaweSweep delegate_;
};

}  // namespace fem::sweep
