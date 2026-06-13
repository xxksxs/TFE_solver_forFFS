#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fastsweep/PadeApproximant.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/sweep/ISweepStrategy.hpp"

#include <complex>
#include <vector>

namespace fem::sweep {

struct AweOptions {
    // Denominator order q for the [q-1/q] scalar Padé fit. The moment
    // recurrence therefore generates 2*q moments.
    int order = 8;

    // Expansion frequency in Hz. 0 means "use band center".
    double expansionFrequencyHz = 0.0;
};

class AweSweep : public ISweepStrategy {
public:
    AweSweep(const ProjectDefinition& project,
             const FEMAssembler& assembler,
             const PortModeSolver& portModeSolver,
             AweOptions options = {});

    SweepResult run(const std::vector<double>& frequencies,
                    const SweepContext& ctx) override;

    const char* name() const override { return "awe"; }

    int buildOffline(double expansionFrequencyHz,
                     linalg::ISparseSolver& solver,
                     const linalg::SolverConfig& solverConfig);
    SParameterPoint evaluate(double frequencyHz) const;

    bool ready() const { return ready_; }
    double expansionFrequency() const { return expansionFrequencyHz_; }
    int requestedOrder() const { return options_.order; }
    int inputDenominatorOrder() const { return inputProjectionPade_.denominatorOrder(); }
    int outputDenominatorOrder() const { return outputProjectionPade_.denominatorOrder(); }

private:
    using Complex = std::complex<double>;

    Complex portProjection(const std::vector<Complex>& edgeDofs, int faceId) const;
    const PortMode& virtualPortMode(int virtualPortIndex) const;

    const ProjectDefinition& project_;
    const FEMAssembler& assembler_;
    const PortModeSolver& portModeSolver_;
    AweOptions options_;

    bool ready_ = false;
    double expansionFrequencyHz_ = 0.0;
    double lambda0_ = 0.0;
    double lambdaScale_ = 1.0;

    FEMAssembler::AffineSystem affine_;
    fastsweep::PadeApproximant inputProjectionPade_;
    fastsweep::PadeApproximant outputProjectionPade_;
};

}  // namespace fem::sweep
