#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fastsweep/GalerkinReducedModel.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"
#include "bpfem/sweep/ISweepStrategy.hpp"

#include <complex>
#include <utility>
#include <vector>

namespace fem::sweep {

struct MgaweOptions {
    // Number of expansion points distributed over the requested sweep band.
    int expansionPointCount = 3;
    bool adaptive = false;
    bool useLastPolynomial = false;
    std::vector<std::pair<double, int>> expansions;
    bool polynomialProjection = false;

    // Local AWE moment vectors generated per expansion point before global
    // orthogonalization/deflation.
    int localOrder = 8;

    // Modified Gram-Schmidt deflation threshold for the unified basis.
    double dropTolerance = 1.0e-10;
};

class MgaweSweep : public ISweepStrategy {
public:
    using Complex = std::complex<double>;

    MgaweSweep(const ProjectDefinition& project,
               const FEMAssembler& assembler,
               const PortModeSolver& portModeSolver,
               MgaweOptions options = {});

    SweepResult run(const std::vector<double>& frequencies,
                    const SweepContext& ctx) override;

    const char* name() const override { return "mgawe"; }

    int buildOffline(const std::vector<double>& expansionFrequencies,
                     linalg::ISparseSolver& solver,
                     const linalg::SolverConfig& solverConfig);
    SParameterPoint evaluate(double frequencyHz) const;
    std::vector<Complex> reconstructField(double frequencyHz) const;

    int dimension() const { return romDim_; }
    int retainedColumns() const { return retainedColumns_; }
    int deflatedColumns() const { return deflatedColumns_; }
    bool ready() const { return ready_; }
    double basisOrthogonalityError() const { return model_.basisOrthogonalityError(); }
    bool reducedSolveSucceeded() const { return model_.lastSolveSucceeded(); }
    double lastOrthogonalizationSec() const { return orthogonalizationSec_; }
    double lastRomProjectionSec() const { return romProjectionSec_; }

private:
    bool appendIfIndependent(std::vector<Complex>&& w);
    void projectReducedModel(const std::vector<std::vector<Complex>>& portVectors);
    int build(std::vector<double>& expansionFrequencies,
              const std::vector<double>& sweepFrequencies,
              linalg::ISparseSolver& solver, const linalg::SolverConfig& solverConfig);
    std::vector<Complex> solveReduced(double frequencyHz) const;
    const PortMode& virtualPortMode(int virtualPortIndex) const;

    const ProjectDefinition& project_;
    const FEMAssembler& assembler_;
    const PortModeSolver& portModeSolver_;
    MgaweOptions options_;

    bool ready_ = false;
    int romDim_ = 0;
    int fullDim_ = 0;
    int numProjectPorts_ = 0;
    int numVirtualPorts_ = 0;
    int retainedColumns_ = 0;
    int deflatedColumns_ = 0;
    double orthogonalizationSec_ = 0.0;
    double romProjectionSec_ = 0.0;

    FEMAssembler::AffineSystem affine_;
    std::vector<std::vector<Complex>> basis_;
    fastsweep::GalerkinReducedModel model_;
    std::vector<fastsweep::PolynomialPortMomentBuilder::PortTaylorSeries> polynomials_;
    std::vector<int> localOrders_;
    double maxRelativeResidual_ = 0.0;
    bool adaptiveConverged_ = false;
};

}  // namespace fem::sweep
