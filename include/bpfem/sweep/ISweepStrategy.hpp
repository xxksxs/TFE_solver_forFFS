#pragma once

#include "bpfem/core/Logger.hpp"
#include "bpfem/core/Types.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/post/ResultExtractor.hpp"

#include <complex>
#include <functional>
#include <vector>

namespace fem::sweep {

// Per-call data passed to a sweep strategy. Holds references only.
//
// `solver` is the linear-system backend used by direct strategies. Adaptive
// strategies (ALPS) build their own internal direct solver because their
// expansion-point factorization needs to outlive a single solve and lives
// on a different sparsity pattern than `assembler.assemble(f)`.
struct SweepContext {
    const ProjectDefinition& project;
    const FEMAssembler& assembler;
    const PortModeSolver& portModeSolver;
    const ResultExtractor& extractor;
    linalg::ISparseSolver& solver;
    Logger& log;

    // Per-frequency callback invoked by direct-style strategies after a
    // successful solve. ALPS-style strategies may ignore it (no per-frequency
    // edge-DOF vector to publish). Default empty closure = no field output.
    using FieldSink = std::function<void(double frequencyHz,
                                         const std::vector<std::complex<double>>& edgeDofs)>;
    FieldSink onFieldSolved;

    // Per-frequency progress callback. DirectSweep invokes this after each
    // assemble + solve with timing breakdown so Application can publish
    // [stat] log rows and append to the JSON sweep array. AlpsSweep does
    // not call it (its online evaluate is too cheap to merit per-point
    // breakdown). Default empty closure = no progress instrumentation.
    using ProgressSink = std::function<void(double frequencyHz,
                                            double assembleSec,
                                            double solveSec,
                                            int iterations,
                                            double residual)>;
    ProgressSink onFrequencyStat;

    // Linear-solver tolerance / iteration budget forwarded by Application.
    // Direct backends ignore both.
    int linearMaxIterations = 400;
    double linearTolerance = 1.0e-7;
};

// Output of a sweep run.
//
// `lastEdgeDofs` is populated by direct-style strategies so Application can
// write `field_last.vtu`. ALPS leaves it empty because its reduced model does
// not retain a full-space basis; this is documented and unchanged from prior
// behavior.
struct SweepResult {
    std::vector<SParameterPoint> points;
    std::vector<std::complex<double>> lastEdgeDofs;
    double lastFrequencyHz = 0.0;
};

// Frequency-sweep strategy interface.
//
// Concrete kinds:
//   - DirectSweep   per-frequency assemble + ISparseSolver::solve, no MOR
//   - AlpsSweep     single-point Krylov MOR, complex-symmetric Galerkin
class ISweepStrategy {
public:
    virtual ~ISweepStrategy() = default;

    // `frequencies` is the user-supplied list. Implementations may evaluate
    // those points exactly (DirectSweep) or use them as a basis for adaptive
    // refinement and return a different set (future AdaptiveSweep). The
    // returned `SweepResult::points` must be sorted by frequency.
    virtual SweepResult run(const std::vector<double>& frequencies,
                            const SweepContext& ctx) = 0;

    virtual const char* name() const = 0;
};

}  // namespace fem::sweep

