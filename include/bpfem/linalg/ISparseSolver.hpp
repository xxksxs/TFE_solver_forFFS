#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"

#include <complex>
#include <vector>

namespace fem::linalg {

// Common configuration object for all sparse solvers. Direct solvers ignore
// max-iteration / tolerance fields; iterative solvers respect them.
struct SolverConfig {
    int maxIterations = 400;
    double tolerance = 1.0e-7;
};

// Abstract sparse complex linear solver: A * x = b.
//
// Concrete backends:
//   - PardisoBackend  (MKL PARDISO direct solve, available iff BPFEM_USE_MKL)
//   - BiCGStabBackend (BiCGSTAB iterative fallback, always available)
//
// Construct a backend through fem::factory::makeSparseSolver(opts) so that the
// caller does not need to know which backend the build supports.
class ISparseSolver {
public:
    virtual ~ISparseSolver() = default;

    // Solve A x = b. Implementations may keep internal symbolic factorization
    // state across calls when patternReuseEnabled() is true.
    virtual SolveResult solve(const SparseMatrix& A,
                              const std::vector<std::complex<double>>& b,
                              const SolverConfig& cfg = {}) = 0;

    // When the caller knows the sparsity pattern stays constant across solves
    // (typical for frequency sweeps with a fixed mesh + DOF map), enabling
    // this lets direct backends reuse the symbolic factorization. Iterative
    // backends ignore the flag. Default is true: matches the existing
    // sweep-loop behavior.
    virtual void rememberPatternForReuse(bool /*enable*/) {}

    // Human-readable backend name for logging.
    virtual const char* name() const = 0;
};

}  // namespace fem::linalg

