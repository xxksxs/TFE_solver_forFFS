#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/linalg/IPreconditioner.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"

#include <complex>
#include <memory>
#include <vector>

namespace fem {

// BiCGSTAB iterative solver for complex non-Hermitian sparse linear systems
// A x = b. Pure C++17, no external dependencies.
//
// Convergence on the H(curl) frequency-domain Maxwell systems this code base
// solves is poor: the operator A = K - k0^2 M + j beta m m^T is indefinite
// (eigenvalues straddle the origin), and BiCGSTAB on indefinite systems is
// numerically unstable. With the available preconditioners (Jacobi, ILU(0))
// the iteration does not converge on the BP filter benchmark in 4000+ iters
// (residual stalls near 1.0). This is *not* a bug in this implementation;
// it matches published behavior for indefinite Maxwell -- AMS / Schwarz /
// shifted Laplacian preconditioners are required for real iterative solves.
//
// Production runs should use PARDISO (--linear-solver direct or --linear-solver auto
// with BPFEM_USE_MKL=ON). BiCGStabSolver is retained as the no-MKL fallback
// and as a plug-in target for future H(curl)-aware preconditioners.
class BiCGStabSolver {
public:
    // Optional right-preconditioner. Default = nullptr (no preconditioning),
    // which preserves the exact iteration sequence the unpreconditioned
    // implementation produces (bit-for-bit equivalent to the pre-Phase 1
    // BiCGSTAB code path). When set, the solver runs right-preconditioned
    // BiCGSTAB: it conceptually solves A M^{-1} y = b for y and recovers
    // x = M^{-1} y, but in code this is implemented in place with two
    // preconditioner applies per iteration -- one on p, one on s.
    void setPreconditioner(std::shared_ptr<linalg::IPreconditioner> precon);

    // Solve A x = b. Returns SolveResult with the final iterate, the
    // iteration count actually performed, and the relative residual
    // ||A x - b|| / ||b||. The iteration stops on first of:
    //   - residual <= tolerance
    //   - iter == maxIterations  (residual is whatever it was at that step)
    //   - numerical breakdown (rho, denom, omega, or t.t hits ~ 1e-30)
    // The function never throws; the caller decides whether the returned
    // residual is acceptable.
    SolveResult solve(const SparseMatrix& a,
                      const std::vector<std::complex<double>>& b,
                      int maxIterations,
                      double tolerance) const;

private:
    static std::complex<double> dot(const std::vector<std::complex<double>>& a, const std::vector<std::complex<double>>& b);
    static double norm(const std::vector<std::complex<double>>& x);

    std::shared_ptr<linalg::IPreconditioner> precon_;
};

}  // namespace fem

