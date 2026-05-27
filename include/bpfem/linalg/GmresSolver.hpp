#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/linalg/IPreconditioner.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"

#include <complex>
#include <memory>
#include <vector>

namespace fem {

// Restarted GMRES(m) for complex non-Hermitian linear systems A x = b.
//
// Uses Arnoldi orthogonalization with modified Gram-Schmidt (one pass of
// re-orthogonalization for stability) and Givens rotations to maintain the
// upper Hessenberg least-squares subproblem in triangular form.
//
// When a preconditioner M is set via setPreconditioner(), runs *right*
// preconditioned GMRES: solves A M^{-1} y = b, then x = M^{-1} y. With
// M = I (no preconditioner), the iteration sequence reproduces unrestarted
// GMRES exactly (within the configured restart length).
class GmresSolver {
public:
    void setPreconditioner(std::shared_ptr<linalg::IPreconditioner> precon);

    // Restart length (Krylov subspace dimension before restart). Defaults to
    // 30, which matches typical PETSc / hypre / SLEPc presets and is a good
    // balance between memory (m * n complex doubles) and convergence.
    void setRestartLength(int m) { restartLength_ = (m > 0) ? m : 30; }

    SolveResult solve(const SparseMatrix& A,
                      const std::vector<std::complex<double>>& b,
                      int maxIterations,
                      double tolerance) const;

private:
    int restartLength_ = 30;
    std::shared_ptr<linalg::IPreconditioner> precon_;
};

}  // namespace fem

