#pragma once

#include "bpfem/linalg/BiCGStabSolver.hpp"
#include "bpfem/linalg/IPreconditioner.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"

#include <memory>

namespace fem::linalg {

// ISparseSolver wrapper around fem::BiCGStabSolver. Composition keeps the
// existing iterative-solver implementation untouched.
//
// SolverConfig::maxIterations and SolverConfig::tolerance flow straight into
// BiCGStabSolver::solve.
//
// Preconditioner: if set via setPreconditioner() the iteration runs in
// right-preconditioned form; if unset, behavior is bit-for-bit identical to
// the original unpreconditioned BiCGSTAB.
class BiCGStabBackend : public ISparseSolver {
public:
    BiCGStabBackend() = default;

    void setPreconditioner(std::shared_ptr<IPreconditioner> precon);

    SolveResult solve(const SparseMatrix& A,
                      const std::vector<std::complex<double>>& b,
                      const SolverConfig& cfg = {}) override;

    const char* name() const override { return "BiCGSTAB"; }

private:
    BiCGStabSolver solver_;
};

}  // namespace fem::linalg

