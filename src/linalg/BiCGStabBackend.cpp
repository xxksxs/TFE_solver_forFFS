#include "bpfem/linalg/BiCGStabBackend.hpp"

#include <utility>

namespace fem::linalg {

void BiCGStabBackend::setPreconditioner(std::shared_ptr<IPreconditioner> precon) {
    solver_.setPreconditioner(std::move(precon));
}

SolveResult BiCGStabBackend::solve(const SparseMatrix& A,
                                   const std::vector<std::complex<double>>& b,
                                   const SolverConfig& cfg) {
    return solver_.solve(A, b, cfg.maxIterations, cfg.tolerance);
}

}  // namespace fem::linalg

