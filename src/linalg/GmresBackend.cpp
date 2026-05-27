#include "bpfem/linalg/GmresBackend.hpp"

#include <utility>

namespace fem::linalg {

void GmresBackend::setPreconditioner(std::shared_ptr<IPreconditioner> precon) {
    solver_.setPreconditioner(std::move(precon));
}

SolveResult GmresBackend::solve(const SparseMatrix& A,
                                const std::vector<std::complex<double>>& b,
                                const SolverConfig& cfg) {
    return solver_.solve(A, b, cfg.maxIterations, cfg.tolerance);
}

}  // namespace fem::linalg

