#include "bpfem/linalg/PardisoBackend.hpp"

#include <stdexcept>

namespace fem::linalg {

PardisoBackend::PardisoBackend() = default;

SolveResult PardisoBackend::solve(const SparseMatrix& A,
                                  const std::vector<std::complex<double>>& b,
                                  const SolverConfig& /*cfg*/) {
#ifdef BPFEM_USE_MKL
    // SolverConfig::maxIterations / tolerance are ignored: PARDISO is direct.
    return solver_.solve(A, b);
#else
    (void)A;
    (void)b;
    throw std::runtime_error(
        "PardisoBackend: this binary was built without BPFEM_USE_MKL, "
        "but the sparse solver factory selected PARDISO. Rebuild with -DBPFEM_USE_MKL=ON "
        "or pass --linear-solver bicgstab.");
#endif
}

void PardisoBackend::rememberPatternForReuse(bool /*enable*/) {
    // PARDISO already reuses the symbolic factorization automatically when the
    // pattern repeats (see MklPardisoSolver::hasSamePattern). Nothing to set.
}

}  // namespace fem::linalg

