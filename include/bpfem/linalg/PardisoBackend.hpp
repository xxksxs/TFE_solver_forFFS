#pragma once

#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/linalg/MklPardisoSolver.hpp"

namespace fem::linalg {

// ISparseSolver wrapper around fem::MklPardisoSolver.
//
// Composition (not inheritance) so the existing MklPardisoSolver class stays
// untouched. The wrapper holds a single MklPardisoSolver instance for its
// whole lifetime, which is what enables PARDISO to skip reorder / analysis
// phases when the sparsity pattern is constant across sweep frequency points.
//
// Available only when the binary is built with BPFEM_USE_MKL. Outside that
// build flavor, factory code must not produce a PardisoBackend; we still
// declare the class so headers compile, but the .cpp guards against use.
class PardisoBackend : public ISparseSolver {
public:
    PardisoBackend();

    SolveResult solve(const SparseMatrix& A,
                      const std::vector<std::complex<double>>& b,
                      const SolverConfig& cfg = {}) override;

    void rememberPatternForReuse(bool enable) override;

    const char* name() const override { return "PARDISO"; }

private:
#ifdef BPFEM_USE_MKL
    MklPardisoSolver solver_;
#endif
};

}  // namespace fem::linalg

