#pragma once

#include "bpfem/linalg/IFactorizedSparseSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/linalg/MklPardisoSolver.hpp"
#include "bpfem/linalg/PreconPardisoSingle.hpp"  // Compatibility for earlier public includes.

namespace fem::linalg {

// ISparseSolver wrapper around fem::MklPardisoSolver.
//
// Composition (not inheritance) so the existing MklPardisoSolver class stays
// untouched. The wrapper holds a single MklPardisoSolver instance for its
// whole lifetime, which is what enables PARDISO to skip reorder / analysis
// phases when the sparsity pattern is constant across sweep frequency points.
//
// Explicit PARDISO selection is accepted in all builds. Factorization and
// solving report an unavailable backend when BPFEM_USE_MKL is disabled.
class PardisoBackend : public ISparseSolver, public IFactorizedSparseSolver {
public:
    explicit PardisoBackend(bool outOfCore = false);

    SolveResult solve(const SparseMatrix& A,
                      const std::vector<std::complex<double>>& b,
                      const SolverConfig& cfg = {}) override;

    void rememberPatternForReuse(bool enable) override;

    void factorize(const SparseMatrix& matrix,
                   const SolverConfig& config = {}) override;

    SolveResult solveFactorized(
        const std::vector<std::complex<double>>& rhs,
        const SolverConfig& config = {}) override;

    std::vector<SolveResult> solveFactorizedBatch(
        const std::vector<std::vector<std::complex<double>>>& rhsList,
        const SolverConfig& config = {}) override;

    void clearFactorization() override;

    FactorizationStatistics factorizationStatistics() const override;

    const char* name() const override { return outOfCore_ ? "PARDISO OOC" : "PARDISO"; }

private:
    bool outOfCore_;
#ifdef BPFEM_USE_MKL
    MklPardisoSolver solver_;
#endif
};

}  // namespace fem::linalg

