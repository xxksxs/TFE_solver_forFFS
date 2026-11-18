#pragma once

#include "bpfem/linalg/IFactorizedSparseSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/linalg/MklPardisoSolver.hpp"
#include "bpfem/linalg/IPreconditioner.hpp"

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

class PreconPardisoSingle final : public IPreconditioner {
public:
    PreconPardisoSingle();
    void setup(const SparseMatrix& matrix) override;
    void clear() noexcept override;
    void apply(const std::vector<std::complex<double>>& input,
               std::vector<std::complex<double>>& output) const override;
    const char* name() const override { return "PARDISO single"; }
    bool requiresFlexibleGmres() const override { return true; }

private:
    mutable MklPardisoSolver solver_{false, MklPardisoSolver::Precision::Single};
};

}  // namespace fem::linalg

