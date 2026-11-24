#pragma once

#include "bpfem/linalg/IFactorizedSparseSolver.hpp"

#include <filesystem>
#include <memory>

namespace fem::linalg {

enum class MumpsMode { InCore, OutOfCore, BlockLowRank };

class MumpsBackend final : public ISparseSolver, public IFactorizedSparseSolver {
public:
    explicit MumpsBackend(MumpsMode mode, double blrTolerance,
                          std::filesystem::path oocDirectory = {});
    ~MumpsBackend() override;

    SolveResult solve(const SparseMatrix& matrix,
                      const std::vector<std::complex<double>>& rhs,
                      const SolverConfig& config = {}) override;
    void factorize(const SparseMatrix& matrix, const SolverConfig& config = {}) override;
    SolveResult solveFactorized(const std::vector<std::complex<double>>& rhs,
                                const SolverConfig& config = {}) override;
    std::vector<SolveResult> solveFactorizedBatch(
        const std::vector<std::vector<std::complex<double>>>& rhsList,
        const SolverConfig& config = {}) override;
    void clearFactorization() override;
    FactorizationStatistics factorizationStatistics() const override { return statistics_; }
    const char* name() const override { return "MUMPS"; }

private:
    struct State;
    MumpsMode mode_;
    double blrTolerance_;
    std::filesystem::path oocDirectory_;
    std::unique_ptr<State> state_;
    FactorizationStatistics statistics_;
};

}  // namespace fem::linalg
