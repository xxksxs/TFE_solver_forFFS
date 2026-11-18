#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/linalg/IFactorizedSparseSolver.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"

#ifdef BPFEM_USE_MKL
#include <mkl.h>
#endif

#include <array>
#include <complex>
#include <vector>

namespace fem {

class MklPardisoSolver {
public:
    enum class Precision { Double, Single };
    explicit MklPardisoSolver(bool outOfCore = false, Precision precision = Precision::Double);
    ~MklPardisoSolver();

    MklPardisoSolver(const MklPardisoSolver&) = delete;
    MklPardisoSolver& operator=(const MklPardisoSolver&) = delete;

    SolveResult solve(const SparseMatrix& matrix, const std::vector<std::complex<double>>& rhs,
                      double tolerance = 1e-7);
    void factorize(const SparseMatrix& matrix);
    SolveResult solveFactorized(const std::vector<std::complex<double>>& rhs);
    std::vector<SolveResult> solveFactorizedBatch(
        const std::vector<std::vector<std::complex<double>>>& rhsList);
    void clearFactorization();
    linalg::FactorizationStatistics factorizationStatistics() const { return statistics_; }

private:
    bool outOfCore_;
    Precision precision_;
#ifdef BPFEM_USE_MKL
    void initializeParameters();
    bool hasSamePattern(MKL_INT n, const std::vector<int>& rowOffsets, const std::vector<int>& columnIndices) const;
    void analyzePattern(MKL_INT n, const std::vector<int>& rowOffsets,
                        const std::vector<int>& columnIndices, void* values);
    template<class Complex>
    void factorizeTyped(const SparseMatrix& matrix, std::vector<Complex>& values);
    template<class Complex>
    std::vector<SolveResult> solveTyped(
        const std::vector<std::vector<std::complex<double>>>& rhsList,
        std::vector<Complex>& matrix, std::vector<Complex>& rhs, std::vector<Complex>& solution);
    void release();

    std::array<void*, 64> pt_{};
    std::array<MKL_INT, 64> iparm_{};
    std::vector<MKL_INT> ia_;
    std::vector<MKL_INT> ja_;
    std::vector<MKL_Complex16> factorValues_;
    std::vector<MKL_Complex16> rhsValues_;
    std::vector<MKL_Complex16> solutionValues_;
    std::vector<MKL_Complex8> singleFactorValues_;
    std::vector<MKL_Complex8> singleRhsValues_;
    std::vector<MKL_Complex8> singleSolutionValues_;
    MKL_INT n_ = 0;
    bool analyzed_ = false;
    bool factorized_ = false;
#endif
    linalg::FactorizationStatistics statistics_;
};

}  // namespace fem
