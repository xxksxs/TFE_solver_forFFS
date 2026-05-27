#pragma once

#include "bpfem/core/Types.hpp"
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
    MklPardisoSolver();
    ~MklPardisoSolver();

    MklPardisoSolver(const MklPardisoSolver&) = delete;
    MklPardisoSolver& operator=(const MklPardisoSolver&) = delete;

    SolveResult solve(const SparseMatrix& matrix, const std::vector<std::complex<double>>& rhs);

private:
#ifdef BPFEM_USE_MKL
    void initializeParameters();
    bool hasSamePattern(MKL_INT n, const std::vector<int>& rowOffsets, const std::vector<int>& columnIndices) const;
    void analyzePattern(MKL_INT n, const std::vector<int>& rowOffsets, const std::vector<int>& columnIndices, std::vector<MKL_Complex16>& values);
    void release();

    std::array<void*, 64> pt_{};
    std::array<MKL_INT, 64> iparm_{};
    std::vector<MKL_INT> ia_;
    std::vector<MKL_INT> ja_;
    MKL_INT n_ = 0;
    bool analyzed_ = false;
#endif
};

}  // namespace fem
