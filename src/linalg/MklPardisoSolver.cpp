#include "bpfem/linalg/MklPardisoSolver.hpp"

#ifdef BPFEM_USE_MKL
#include <mkl.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem {
namespace {

double relativeResidual(const SparseMatrix& matrix, const std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& b) {
    const auto ax = matrix.multiply(x);
    double r2 = 0.0;
    double b2 = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) {
        r2 += std::norm(ax[i] - b[i]);
        b2 += std::norm(b[i]);
    }
    return std::sqrt(r2) / std::max(std::sqrt(b2), 1.0e-300);
}

}  // namespace

MklPardisoSolver::MklPardisoSolver() {
#ifdef BPFEM_USE_MKL
    initializeParameters();
#endif
}

MklPardisoSolver::~MklPardisoSolver() {
#ifdef BPFEM_USE_MKL
    release();
#endif
}

#ifdef BPFEM_USE_MKL
void MklPardisoSolver::initializeParameters() {
    pt_.fill(nullptr);
    iparm_.fill(0);
    iparm_[0] = 1;
    iparm_[1] = 2;
    iparm_[7] = 2;
    iparm_[9] = 13;
    iparm_[10] = 1;
    iparm_[12] = 1;
    iparm_[17] = -1;
    iparm_[18] = -1;
    iparm_[26] = 1;
    iparm_[34] = 1;
}

bool MklPardisoSolver::hasSamePattern(MKL_INT n, const std::vector<int>& rowOffsets, const std::vector<int>& columnIndices) const {
    if (!analyzed_ || n != n_ || rowOffsets.size() != ia_.size() || columnIndices.size() != ja_.size()) {
        return false;
    }
    return std::equal(rowOffsets.begin(), rowOffsets.end(), ia_.begin())
        && std::equal(columnIndices.begin(), columnIndices.end(), ja_.begin());
}

void MklPardisoSolver::analyzePattern(MKL_INT n, const std::vector<int>& rowOffsets, const std::vector<int>& columnIndices, std::vector<MKL_Complex16>& values) {
    if (analyzed_) {
        release();
    }
    n_ = n;
    ia_.assign(rowOffsets.begin(), rowOffsets.end());
    ja_.assign(columnIndices.begin(), columnIndices.end());

    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT mtype = 6;
    MKL_INT phase = 11;
    MKL_INT nrhs = 1;
    MKL_INT msglvl = 0;
    MKL_INT error = 0;
    MKL_INT idum = 0;
    MKL_Complex16 ddum{0.0, 0.0};

    pardiso(pt_.data(), &maxfct, &mnum, &mtype, &phase, &n_, values.data(), ia_.data(), ja_.data(), &idum, &nrhs, iparm_.data(), &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        release();
        throw std::runtime_error("MKL PARDISO analysis failed with error " + std::to_string(error));
    }
    analyzed_ = true;
}

void MklPardisoSolver::release() {
    if (!analyzed_) {
        return;
    }
    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT mtype = 6;
    MKL_INT phase = -1;
    MKL_INT nrhs = 1;
    MKL_INT msglvl = 0;
    MKL_INT error = 0;
    MKL_INT idum = 0;
    MKL_Complex16 ddum{0.0, 0.0};
    pardiso(pt_.data(), &maxfct, &mnum, &mtype, &phase, &n_, &ddum, ia_.data(), ja_.data(), &idum, &nrhs, iparm_.data(), &msglvl, &ddum, &ddum, &error);
    pt_.fill(nullptr);
    ia_.clear();
    ja_.clear();
    n_ = 0;
    analyzed_ = false;
}
#endif

SolveResult MklPardisoSolver::solve(const SparseMatrix& matrix, const std::vector<std::complex<double>>& rhs) {
#ifndef BPFEM_USE_MKL
    (void)matrix;
    (void)rhs;
    throw std::runtime_error("MKL PARDISO solver was requested, but this binary was built without BPFEM_USE_MKL.");
#else
    const MKL_INT n = static_cast<MKL_INT>(matrix.size());
    const auto& rowOffsets = matrix.rowOffsets();
    const auto& columnIndices = matrix.columnIndices();
    const auto& csrValues = matrix.values();
    std::vector<MKL_Complex16> a;
    a.reserve(csrValues.size());
    for (const auto& value : csrValues) {
        a.push_back(MKL_Complex16{value.real(), value.imag()});
    }

    std::vector<MKL_Complex16> b(rhs.size());
    std::vector<MKL_Complex16> x(rhs.size(), MKL_Complex16{0.0, 0.0});
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        b[i] = MKL_Complex16{rhs[i].real(), rhs[i].imag()};
    }

    if (!hasSamePattern(n, rowOffsets, columnIndices)) {
        analyzePattern(n, rowOffsets, columnIndices, a);
    }

    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT mtype = 6;
    MKL_INT phase = 23;
    MKL_INT nrhs = 1;
    MKL_INT msglvl = 0;
    MKL_INT error = 0;
    MKL_INT idum = 0;

    iparm_[17] = -1;
    iparm_[18] = -1;
    pardiso(pt_.data(), &maxfct, &mnum, &mtype, &phase, &n_, a.data(), ia_.data(), ja_.data(), &idum, &nrhs, iparm_.data(), &msglvl, b.data(), x.data(), &error);
    if (error != 0) {
        release();
        throw std::runtime_error("MKL PARDISO factor/solve failed with error " + std::to_string(error));
    }

    std::vector<std::complex<double>> field(rhs.size(), 0.0);
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        field[i] = std::complex<double>(x[i].real, x[i].imag);
    }
    return SolveResult{field, 1, relativeResidual(matrix, field, rhs)};
#endif
}

}  // namespace fem
