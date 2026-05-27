#include "bpfem/linalg/PreconJacobi.hpp"

#include <cmath>
#include <complex>

namespace fem::linalg {

namespace {

constexpr double kSmallDiagonalGuard = 1.0e-30;

}  // namespace

void PreconJacobi::setup(const SparseMatrix& A) {
    const std::size_t n = A.size();
    invDiag_.assign(n, std::complex<double>(0.0, 0.0));

    // SparseMatrix uses upper-triangular CSR: row r stores columns c with
    // c >= r in ascending order, so the diagonal entry of row r (if present)
    // is the first column index in that row's slice.
    const auto& rowOffsets = A.rowOffsets();
    const auto& columnIndices = A.columnIndices();
    const auto& values = A.values();
    if (rowOffsets.size() < n + 1 || values.empty()) {
        return;
    }
    for (std::size_t r = 0; r < n; ++r) {
        const int start = rowOffsets[r];
        const int end = rowOffsets[r + 1];
        if (start >= end) {
            continue;
        }
        if (columnIndices[static_cast<std::size_t>(start)] != static_cast<int>(r)) {
            continue;
        }
        const std::complex<double> diag = values[static_cast<std::size_t>(start)];
        if (std::abs(diag) > kSmallDiagonalGuard) {
            invDiag_[r] = std::complex<double>(1.0, 0.0) / diag;
        }
    }
}

void PreconJacobi::apply(const std::vector<std::complex<double>>& x,
                         std::vector<std::complex<double>>& y) const {
    if (y.size() != x.size()) {
        y.assign(x.size(), std::complex<double>(0.0, 0.0));
    }
    const std::size_t n = std::min(x.size(), invDiag_.size());
    for (std::size_t i = 0; i < n; ++i) {
        y[i] = invDiag_[i] * x[i];
    }
    for (std::size_t i = n; i < x.size(); ++i) {
        y[i] = x[i];  // no diagonal data => identity
    }
}

}  // namespace fem::linalg

