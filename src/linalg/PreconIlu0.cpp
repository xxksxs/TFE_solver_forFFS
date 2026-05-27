#include "bpfem/linalg/PreconIlu0.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

namespace fem::linalg {

namespace {

// Expand SparseMatrix's upper-triangular CSR (only entries with col >= row,
// representing a complex-symmetric operator) into the full CSR pattern that
// stores both A_ij and A_ji. The resulting layout is row-major with each
// row sorted by column ascending; this is what scalar ILU(0) needs to work
// in place.
struct FullCsr {
    std::size_t n = 0;
    std::vector<int> rowOffsets;
    std::vector<int> columnIndices;
    std::vector<std::complex<double>> values;
    std::vector<int> diagPos;  // index in values/columnIndices of the diagonal entry of each row
};

FullCsr expandSymmetricUpperToFull(const SparseMatrix& A) {
    const std::size_t n = A.size();
    const auto& upperRowOffsets = A.rowOffsets();
    const auto& upperCols = A.columnIndices();
    const auto& upperVals = A.values();

    // First pass: count entries per row (full matrix).
    std::vector<int> rowCount(n, 0);
    for (std::size_t r = 0; r < n; ++r) {
        const int begin = upperRowOffsets[r];
        const int end = upperRowOffsets[r + 1];
        for (int k = begin; k < end; ++k) {
            const int c = upperCols[static_cast<std::size_t>(k)];
            ++rowCount[r];
            if (c != static_cast<int>(r)) {
                ++rowCount[static_cast<std::size_t>(c)];
            }
        }
    }

    FullCsr out;
    out.n = n;
    out.rowOffsets.assign(n + 1, 0);
    for (std::size_t r = 0; r < n; ++r) {
        out.rowOffsets[r + 1] = out.rowOffsets[r] + rowCount[r];
    }
    const std::size_t totalNnz = static_cast<std::size_t>(out.rowOffsets[n]);
    out.columnIndices.assign(totalNnz, 0);
    out.values.assign(totalNnz, std::complex<double>(0.0, 0.0));
    out.diagPos.assign(n, -1);

    // Second pass: emit unsorted entries into per-row scratch, then sort.
    std::vector<int> writeCursor(n, 0);
    for (std::size_t r = 0; r < n; ++r) {
        const int begin = upperRowOffsets[r];
        const int end = upperRowOffsets[r + 1];
        for (int k = begin; k < end; ++k) {
            const int c = upperCols[static_cast<std::size_t>(k)];
            const auto v = upperVals[static_cast<std::size_t>(k)];
            const int slotR = out.rowOffsets[r] + writeCursor[r]++;
            out.columnIndices[static_cast<std::size_t>(slotR)] = c;
            out.values[static_cast<std::size_t>(slotR)] = v;
            if (c != static_cast<int>(r)) {
                const std::size_t cc = static_cast<std::size_t>(c);
                const int slotC = out.rowOffsets[cc] + writeCursor[cc]++;
                out.columnIndices[static_cast<std::size_t>(slotC)] = static_cast<int>(r);
                out.values[static_cast<std::size_t>(slotC)] = v;  // symmetric
            }
        }
    }

    // Sort each row by column ascending and find the diagonal slot.
    for (std::size_t r = 0; r < n; ++r) {
        const int begin = out.rowOffsets[r];
        const int end = out.rowOffsets[r + 1];
        const std::size_t len = static_cast<std::size_t>(end - begin);
        // Build (col, value) pairs, sort, write back.
        std::vector<std::pair<int, std::complex<double>>> entries(len);
        for (int k = 0; k < end - begin; ++k) {
            entries[static_cast<std::size_t>(k)] = {
                out.columnIndices[static_cast<std::size_t>(begin + k)],
                out.values[static_cast<std::size_t>(begin + k)]
            };
        }
        std::sort(entries.begin(), entries.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        for (std::size_t k = 0; k < len; ++k) {
            out.columnIndices[static_cast<std::size_t>(begin + static_cast<int>(k))] = entries[k].first;
            out.values[static_cast<std::size_t>(begin + static_cast<int>(k))] = entries[k].second;
            if (entries[k].first == static_cast<int>(r)) {
                out.diagPos[r] = begin + static_cast<int>(k);
            }
        }
    }

    return out;
}

}  // namespace

void PreconIlu0::setup(const SparseMatrix& A) {
    FullCsr full = expandSymmetricUpperToFull(A);
    n_ = full.n;
    rowOffsets_ = std::move(full.rowOffsets);
    columnIndices_ = std::move(full.columnIndices);
    values_ = std::move(full.values);
    diagPos_ = std::move(full.diagPos);

    // Optional diagonal shift A_ii <- A_ii + alpha * |A_ii| (sign of alpha
    // chosen along the original A_ii direction; keeps the operator on the
    // same complex-plane half).
    if (diagonalShift_ > 0.0) {
        for (std::size_t r = 0; r < n_; ++r) {
            const int dp = diagPos_[r];
            if (dp < 0) continue;
            const std::complex<double> d = values_[static_cast<std::size_t>(dp)];
            const double mag = std::abs(d);
            if (mag > 0.0) {
                values_[static_cast<std::size_t>(dp)] = d * (1.0 + diagonalShift_);
            } else {
                values_[static_cast<std::size_t>(dp)] = std::complex<double>(diagonalShift_, 0.0);
            }
        }
    }

    // Standard scalar ILU(0) (Saad 1996, Alg. 10.3) operating on the full CSR.
    // For each row i (in order):
    //   for each k in row i with k < i:
    //     a_{ik} <- a_{ik} / a_{kk}
    //     for each j in row i with j > k:
    //       if (k, j) is a non-zero in the original pattern (i.e. exists in
    //       row k at column j), then a_{ij} <- a_{ij} - a_{ik} * a_{kj}
    //   pivot guard on a_{ii}.
    //
    // We use a column-position map per inner step to make the "exists in row k"
    // test O(1).
    std::vector<int> rowKColPos(n_, -1);  // col -> position in row k, scratch
    for (std::size_t i = 0; i < n_; ++i) {
        const int rowStart = rowOffsets_[i];
        const int rowEnd = rowOffsets_[i + 1];

        for (int idxIK = rowStart; idxIK < rowEnd; ++idxIK) {
            const int k = columnIndices_[static_cast<std::size_t>(idxIK)];
            if (k >= static_cast<int>(i)) {
                continue;  // we only iterate the strictly-lower (L) entries
            }
            // a_{ik} <- a_{ik} / a_{kk}
            const int dpK = diagPos_[static_cast<std::size_t>(k)];
            if (dpK < 0) continue;
            std::complex<double> diagK = values_[static_cast<std::size_t>(dpK)];
            if (std::abs(diagK) < pivotEpsilon_) {
                // Replace with epsilon, sign-preserving when possible.
                if (std::abs(diagK) > 0.0) {
                    diagK = diagK / std::abs(diagK) * pivotEpsilon_;
                } else {
                    diagK = std::complex<double>(pivotEpsilon_, 0.0);
                }
                values_[static_cast<std::size_t>(dpK)] = diagK;
            }
            const std::complex<double> aIK = values_[static_cast<std::size_t>(idxIK)] / diagK;
            values_[static_cast<std::size_t>(idxIK)] = aIK;

            // Build column-position map of row k (only on first inner k of
            // this i, but we rebuild per k since the cost is just O(|row k|)).
            const int kRowStart = rowOffsets_[static_cast<std::size_t>(k)];
            const int kRowEnd = rowOffsets_[static_cast<std::size_t>(k) + 1];
            for (int p = kRowStart; p < kRowEnd; ++p) {
                rowKColPos[static_cast<std::size_t>(columnIndices_[static_cast<std::size_t>(p)])] = p;
            }

            for (int idxIJ = idxIK + 1; idxIJ < rowEnd; ++idxIJ) {
                const int j = columnIndices_[static_cast<std::size_t>(idxIJ)];
                const int posKJ = rowKColPos[static_cast<std::size_t>(j)];
                if (posKJ < 0) continue;  // (k, j) not in pattern, drop
                values_[static_cast<std::size_t>(idxIJ)] -=
                    aIK * values_[static_cast<std::size_t>(posKJ)];
            }

            // Reset scratch.
            for (int p = kRowStart; p < kRowEnd; ++p) {
                rowKColPos[static_cast<std::size_t>(columnIndices_[static_cast<std::size_t>(p)])] = -1;
            }
        }

        // Pivot guard on the new U_ii.
        const int dpI = diagPos_[i];
        if (dpI >= 0) {
            std::complex<double> d = values_[static_cast<std::size_t>(dpI)];
            if (std::abs(d) < pivotEpsilon_) {
                if (std::abs(d) > 0.0) {
                    d = d / std::abs(d) * pivotEpsilon_;
                } else {
                    d = std::complex<double>(pivotEpsilon_, 0.0);
                }
                values_[static_cast<std::size_t>(dpI)] = d;
            }
        }
    }
}

void PreconIlu0::apply(const std::vector<std::complex<double>>& x,
                       std::vector<std::complex<double>>& y) const {
    if (n_ == 0) {
        y = x;
        return;
    }
    if (y.size() != x.size()) {
        y.assign(x.size(), std::complex<double>(0.0, 0.0));
    }
    // Forward substitution: solve (L + I) z = x, with L stored in the
    // strictly-lower entries of `values_` (implicit unit diagonal of L).
    // Reusing `y` as scratch for z, then overwriting with the upper solve.
    std::vector<std::complex<double>> z(n_, std::complex<double>(0.0, 0.0));
    for (std::size_t i = 0; i < n_; ++i) {
        const int rowStart = rowOffsets_[i];
        const int rowEnd = rowOffsets_[i + 1];
        std::complex<double> sum = x[i];
        for (int idx = rowStart; idx < rowEnd; ++idx) {
            const int j = columnIndices_[static_cast<std::size_t>(idx)];
            if (j >= static_cast<int>(i)) break;  // sorted ascending; remaining are upper
            sum -= values_[static_cast<std::size_t>(idx)] * z[static_cast<std::size_t>(j)];
        }
        z[i] = sum;
    }
    // Backward substitution: solve (D + U) y = z, with U_ii on diagPos and
    // strictly-upper entries holding U_ij for j > i.
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n_) - 1; i >= 0; --i) {
        const int rowStart = rowOffsets_[static_cast<std::size_t>(i)];
        const int rowEnd = rowOffsets_[static_cast<std::size_t>(i) + 1];
        std::complex<double> sum = z[static_cast<std::size_t>(i)];
        for (int idx = rowStart; idx < rowEnd; ++idx) {
            const int j = columnIndices_[static_cast<std::size_t>(idx)];
            if (j <= i) continue;
            sum -= values_[static_cast<std::size_t>(idx)] * y[static_cast<std::size_t>(j)];
        }
        const int dp = diagPos_[static_cast<std::size_t>(i)];
        if (dp < 0) {
            y[static_cast<std::size_t>(i)] = sum;  // degenerate row, leave unscaled
        } else {
            y[static_cast<std::size_t>(i)] = sum / values_[static_cast<std::size_t>(dp)];
        }
    }
}

}  // namespace fem::linalg

