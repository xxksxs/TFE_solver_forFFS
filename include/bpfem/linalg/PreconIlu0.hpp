#pragma once

#include "bpfem/linalg/IPreconditioner.hpp"

#include <complex>
#include <cstddef>
#include <vector>

namespace fem::linalg {

// ILU(0) preconditioner: factorize A = (L + I)(D + U) - R where (L, D, U)
// share the sparsity of A and R contains the entries that would have filled
// in but were dropped. apply() does
//   y = M^{-1} x   with   M = (L + I) (D + U).
//
// Storage:
// We expand the upper-triangular CSR storage of SparseMatrix into a full CSR
// view (every off-diagonal stored both as A_ij and A_ji), then run the
// standard scalar ILU(0) algorithm in place on that copy. The result is
// stored as a single CSR array `lu_` with rowOffsets_ / columnIndices_ that
// mirrors the full matrix pattern; the diagonal entries hold U_ii, the
// strictly lower portion holds L_ij (L has implicit unit diagonal), and the
// strictly upper portion holds U_ij.
//
// Pivot guarding: when |U_ii| < pivotEpsilon, the entry is replaced by
// pivotEpsilon with its sign chosen along the original diagonal direction.
// This avoids divisions by zero when the operator is indefinite (typical
// for frequency-domain curl-curl systems near a resonance).
class PreconIlu0 : public IPreconditioner {
public:
    PreconIlu0() = default;

    // Configuration knobs. Defaults are chosen to match standard ILU(0)
    // textbook recipes; for stronger regularization on highly indefinite
    // frequency-domain systems set diagonalShift > 0 (typical 1e-4 .. 1e-2).
    void setDiagonalShift(double alpha) { diagonalShift_ = alpha; }
    void setPivotEpsilon(double eps) { pivotEpsilon_ = eps; }

    void setup(const SparseMatrix& A) override;

    void apply(const std::vector<std::complex<double>>& x,
               std::vector<std::complex<double>>& y) const override;

    const char* name() const override { return "ILU(0)"; }

private:
    std::size_t n_ = 0;
    // Full CSR storage of LU factor. Each row r has its entries sorted by
    // column index ascending; the entry with column == r holds U_ii.
    std::vector<int> rowOffsets_;
    std::vector<int> columnIndices_;
    std::vector<std::complex<double>> values_;
    // Index of the diagonal entry within each row in `values_` /
    // `columnIndices_`. -1 if this row has no diagonal entry (degenerate).
    std::vector<int> diagPos_;

    double diagonalShift_ = 0.0;
    double pivotEpsilon_ = 1.0e-30;
};

}  // namespace fem::linalg

