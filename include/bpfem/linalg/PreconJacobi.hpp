#pragma once

#include "bpfem/linalg/IPreconditioner.hpp"

#include <complex>
#include <vector>

namespace fem::linalg {

// Diagonal (Jacobi) preconditioner: M = diag(A), so apply() does
//   y_i = x_i / A_ii.
//
// Cheap and parameter-free. For the curl-curl H(curl) systems this code base
// solves it is *not* a strong preconditioner (the dominant block is
// elliptic-like with high anisotropy; AMS / SAI / Schwarz are needed for
// real performance). It is included as a sanity-check baseline and to
// exercise the IPreconditioner abstraction on a real problem.
class PreconJacobi : public IPreconditioner {
public:
    PreconJacobi() = default;

    void setup(const SparseMatrix& A) override;

    void apply(const std::vector<std::complex<double>>& x,
               std::vector<std::complex<double>>& y) const override;

    const char* name() const override { return "Jacobi"; }

private:
    // Pre-computed inverse of diag(A). Entries with |A_ii| < 1e-30 are stored
    // as 0 to avoid blowups (rare in healthy assemblies but can happen if a
    // boundary condition zeroed a row).
    std::vector<std::complex<double>> invDiag_;
};

}  // namespace fem::linalg

