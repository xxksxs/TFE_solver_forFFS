// ImpedanceBC: surface impedance Z_s boundary placeholder.
//
// Status: skeleton. apply() throws so any accidental registration is loud
// rather than producing silent-zero contributions. declareSparsity() is a
// no-op for the same reason: nothing will be written to the matrix, so no
// pattern entries need to be reserved.
//
// To finish this BC:
//   1. Have the constructor / a setter capture the surface face IDs that
//      should carry the Z_s = R + jX boundary, plus the impedance value.
//   2. In declareSparsity, walk all surface triangles whose facet maps to
//      one of those face IDs, and for each pair of edge-element DOFs on the
//      triangle add the (row, col) pair via builder.add(...).
//   3. In apply, assemble the surface mass-like contribution
//          (j omega mu_0 / Z_s) * integral_T (n x N_i) . (n x N_j) dS
//      and matrix.add it. This is purely an LHS contribution; rhs is left
//      untouched.
// See docs/optimization/strategy-interfaces/example.md for a full worked
// example using exactly this BC class.

#include "bpfem/bc/ImpedanceBC.hpp"

#include <stdexcept>

namespace fem::bc {

void ImpedanceBC::apply(SparseMatrix& /*matrix*/,
                        std::vector<std::complex<double>>& /*rhs*/,
                        const AssemblyContext& /*ctx*/) const {
    throw std::runtime_error("ImpedanceBC is not implemented yet");
}

void ImpedanceBC::declareSparsity(SparsePatternBuilder& /*builder*/,
                                  const AssemblyContext& /*ctx*/) const {
    // Intentionally empty: apply() throws, so no matrix entries will ever
    // be written by this BC. Nothing to reserve in the sparsity pattern.
}

}  // namespace fem::bc
