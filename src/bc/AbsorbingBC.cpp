// AbsorbingBC: 1st-order Engquist-Majda / Sommerfeld absorbing BC, placeholder.
//
// Status: skeleton. The boundary acts on a single open-air face and adds the
// surface contribution
//     n x curl(E) + j k_0 (n x (n x E)) = 0
// which translates to a complex surface mass-like operator:
//     +j k_0 * integral_T (n x N_i) . (n x N_j) dS
// added to the LHS for every pair of edge-element DOFs on the absorbing
// face. Implementation pattern is identical to ImpedanceBC; only the scalar
// coefficient is different (j*k0 vs j*omega*mu_0 / Z_s).
//
// Finish-here checklist mirrors ImpedanceBC.cpp; see that file's leading
// comment for the step-by-step.

#include "bpfem/bc/AbsorbingBC.hpp"

#include <stdexcept>

namespace fem::bc {

void AbsorbingBC::apply(SparseMatrix& /*matrix*/,
                        std::vector<std::complex<double>>& /*rhs*/,
                        const AssemblyContext& /*ctx*/) const {
    throw std::runtime_error("AbsorbingBC is not implemented yet");
}

void AbsorbingBC::declareSparsity(SparsePatternBuilder& /*builder*/,
                                  const AssemblyContext& /*ctx*/) const {
    // No-op while unimplemented; safe because apply() throws before any
    // matrix entry is ever needed.
}

}  // namespace fem::bc
