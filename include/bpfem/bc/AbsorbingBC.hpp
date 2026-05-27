#pragma once

#include "bpfem/bc/IBoundaryCondition.hpp"

namespace fem::bc {

// First-order absorbing boundary condition (1st-order Engquist-Majda /
// Sommerfeld), placeholder skeleton.
//
// The mathematical form is
//   n x curl(E) + j k_0 (n x (n x E)) = 0  on faces tagged absorbing
// which adds a surface mass-like rank-N operator to the system.
//
// Status: NOT IMPLEMENTED. Constructing it is a no-op; calling apply()
// throws. Filed in this PR so future work has a concrete file to extend
// without reshuffling the BC factory.
class AbsorbingBC : public IBoundaryCondition {
public:
    AbsorbingBC() = default;

    void apply(SparseMatrix& matrix,
               std::vector<std::complex<double>>& rhs,
               const AssemblyContext& ctx) const override;

    void declareSparsity(SparsePatternBuilder& builder,
                         const AssemblyContext& ctx) const override;

    const char* name() const override { return "AbsorbingBC"; }
};

}  // namespace fem::bc

