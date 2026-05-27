#pragma once

#include "bpfem/bc/IBoundaryCondition.hpp"

namespace fem::bc {

// Wave-port boundary condition.
//
// Reads the project's port list from AssemblyContext::project and consumes
// either the multi-mode (TFE) or single-mode (NPM, APM) entry that the
// caller has populated in PortModeSolver. For each retained mode it adds
//
//   matrix  +=  +j * beta * m * m^T              // rank-1 absorber
//   rhs     +=  2j * beta * a_inc * m            // 2-way Robin source
//
// where m = couplingWeights, beta = sqrt(max(0, k0^2 - k_c^2)), and a_inc
// is built from PortDefinition::magnitudeW / phaseDeg combined with the
// per-frequency Poynting normalization factor `s = 1 / sqrt(P_norm)` so
// that S parameters come out in sqrt(W).
class WavePortBC : public IBoundaryCondition {
public:
    WavePortBC() = default;

    void apply(SparseMatrix& matrix,
               std::vector<std::complex<double>>& rhs,
               const AssemblyContext& ctx) const override;

    void declareSparsity(SparsePatternBuilder& builder,
                         const AssemblyContext& ctx) const override;

    std::vector<AffinePortContribution> affineContributions(
        const AssemblyContext& ctx) const override;

    const char* name() const override { return "WavePort"; }
};

}  // namespace fem::bc

