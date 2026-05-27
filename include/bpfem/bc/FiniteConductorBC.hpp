#pragma once

#include "bpfem/bc/IBoundaryCondition.hpp"

namespace fem::bc {

// Leontovich surface impedance boundary condition for a good but finite
// conductor (skin-effect approximation), placeholder skeleton.
//
// Form:
//   n x curl(E) - (j omega mu_0 / Z_s) (n x (n x E)) = 0
// with the metal's intrinsic impedance
//   Z_s = (1 + j) * sqrt(omega mu_0 / (2 sigma_metal))
// derived from the conductor's conductivity sigma_metal at runtime, hence
// frequency-dependent. Used to model lossy waveguide walls without meshing
// the metal volume.
//
// Status: NOT IMPLEMENTED.
class FiniteConductorBC : public IBoundaryCondition {
public:
    FiniteConductorBC() = default;
    explicit FiniteConductorBC(double conductivitySperM)
        : conductivity_(conductivitySperM) {}

    void apply(SparseMatrix& matrix,
               std::vector<std::complex<double>>& rhs,
               const AssemblyContext& ctx) const override;

    void declareSparsity(SparsePatternBuilder& builder,
                         const AssemblyContext& ctx) const override;

    const char* name() const override { return "FiniteConductorBC"; }

    double conductivity() const { return conductivity_; }

private:
    double conductivity_ = 0.0;
};

}  // namespace fem::bc

