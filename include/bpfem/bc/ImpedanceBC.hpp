#pragma once

#include "bpfem/bc/IBoundaryCondition.hpp"

#include <complex>

namespace fem::bc {

// Surface-impedance (Z_s) boundary condition, placeholder skeleton.
//
// Form (lossy SIBC):
//   n x curl(E) - (j omega mu_0 / Z_s) (n x (n x E)) = 0
// Distinct from FiniteConductorBC in that Z_s is user-specified (e.g. for
// impedance walls in transmission-line problems) rather than derived from
// (sigma, mu_r, freq).
//
// Status: NOT IMPLEMENTED.
class ImpedanceBC : public IBoundaryCondition {
public:
    ImpedanceBC() = default;
    explicit ImpedanceBC(std::complex<double> surfaceImpedance)
        : surfaceImpedance_(surfaceImpedance) {}

    void apply(SparseMatrix& matrix,
               std::vector<std::complex<double>>& rhs,
               const AssemblyContext& ctx) const override;

    void declareSparsity(SparsePatternBuilder& builder,
                         const AssemblyContext& ctx) const override;

    const char* name() const override { return "ImpedanceBC"; }

    std::complex<double> surfaceImpedance() const { return surfaceImpedance_; }

private:
    std::complex<double> surfaceImpedance_{0.0, 0.0};
};

}  // namespace fem::bc

