// FiniteConductorBC: Leontovich surface impedance for a good but finite
// conductor (skin-effect approximation), placeholder skeleton.
//
// Status: same as ImpedanceBC -- apply() throws, declareSparsity() is a
// no-op. The mathematics differs only in that Z_s is *frequency-dependent*
// and computed from the conductor's intrinsic impedance:
//     Z_s = (1 + j) sqrt(omega mu_0 / (2 sigma_metal))
// so the implementation should compute Z_s at ctx.frequencyHz before
// scaling the surface mass contribution. Otherwise the assembly pattern
// matches ImpedanceBC: walk the tagged faces, write a complex surface mass
// integral into the LHS, leave rhs alone.

#include "bpfem/bc/FiniteConductorBC.hpp"

#include <stdexcept>

namespace fem::bc {

void FiniteConductorBC::apply(SparseMatrix& /*matrix*/,
                              std::vector<std::complex<double>>& /*rhs*/,
                              const AssemblyContext& /*ctx*/) const {
    throw std::runtime_error("FiniteConductorBC is not implemented yet");
}

void FiniteConductorBC::declareSparsity(SparsePatternBuilder& /*builder*/,
                                        const AssemblyContext& /*ctx*/) const {
    // Intentionally empty: apply() throws.
}

}  // namespace fem::bc
