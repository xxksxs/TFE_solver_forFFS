#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/linalg/SparseMatrix.hpp"
#include "bpfem/sweep/ISweepStrategy.hpp"

#include <complex>
#include <cstddef>
#include <vector>

namespace fem::sweep {

// Adaptive Lanczos-Pade Sweep (ALPS) MVP.
//
// First iteration: single expansion point + block Krylov (shift-and-invert) +
// complex-symmetric Galerkin projection. Multi-point expansion, residual-driven
// adaptivity, and explicit pole-residue extraction are deferred.
//
// Mathematics matches docs/optimization/alps-sweep/theory-cn.tex section 5
// (single-point form) restricted to lossless materials. With sigma == 0 the
// FEM system A(omega) factors as
//     A(omega) = K - k0(omega)^2 * M + j sum_p beta_p(omega) m_p m_p^T
// where K and M are real symmetric and beta_p(omega) = sqrt(k0^2 - k_{c,p}^2).
// We pick an expansion point omega_0 inside the sweep band, factorize
//     A_0 = A(omega_0)
// once with PARDISO, build the right block-Krylov subspace
//     V = orth([ A_0^{-1} m_1, ..., A_0^{-1} m_{Np}, A_0^{-1} M V_1, ... ])
// of dimension q, and project to obtain the small complex-symmetric ROM
//     A_tilde(omega) = V^T K V - k0(omega)^2 V^T M V + j sum_p beta_p V^T m_p m_p^T V.
// Each online frequency only solves a q x q dense complex linear system.
struct AlpsOptions {
    int krylovOrder = 30;            // q (max columns) per port; total ROM size ~ Np * q
    int maxRestarts = 1;             // currently unused (single-pass MGS)
    double dropTolerance = 1.0e-12;  // deflation threshold for orthogonalization

    // Expansion frequency in Hz. 0 means "use band center" and `run()` will
    // default to 0.5 * (frequencies.front() + frequencies.back()).
    double expansionFrequencyHz = 0.0;
};

class AlpsSweep : public ISweepStrategy {
public:
    AlpsSweep(const ProjectDefinition& project,
              const FEMAssembler& assembler,
              const PortModeSolver& portModeSolver,
              AlpsOptions options = {});

    // ISweepStrategy entry point. Builds the offline ROM at the configured
    // expansion frequency (or band center if 0), then evaluates every
    // requested frequency. SweepResult::lastEdgeDofs is left empty: ALPS
    // does not retain a full-space basis in this MVP, which is why
    // `--sweep alps` cannot write per-frequency VTU. ctx.onFieldSolved is
    // ignored.
    SweepResult run(const std::vector<double>& frequencies,
                    const SweepContext& ctx) override;

    const char* name() const override { return "alps"; }

    // Exposed for unit tests and offline diagnostics.
    int buildOffline(double expansionFrequencyHz);
    SParameterPoint evaluate(double frequencyHz) const;

    int dimension() const { return romDim_; }
    double expansionFrequency() const { return expansionFrequencyHz_; }
    bool ready() const { return ready_; }

private:
    using Complex = std::complex<double>;

    const ProjectDefinition& project_;
    const FEMAssembler& assembler_;
    const PortModeSolver& portModeSolver_;
    AlpsOptions options_;

    bool ready_ = false;
    int romDim_ = 0;
    double expansionFrequencyHz_ = 0.0;
    int fullDim_ = 0;
    int numProjectPorts_ = 0;
    int numVirtualPorts_ = 0;

    FEMAssembler::AffineSystem affine_;

    std::vector<Complex> Ktilde_;
    std::vector<Complex> Mtilde_;
    std::vector<std::vector<Complex>> portModeReduced_;
};

}  // namespace fem::sweep

