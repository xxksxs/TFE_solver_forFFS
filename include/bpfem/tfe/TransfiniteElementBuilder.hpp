#pragma once

#include "bpfem/apm/AnalyticPortBuilder.hpp"   // RectPortGeometry, detector (informational only)
#include "bpfem/core/Types.hpp"
#include "bpfem/fem/EdgeTopology.hpp"
#include "bpfem/fem/PortModeSolver.hpp"

namespace fem::tfe {

// Configuration for the true Transfinite Element (TFE) port-region treatment.
//
// TFE keeps the same per-mode FE-projection machinery the rest of the
// pipeline expects, but admits *multiple* transverse modes per port. Each
// retained mode contributes one rank-1 port operator
//   +j * beta_p^{(k)} * m_p^{(k)} (m_p^{(k)})^T
// to the global system after Schur-eliminating the per-mode coefficients
// alpha (which are diagonal because port modes are L2-orthogonal on the
// cross section).
//
// The mode shapes are computed *numerically* by solving the 2D H(curl)
// generalized eigenproblem on the port surface mesh:
//   K_port v = k_c^2 M_port v.
// This is what HFSS and other commercial codes call the "port mode solver",
// and is what makes TFE physically distinct from the analytic-port-mode
// (APM) path: rather than projecting a closed-form rectangular waveguide
// shape, TFE picks up whatever transverse modes the actual port mesh
// supports (rectangular, ridged, circular, off-axis, etc.).
struct TfeOptions {
    // Maximum number of physical (non-spurious) modes retained per port,
    // sorted ascending by k_c^2. The dominant mode is at index 0 and is
    // used as the excitation mode.
    int modesPerPort = 1;
};

class TransfiniteElementBuilder {
public:
    TransfiniteElementBuilder(const Mesh& mesh,
                              const EdgeTopology& topology,
                              const PortModeSolver& portModeSolver,
                              TfeOptions options = {});

    // Build the multi-mode port description for a single port face by
    // running the numerical 2D H(curl) eigensolve on the port mesh and
    // keeping the lowest options_.modesPerPort non-spurious eigenmodes.
    MultiPortMode build(int faceId) const;

    // Optional rectangular-cross-section detection. Used only for human-
    // readable logging in Application; TFE itself does not require the port
    // to be rectangular.
    apm::RectPortGeometry detectRectangular(int faceId) const;

private:
    const Mesh& mesh_;
    const EdgeTopology& topology_;
    const PortModeSolver& portModeSolver_;
    TfeOptions options_;
};

}  // namespace fem::tfe

