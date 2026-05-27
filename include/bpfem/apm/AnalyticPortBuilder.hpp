#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fem/EdgeTopology.hpp"
#include "bpfem/fem/PortModeSolver.hpp"

namespace fem::apm {

// Detected rectangular cross-section of a waveguide port. Built once per port
// face by projecting surface triangle vertices onto the port plane and fitting
// an axis-aligned bounding box in PCA-aligned local coordinates.
//
// This type was previously called fem::tfe::RectPortGeometry and lived in the
// (incorrectly named) tfe module. The "real" TFE method now occupies that
// namespace; analytic-port-mode (APM) construction lives here.
struct RectPortGeometry {
    Vec3 origin;     // (u, v) = (0, 0) in 3D
    Vec3 axisU;      // unit vector along width 'a' (UV plane in port)
    Vec3 axisV;      // unit vector along height 'b'
    Vec3 normal;     // outward port normal
    double a = 0.0;  // longer side (m)
    double b = 0.0;  // shorter side (m)
    bool valid = false;
};

// Analytic Port Mode (APM) builder for rectangular waveguides.
//
// MVP scope: dominant TE10 mode only. The analytic mode shape e_t^analytic is
// L2-projected onto the H(curl) port-face DOFs and renormalized so that
// ||e_t^FE||_{L2} = 1. The product is a fem::PortMode object with the same
// contract used by FEMAssembler::applyWavePorts and ResultExtractor, so the
// pipeline is unaware of which path produced the mode. Numerical port-mode
// solver (NPM) remains available internally; --port-method analytic and tfe
// both ultimately use the same FEM assembly contract.
class AnalyticPortBuilder {
public:
    AnalyticPortBuilder(const Mesh& mesh, const EdgeTopology& topology);

    PortMode build(int faceId) const;
    RectPortGeometry detectRectangular(int faceId) const;

private:
    const Mesh& mesh_;
    const EdgeTopology& topology_;
};

}  // namespace fem::apm
