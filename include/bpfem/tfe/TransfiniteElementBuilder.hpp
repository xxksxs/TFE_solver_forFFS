#pragma once

#include "bpfem/apm/AnalyticPortBuilder.hpp"   // RectPortGeometry, detector (informational only)
#include "bpfem/core/Types.hpp"
#include "bpfem/fem/EdgeTopology.hpp"
#include "bpfem/fem/PortModeSolver.hpp"

namespace fem::tfe {

// Compatibility adapter for the former numerical-multimode port API.
// The active full-edge port path uses PortModeSolver and WavePortBC.
// The modal trace replacement is implemented by TfeSystem.
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

    // Optional rectangular-cross-section detection for geometry reporting.
    apm::RectPortGeometry detectRectangular(int faceId) const;

private:
    const Mesh& mesh_;
    const EdgeTopology& topology_;
    const PortModeSolver& portModeSolver_;
    TfeOptions options_;
};

}  // namespace fem::tfe

