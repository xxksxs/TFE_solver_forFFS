#pragma once

#include "bpfem/core/Types.hpp"
#include "bpfem/fem/EdgeTopology.hpp"

#include <cstdint>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fem {

struct PortQuadraturePoint {
    Vec3 modeFieldValue;
    Vec3 normal;
    double weight = 0.0;
};

struct PortMode {
    int faceId = -1;
    double cutoffWavenumberSquared = 0.0;
    std::vector<std::pair<int, double>> edgeDofs;
    std::vector<std::pair<int, double>> couplingWeights;
    std::vector<PortQuadraturePoint> quadrature;
};

// Multi-mode port description used by the transfinite-element (TFE) path.
// Carries one PortMode per analytic transverse mode; the FEM assembler treats
// each entry as an independent rank-1 port operator. The single-mode path
// (NPM, APM) still uses the plain PortMode contract above.
struct MultiPortMode {
    int faceId = -1;
    std::vector<PortMode> modes;       // sorted by cutoff k_c^2 ascending
    int excitationModeIndex = 0;       // index into modes; usually TE10 (0)
};

double poyntingPowerIntegral(const PortMode& mode, double frequencyHz);
double powerNormalizationFactor(const PortMode& mode, double frequencyHz);

class PortModeSolver {
public:
    PortModeSolver(const Mesh& mesh, const EdgeTopology& topology);

    const PortMode& solve(int faceId) const;

    // Numerical multi-mode eigensolve on the port-face mesh. Returns the
    // lowest `count` non-spurious modes (ascending k_c^2) computed from the
    // 2D H(curl) generalized eigenproblem K_port v = k_c^2 M_port v on the
    // port triangulation. This is the same solver used internally by
    // computeMode, generalized to multiple modes; it is purely numerical and
    // makes no assumption about the port cross-section shape.
    //
    // Spurious near-zero eigenvalues (the gradient null space) are filtered
    // out with the same threshold rule as solve(faceId). The MultiPortMode
    // returned has its excitationModeIndex set to 0, i.e. the dominant
    // (lowest physical k_c^2) mode.
    MultiPortMode computeMultiMode(int faceId, int count) const;

    // Single-mode injection (used by APM / numerical paths).
    void setPrecomputed(int faceId, PortMode mode) const;

    // Multi-mode injection (used by the TFE path). When this is set, the
    // FEM assembler will iterate over MultiPortMode::modes and add one
    // rank-1 port term per mode. The legacy single-mode cache is also
    // populated with the dominant (excitationModeIndex) entry so existing
    // log / ResultExtractor code paths still see a valid PortMode.
    void setMultiMode(int faceId, MultiPortMode multi) const;

    // Returns nullptr when no multi-mode entry has been registered for faceId.
    const MultiPortMode* multiMode(int faceId) const;

private:
    const Mesh& mesh_;
    const EdgeTopology& topology_;
    std::unordered_map<std::uint64_t, std::pair<int, int>> edgeDofsByNodePair_;
    mutable std::map<int, PortMode> cache_;
    mutable std::map<int, MultiPortMode> multiCache_;

    PortMode computeMode(int faceId) const;
    std::pair<int, int> edgeDofsForNodes(int aNode, int bNode) const;
    static std::uint64_t edgeKey(int aNode, int bNode);
};

}  // namespace fem
