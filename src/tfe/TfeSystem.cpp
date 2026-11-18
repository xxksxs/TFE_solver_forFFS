#include "bpfem/tfe/TfeSystem.hpp"

#include "bpfem/core/Constants.hpp"
#include "bpfem/core/Math.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace fem::tfe {
namespace {

using Face = std::array<int, 3>;

Face sortedFace(Face face) {
    std::sort(face.begin(), face.end());
    return face;
}

struct PortTriangle {
    int faceId;
    const SurfaceTriangle* triangle;
    const Tetrahedron* adjacent = nullptr;
    int opposite = -1;
};

std::map<int, Material> portMaterials(const Mesh& mesh, const ProjectDefinition& project,
                                      const FEMAssembler& fem) {
    std::map<int, Material> materials;
    std::map<int, Vec3> origins;
    std::map<int, Vec3> normals;
    for (const auto& port : project.ports) {
        if (!materials.try_emplace(port.faceId).second) {
            throw std::runtime_error("TFE ports must have distinct mesh faces");
        }
    }
    std::map<Face, PortTriangle> triangles;
    for (const auto& tri : mesh.surfaceTriangles) {
        const auto facet = mesh.facets.find(tri.facetId);
        if (facet == mesh.facets.end() || !materials.count(facet->second.faceId)) continue;
        if (!triangles.emplace(sortedFace(tri.vertexIds),
                               PortTriangle{facet->second.faceId, &tri}).second) {
            throw std::runtime_error("Repeated TFE port triangle");
        }
    }
    for (const auto& tet : mesh.tetrahedra) {
        for (int skip = 0; skip < 4; ++skip) {
            const auto it = triangles.find(sortedFace({tet.vertexIds[(skip + 1) % 4],
                                                       tet.vertexIds[(skip + 2) % 4],
                                                       tet.vertexIds[(skip + 3) % 4]}));
            if (it == triangles.end()) continue;
            if (it->second.adjacent) throw std::runtime_error("TFE port lies on an interior mesh face");
            it->second.adjacent = &tet;
            it->second.opposite = tet.vertexIds[skip];
        }
    }
    for (const auto& entry : triangles) {
        const auto& port = entry.second;
        if (!port.adjacent) throw std::runtime_error("TFE port triangle has no adjacent volume element");
        const auto& a = mesh.pointsById.at(port.triangle->vertexIds[0]);
        const auto& b = mesh.pointsById.at(port.triangle->vertexIds[1]);
        const auto& c = mesh.pointsById.at(port.triangle->vertexIds[2]);
        Vec3 normal = cross(b - a, c - a);
        const double area = norm(normal);
        if (!(area > 0.0)) throw std::runtime_error("Degenerate TFE port triangle");
        normal = normal / area;
        if (dot(normal, mesh.pointsById.at(port.opposite) - a) > 0.0) normal = -1.0 * normal;

        const Material material = fem.materialForBody(port.adjacent->bodyId);
        if (!std::isfinite(material.relativePermittivity) || !(material.relativePermittivity > 0.0)
            || !std::isfinite(material.relativePermeability) || !(material.relativePermeability > 0.0)
            || material.conductivity != 0.0) {
            throw std::runtime_error("TFE requires a homogeneous, lossless waveguide port");
        }
        if (normals.try_emplace(port.faceId, normal).second) {
            materials.at(port.faceId) = material;
            origins.emplace(port.faceId, a);
            continue;
        }
        const auto& referenceNormal = normals.at(port.faceId);
        const auto& origin = origins.at(port.faceId);
        const auto& referenceMaterial = materials.at(port.faceId);
        const double scale = std::max({norm(b - a), norm(c - a), norm(a - origin)});
        if (norm(normal - referenceNormal) > 1.0e-6 ||
            std::abs(dot(a - origin, referenceNormal)) > 1.0e-6 * scale ||
            material.relativePermittivity != referenceMaterial.relativePermittivity ||
            material.relativePermeability != referenceMaterial.relativePermeability) {
            throw std::runtime_error("TFE port is not a planar homogeneous cross-section on face "
                                     + std::to_string(port.faceId));
        }
    }
    if (normals.size() != materials.size() || materials.empty()) {
        throw std::runtime_error("TFE requires resolved, nonempty waveguide port faces");
    }
    return materials;
}

template <class Visitor>
void projectedEntries(const SparseMatrix& K, const TfeSystem::Trace& trace, Visitor visit) {
    const auto& rows = K.rowOffsets();
    const auto& columns = K.columnIndices();
    for (std::size_t row = 0; row < K.size(); ++row) {
        for (int i = rows[row]; i < rows[row + 1]; ++i) {
            const auto column = static_cast<std::size_t>(columns[i]);
            for (const auto& [r, wr] : trace[row]) {
                for (const auto& [c, wc] : trace[column]) {
                    if (column == row && r > c) continue;
                    // An off-diagonal entry of D represents both D_rc and D_cr.
                    // Both land on the same diagonal only if the new indices coincide.
                    const double factor = column != row && r == c ? 2.0 : 1.0;
                    visit(std::min(r, c), std::max(r, c), i, factor * wr * wc);
                }
            }
        }
    }
}

}  // namespace

TfeSystem::TfeSystem(const Mesh& mesh, const ProjectDefinition& project, const EdgeTopology& topology,
                     const FEMAssembler& fem, const PortModeSolver& modes, int modeOverride)
    : trace_(topology.edgeCount()) {
    const auto materials = portMaterials(mesh, project, fem);
    const auto& pec = fem.pecDofs();
    std::vector<int> portByDof(trace_.size(), -1);
    for (const auto& port : project.ports) {
        for (const int dof : topology.edgesOnFace(port.faceId)) {
            if (pec.count(dof)) continue;
            if (portByDof.at(dof) != -1) throw std::runtime_error("Overlapping effective TFE port traces");
            portByDof[dof] = port.faceId;
            ++traceCount_;
        }
    }
    for (std::size_t i = 0; i < trace_.size(); ++i) {
        if (!pec.count(static_cast<int>(i)) && portByDof[i] == -1) {
            trace_[i].emplace_back(static_cast<int>(interiorCount_++), 1.0);
        }
    }
    std::unordered_set<int> portIds;
    for (const auto& port : project.ports) {
        if (!portIds.insert(port.id).second) throw std::runtime_error("Duplicate TFE physical port ID");
        const auto& medium = materials.at(port.faceId);
        const auto available = modes.computeWaveguideModes(
            port.faceId, modeOverride > 0 ? modeOverride : std::max(1, port.modes), pec);
        for (std::size_t i = 0; i < available.size(); ++i) {
            const auto& mode = available[i];
            const std::map<int, double> coefficients(mode.trace.edgeDofs.begin(),
                                                      mode.trace.edgeDofs.end());
            for (std::size_t j = 0; j <= i; ++j) {
                const auto& weights = available[j].trace.couplingWeights;
                double gram = 0.0;
                for (const auto& [dof, value] : weights) {
                    const auto found = coefficients.find(dof);
                    if (found != coefficients.end()) gram += value * found->second;
                }
                if (std::abs(gram - (i == j ? 1.0 : 0.0)) > 1.0e-7) {
                    throw std::runtime_error("TFE modal trace Gram is not orthonormal on face "
                                             + std::to_string(port.faceId));
                }
            }
            const int index = static_cast<int>(interiorCount_ + channels_.size());
            const double kc2 = mode.trace.cutoffWavenumberSquared;
            channels_.push_back({port.id, static_cast<int>(i + 1), mode.transverseMagnetic,
                c0 * std::sqrt(kc2 / (medium.relativePermeability * medium.relativePermittivity)) / (2.0 * pi),
                medium.relativePermittivity, medium.relativePermeability, kc2});
            for (const auto& [dof, value] : mode.trace.edgeDofs) {
                if (portByDof.at(dof) != port.faceId) throw std::runtime_error("TFE mode outside effective trace");
                if (value != 0.0) trace_[dof].emplace_back(index, value);
            }
        }
    }

    const auto volume = fem.buildVolumeSystem();
    const auto n = interiorCount_ + channels_.size();
    std::vector<std::vector<int>> rows(n);
    for (std::size_t r = 0; r < n; ++r) rows[r].push_back(static_cast<int>(r));
    projectedEntries(volume.K, trace_, [&](int row, int col, int, double) {
        rows[row].push_back(col);
    });
    for (auto& row : rows) {
        std::sort(row.begin(), row.end());
        row.erase(std::unique(row.begin(), row.end()), row.end());
    }
    pattern_ = SparseMatrix::createPattern(n, rows);
    stiffness_ = SparseMatrix(pattern_);
    mass_ = SparseMatrix(pattern_);
    linearFrequency_ = SparseMatrix(pattern_);
    projectedEntries(volume.K, trace_, [&](int row, int col, int entry, double weight) {
        stiffness_.add(row, col, weight * volume.K.values()[entry]);
        mass_.add(row, col, weight * volume.M.values()[entry]);
        linearFrequency_.add(row, col, weight *
            (volume.conductivity.values()[entry] / (c0 * c0 * epsilon0) + volume.radiation.values()[entry]));
    });
}

double TfeSystem::cutoffGuard(std::size_t channel) const {
    return 64.0 * std::numeric_limits<double>::epsilon() * channels_.at(channel).cutoffHz;
}

ModeStatus TfeSystem::status(std::size_t channel, double frequencyHz) const {
    const double offset = frequencyHz - channels_.at(channel).cutoffHz;
    if (std::abs(offset) <= cutoffGuard(channel)) return ModeStatus::Cutoff;
    return offset > 0.0 ? ModeStatus::Propagating : ModeStatus::Evanescent;
}

std::complex<double> TfeSystem::admittance(std::size_t index, double frequencyHz) const {
    if (!std::isfinite(frequencyHz) || !(frequencyHz > 0.0)) {
        throw std::invalid_argument("TFE frequency must be positive and finite");
    }
    const auto& channel = channels_.at(index);
    const auto propagation = status(index, frequencyHz);
    if (propagation == ModeStatus::Cutoff) {
        throw std::runtime_error("TFE normalization is singular at cutoff for port "
                                 + std::to_string(channel.portId) + ", mode " + std::to_string(channel.mode));
    }
    const double omega = 2.0 * pi * frequencyHz;
    const double relative = frequencyHz / channel.cutoffHz;
    const double betaSquared = channel.cutoffSquared * (relative - 1.0) * (relative + 1.0);
    const std::complex<double> beta = propagation == ModeStatus::Propagating
        ? std::complex<double>(std::sqrt(betaSquared), 0.0)
        : std::complex<double>(0.0, -std::sqrt(-betaSquared));
    if (channel.transverseMagnetic) return omega * epsilon0 * channel.relativePermittivity / beta;
    return beta / (omega * mu0 * channel.relativePermeability);
}

TfeSystem::Vector TfeSystem::modalScale(double frequencyHz) const {
    Vector scale(pattern_->n, 1.0);
    for (std::size_t i = 0; i < channels_.size(); ++i) {
        scale[interiorCount_ + i] = 1.0 / std::sqrt(admittance(i, frequencyHz));
    }
    return scale;
}

SparseMatrix TfeSystem::assemble(double frequencyHz) const {
    const auto scale = modalScale(frequencyHz);
    const double omega = 2.0 * pi * frequencyHz;
    const double k0 = omega / c0;
    const std::complex<double> jOmegaMu(0.0, omega * mu0);
    SparseMatrix matrix(pattern_);
    for (std::size_t r = 0; r < pattern_->n; ++r) {
        for (int i = pattern_->rowOffsets[r]; i < pattern_->rowOffsets[r + 1]; ++i) {
            const int c = pattern_->columnIndices[i];
            matrix.add(static_cast<int>(r), c, scale[r] * scale[c] *
                (stiffness_.values()[i] - k0 * k0 * mass_.values()[i]
                 + std::complex<double>(0.0, omega) * linearFrequency_.values()[i]));
        }
        if (r >= interiorCount_) matrix.add(static_cast<int>(r), static_cast<int>(r), jOmegaMu);
    }
    return matrix;
}

TfeSystem::Vector TfeSystem::rhs(double frequencyHz, const Vector& incidence) const {
    if (incidence.size() != channels_.size() || !std::isfinite(frequencyHz) || !(frequencyHz > 0.0)) {
        throw std::invalid_argument("Invalid TFE incidence or frequency");
    }
    Vector b(pattern_->n);
    const std::complex<double> gain(0.0, 4.0 * pi * frequencyHz * mu0);
    for (std::size_t i = 0; i < incidence.size(); ++i) {
        if (incidence[i] != 0.0 && status(i, frequencyHz) != ModeStatus::Propagating) {
            throw std::runtime_error("TFE excitation requires a propagating mode");
        }
        b[interiorCount_ + i] = gain * incidence[i];
    }
    return b;
}

TfeSystem::Vector TfeSystem::restoreEdgeField(double frequencyHz, const Vector& state) const {
    if (state.size() != pattern_->n) throw std::invalid_argument("Invalid TFE state size");
    const auto scale = modalScale(frequencyHz);
    Vector edge(trace_.size());
    for (std::size_t i = 0; i < edge.size(); ++i) {
        for (const auto& [j, weight] : trace_[i]) edge[i] += weight * scale[j] * state[j];
    }
    return edge;
}

TfeSystem::Vector TfeSystem::outgoing(const Vector& state, const Vector& incidence) const {
    if (state.size() != pattern_->n || incidence.size() != channels_.size()) {
        throw std::invalid_argument("Invalid TFE state or incidence size");
    }
    Vector result(channels_.size());
    for (std::size_t i = 0; i < result.size(); ++i) result[i] = state[interiorCount_ + i] - incidence[i];
    return result;
}

}  // namespace fem::tfe
