#include "bpfem/fem/EdgeTopology.hpp"

#include "bpfem/core/Math.hpp"

#include <algorithm>
#include <stdexcept>

namespace fem {

EdgeTopology::EdgeTopology(const Mesh& mesh, int basisOrder) : mesh_(mesh), basisOrder_(basisOrder) {
    if (basisOrder_ < 0 || basisOrder_ > 1) {
        throw std::runtime_error("Unsupported basis order: " + std::to_string(basisOrder_));
    }
    elementDofs_.reserve(mesh_.tetrahedra.size());
    for (const auto& tet : mesh_.tetrahedra) {
        std::vector<LocalDofRef> refs;
        refs.reserve(localDofCount());
        for (const auto& pair : localPairs) {
            addLocalEdgeDofs(refs, tet, pair[0], pair[1]);
        }
        if (basisOrder_ >= 1) {
            for (const auto& face : localFaces) {
                addLocalFaceDofs(refs, tet, face);
            }
        }
        elementDofs_.push_back(std::move(refs));
    }
    for (const auto& tri : mesh_.surfaceTriangles) {
        const auto fit = mesh_.facets.find(tri.facetId);
        if (fit == mesh_.facets.end()) {
            continue;
        }
        addFaceTraceDofs(fit->second.faceId, tri.vertexIds);
    }
}

int EdgeTopology::basisOrder() const {
    return basisOrder_;
}

std::size_t EdgeTopology::edgeCount() const {
    return static_cast<std::size_t>(nextDof_);
}

std::size_t EdgeTopology::geometricEdgeCount() const {
    return edges_.size();
}

std::size_t EdgeTopology::localDofCount() const {
    return basisOrder_ == 0 ? 6 : 20;
}

const std::vector<EdgeInfo>& EdgeTopology::edges() const {
    return edges_;
}

const std::vector<std::vector<LocalDofRef>>& EdgeTopology::elementDofs() const {
    return elementDofs_;
}

std::unordered_set<int> EdgeTopology::edgesOnFace(int faceId) const {
    const auto it = faceDofs_.find(faceId);
    if (it == faceDofs_.end()) {
        return {};
    }
    return it->second;
}

int EdgeTopology::edgeZeroDofForNodes(int aNode, int bNode) const {
    const auto it = edgeIndex_.find(edgeKey(aNode, bNode));
    if (it == edgeIndex_.end()) {
        return -1;
    }
    return edges_[static_cast<std::size_t>(it->second)].zeroDof;
}

std::pair<int, int> EdgeTopology::edgeDofsForNodes(int aNode, int bNode) const {
    const auto it = edgeIndex_.find(edgeKey(aNode, bNode));
    if (it == edgeIndex_.end()) {
        return {-1, -1};
    }
    const auto& edge = edges_[static_cast<std::size_t>(it->second)];
    return {edge.zeroDof, edge.firstDof};
}

std::array<int, 2> EdgeTopology::faceFirstDofsForNodes(int aNode, int bNode, int cNode) const {
    const auto it = faceIndex_.find(faceKey(aNode, bNode, cNode));
    if (it == faceIndex_.end()) {
        return {{-1, -1}};
    }
    return faces_[static_cast<std::size_t>(it->second)].firstDofs;
}

void EdgeTopology::addLocalEdgeDofs(std::vector<LocalDofRef>& refs, const Tetrahedron& tet, int localA, int localB) {
    const int aNode = tet.vertexIds[static_cast<std::size_t>(localA)];
    const int bNode = tet.vertexIds[static_cast<std::size_t>(localB)];
    const int edgeIndex = addEdge(aNode, bNode);
    const EdgeInfo& edge = edges_[static_cast<std::size_t>(edgeIndex)];
    const int orientedA = aNode < bNode ? localA : localB;
    const int orientedB = aNode < bNode ? localB : localA;
    refs.push_back({edge.zeroDof, LocalDofKind::EdgeZero, {{orientedA, orientedB, -1}}, 1.0});
    if (basisOrder_ >= 1) {
        refs.push_back({edge.firstDof, LocalDofKind::EdgeFirst, {{orientedA, orientedB, -1}}, 1.0});
    }
}

void EdgeTopology::addLocalFaceDofs(std::vector<LocalDofRef>& refs, const Tetrahedron& tet, const std::array<int, 3>& localFace) {
    std::array<std::pair<int, int>, 3> nodes{{
        {tet.vertexIds[static_cast<std::size_t>(localFace[0])], localFace[0]},
        {tet.vertexIds[static_cast<std::size_t>(localFace[1])], localFace[1]},
        {tet.vertexIds[static_cast<std::size_t>(localFace[2])], localFace[2]}
    }};
    std::sort(nodes.begin(), nodes.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });
    const int faceIndex = addFace(nodes[0].first, nodes[1].first, nodes[2].first);
    const FaceInfo& face = faces_[static_cast<std::size_t>(faceIndex)];
    const std::array<int, 3> sortedLocalFace{{nodes[0].second, nodes[1].second, nodes[2].second}};
    refs.push_back({face.firstDofs[0], LocalDofKind::FaceFirst0, sortedLocalFace, 1.0});
    refs.push_back({face.firstDofs[1], LocalDofKind::FaceFirst1, sortedLocalFace, 1.0});
}

void EdgeTopology::addFaceTraceDofs(int faceId, const std::array<int, 3>& vertexIds) {
    addFaceEdge(faceId, vertexIds[0], vertexIds[1]);
    addFaceEdge(faceId, vertexIds[1], vertexIds[2]);
    addFaceEdge(faceId, vertexIds[2], vertexIds[0]);
    if (basisOrder_ >= 1) {
        const int faceIndex = addFace(vertexIds[0], vertexIds[1], vertexIds[2]);
        const FaceInfo& face = faces_[static_cast<std::size_t>(faceIndex)];
        faceDofs_[faceId].insert(face.firstDofs[0]);
        faceDofs_[faceId].insert(face.firstDofs[1]);
    }
}

void EdgeTopology::addFaceEdge(int faceId, int aNode, int bNode) {
    const int edgeIndex = addEdge(aNode, bNode);
    const EdgeInfo& edge = edges_[static_cast<std::size_t>(edgeIndex)];
    faceDofs_[faceId].insert(edge.zeroDof);
    if (basisOrder_ >= 1) {
        faceDofs_[faceId].insert(edge.firstDof);
    }
}

int EdgeTopology::addEdge(int aNode, int bNode) {
    const auto key = edgeKey(aNode, bNode);
    const auto it = edgeIndex_.find(key);
    if (it != edgeIndex_.end()) {
        return it->second;
    }
    const int lo = std::min(aNode, bNode);
    const int hi = std::max(aNode, bNode);
    const Vec3& p0 = mesh_.pointsById.at(lo);
    const Vec3& p1 = mesh_.pointsById.at(hi);
    const int index = static_cast<int>(edges_.size());
    EdgeInfo edge;
    edge.v0 = lo;
    edge.v1 = hi;
    edge.length = norm(p1 - p0);
    edge.zeroDof = nextDof_++;
    if (basisOrder_ >= 1) {
        edge.firstDof = nextDof_++;
    }
    edges_.push_back(edge);
    edgeIndex_[key] = index;
    return index;
}

int EdgeTopology::addFace(int aNode, int bNode, int cNode) {
    const auto key = faceKey(aNode, bNode, cNode);
    const auto it = faceIndex_.find(key);
    if (it != faceIndex_.end()) {
        return it->second;
    }
    std::array<int, 3> nodes{{aNode, bNode, cNode}};
    std::sort(nodes.begin(), nodes.end());
    FaceInfo face;
    face.vertices = nodes;
    face.firstDofs[0] = nextDof_++;
    face.firstDofs[1] = nextDof_++;
    const int index = static_cast<int>(faces_.size());
    faces_.push_back(face);
    faceIndex_[key] = index;
    return index;
}

std::uint64_t EdgeTopology::edgeKey(int aNode, int bNode) {
    const int lo = std::min(aNode, bNode);
    const int hi = std::max(aNode, bNode);
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lo)) << 32U) | static_cast<std::uint32_t>(hi);
}

std::uint64_t EdgeTopology::faceKey(int aNode, int bNode, int cNode) {
    std::array<int, 3> nodes{{aNode, bNode, cNode}};
    std::sort(nodes.begin(), nodes.end());
    std::uint64_t seed = 1469598103934665603ULL;
    for (int node : nodes) {
        seed ^= static_cast<std::uint32_t>(node);
        seed *= 1099511628211ULL;
    }
    return seed;
}

}  // namespace fem
