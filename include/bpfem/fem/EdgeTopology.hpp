#pragma once

#include "bpfem/core/Types.hpp"

#include <array>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fem {

enum class LocalDofKind {
    EdgeZero,
    EdgeFirst,
    FaceFirst0,
    FaceFirst1
};

struct EdgeInfo {
    int v0 = -1;
    int v1 = -1;
    double length = 0.0;
    int zeroDof = -1;
    int firstDof = -1;
};

struct FaceInfo {
    std::array<int, 3> vertices{};
    std::array<int, 2> firstDofs{{-1, -1}};
};

struct LocalEdgeRef {
    int globalIndex = -1;
    int localA = -1;
    int localB = -1;
    double sign = 1.0;
};

struct LocalDofRef {
    int globalIndex = -1;
    LocalDofKind kind = LocalDofKind::EdgeZero;
    std::array<int, 3> localNodes{{-1, -1, -1}};
    double sign = 1.0;
};

class EdgeTopology {
public:
    static constexpr std::array<std::array<int, 2>, 6> localPairs = {{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
    static constexpr std::array<std::array<int, 3>, 4> localFaces = {{{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}}};

    explicit EdgeTopology(const Mesh& mesh, int basisOrder = 0);

    int basisOrder() const;
    std::size_t edgeCount() const;
    std::size_t geometricEdgeCount() const;
    std::size_t localDofCount() const;
    const std::vector<EdgeInfo>& edges() const;
    const std::vector<std::vector<LocalDofRef>>& elementDofs() const;
    std::unordered_set<int> edgesOnFace(int faceId) const;
    int edgeZeroDofForNodes(int aNode, int bNode) const;
    std::pair<int, int> edgeDofsForNodes(int aNode, int bNode) const;
    std::array<int, 2> faceFirstDofsForNodes(int aNode, int bNode, int cNode) const;

private:
    const Mesh& mesh_;
    int basisOrder_ = 0;
    int nextDof_ = 0;
    std::vector<EdgeInfo> edges_;
    std::unordered_map<std::uint64_t, int> edgeIndex_;
    std::vector<FaceInfo> faces_;
    std::unordered_map<std::uint64_t, int> faceIndex_;
    std::vector<std::vector<LocalDofRef>> elementDofs_;
    std::map<int, std::unordered_set<int>> faceDofs_;

    void addLocalEdgeDofs(std::vector<LocalDofRef>& refs, const Tetrahedron& tet, int localA, int localB);
    void addLocalFaceDofs(std::vector<LocalDofRef>& refs, const Tetrahedron& tet, const std::array<int, 3>& localFace);
    void addFaceTraceDofs(int faceId, const std::array<int, 3>& vertexIds);
    void addFaceEdge(int faceId, int aNode, int bNode);
    int addEdge(int aNode, int bNode);
    int addFace(int aNode, int bNode, int cNode);
    static std::uint64_t edgeKey(int aNode, int bNode);
    static std::uint64_t faceKey(int aNode, int bNode, int cNode);
};

}  // namespace fem
