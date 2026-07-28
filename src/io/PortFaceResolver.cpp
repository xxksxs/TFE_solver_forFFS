#include "bpfem/io/PortFaceResolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace fem {

namespace {

struct Bounds {
    Vec3 min;
    Vec3 max;
    bool initialized = false;
};

// 把一个节点并入包围盒。
void includePoint(Bounds& bounds, const Vec3& point) {
    if (!bounds.initialized) {
        bounds.min = point;
        bounds.max = point;
        bounds.initialized = true;
        return;
    }
    bounds.min.x = std::min(bounds.min.x, point.x);
    bounds.min.y = std::min(bounds.min.y, point.y);
    bounds.min.z = std::min(bounds.min.z, point.z);
    bounds.max.x = std::max(bounds.max.x, point.x);
    bounds.max.y = std::max(bounds.max.y, point.y);
    bounds.max.z = std::max(bounds.max.z, point.z);
}

// 返回两个三维包围盒六个边界坐标中的最大绝对差。
double boundsError(const Bounds& candidate, const BodyInfo& body) {
    return std::max({
        std::abs(candidate.min.x - body.boundsMin.x),
        std::abs(candidate.min.y - body.boundsMin.y),
        std::abs(candidate.min.z - body.boundsMin.z),
        std::abs(candidate.max.x - body.boundsMax.x),
        std::abs(candidate.max.y - body.boundsMax.y),
        std::abs(candidate.max.z - body.boundsMax.z)
    });
}

// 从保留下来的表面三角形计算每个 AEDT faceId 的几何包围盒。
std::map<int, Bounds> buildFaceBounds(const Mesh& mesh) {
    std::map<int, Bounds> result;
    for (const auto& triangle : mesh.surfaceTriangles) {
        const auto facetIt = mesh.facets.find(triangle.facetId);
        if (facetIt == mesh.facets.end()) {
            continue;
        }
        Bounds& bounds = result[facetIt->second.faceId];
        for (int vertexId : triangle.vertexIds) {
            const auto pointIt = mesh.pointsById.find(vertexId);
            if (pointIt != mesh.pointsById.end()) {
                includePoint(bounds, pointIt->second);
            }
        }
    }
    return result;
}

// 根据模型整体尺寸生成严格但可容纳文本舍入误差的几何比较容差。
double matchingTolerance(const std::map<int, Bounds>& faceBounds) {
    Bounds global;
    for (const auto& [faceId, bounds] : faceBounds) {
        (void)faceId;
        if (!bounds.initialized) {
            continue;
        }
        includePoint(global, bounds.min);
        includePoint(global, bounds.max);
    }
    if (!global.initialized) {
        return 1.0e-12;
    }
    const double dx = global.max.x - global.min.x;
    const double dy = global.max.y - global.min.y;
    const double dz = global.max.z - global.min.z;
    const double scale = std::sqrt(dx * dx + dy * dy + dz * dz);
    return std::max(1.0e-12, scale * 1.0e-8);
}

// 构造包含最近候选面的诊断信息，避免端口匹配失败时只能看到 faceId=-1。
std::string nearestFaceMessage(const std::map<int, Bounds>& faceBounds, const BodyInfo& body) {
    int nearestFaceId = -1;
    double nearestError = std::numeric_limits<double>::infinity();
    for (const auto& [faceId, bounds] : faceBounds) {
        if (!bounds.initialized) {
            continue;
        }
        const double error = boundsError(bounds, body);
        if (error < nearestError) {
            nearestError = error;
            nearestFaceId = faceId;
        }
    }

    std::ostringstream out;
    if (nearestFaceId >= 0) {
        out << "; nearest face=" << nearestFaceId << ", bounds error=" << nearestError;
    }
    return out.str();
}

}  // namespace

void PortFaceResolver::resolve(ProjectDefinition& project, const Mesh& mesh) {
    const auto faceBounds = buildFaceBounds(mesh);
    const double tolerance = matchingTolerance(faceBounds);

    for (auto& port : project.ports) {
        if (port.faceId >= 0) {
            if (faceBounds.find(port.faceId) == faceBounds.end()) {
                throw std::runtime_error(
                    "PortFaceResolver: AEDT port " + std::to_string(port.id)
                    + " references face " + std::to_string(port.faceId)
                    + ", but that face has no retained NGMesh surface triangles");
            }
            continue;
        }
        if (port.objectId < 0) {
            throw std::runtime_error(
                "PortFaceResolver: AEDT port " + std::to_string(port.id)
                + " contains neither Faces(...) nor Objects(...)");
        }

        const auto bodyIt = mesh.bodies.find(port.objectId);
        if (bodyIt == mesh.bodies.end() || !bodyIt->second.hasBounds) {
            throw std::runtime_error(
                "PortFaceResolver: sheet object " + std::to_string(port.objectId)
                + " for AEDT port " + std::to_string(port.id)
                + " has no NGMesh body bounds");
        }

        std::vector<int> matches;
        for (const auto& [faceId, bounds] : faceBounds) {
            if (bounds.initialized && boundsError(bounds, bodyIt->second) <= tolerance) {
                matches.push_back(faceId);
            }
        }
        if (matches.size() != 1) {
            throw std::runtime_error(
                "PortFaceResolver: sheet object " + std::to_string(port.objectId)
                + " for AEDT port " + std::to_string(port.id)
                + " matched " + std::to_string(matches.size())
                + " NGMesh faces; expected exactly one"
                + nearestFaceMessage(faceBounds, bodyIt->second));
        }
        port.faceId = matches.front();
    }
}

}  // namespace fem
