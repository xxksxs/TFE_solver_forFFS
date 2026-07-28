#pragma once

#include "bpfem/core/Types.hpp"

#include <filesystem>
#include <string>

namespace fem {

class NGMeshParser {
public:
    // 读取 Ansys NGMesh，并移除不参与 FEM 求解的 background 体。
    Mesh parse(const std::filesystem::path& path) const;

private:
    // 读取 body 的编号、名称和体单元数。
    static void parseBody(const std::string& line, Mesh& mesh);
    // 读取紧随 body 记录之后、单独成行的 face_ids。
    static void parseBodyFaceIds(const std::string& line, int bodyId, Mesh& mesh);
    // 读取 body 包围盒的一项，并返回对应的六位完成掩码。
    static unsigned parseBodyBound(const std::string& line, int bodyId, Mesh& mesh);
    // 读取网格节点。
    static void parsePoint(const std::string& line, Mesh& mesh);
    // 读取 facet 及其对应的 AEDT faceId。
    static int parseFacet(const std::string& line, Mesh& mesh);
    // 读取当前 facet 下的表面三角形。
    static void parseSurfaceTriangle(const std::string& line, int currentFacetId, Mesh& mesh);
    // 读取四面体体单元。
    static void parseTetrahedron(const std::string& line, Mesh& mesh);
};

}  // namespace fem
