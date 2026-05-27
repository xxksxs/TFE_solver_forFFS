#include "bpfem/io/NGMeshParser.hpp"

#include "bpfem/core/Utilities.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace fem {

namespace {

std::string lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

bool isBackgroundBody(const BodyInfo& body) {
    return lowercase(body.name) == "background";
}

void removeBackgroundMesh(Mesh& mesh) {
    std::unordered_set<int> backgroundBodyIds;
    std::unordered_set<int> backgroundFaceIds;
    for (const auto& [bodyId, body] : mesh.bodies) {
        if (!isBackgroundBody(body)) {
            continue;
        }
        backgroundBodyIds.insert(bodyId);
        backgroundFaceIds.insert(body.faceIds.begin(), body.faceIds.end());
    }
    if (backgroundBodyIds.empty()) {
        return;
    }

    mesh.tetrahedra.erase(std::remove_if(mesh.tetrahedra.begin(), mesh.tetrahedra.end(), [&](const Tetrahedron& tet) {
        return backgroundBodyIds.count(tet.bodyId) != 0;
    }), mesh.tetrahedra.end());

    mesh.surfaceTriangles.erase(std::remove_if(mesh.surfaceTriangles.begin(), mesh.surfaceTriangles.end(), [&](const SurfaceTriangle& tri) {
        const auto fit = mesh.facets.find(tri.facetId);
        return fit != mesh.facets.end() && backgroundFaceIds.count(fit->second.faceId) != 0;
    }), mesh.surfaceTriangles.end());

    for (auto it = mesh.facets.begin(); it != mesh.facets.end();) {
        if (backgroundFaceIds.count(it->second.faceId) != 0) {
            it = mesh.facets.erase(it);
        } else {
            it->second.surfaceTriangleIds.clear();
            ++it;
        }
    }
    for (const auto& tri : mesh.surfaceTriangles) {
        const auto fit = mesh.facets.find(tri.facetId);
        if (fit != mesh.facets.end()) {
            fit->second.surfaceTriangleIds.push_back(tri.id);
        }
    }

    for (auto it = mesh.bodies.begin(); it != mesh.bodies.end();) {
        if (backgroundBodyIds.count(it->first) != 0) {
            it = mesh.bodies.erase(it);
        } else {
            ++it;
        }
    }

    std::unordered_set<int> usedPointIds;
    for (const auto& tet : mesh.tetrahedra) {
        usedPointIds.insert(tet.vertexIds.begin(), tet.vertexIds.end());
    }
    for (const auto& tri : mesh.surfaceTriangles) {
        usedPointIds.insert(tri.vertexIds.begin(), tri.vertexIds.end());
    }

    mesh.pointIds.erase(std::remove_if(mesh.pointIds.begin(), mesh.pointIds.end(), [&](int pointId) {
        return usedPointIds.count(pointId) == 0;
    }), mesh.pointIds.end());
    for (auto it = mesh.pointsById.begin(); it != mesh.pointsById.end();) {
        if (usedPointIds.count(it->first) == 0) {
            it = mesh.pointsById.erase(it);
        } else {
            ++it;
        }
    }
}

}  // namespace

Mesh NGMeshParser::parse(const std::filesystem::path& path) const {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open mesh file: " + path.string());
    }

    Mesh mesh;
    std::string line;
    int currentFacetId = -1;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) {
            continue;
        }

        if (line.rfind("user_unit_name", 0) == 0) {
            std::istringstream ss(line);
            std::string key;
            ss >> key >> mesh.unitName;
        } else if (line.rfind("user_units_per_one_meter", 0) == 0) {
            std::istringstream ss(line);
            std::string key;
            ss >> key >> mesh.userUnitsPerMeter;
        } else if (line.rfind("body_id", 0) == 0) {
            parseBody(line, mesh);
        } else if (line.rfind("pid", 0) == 0) {
            parsePoint(line, mesh);
        } else if (line.rfind("facet_id", 0) == 0) {
            currentFacetId = parseFacet(line, mesh);
        } else if (line.rfind("seid", 0) == 0) {
            parseSurfaceTriangle(line, currentFacetId, mesh);
        } else if (line.rfind("veid", 0) == 0) {
            parseTetrahedron(line, mesh);
        }
    }
    std::sort(mesh.pointIds.begin(), mesh.pointIds.end());
    removeBackgroundMesh(mesh);
    return mesh;
}

void NGMeshParser::parseBody(const std::string& line, Mesh& mesh) {
    std::istringstream ss(line);
    BodyInfo body;
    std::string token;
    ss >> token >> body.id;
    while (ss >> token) {
        if (token == "body_name") {
            ss >> body.name;
        } else if (token == "nvelems_on_body") {
            ss >> body.volumeElementCount;
        } else if (token == "face_ids") {
            int faceId = 0;
            while (ss >> faceId) {
                body.faceIds.push_back(faceId);
            }
        }
    }
    if (body.id >= 0) {
        mesh.bodies[body.id] = body;
    }
}

void NGMeshParser::parsePoint(const std::string& line, Mesh& mesh) {
    std::istringstream ss(line);
    std::string token;
    int id = -1;
    Vec3 p;
    ss >> token >> id >> token >> p.x >> p.y >> p.z;
    if (id > 0) {
        mesh.pointsById[id] = p;
        mesh.pointIds.push_back(id);
    }
}

int NGMeshParser::parseFacet(const std::string& line, Mesh& mesh) {
    std::istringstream ss(line);
    std::string token;
    FacetInfo facet;
    ss >> token >> facet.id;
    while (ss >> token) {
        if (token == "face_ids") {
            ss >> facet.faceId;
            break;
        }
    }
    if (facet.id >= 0) {
        mesh.facets[facet.id] = facet;
    }
    return facet.id;
}

void NGMeshParser::parseSurfaceTriangle(const std::string& line, int currentFacetId, Mesh& mesh) {
    std::istringstream ss(line);
    std::string token;
    SurfaceTriangle tri;
    tri.facetId = currentFacetId;
    ss >> token >> tri.id;
    while (ss >> token) {
        if (token == "facet_id") {
            ss >> tri.facetId;
        } else if (token == "vert_ids") {
            ss >> tri.vertexIds[0] >> tri.vertexIds[1] >> tri.vertexIds[2];
            break;
        }
    }
    mesh.surfaceTriangles.push_back(tri);
    mesh.facets[tri.facetId].surfaceTriangleIds.push_back(tri.id);
}

void NGMeshParser::parseTetrahedron(const std::string& line, Mesh& mesh) {
    std::istringstream ss(line);
    std::string token;
    Tetrahedron tet;
    ss >> token >> tet.id;
    while (ss >> token) {
        if (token == "body_id") {
            ss >> tet.bodyId;
        } else if (token == "vert_ids") {
            ss >> tet.vertexIds[0] >> tet.vertexIds[1] >> tet.vertexIds[2] >> tet.vertexIds[3];
            break;
        }
    }
    mesh.tetrahedra.push_back(tet);
}

}  // namespace fem
