#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace fem {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Material {
    std::string name;
    double relativePermittivity = 1.0;
    double relativePermeability = 1.0;
    double conductivity = 0.0;
};

struct PortDefinition {
    int id = -1;
    int faceId = -1;
    int modes = 1;
    bool excited = false;
    double magnitudeW = 0.0;
    double phaseDeg = 0.0;
    // HFSS 允许把波端口绑定到 sheet object；解析网格后再由其几何范围解析 faceId。
    int objectId = -1;
};

struct SweepDefinition {
    double setupFrequencyHz = 0.0;
    double startHz = 0.0;
    double endHz = 0.0;
    int count = 1;
};

struct ProjectDefinition {
    std::map<std::string, Material> materials;
    std::map<std::string, double> variablesSI;
    SweepDefinition sweep;
    std::vector<PortDefinition> ports;
    std::string backgroundMaterial = "vacuum";
};

struct BodyInfo {
    int id = -1;
    std::string name;
    int volumeElementCount = 0;
    std::vector<int> faceIds;
    Vec3 boundsMin;
    Vec3 boundsMax;
    bool hasBounds = false;
};

struct SurfaceTriangle {
    int id = -1;
    int facetId = -1;
    std::array<int, 3> vertexIds{};
};

struct FacetInfo {
    int id = -1;
    int faceId = -1;
    std::vector<int> surfaceTriangleIds;
};

struct Tetrahedron {
    int id = -1;
    int bodyId = -1;
    std::array<int, 4> vertexIds{};
};

struct Mesh {
    std::string unitName;
    double userUnitsPerMeter = 1.0;
    std::map<int, BodyInfo> bodies;
    std::map<int, Vec3> pointsById;
    std::vector<int> pointIds;
    std::map<int, FacetInfo> facets;
    std::vector<SurfaceTriangle> surfaceTriangles;
    std::vector<Tetrahedron> tetrahedra;

    std::unordered_map<int, std::size_t> pointIndex() const;
};

struct SolveResult {
    std::vector<std::complex<double>> field;
    int iterations = 0;
    double residual = 0.0;
};

struct SParameterPoint {
    double frequencyHz = 0.0;
    std::complex<double> s11;
    std::complex<double> s21;
};

}  // namespace fem
