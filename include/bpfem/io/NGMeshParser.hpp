#pragma once

#include "bpfem/core/Types.hpp"

#include <filesystem>
#include <string>

namespace fem {

class NGMeshParser {
public:
    Mesh parse(const std::filesystem::path& path) const;

private:
    static void parseBody(const std::string& line, Mesh& mesh);
    static void parsePoint(const std::string& line, Mesh& mesh);
    static int parseFacet(const std::string& line, Mesh& mesh);
    static void parseSurfaceTriangle(const std::string& line, int currentFacetId, Mesh& mesh);
    static void parseTetrahedron(const std::string& line, Mesh& mesh);
};

}  // namespace fem
