#pragma once

#include "bpfem/core/Types.hpp"

#include <filesystem>
#include <optional>
#include <string>

namespace fem {

class AEDTParser {
public:
    ProjectDefinition parse(const std::filesystem::path& path) const;

private:
    static void parseVariables(const std::string& text, ProjectDefinition& project);
    static void parseMaterials(const std::string& text, ProjectDefinition& project);
    static void parseSweep(const std::string& text, ProjectDefinition& project);
    static void parsePorts(const std::string& text, ProjectDefinition& project);
    static std::optional<std::string> captureValue(const std::string& text, const std::string& pattern);
    static double valueOr(const std::string& text, const std::string& pattern, double fallback);
};

}  // namespace fem
