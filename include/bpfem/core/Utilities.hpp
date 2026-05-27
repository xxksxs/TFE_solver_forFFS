#pragma once

#include "bpfem/core/Types.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace fem {

std::string readTextFile(const std::filesystem::path& path);
std::string trim(std::string s);
double parseEngineeringValue(std::string value);
std::vector<double> buildFrequencies(const SweepDefinition& sweep, int maxPoints);

}  // namespace fem
