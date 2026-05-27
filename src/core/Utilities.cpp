#include "bpfem/core/Utilities.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace fem {

std::string readTextFile(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot open file: " + path.string());
    }
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

std::string trim(std::string s) {
    auto notSpace = [](unsigned char ch) { return !std::isspace(ch); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
    return s;
}

double parseEngineeringValue(std::string value) {
    value = trim(value);
    value.erase(std::remove(value.begin(), value.end(), '\''), value.end());
    value.erase(std::remove(value.begin(), value.end(), '"'), value.end());

    static const std::regex pattern(R"(^\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*([A-Za-z]*)\s*$)");
    std::smatch match;
    if (!std::regex_match(value, match, pattern)) {
        throw std::runtime_error("Unable to parse engineering value: " + value);
    }

    const double numeric = std::stod(match[1].str());
    const std::string unit = match[2].str();

    static const std::unordered_map<std::string, double> multipliers = {
        {"", 1.0}, {"m", 1.0}, {"mm", 1.0e-3}, {"um", 1.0e-6}, {"nm", 1.0e-9},
        {"Hz", 1.0}, {"kHz", 1.0e3}, {"MHz", 1.0e6}, {"GHz", 1.0e9},
        {"THz", 1.0e12}, {"deg", 1.0}, {"W", 1.0}, {"cel", 1.0}, {"ps", 1.0e-12}
    };

    auto it = multipliers.find(unit);
    return numeric * (it == multipliers.end() ? 1.0 : it->second);
}

std::vector<double> buildFrequencies(const SweepDefinition& sweep, int maxPoints) {
    const int count = std::max(1, maxPoints > 0 ? std::min(maxPoints, sweep.count) : sweep.count);
    std::vector<double> frequencies;
    frequencies.reserve(static_cast<std::size_t>(count));
    if (count == 1) {
        frequencies.push_back(sweep.setupFrequencyHz > 0.0 ? sweep.setupFrequencyHz : sweep.startHz);
        return frequencies;
    }
    for (int i = 0; i < count; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(count - 1);
        frequencies.push_back((1.0 - t) * sweep.startHz + t * sweep.endHz);
    }
    return frequencies;
}

}  // namespace fem
