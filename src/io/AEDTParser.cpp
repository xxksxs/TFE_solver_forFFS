#include "bpfem/io/AEDTParser.hpp"

#include "bpfem/core/Utilities.hpp"

#include <algorithm>
#include <regex>

namespace fem {

ProjectDefinition AEDTParser::parse(const std::filesystem::path& path) const {
    const std::string text = readTextFile(path);
    ProjectDefinition project;
    parseVariables(text, project);
    parseMaterials(text, project);
    parseSweep(text, project);
    parsePorts(text, project);
    return project;
}

void AEDTParser::parseVariables(const std::string& text, ProjectDefinition& project) {
    const std::regex variablePattern(R"(VariableProp\('([^']+)'\s*,\s*'[^']*'\s*,\s*'[^']*'\s*,\s*'([^']+)'\))");
    for (std::sregex_iterator it(text.begin(), text.end(), variablePattern), end; it != end; ++it) {
        const std::string name = (*it)[1].str();
        try {
            project.variablesSI[name] = parseEngineeringValue((*it)[2].str());
        } catch (...) {
            project.variablesSI[name] = 0.0;
        }
    }
}

void AEDTParser::parseMaterials(const std::string& text, ProjectDefinition& project) {
    const std::string beginMarker = "$begin 'Materials'";
    const std::string endMarker = "$end 'Materials'";
    const std::size_t begin = text.find(beginMarker);
    const std::size_t end = text.find(endMarker, begin == std::string::npos ? 0 : begin);
    const std::string materialText = begin == std::string::npos || end == std::string::npos
                                   ? text
                                   : text.substr(begin + beginMarker.size(), end - begin - beginMarker.size());
    const std::regex blockPattern(R"(\$begin '([^']+)'([\s\S]*?)\$end '\1')");
    for (std::sregex_iterator it(materialText.begin(), materialText.end(), blockPattern), last; it != last; ++it) {
        Material material;
        material.name = (*it)[1].str();
        const std::string body = (*it)[2].str();
        if (body.find("permittivity=") == std::string::npos && body.find("conductivity=") == std::string::npos) {
            continue;
        }
        const auto eps = captureValue(body, R"(permittivity='([^']+)')");
        const auto mu = captureValue(body, R"(permeability='([^']+)')");
        const auto sigma = captureValue(body, R"(conductivity='([^']+)')");
        if (eps) material.relativePermittivity = parseEngineeringValue(*eps);
        if (mu) material.relativePermeability = parseEngineeringValue(*mu);
        if (sigma) material.conductivity = parseEngineeringValue(*sigma);
        project.materials[material.name] = material;
    }
    if (project.materials.empty()) {
        project.materials["vacuum"] = Material{"vacuum", 1.0, 1.0, 0.0};
    }
}

void AEDTParser::parseSweep(const std::string& text, ProjectDefinition& project) {
    project.sweep.setupFrequencyHz = valueOr(text, R"(Frequency='([^']+)')", 41.5e9);
    project.sweep.startHz = valueOr(text, R"(RangeStart='([^']+)')", project.sweep.setupFrequencyHz);
    project.sweep.endHz = valueOr(text, R"(RangeEnd='([^']+)')", project.sweep.setupFrequencyHz);
    project.sweep.count = static_cast<int>(valueOr(text, R"(RangeCount=([0-9]+))", 1.0));
    project.sweep.count = std::max(1, project.sweep.count);
}

void AEDTParser::parsePorts(const std::string& text, ProjectDefinition& project) {
    const std::regex portPattern(R"(\$begin '([0-9]+)'([\s\S]*?)\$end '\1')");
    for (std::sregex_iterator it(text.begin(), text.end(), portPattern), end; it != end; ++it) {
        const std::string body = (*it)[2].str();
        if (body.find("BoundType='Wave Port'") == std::string::npos) {
            continue;
        }
        PortDefinition port;
        port.id = static_cast<int>(valueOr(body, R"(ID=([0-9]+))", -1.0));
        port.faceId = static_cast<int>(valueOr(body, R"(Faces\(([0-9]+)\))", -1.0));
        port.modes = static_cast<int>(valueOr(body, R"(NumModes=([0-9]+))", 1.0));
        project.ports.push_back(port);
    }

    const std::regex sourcePattern(R"(SourceEntry\(ID=([0-9]+),[^\)]*Magnitude='([^']+)'\s*,\s*Phase='([^']+)')");
    for (std::sregex_iterator it(text.begin(), text.end(), sourcePattern), end; it != end; ++it) {
        const int id = std::stoi((*it)[1].str());
        const double mag = parseEngineeringValue((*it)[2].str());
        const double phase = parseEngineeringValue((*it)[3].str());
        for (auto& port : project.ports) {
            if (port.id == id) {
                port.magnitudeW = mag;
                port.phaseDeg = phase;
                port.excited = mag > 0.0;
            }
        }
    }

    std::sort(project.ports.begin(), project.ports.end(), [](const auto& a, const auto& b) { return a.id < b.id; });
}

std::optional<std::string> AEDTParser::captureValue(const std::string& text, const std::string& pattern) {
    std::smatch match;
    if (std::regex_search(text, match, std::regex(pattern))) {
        return match[1].str();
    }
    return std::nullopt;
}

double AEDTParser::valueOr(const std::string& text, const std::string& pattern, double fallback) {
    const auto value = captureValue(text, pattern);
    if (!value) {
        return fallback;
    }
    try {
        return parseEngineeringValue(*value);
    } catch (...) {
        return fallback;
    }
}

}  // namespace fem
