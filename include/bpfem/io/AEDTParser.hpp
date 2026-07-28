#pragma once

#include "bpfem/core/Types.hpp"

#include <filesystem>
#include <optional>
#include <string>

namespace fem {

class AEDTParser {
public:
    // 读取 AEDT 文本工程中的材料、扫频和波端口定义。
    ProjectDefinition parse(const std::filesystem::path& path) const;

private:
    // 读取工程变量并统一转换为 SI 数值。
    static void parseVariables(const std::string& text, ProjectDefinition& project);
    // 读取各向同性材料的介电、磁导和电导参数。
    static void parseMaterials(const std::string& text, ProjectDefinition& project);
    // 读取求解频率和线性扫频区间。
    static void parseSweep(const std::string& text, ProjectDefinition& project);
    // 同时支持 Faces(...) 和 Objects(...) 两种 HFSS 波端口写法。
    static void parsePorts(const std::string& text, ProjectDefinition& project);
    // 用正则表达式捕获第一个字符串值。
    static std::optional<std::string> captureValue(const std::string& text, const std::string& pattern);
    // 捕获并解析工程量；缺失或格式错误时返回给定默认值。
    static double valueOr(const std::string& text, const std::string& pattern, double fallback);
};

}  // namespace fem
