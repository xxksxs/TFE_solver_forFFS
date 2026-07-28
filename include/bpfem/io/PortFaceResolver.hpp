#pragma once

#include "bpfem/core/Types.hpp"

namespace fem {

class PortFaceResolver {
public:
    // 校验直接 face 端口，并把 sheet object 端口映射到唯一重合的网格表面。
    static void resolve(ProjectDefinition& project, const Mesh& mesh);
};

}  // namespace fem
