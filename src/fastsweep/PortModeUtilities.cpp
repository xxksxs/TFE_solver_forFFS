#include "bpfem/fastsweep/PortModeUtilities.hpp"

#include <cstddef>

namespace fem::fastsweep {

// 根据虚拟端口索引返回对应端口模式；多模端口按同一 face 的出现顺序定位模式。
const PortMode& virtualPortMode(const PortModeSolver& portModeSolver,
                                const FEMAssembler::AffineSystem& affine,
                                int virtualPortIndex) {
    const int faceId = affine.portFaceIds[static_cast<std::size_t>(virtualPortIndex)];
    if (const auto* multi = portModeSolver.multiMode(faceId); multi != nullptr) {
        int modeIndex = 0;
        for (int i = 0; i < virtualPortIndex; ++i) {
            if (affine.portFaceIds[static_cast<std::size_t>(i)] == faceId) {
                ++modeIndex;
            }
        }
        return multi->modes[static_cast<std::size_t>(modeIndex)];
    }
    return portModeSolver.solve(faceId);
}

// 为每个工程端口找出用于 S 参数提取的主导虚拟端口。
std::vector<int> dominantVirtualPortsByProject(const PortModeSolver& portModeSolver,
                                               const FEMAssembler::AffineSystem& affine,
                                               int projectPortCount) {
    std::vector<int> dominant(static_cast<std::size_t>(projectPortCount), -1);
    const int virtualCount = static_cast<int>(affine.portFaceIds.size());
    for (int v = 0; v < virtualCount; ++v) {
        const int faceId = affine.portFaceIds[static_cast<std::size_t>(v)];
        const int projectIndex = affine.projectPortIndex[static_cast<std::size_t>(v)];
        if (projectIndex < 0 || projectIndex >= projectPortCount) {
            continue;
        }
        if (const auto* multi = portModeSolver.multiMode(faceId); multi != nullptr) {
            int modeIndex = 0;
            for (int u = 0; u < v; ++u) {
                if (affine.portFaceIds[static_cast<std::size_t>(u)] == faceId) {
                    ++modeIndex;
                }
            }
            if (modeIndex == multi->excitationModeIndex) {
                dominant[static_cast<std::size_t>(projectIndex)] = v;
            }
        } else {
            dominant[static_cast<std::size_t>(projectIndex)] = v;
        }
    }
    return dominant;
}

}  // namespace fem::fastsweep
