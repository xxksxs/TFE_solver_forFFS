#include "bpfem/fastsweep/PortModeUtilities.hpp"

#include <cstddef>

namespace fem::fastsweep {

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
