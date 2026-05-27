#include "bpfem/tfe/TransfiniteElementBuilder.hpp"

#include <stdexcept>
#include <string>

namespace fem::tfe {

TransfiniteElementBuilder::TransfiniteElementBuilder(const Mesh& mesh,
                                                     const EdgeTopology& topology,
                                                     const PortModeSolver& portModeSolver,
                                                     TfeOptions options)
    : mesh_(mesh), topology_(topology), portModeSolver_(portModeSolver), options_(options) {}

apm::RectPortGeometry TransfiniteElementBuilder::detectRectangular(int faceId) const {
    apm::AnalyticPortBuilder helper(mesh_, topology_);
    return helper.detectRectangular(faceId);
}

MultiPortMode TransfiniteElementBuilder::build(int faceId) const {
    // Numerical multi-mode port eigensolve. The port mesh itself decides
    // what modes are supported; no rectangular-cross-section assumption is
    // made. The dominant (lowest k_c^2) mode is kept at index 0 so the
    // existing single-mode S-parameter extractor stays consistent.
    const int requested = std::max(1, options_.modesPerPort);
    MultiPortMode multi = portModeSolver_.computeMultiMode(faceId, requested);
    if (multi.modes.empty()) {
        throw std::runtime_error("TransfiniteElementBuilder: numerical port eigensolve produced no modes for face "
                                 + std::to_string(faceId));
    }
    multi.faceId = faceId;
    multi.excitationModeIndex = 0;
    return multi;
}

}  // namespace fem::tfe

