#include "bpfem/core/Types.hpp"

namespace fem {

std::unordered_map<int, std::size_t> Mesh::pointIndex() const {
    std::unordered_map<int, std::size_t> index;
    index.reserve(pointIds.size());
    for (std::size_t i = 0; i < pointIds.size(); ++i) {
        index[pointIds[i]] = i;
    }
    return index;
}

}  // namespace fem
