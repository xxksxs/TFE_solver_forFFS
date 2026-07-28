#include "bpfem/fastsweep/UpperTriangularBlockAction.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace fem::fastsweep {

namespace {

// 将固定步长上三角矩阵的二维索引转换为一维位置。
std::size_t matrixIndex(int row, int col, int stride) {
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(stride)
        + static_cast<std::size_t>(col);
}

// 校验固定步长矩阵及当前已生成阶数，防止递推读取尚未形成的 U 区域。
void validateStorage(const std::vector<Complex>& upper,
                     int stride,
                     int availableOrder,
                     double singularTolerance) {
    if (stride < 1 || availableOrder < 1 || availableOrder > stride) {
        throw std::invalid_argument(
            "UpperTriangularBlockAction: invalid stride or available order");
    }
    const std::size_t required =
        static_cast<std::size_t>(stride) * static_cast<std::size_t>(stride);
    if (upper.size() != required) {
        throw std::invalid_argument(
            "UpperTriangularBlockAction: U storage must have stride squared entries");
    }
    if (!std::isfinite(singularTolerance) || singularTolerance <= 0.0) {
        throw std::invalid_argument(
            "UpperTriangularBlockAction: singular tolerance must be finite and positive");
    }
}

}  // namespace

// 对 U 的连续主子块执行上三角回代，不显式形成任何逆矩阵。
std::vector<Complex> UpperTriangularBlockAction::solveSubblock(
    const std::vector<Complex>& upper,
    int stride,
    int availableOrder,
    int start,
    const std::vector<Complex>& rhs,
    double singularTolerance) {
    validateStorage(upper, stride, availableOrder, singularTolerance);
    const int size = static_cast<int>(rhs.size());
    if (size < 1 || start < 0 || start + size > availableOrder) {
        throw std::invalid_argument(
            "UpperTriangularBlockAction: subblock is outside the available U region");
    }

    std::vector<Complex> solution(rhs);
    double blockScale = 0.0;
    for (int row = 0; row < size; ++row) {
        for (int col = row; col < size; ++col) {
            blockScale = std::max(
                blockScale,
                std::abs(upper[matrixIndex(start + row, start + col, stride)]));
        }
    }
    const double diagonalFloor = singularTolerance * std::max(blockScale, 1.0);
    for (int row = size - 1; row >= 0; --row) {
        Complex value = solution[static_cast<std::size_t>(row)];
        for (int col = row + 1; col < size; ++col) {
            value -= upper[matrixIndex(start + row, start + col, stride)]
                * solution[static_cast<std::size_t>(col)];
        }
        const Complex diagonal = upper[matrixIndex(start + row, start + row, stride)];
        if (!std::isfinite(diagonal.real()) || !std::isfinite(diagonal.imag())
            || std::abs(diagonal) <= diagonalFloor) {
            throw UpperTriangularBreakdown(
                "UpperTriangularBlockAction: near-singular U subblock");
        }
        solution[static_cast<std::size_t>(row)] = value / diagonal;
    }
    return solution;
}

// 按 B_w^-1...B_m^-1 的右作用顺序，从 t=m 到 t=w 连续执行三角回代。
std::vector<Complex> UpperTriangularBlockAction::applyToLastUnitVector(
    const std::vector<Complex>& upper,
    int stride,
    int availableOrder,
    int n,
    int m,
    int w,
    double singularTolerance,
    std::uint64_t* triangularSolveCount) {
    validateStorage(upper, stride, availableOrder, singularTolerance);
    if (n < 2 || w < 1 || m < w || m > n - 1 || availableOrder < n - 1) {
        throw std::invalid_argument(
            "UpperTriangularBlockAction: invalid one-based WCAWE indices");
    }

    const int blockSize = n - m;
    std::vector<Complex> result(
        static_cast<std::size_t>(blockSize), Complex(0.0, 0.0));
    result.back() = Complex(1.0, 0.0);
    for (int t = m; t >= w; --t) {
        result = solveSubblock(
            upper, stride, availableOrder, t - 1, result, singularTolerance);
        if (triangularSolveCount != nullptr) {
            ++(*triangularSolveCount);
        }
    }
    return result;
}

}  // namespace fem::fastsweep
