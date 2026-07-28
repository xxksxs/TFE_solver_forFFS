#pragma once

#include "bpfem/fastsweep/LinearAlgebra.hpp"

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace fem::fastsweep {

// 表示式 (9) 所需 U 子块已经数值近奇异，调用方必须终止当前单点递推链。
class UpperTriangularBreakdown : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// 计算 Slone 2003 式 (9) 中连续上三角主子块逆乘积对单位向量的作用。
class UpperTriangularBlockAction {
public:
    // 返回 P_Uw(n,m)e_{n-m}。n、m、w 使用论文的一基索引，U 为固定 stride 的行主序数组。
    static std::vector<Complex> applyToLastUnitVector(
        const std::vector<Complex>& upper,
        int stride,
        int availableOrder,
        int n,
        int m,
        int w,
        double singularTolerance = 1.0e-14,
        std::uint64_t* triangularSolveCount = nullptr);

    // 求解 U[start:start+size-1,start:start+size-1] x=rhs，start 使用零基索引。
    static std::vector<Complex> solveSubblock(
        const std::vector<Complex>& upper,
        int stride,
        int availableOrder,
        int start,
        const std::vector<Complex>& rhs,
        double singularTolerance = 1.0e-14);
};

}  // namespace fem::fastsweep
