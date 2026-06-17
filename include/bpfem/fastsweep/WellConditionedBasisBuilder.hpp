#pragma once

#include "bpfem/fastsweep/LinearAlgebra.hpp"

#include <vector>

namespace fem::fastsweep {

struct BasisConditionRecord {
    int order = 0;
    double aweConditionProxy = 0.0;
    double wcaweConditionProxy = 1.0;
    double rDiagonalAbs = 0.0;
    double orthogonalityError = 0.0;
    double momentReconstructionError = 0.0;
};

class WellConditionedBasisBuilder {
public:
    // 创建良条件基构造器，dropTolerance 控制近线性相关方向的丢弃。
    explicit WellConditionedBasisBuilder(double dropTolerance = 1.0e-12);

    // 添加一个 AWE 矩向量；成功时扩展正交基，失败时计入 deflation。
    bool append(const std::vector<Complex>& moment);

    // 清空基、上三角 R 和诊断记录。
    void clear();

    // 返回当前 WCAWE 正交基。
    const std::vector<std::vector<Complex>>& basis() const { return basis_; }

    // 返回满足 X≈V R 的上三角系数矩阵（行主序存储）。
    const std::vector<Complex>& triangularR() const { return triangularR_; }

    // 返回每一阶的条件数和重构误差诊断记录。
    const std::vector<BasisConditionRecord>& conditionRecords() const { return records_; }

    // 返回当前保留的基向量数量。
    int dimension() const { return static_cast<int>(basis_.size()); }

    // 返回因近线性相关而被丢弃的候选列数。
    int deflatedColumns() const { return deflatedColumns_; }

    // 返回 X≈V R 的最大相对重构误差。
    double maxMomentReconstructionError() const { return maxMomentReconstructionError_; }

private:
    double dropTolerance_ = 1.0e-12;
    int deflatedColumns_ = 0;
    double maxMomentReconstructionError_ = 0.0;
    std::vector<std::vector<Complex>> basis_;
    std::vector<Complex> triangularR_;
    std::vector<double> rDiagonalAbs_;
    std::vector<BasisConditionRecord> records_;
};

}  // namespace fem::fastsweep
