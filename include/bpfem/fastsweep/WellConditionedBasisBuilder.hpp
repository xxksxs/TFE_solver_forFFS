#pragma once

#include "bpfem/fastsweep/LinearAlgebra.hpp"

#include <utility>
#include <vector>

namespace fem::fastsweep {

struct WcaweAppendResult {
    bool accepted = false;
    int column = 0;
    double candidateNormBeforeMgs = 0.0;
    double diagonalAbs = 0.0;
    double basisRelationResidual = 0.0;
    double orthogonalityError = 0.0;
    int reorthogonalizationPasses = 0;
};

class WellConditionedBasisBuilder {
public:
    // 创建空构造器；开始一次递推前必须调用 reset 固定 U 的步长。
    WellConditionedBasisBuilder() = default;

    // 为 requestedOrder 阶递推分配固定步长 U，并设置 breakdown 与再正交参数。
    void reset(int requestedOrder,
               double breakdownTolerance = 1.0e-12,
               int reorthogonalizationPasses = 2);

    // 清空基和上三角 U，结束当前递推会话。
    void clear();

    // 正交化已经包含论文校正项的候选向量，并把 MGS 系数立即写入 U 当前列。
    WcaweAppendResult appendCandidate(const std::vector<Complex>& candidate);

    // 返回当前 WCAWE Hermitian 正交基。
    const std::vector<std::vector<Complex>>& basis() const { return basis_; }

    // 将正交基的所有权移交给调用方，避免同时保留两份全尺寸基。
    std::vector<std::vector<Complex>> takeBasis() { return std::move(basis_); }

    // 返回满足 V_tilde=V*U 的固定步长上三角矩阵。
    const std::vector<Complex>& upperTriangularU() const { return upperTriangularU_; }

    // 返回当前保留的基向量数量。
    int dimension() const { return static_cast<int>(basis_.size()); }

    // 返回 U 的固定行步长，即本次请求的目标阶数。
    int stride() const { return requestedOrder_; }

private:
    int requestedOrder_ = 0;
    int reorthogonalizationPasses_ = 2;
    double breakdownTolerance_ = 1.0e-12;
    std::vector<std::vector<Complex>> basis_;
    std::vector<Complex> upperTriangularU_;
};

}  // namespace fem::fastsweep
