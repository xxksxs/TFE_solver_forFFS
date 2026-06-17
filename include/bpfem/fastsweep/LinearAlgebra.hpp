#pragma once

#include <complex>
#include <utility>
#include <vector>

namespace fem::fastsweep {

using Complex = std::complex<double>;

// 将稀疏实向量提升为长度为 n 的全阶复向量。
std::vector<Complex> liftRealSparseVector(const std::vector<std::pair<int, double>>& sparse,
                                          std::size_t n);

// 复对称 Galerkin 投影使用的非共轭双线性内积。
Complex bilinear(const std::vector<Complex>& a, const std::vector<Complex>& b);

// 标准 Hermitian 内积，用于正交化和范数计算。
Complex hdot(const std::vector<Complex>& a, const std::vector<Complex>& b);

// 计算复向量的二范数。
double norm2(const std::vector<Complex>& a);

// 对候选向量做改进 Gram-Schmidt 正交化，并返回正交化后的范数。
double modifiedGramSchmidt(std::vector<Complex>& w,
                           const std::vector<std::vector<Complex>>& basis,
                           int reorthogonalizationPasses = 2);

// 计算一组基向量的最大正交性误差。
double orthogonalityError(const std::vector<std::vector<Complex>>& basis);

// 解小规模稠密复线性系统；a 和 rhs 会被原地覆盖为 LU/解。
bool denseSolve(std::vector<Complex>& a,
                std::vector<Complex>& rhs,
                int n,
                double* pivotRatio = nullptr);

}  // namespace fem::fastsweep
