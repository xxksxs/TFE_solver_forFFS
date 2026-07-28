#pragma once

#include "bpfem/linalg/SparseMatrix.hpp"

#include <complex>
#include <functional>
#include <vector>

namespace fem::fastsweep {

class PolynomialMomentRecurrence {
public:
    using Complex = std::complex<double>;
    using SolveFunction = std::function<std::vector<Complex>(const std::vector<Complex>& rhs)>;
    using LinearOperator = std::function<std::vector<Complex>(const std::vector<Complex>& x)>;
    using RhsCoefficientFunction = std::function<std::vector<Complex>(std::size_t order)>;

    // 为一般多项式矩阵方程生成矩向量：
    //
    //   A(t) = A0 + t A1 + t^2 A2 + ...
    //   b(t) = b0 + t b1 + t^2 b2 + ...
    //
    // 系数使用归一化展开变量 t；matrixCoefficientOperators[r - 1] 表示 Ar 的作用。
    static std::vector<std::vector<Complex>> generatePolynomialMoments(
        const std::vector<std::vector<Complex>>& rhsCoefficients,
        const std::vector<LinearOperator>& matrixCoefficientOperators,
        const SolveFunction& solveAtExpansion);

    // 按阶即时生成 RHS 系数，避免同时保存 order 个全尺寸复向量。
    static std::vector<std::vector<Complex>> generatePolynomialMoments(
        std::size_t momentCount,
        const RhsCoefficientFunction& rhsCoefficientAt,
        const std::vector<LinearOperator>& matrixCoefficientOperators,
        const SolveFunction& solveAtExpansion);

    // 为一阶质量矩阵模型生成矩向量：
    //
    //   A(lambda) = A0 - (lambda - lambda0) M, b(lambda) = b0.
    //
    // 递推关系为 x0 = A0^-1 b0，xj = A0^-1 M x{j-1}。
    static std::vector<std::vector<Complex>> generateFirstOrderMassMoments(
        const SparseMatrix& massMatrix,
        const std::vector<Complex>& rhs0,
        int momentCount,
        const SolveFunction& solveAtExpansion);
};

}  // namespace fem::fastsweep
