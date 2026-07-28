#include "bpfem/fastsweep/PolynomialMomentRecurrence.hpp"

#include <stdexcept>
#include <string>

namespace fem::fastsweep {

// 对一般多项式矩阵方程执行矩递推：每一阶矩都通过展开点矩阵求解得到。
std::vector<std::vector<PolynomialMomentRecurrence::Complex>>
PolynomialMomentRecurrence::generatePolynomialMoments(
    const std::vector<std::vector<Complex>>& rhsCoefficients,
    const std::vector<LinearOperator>& matrixCoefficientOperators,
    const SolveFunction& solveAtExpansion) {
    if (rhsCoefficients.empty()) {
        throw std::runtime_error(
            "PolynomialMomentRecurrence: at least one RHS coefficient is required");
    }
    return generatePolynomialMoments(
        rhsCoefficients.size(),
        [&rhsCoefficients](std::size_t order) { return rhsCoefficients[order]; },
        matrixCoefficientOperators,
        solveAtExpansion);
}

// 按阶生成 RHS，并对一般多项式矩阵方程执行矩递推。
std::vector<std::vector<PolynomialMomentRecurrence::Complex>>
PolynomialMomentRecurrence::generatePolynomialMoments(
    std::size_t momentCount,
    const RhsCoefficientFunction& rhsCoefficientAt,
    const std::vector<LinearOperator>& matrixCoefficientOperators,
    const SolveFunction& solveAtExpansion) {
    if (momentCount == 0) {
        throw std::runtime_error(
            "PolynomialMomentRecurrence: at least one moment is required");
    }

    std::vector<std::vector<Complex>> moments;
    moments.reserve(momentCount);
    std::size_t dimension = 0;
    for (std::size_t j = 0; j < momentCount; ++j) {
        std::vector<Complex> rhs = rhsCoefficientAt(j);
        if (j == 0) {
            dimension = rhs.size();
            if (dimension == 0) {
                throw std::runtime_error(
                    "PolynomialMomentRecurrence: RHS coefficient vectors must be non-empty");
            }
        } else if (rhs.size() != dimension) {
            throw std::runtime_error(
                "PolynomialMomentRecurrence: inconsistent RHS coefficient size");
        }

        const std::size_t maxCoeff = std::min(j, matrixCoefficientOperators.size());
        for (std::size_t r = 1; r <= maxCoeff; ++r) {
            const std::vector<Complex> applied =
                matrixCoefficientOperators[r - 1](moments[j - r]);
            if (applied.size() != dimension) {
                throw std::runtime_error(
                    "PolynomialMomentRecurrence: matrix coefficient returned the wrong size");
            }
            for (std::size_t i = 0; i < dimension; ++i) {
                rhs[i] -= applied[i];
            }
        }

        std::vector<Complex> x = solveAtExpansion(rhs);
        if (x.size() != dimension) {
            throw std::runtime_error(
                "PolynomialMomentRecurrence: solver returned a vector with the wrong size");
        }
        moments.push_back(std::move(x));
    }
    return moments;
}

// 针对 A(lambda)=A0-(lambda-lambda0)M 的一阶质量矩阵模型生成局部矩。
std::vector<std::vector<PolynomialMomentRecurrence::Complex>>
PolynomialMomentRecurrence::generateFirstOrderMassMoments(
    const SparseMatrix& massMatrix,
    const std::vector<Complex>& rhs0,
    int momentCount,
    const SolveFunction& solveAtExpansion) {
    if (momentCount < 1) {
        throw std::runtime_error("PolynomialMomentRecurrence: momentCount must be >= 1");
    }
    if (rhs0.size() != massMatrix.size()) {
        throw std::runtime_error("PolynomialMomentRecurrence: rhs size does not match mass matrix");
    }

    std::vector<std::vector<Complex>> moments;
    moments.reserve(static_cast<std::size_t>(momentCount));
    moments.push_back(solveAtExpansion(rhs0));
    if (moments.front().size() != rhs0.size()) {
        throw std::runtime_error("PolynomialMomentRecurrence: solver returned a vector with the wrong size");
    }

    for (int j = 1; j < momentCount; ++j) {
        const std::vector<Complex> rhs = massMatrix.multiply(moments.back());
        std::vector<Complex> x = solveAtExpansion(rhs);
        if (x.size() != rhs0.size()) {
            throw std::runtime_error("PolynomialMomentRecurrence: solver returned a vector with the wrong size");
        }
        moments.push_back(std::move(x));
    }
    return moments;
}

}  // namespace fem::fastsweep
