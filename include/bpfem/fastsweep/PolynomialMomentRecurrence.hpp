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

    // Generates moments for
    //
    //   A(t) = A0 + t A1 + t^2 A2 + ...
    //   b(t) = b0 + t b1 + t^2 b2 + ...
    //
    // Coefficients are supplied in the normalized expansion variable `t`.
    // matrixCoefficientOperators[r - 1] applies Ar for r >= 1.
    static std::vector<std::vector<Complex>> generatePolynomialMoments(
        const std::vector<std::vector<Complex>>& rhsCoefficients,
        const std::vector<LinearOperator>& matrixCoefficientOperators,
        const SolveFunction& solveAtExpansion);

    // Generates moments for the local first-order model
    //
    //   A(lambda) = A0 - (lambda - lambda0) M, b(lambda) = b0.
    //
    // The recurrence is x0 = A0^-1 b0 and xj = A0^-1 M x{j-1}.
    static std::vector<std::vector<Complex>> generateFirstOrderMassMoments(
        const SparseMatrix& massMatrix,
        const std::vector<Complex>& rhs0,
        int momentCount,
        const SolveFunction& solveAtExpansion);
};

}  // namespace fem::fastsweep
