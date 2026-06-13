#include "bpfem/fastsweep/PadeApproximant.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace fem::fastsweep {

namespace {

using Complex = PadeApproximant::Complex;

std::size_t idx(int row, int col, int n) {
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(n)
        + static_cast<std::size_t>(col);
}

bool denseSolve(std::vector<Complex>& A, std::vector<Complex>& rhs, int n, double& pivotRatio) {
    double minPivot = std::numeric_limits<double>::infinity();
    double maxPivot = 0.0;

    for (int k = 0; k < n; ++k) {
        int pivotRow = k;
        double pivotMag = std::abs(A[idx(k, k, n)]);
        for (int r = k + 1; r < n; ++r) {
            const double mag = std::abs(A[idx(r, k, n)]);
            if (mag > pivotMag) {
                pivotMag = mag;
                pivotRow = r;
            }
        }
        if (pivotMag < 1.0e-30) {
            return false;
        }
        minPivot = std::min(minPivot, pivotMag);
        maxPivot = std::max(maxPivot, pivotMag);
        if (pivotRow != k) {
            for (int c = 0; c < n; ++c) {
                std::swap(A[idx(k, c, n)], A[idx(pivotRow, c, n)]);
            }
            std::swap(rhs[static_cast<std::size_t>(k)], rhs[static_cast<std::size_t>(pivotRow)]);
        }

        const Complex pivot = A[idx(k, k, n)];
        for (int r = k + 1; r < n; ++r) {
            const Complex factor = A[idx(r, k, n)] / pivot;
            A[idx(r, k, n)] = factor;
            for (int c = k + 1; c < n; ++c) {
                A[idx(r, c, n)] -= factor * A[idx(k, c, n)];
            }
            rhs[static_cast<std::size_t>(r)] -= factor * rhs[static_cast<std::size_t>(k)];
        }
    }

    for (int r = n - 1; r >= 0; --r) {
        Complex sum = rhs[static_cast<std::size_t>(r)];
        for (int c = r + 1; c < n; ++c) {
            sum -= A[idx(r, c, n)] * rhs[static_cast<std::size_t>(c)];
        }
        rhs[static_cast<std::size_t>(r)] = sum / A[idx(r, r, n)];
    }

    pivotRatio = maxPivot / std::max(minPivot, 1.0e-300);
    return true;
}

Complex horner(const std::vector<Complex>& coeffs, Complex t) {
    Complex out(0.0, 0.0);
    for (auto it = coeffs.rbegin(); it != coeffs.rend(); ++it) {
        out = out * t + *it;
    }
    return out;
}

}  // namespace

PadeApproximant PadeApproximant::build(const std::vector<Complex>& moments,
                                       int numeratorOrder,
                                       int denominatorOrder) {
    if (numeratorOrder < 0 || denominatorOrder < 0) {
        throw std::runtime_error("PadeApproximant: orders must be non-negative");
    }
    const int required = numeratorOrder + denominatorOrder + 1;
    if (static_cast<int>(moments.size()) < required) {
        throw std::runtime_error("PadeApproximant: not enough moments for requested order");
    }

    PadeApproximant out;
    out.denominator_.assign(static_cast<std::size_t>(denominatorOrder + 1), Complex(0.0, 0.0));
    out.denominator_[0] = Complex(1.0, 0.0);

    if (denominatorOrder > 0) {
        std::vector<Complex> system(static_cast<std::size_t>(denominatorOrder)
                                    * static_cast<std::size_t>(denominatorOrder),
                                    Complex(0.0, 0.0));
        std::vector<Complex> rhs(static_cast<std::size_t>(denominatorOrder), Complex(0.0, 0.0));
        for (int row = 0; row < denominatorOrder; ++row) {
            const int k = numeratorOrder + 1 + row;
            rhs[static_cast<std::size_t>(row)] = -moments[static_cast<std::size_t>(k)];
            for (int col = 0; col < denominatorOrder; ++col) {
                const int j = col + 1;
                system[idx(row, col, denominatorOrder)] = moments[static_cast<std::size_t>(k - j)];
            }
        }

        double pivotRatio = 1.0;
        if (!denseSolve(system, rhs, denominatorOrder, pivotRatio)) {
            throw std::runtime_error("PadeApproximant: singular moment system");
        }
        out.pivotRatio_ = pivotRatio;
        for (int j = 1; j <= denominatorOrder; ++j) {
            out.denominator_[static_cast<std::size_t>(j)] = rhs[static_cast<std::size_t>(j - 1)];
        }
    }

    out.numerator_.assign(static_cast<std::size_t>(numeratorOrder + 1), Complex(0.0, 0.0));
    for (int i = 0; i <= numeratorOrder; ++i) {
        Complex coeff(0.0, 0.0);
        const int maxJ = std::min(i, denominatorOrder);
        for (int j = 0; j <= maxJ; ++j) {
            coeff += out.denominator_[static_cast<std::size_t>(j)]
                * moments[static_cast<std::size_t>(i - j)];
        }
        out.numerator_[static_cast<std::size_t>(i)] = coeff;
    }
    return out;
}

PadeApproximant PadeApproximant::buildBestEffort(const std::vector<Complex>& moments,
                                                 int requestedNumeratorOrder,
                                                 int requestedDenominatorOrder) {
    for (int den = requestedDenominatorOrder; den >= 0; --den) {
        const int maxNum = static_cast<int>(moments.size()) - den - 1;
        const int num = std::min(requestedNumeratorOrder, maxNum);
        if (num < 0) {
            continue;
        }
        try {
            return build(moments, num, den);
        } catch (const std::runtime_error&) {
            if (den == 0) {
                throw;
            }
        }
    }
    throw std::runtime_error("PadeApproximant: no stable Padé order could be built");
}

PadeApproximant::Complex PadeApproximant::evaluate(Complex t) const {
    const Complex den = horner(denominator_, t);
    if (std::abs(den) < 1.0e-30) {
        throw std::runtime_error("PadeApproximant: denominator is near zero");
    }
    return horner(numerator_, t) / den;
}

}  // namespace fem::fastsweep
