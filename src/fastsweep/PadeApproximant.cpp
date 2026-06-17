#include "bpfem/fastsweep/PadeApproximant.hpp"

#include "bpfem/fastsweep/LinearAlgebra.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace fem::fastsweep {

namespace {

using Complex = PadeApproximant::Complex;

// 将 Padé 小系统的二维索引映射到一维行主序数组。
std::size_t idx(int row, int col, int n) {
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(n)
        + static_cast<std::size_t>(col);
}

// 用 Horner 形式稳定评估多项式系数序列。
Complex horner(const std::vector<Complex>& coeffs, Complex t) {
    Complex out(0.0, 0.0);
    for (auto it = coeffs.rbegin(); it != coeffs.rend(); ++it) {
        out = out * t + *it;
    }
    return out;
}

}  // namespace

// 由 Taylor 矩构造指定 [L/M] Padé 近似，并记录分母线性系统的主元比。
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
        if (!denseSolve(system, rhs, denominatorOrder, &pivotRatio)) {
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

// 在高阶 Padé 小系统奇异时逐步降低分母阶数，构造一个可用的近似。
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

// 在给定局部变量 t 上评估 Padé 有理函数。
PadeApproximant::Complex PadeApproximant::evaluate(Complex t) const {
    const Complex den = horner(denominator_, t);
    if (std::abs(den) < 1.0e-30) {
        throw std::runtime_error("PadeApproximant: denominator is near zero");
    }
    return horner(numerator_, t) / den;
}

}  // namespace fem::fastsweep
