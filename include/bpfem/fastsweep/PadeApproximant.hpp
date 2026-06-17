#pragma once

#include <complex>
#include <vector>

namespace fem::fastsweep {

class PadeApproximant {
public:
    using Complex = std::complex<double>;

    // 用 Taylor 矩构造指定阶数的 Padé 有理近似。
    static PadeApproximant build(const std::vector<Complex>& moments,
                                 int numeratorOrder,
                                 int denominatorOrder);

    // 在高阶 Padé 小系统不稳定时自动降阶，返回可构造的近似。
    static PadeApproximant buildBestEffort(const std::vector<Complex>& moments,
                                           int requestedNumeratorOrder,
                                           int requestedDenominatorOrder);

    // 在局部变量 t 上评估 Padé 有理函数。
    Complex evaluate(Complex t) const;

    // 返回分子多项式阶数。
    int numeratorOrder() const { return static_cast<int>(numerator_.size()) - 1; }

    // 返回分母多项式阶数。
    int denominatorOrder() const { return static_cast<int>(denominator_.size()) - 1; }

    // 返回构造 Padé 分母时小系统的主元比，用于诊断病态程度。
    double pivotRatio() const { return pivotRatio_; }

    // 返回分子多项式系数。
    const std::vector<Complex>& numerator() const { return numerator_; }

    // 返回分母多项式系数，首项按约定为 1。
    const std::vector<Complex>& denominator() const { return denominator_; }

private:
    std::vector<Complex> numerator_;
    std::vector<Complex> denominator_;
    double pivotRatio_ = 1.0;
};

}  // namespace fem::fastsweep
