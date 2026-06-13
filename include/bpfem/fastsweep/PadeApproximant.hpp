#pragma once

#include <complex>
#include <vector>

namespace fem::fastsweep {

class PadeApproximant {
public:
    using Complex = std::complex<double>;

    static PadeApproximant build(const std::vector<Complex>& moments,
                                 int numeratorOrder,
                                 int denominatorOrder);

    static PadeApproximant buildBestEffort(const std::vector<Complex>& moments,
                                           int requestedNumeratorOrder,
                                           int requestedDenominatorOrder);

    Complex evaluate(Complex t) const;

    int numeratorOrder() const { return static_cast<int>(numerator_.size()) - 1; }
    int denominatorOrder() const { return static_cast<int>(denominator_.size()) - 1; }
    double pivotRatio() const { return pivotRatio_; }

    const std::vector<Complex>& numerator() const { return numerator_; }
    const std::vector<Complex>& denominator() const { return denominator_; }

private:
    std::vector<Complex> numerator_;
    std::vector<Complex> denominator_;
    double pivotRatio_ = 1.0;
};

}  // namespace fem::fastsweep
