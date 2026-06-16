#include "bpfem/fastsweep/WellConditionedBasisBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace fem::fastsweep {

WellConditionedBasisBuilder::WellConditionedBasisBuilder(double dropTolerance)
    : dropTolerance_(dropTolerance > 0.0 ? dropTolerance : 1.0e-12) {}

void WellConditionedBasisBuilder::clear() {
    deflatedColumns_ = 0;
    maxMomentReconstructionError_ = 0.0;
    basis_.clear();
    triangularR_.clear();
    rDiagonalAbs_.clear();
    records_.clear();
}

bool WellConditionedBasisBuilder::append(const std::vector<Complex>& moment) {
    std::vector<Complex> w = moment;
    const double pre = norm2(w);
    if (pre < dropTolerance_) {
        ++deflatedColumns_;
        return false;
    }

    const int oldDim = static_cast<int>(basis_.size());
    std::vector<Complex> coeffs(static_cast<std::size_t>(oldDim), Complex(0.0, 0.0));
    for (int pass = 0; pass < 2; ++pass) {
        for (int j = 0; j < oldDim; ++j) {
            const Complex h = hdot(basis_[static_cast<std::size_t>(j)], w);
            coeffs[static_cast<std::size_t>(j)] += h;
            for (std::size_t i = 0; i < w.size(); ++i) {
                w[i] -= h * basis_[static_cast<std::size_t>(j)][i];
            }
        }
    }

    const double post = norm2(w);
    if (post < dropTolerance_ * std::max(pre, 1.0)) {
        ++deflatedColumns_;
        return false;
    }

    for (auto& z : w) {
        z /= post;
    }
    basis_.push_back(std::move(w));
    rDiagonalAbs_.push_back(post);

    const int newDim = oldDim + 1;
    std::vector<Complex> newR(static_cast<std::size_t>(newDim) * static_cast<std::size_t>(newDim),
                              Complex(0.0, 0.0));
    for (int r = 0; r < oldDim; ++r) {
        for (int c = 0; c < oldDim; ++c) {
            newR[static_cast<std::size_t>(r) * static_cast<std::size_t>(newDim)
                 + static_cast<std::size_t>(c)] =
                triangularR_[static_cast<std::size_t>(r) * static_cast<std::size_t>(oldDim)
                             + static_cast<std::size_t>(c)];
        }
    }
    for (int r = 0; r < oldDim; ++r) {
        newR[static_cast<std::size_t>(r) * static_cast<std::size_t>(newDim)
             + static_cast<std::size_t>(oldDim)] = coeffs[static_cast<std::size_t>(r)];
    }
    newR[static_cast<std::size_t>(oldDim) * static_cast<std::size_t>(newDim)
         + static_cast<std::size_t>(oldDim)] = Complex(post, 0.0);
    triangularR_ = std::move(newR);

    std::vector<Complex> reconstructed(moment.size(), Complex(0.0, 0.0));
    for (int j = 0; j < newDim; ++j) {
        const Complex coeff = triangularR_[static_cast<std::size_t>(j)
                                           * static_cast<std::size_t>(newDim)
                                           + static_cast<std::size_t>(oldDim)];
        const auto& v = basis_[static_cast<std::size_t>(j)];
        for (std::size_t i = 0; i < reconstructed.size(); ++i) {
            reconstructed[i] += coeff * v[i];
        }
    }
    for (std::size_t i = 0; i < reconstructed.size(); ++i) {
        reconstructed[i] -= moment[i];
    }
    const double reconstructionError = norm2(reconstructed) / std::max(pre, 1.0e-300);
    maxMomentReconstructionError_ = std::max(maxMomentReconstructionError_, reconstructionError);

    double minDiag = rDiagonalAbs_.front();
    double maxDiag = rDiagonalAbs_.front();
    for (double d : rDiagonalAbs_) {
        minDiag = std::min(minDiag, d);
        maxDiag = std::max(maxDiag, d);
    }

    BasisConditionRecord rec;
    rec.order = newDim;
    rec.aweConditionProxy = maxDiag / std::max(minDiag, 1.0e-300);
    rec.wcaweConditionProxy = 1.0 + orthogonalityError(basis_);
    rec.rDiagonalAbs = post;
    rec.orthogonalityError = rec.wcaweConditionProxy - 1.0;
    rec.momentReconstructionError = reconstructionError;
    records_.push_back(rec);
    return true;
}

}  // namespace fem::fastsweep
