#include "bpfem/fastsweep/TfePolynomialModel.hpp"

#include <cmath>
#include <stdexcept>

namespace fem::fastsweep {

TfePolynomialModel::TfePolynomialModel(const tfe::TfeSystem& system, double leftHz, double rightHz)
    : midpoint_(leftHz + 0.5 * (rightHz - leftHz)), halfWidth_(0.5 * (rightHz - leftHz)),
      nodes_{-std::sqrt(3.0) / 2.0, 0.0, std::sqrt(3.0) / 2.0},
      frequencies_{midpoint_ + halfWidth_ * nodes_[0], midpoint_, midpoint_ + halfWidth_ * nodes_[2]},
      samples_{system.assemble(frequencies_[0]), system.assemble(frequencies_[1]),
               system.assemble(frequencies_[2])},
      pattern_(system.pattern()) {
    if (!(halfWidth_ > 0.0) || !(frequencies_[0] < frequencies_[1] && frequencies_[1] < frequencies_[2])) {
        throw std::invalid_argument("TFE Chebyshev interpolation needs three distinct sample frequencies");
    }
}

double TfePolynomialModel::coordinate(double frequencyHz) const {
    return (frequencyHz - midpoint_) / halfWidth_;
}

TfePolynomialModel::Weights TfePolynomialModel::weights(double t, unsigned order) const {
    Weights result{};
    for (std::size_t i = 0; i < 3; ++i) {
        const double a = nodes_[(i + 1) % 3];
        const double b = nodes_[(i + 2) % 3];
        const double denominator = (nodes_[i] - a) * (nodes_[i] - b);
        if (order == 0) result[i] = (t - a) * (t - b) / denominator;
        if (order == 1) result[i] = (2.0 * t - a - b) / denominator;
        if (order == 2) result[i] = 1.0 / denominator;
    }
    return result;
}

SparseMatrix TfePolynomialModel::matrix(const Weights& weights) const {
    SparseMatrix out(pattern_);
    for (std::size_t r = 0; r < pattern_->n; ++r) {
        for (int i = pattern_->rowOffsets[r]; i < pattern_->rowOffsets[r + 1]; ++i) {
            std::complex<double> value{};
            for (std::size_t j = 0; j < samples_.size(); ++j) value += weights[j] * samples_[j].values()[i];
            out.add(static_cast<int>(r), pattern_->columnIndices[i], value);
        }
    }
    return out;
}

TfePolynomialModel::Vector TfePolynomialModel::apply(const Weights& weights, const Vector& state) const {
    Vector result(pattern_->n);
    for (std::size_t j = 0; j < samples_.size(); ++j) {
        const auto applied = samples_[j].multiply(state);
        for (std::size_t i = 0; i < result.size(); ++i) result[i] += weights[j] * applied[i];
    }
    return result;
}

}  // namespace fem::fastsweep
