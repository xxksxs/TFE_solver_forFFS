#pragma once

#include "bpfem/tfe/TfeSystem.hpp"

#include <array>

namespace fem::fastsweep {

class TfePolynomialModel {
public:
    using Vector = tfe::TfeSystem::Vector;
    using Weights = std::array<double, 3>;

    TfePolynomialModel(const tfe::TfeSystem& system, double leftHz, double rightHz);

    double coordinate(double frequencyHz) const;
    const Weights& frequencies() const { return frequencies_; }
    const std::array<SparseMatrix, 3>& samples() const { return samples_; }
    Weights weights(double t, unsigned order = 0) const;
    SparseMatrix matrix(const Weights& weights) const;
    Vector apply(const Weights& weights, const Vector& state) const;

private:
    double midpoint_;
    double halfWidth_;
    Weights nodes_;
    Weights frequencies_;
    std::array<SparseMatrix, 3> samples_;
    std::shared_ptr<const SparsePattern> pattern_;
};

}  // namespace fem::fastsweep
