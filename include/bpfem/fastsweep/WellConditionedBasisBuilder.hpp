#pragma once

#include "bpfem/fastsweep/LinearAlgebra.hpp"

#include <vector>

namespace fem::fastsweep {

struct BasisConditionRecord {
    int order = 0;
    double aweConditionProxy = 0.0;
    double wcaweConditionProxy = 1.0;
    double rDiagonalAbs = 0.0;
    double orthogonalityError = 0.0;
    double momentReconstructionError = 0.0;
};

class WellConditionedBasisBuilder {
public:
    explicit WellConditionedBasisBuilder(double dropTolerance = 1.0e-12);

    bool append(const std::vector<Complex>& moment);
    void clear();

    const std::vector<std::vector<Complex>>& basis() const { return basis_; }
    const std::vector<Complex>& triangularR() const { return triangularR_; }
    const std::vector<BasisConditionRecord>& conditionRecords() const { return records_; }

    int dimension() const { return static_cast<int>(basis_.size()); }
    int deflatedColumns() const { return deflatedColumns_; }
    double maxMomentReconstructionError() const { return maxMomentReconstructionError_; }

private:
    double dropTolerance_ = 1.0e-12;
    int deflatedColumns_ = 0;
    double maxMomentReconstructionError_ = 0.0;
    std::vector<std::vector<Complex>> basis_;
    std::vector<Complex> triangularR_;
    std::vector<double> rDiagonalAbs_;
    std::vector<BasisConditionRecord> records_;
};

}  // namespace fem::fastsweep
