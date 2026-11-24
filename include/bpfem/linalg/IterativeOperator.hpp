#pragma once

#include "bpfem/linalg/SparseMatrix.hpp"

namespace fem::linalg {

class IterativeOperator {
public:
    IterativeOperator(const SparseMatrix& matrix, bool useMkl);
    ~IterativeOperator();

    std::vector<std::complex<double>> multiply(
        const std::vector<std::complex<double>>& input) const;

private:
    struct State;
    const SparseMatrix& matrix_;
    std::unique_ptr<State> state_;
};

}  // namespace fem::linalg
