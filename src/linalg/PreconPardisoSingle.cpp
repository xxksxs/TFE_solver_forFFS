#include "bpfem/linalg/PreconPardisoSingle.hpp"

#include <stdexcept>

namespace fem::linalg {

PreconPardisoSingle::PreconPardisoSingle() {
#ifndef BPFEM_USE_MKL
    throw std::invalid_argument("--precon pardiso-single requires Intel oneMKL");
#endif
}

void PreconPardisoSingle::setup(const SparseMatrix& matrix) {
    solver_.factorize(matrix);
}

void PreconPardisoSingle::clear() noexcept {
    solver_.clearFactorization();
}

void PreconPardisoSingle::apply(const std::vector<std::complex<double>>& input,
                               std::vector<std::complex<double>>& output) const {
    output = solver_.solveFactorized(input).field;
}

}  // namespace fem::linalg
