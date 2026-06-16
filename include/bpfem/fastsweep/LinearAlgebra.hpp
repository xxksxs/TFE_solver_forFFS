#pragma once

#include <complex>
#include <utility>
#include <vector>

namespace fem::fastsweep {

using Complex = std::complex<double>;

std::vector<Complex> liftRealSparseVector(const std::vector<std::pair<int, double>>& sparse,
                                          std::size_t n);

Complex bilinear(const std::vector<Complex>& a, const std::vector<Complex>& b);
Complex hdot(const std::vector<Complex>& a, const std::vector<Complex>& b);
double norm2(const std::vector<Complex>& a);

double modifiedGramSchmidt(std::vector<Complex>& w,
                           const std::vector<std::vector<Complex>>& basis,
                           int reorthogonalizationPasses = 2);

double orthogonalityError(const std::vector<std::vector<Complex>>& basis);

bool denseSolve(std::vector<Complex>& a,
                std::vector<Complex>& rhs,
                int n,
                double* pivotRatio = nullptr);

}  // namespace fem::fastsweep
