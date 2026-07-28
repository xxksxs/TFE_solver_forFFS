#include "bpfem/fastsweep/WellConditionedBasisBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace fem::fastsweep {

namespace {

// 将固定步长 U 的二维索引转换为行主序位置。
std::size_t matrixIndex(int row, int col, int stride) {
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(stride)
        + static_cast<std::size_t>(col);
}

}  // namespace

// 初始化一次 WCAWE 基构造会话，并固定 U 的存储步长。
void WellConditionedBasisBuilder::reset(int requestedOrder,
                                        double breakdownTolerance,
                                        int reorthogonalizationPasses) {
    if (requestedOrder < 1) {
        throw std::invalid_argument(
            "WellConditionedBasisBuilder: requested order must be positive");
    }
    if (!std::isfinite(breakdownTolerance) || breakdownTolerance <= 0.0) {
        throw std::invalid_argument(
            "WellConditionedBasisBuilder: breakdown tolerance must be finite and positive");
    }
    if (reorthogonalizationPasses < 1 || reorthogonalizationPasses > 2) {
        throw std::invalid_argument(
            "WellConditionedBasisBuilder: MGS passes must be one or two");
    }
    requestedOrder_ = requestedOrder;
    breakdownTolerance_ = breakdownTolerance;
    reorthogonalizationPasses_ = reorthogonalizationPasses;
    basis_.clear();
    basis_.reserve(static_cast<std::size_t>(requestedOrder));
    upperTriangularU_.assign(
        static_cast<std::size_t>(requestedOrder) * static_cast<std::size_t>(requestedOrder),
        Complex(0.0, 0.0));
}

// 清空已有基和 U，防止下一次离线阶段读取旧递推状态。
void WellConditionedBasisBuilder::clear() {
    requestedOrder_ = 0;
    basis_.clear();
    upperTriangularU_.clear();
}

// 对论文递推候选做一至两遍 MGS，并保持 V_tilde=V*U 的当前列关系。
WcaweAppendResult WellConditionedBasisBuilder::appendCandidate(
    const std::vector<Complex>& candidate) {
    if (requestedOrder_ < 1
        || upperTriangularU_.size()
            != static_cast<std::size_t>(requestedOrder_)
                * static_cast<std::size_t>(requestedOrder_)) {
        throw std::logic_error(
            "WellConditionedBasisBuilder: reset must be called before appendCandidate");
    }
    if (basis_.size() >= static_cast<std::size_t>(requestedOrder_)) {
        throw std::logic_error(
            "WellConditionedBasisBuilder: requested order has already been reached");
    }
    if (candidate.empty()) {
        throw std::invalid_argument(
            "WellConditionedBasisBuilder: candidate must be non-empty");
    }
    if (!basis_.empty() && candidate.size() != basis_.front().size()) {
        throw std::invalid_argument(
            "WellConditionedBasisBuilder: candidate size does not match the basis");
    }

    WcaweAppendResult result;
    const int oldDim = static_cast<int>(basis_.size());
    result.column = oldDim + 1;
    result.reorthogonalizationPasses = reorthogonalizationPasses_;
    result.candidateNormBeforeMgs = norm2(candidate);
    if (!std::isfinite(result.candidateNormBeforeMgs)
        || result.candidateNormBeforeMgs <= breakdownTolerance_) {
        return result;
    }

    std::vector<Complex> work(candidate);
    std::vector<Complex> coefficients(
        static_cast<std::size_t>(oldDim), Complex(0.0, 0.0));
    for (int pass = 0; pass < reorthogonalizationPasses_; ++pass) {
        for (int column = 0; column < oldDim; ++column) {
            const Complex coefficient =
                hdot(basis_[static_cast<std::size_t>(column)], work);
            coefficients[static_cast<std::size_t>(column)] += coefficient;
            for (std::size_t row = 0; row < work.size(); ++row) {
                work[row] -= coefficient
                    * basis_[static_cast<std::size_t>(column)][row];
            }
        }
    }

    const double diagonal = norm2(work);
    result.diagonalAbs = diagonal;
    if (!std::isfinite(diagonal)
        || diagonal
            <= breakdownTolerance_ * std::max(result.candidateNormBeforeMgs, 1.0)) {
        return result;
    }

    for (Complex& value : work) {
        value /= diagonal;
    }
    basis_.push_back(std::move(work));
    for (int row = 0; row < oldDim; ++row) {
        upperTriangularU_[matrixIndex(row, oldDim, requestedOrder_)] =
            coefficients[static_cast<std::size_t>(row)];
    }
    upperTriangularU_[matrixIndex(oldDim, oldDim, requestedOrder_)] =
        Complex(diagonal, 0.0);

    std::vector<Complex> reconstructed(candidate.size(), Complex(0.0, 0.0));
    for (int column = 0; column <= oldDim; ++column) {
        const Complex coefficient =
            upperTriangularU_[matrixIndex(column, oldDim, requestedOrder_)];
        const auto& basisVector = basis_[static_cast<std::size_t>(column)];
        for (std::size_t row = 0; row < reconstructed.size(); ++row) {
            reconstructed[row] += coefficient * basisVector[row];
        }
    }
    for (std::size_t row = 0; row < reconstructed.size(); ++row) {
        reconstructed[row] -= candidate[row];
    }

    result.basisRelationResidual =
        norm2(reconstructed) / std::max(result.candidateNormBeforeMgs, 1.0e-300);
    result.orthogonalityError = orthogonalityError(basis_);
    result.accepted = true;
    return result;
}

}  // namespace fem::fastsweep
