#pragma once

#include "bpfem/linalg/IFactorizedSparseSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/linalg/SolveResidual.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace fem::linalg {

// 在一个固定矩阵上统一管理“分解一次、多 RHS 求解”和迭代后端回退路径。
// 同一个后端同一时间仅允许一个会话，并发构造时只有一个会话能占用后端。
// 会话存活期间按顺序调用 solve / solveBatch，且不得直接操作或移动后端。
class FactorizedSolveSession {
public:
    // 若后端支持可复用分解，构造时立即完成该矩阵的数值分解。
    FactorizedSolveSession(ISparseSolver& solver,
                           const SparseMatrix& matrix,
                           const SolverConfig& config = {},
                           bool computeResidual = true)
        : solver_(solver), matrix_(matrix), config_(config),
          computeResidual_(computeResidual),
          factorized_(dynamic_cast<IFactorizedSparseSolver*>(&solver)) {
        if (solver_.fixedMatrixSessionActive_.test_and_set(std::memory_order_acquire)) {
            throw std::logic_error("Sparse solver already has a fixed-matrix session");
        }
        try {
            if (factorized_ != nullptr) {
                factorized_->factorize(matrix_, config_);
            } else {
                solver_.beginMatrix(matrix_);
            }
        } catch (...) {
            solver_.fixedMatrixSessionActive_.clear(std::memory_order_release);
            throw;
        }
    }

    FactorizedSolveSession(const FactorizedSolveSession&) = delete;
    FactorizedSolveSession& operator=(const FactorizedSolveSession&) = delete;

    ~FactorizedSolveSession() {
        if (factorized_ == nullptr) {
            solver_.endMatrix();
        }
        solver_.fixedMatrixSessionActive_.clear(std::memory_order_release);
    }

    // 求解一个 RHS；直接后端使用已有数值因子，迭代后端保持旧行为。
    SolveResult solve(const std::vector<std::complex<double>>& rhs) {
        if (factorized_ == nullptr) {
            SolveResult result = solver_.solve(matrix_, rhs, config_);
            updateResidual(rhs, result);
            return result;
        }
        SolveResult result = factorized_->solveFactorized(rhs, config_);
        if (computeResidual_) {
            updateResidual(rhs, result);
        }
        return result;
    }

    // 对同一矩阵批量求解多个 RHS，并统一补算全阶相对残差。
    std::vector<SolveResult> solveBatch(
        const std::vector<std::vector<std::complex<double>>>& rhsList) {
        if (rhsList.empty()) {
            return {};
        }
        std::vector<SolveResult> results;
        if (factorized_ == nullptr) {
            results.reserve(rhsList.size());
            for (const auto& rhs : rhsList) {
                results.push_back(solver_.solve(matrix_, rhs, config_));
            }
        } else {
            results = factorized_->solveFactorizedBatch(rhsList, config_);
        }
        if (results.size() != rhsList.size()) {
            throw std::runtime_error(
                "FactorizedSolveSession: batch solver returned the wrong result count");
        }
        for (std::size_t i = 0; i < rhsList.size(); ++i) {
            if (computeResidual_) {
                updateResidual(rhsList[i], results[i]);
            }
        }
        return results;
    }

    // 指示当前会话是否使用一次分解、多次回代的能力接口。
    bool reusesFactorization() const { return factorized_ != nullptr; }

private:
    // 由原始稀疏矩阵计算相对残差，避免依赖求解器内部工作数组。
    void updateResidual(const std::vector<std::complex<double>>& rhs,
                        SolveResult& result) const {
        if (result.status == SolveStatus::Unchecked) {
            checkResidual(matrix_, rhs, config_.tolerance, result);
        }
    }

    ISparseSolver& solver_;
    const SparseMatrix& matrix_;
    SolverConfig config_;
    bool computeResidual_ = true;
    IFactorizedSparseSolver* factorized_ = nullptr;
};

}  // namespace fem::linalg