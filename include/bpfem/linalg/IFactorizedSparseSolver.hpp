#pragma once

#include "bpfem/linalg/ISparseSolver.hpp"

#include <cstdint>
#include <vector>

namespace fem::linalg {

// 累计记录直接求解器的符号分析、数值分解和回代开销。
struct FactorizationStatistics {
    std::uint64_t symbolicAnalysisCount = 0;
    std::uint64_t numericFactorizationCount = 0;
    // 实际求解的 RHS 向量数；批量回代时一次调用可包含多个 RHS。
    std::uint64_t factorizedRhsSolveCount = 0;
    // phase=33 或等价回代入口的调用次数。
    std::uint64_t factorizedSolveCallCount = 0;
    std::uint64_t batchRhsMax = 0;
    double symbolicAnalysisSec = 0.0;
    double numericFactorizationSec = 0.0;
    double factorizedRhsSolveSec = 0.0;
};

// 可复用数值分解的直接求解器能力接口。
//
// 调用方先对一个矩阵执行 factorize()，随后可针对多个右端项调用
// solveFactorized()。新的 factorize() 会替换已有数值因子。
class IFactorizedSparseSolver {
public:
    virtual ~IFactorizedSparseSolver() = default;

    virtual void factorize(const SparseMatrix& matrix,
                           const SolverConfig& config = {}) = 0;

    virtual SolveResult solveFactorized(
        const std::vector<std::complex<double>>& rhs,
        const SolverConfig& config = {}) = 0;

    // 批量求解共享同一数值因子的 RHS。默认实现保证旧后端仍可工作，
    // PARDISO 后端会覆盖该函数并用一次 phase=33 完成整批回代。
    virtual std::vector<SolveResult> solveFactorizedBatch(
        const std::vector<std::vector<std::complex<double>>>& rhsList,
        const SolverConfig& config = {}) {
        std::vector<SolveResult> results;
        results.reserve(rhsList.size());
        for (const auto& rhs : rhsList) {
            results.push_back(solveFactorized(rhs, config));
        }
        return results;
    }

    virtual void clearFactorization() = 0;

    virtual FactorizationStatistics factorizationStatistics() const = 0;
};

}  // namespace fem::linalg
