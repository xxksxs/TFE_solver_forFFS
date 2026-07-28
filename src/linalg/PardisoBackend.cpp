#include "bpfem/linalg/PardisoBackend.hpp"

#include <stdexcept>

namespace fem::linalg {

// 创建持有单个 MKL PARDISO 状态对象的后端。
PardisoBackend::PardisoBackend() = default;

// 保留一次性求解语义，内部由 PARDISO 分解与回代两阶段完成。
SolveResult PardisoBackend::solve(const SparseMatrix& matrix,
                                  const std::vector<std::complex<double>>& rhs,
                                  const SolverConfig& /*config*/) {
#ifdef BPFEM_USE_MKL
    return solver_.solve(matrix, rhs);
#else
    (void)matrix;
    (void)rhs;
    throw std::runtime_error(
        "PardisoBackend: this binary was built without BPFEM_USE_MKL, "
        "but the sparse solver factory selected PARDISO. Rebuild with -DBPFEM_USE_MKL=ON "
        "or pass --linear-solver bicgstab.");
#endif
}

// 兼容原有符号模式复用开关；PARDISO 会自动按 CSR 模式判断是否复用。
void PardisoBackend::rememberPatternForReuse(bool /*enable*/) {
}

// 对一个展开点矩阵执行符号分析（必要时）和数值分解。
void PardisoBackend::factorize(const SparseMatrix& matrix,
                               const SolverConfig& /*config*/) {
#ifdef BPFEM_USE_MKL
    solver_.factorize(matrix);
#else
    (void)matrix;
    throw std::runtime_error(
        "PardisoBackend::factorize requires a build with BPFEM_USE_MKL.");
#endif
}

// 对当前数值因子执行单 RHS 回代。
SolveResult PardisoBackend::solveFactorized(
    const std::vector<std::complex<double>>& rhs,
    const SolverConfig& /*config*/) {
#ifdef BPFEM_USE_MKL
    return solver_.solveFactorized(rhs);
#else
    (void)rhs;
    throw std::runtime_error(
        "PardisoBackend::solveFactorized requires a build with BPFEM_USE_MKL.");
#endif
}

// 批量回代共享同一 PARDISO phase=33 调用。
std::vector<SolveResult> PardisoBackend::solveFactorizedBatch(
    const std::vector<std::vector<std::complex<double>>>& rhsList,
    const SolverConfig& /*config*/) {
#ifdef BPFEM_USE_MKL
    return solver_.solveFactorizedBatch(rhsList);
#else
    (void)rhsList;
    throw std::runtime_error(
        "PardisoBackend::solveFactorizedBatch requires a build with BPFEM_USE_MKL.");
#endif
}

// 清除当前 PARDISO 分解状态。
void PardisoBackend::clearFactorization() {
#ifdef BPFEM_USE_MKL
    solver_.clearFactorization();
#endif
}

// 返回本次运行累计的分解与回代统计。
FactorizationStatistics PardisoBackend::factorizationStatistics() const {
#ifdef BPFEM_USE_MKL
    return solver_.factorizationStatistics();
#else
    return {};
#endif
}

}  // namespace fem::linalg
