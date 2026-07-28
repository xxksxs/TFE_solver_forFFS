#include "bpfem/linalg/MklPardisoSolver.hpp"

#ifdef BPFEM_USE_MKL
#include <mkl.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem {
namespace {

// 使用原始稀疏矩阵计算相对残差，避免读取 PARDISO 的数值工作数组。
double relativeResidual(const SparseMatrix& matrix,
                        const std::vector<std::complex<double>>& x,
                        const std::vector<std::complex<double>>& b) {
    const auto ax = matrix.multiply(x);
    double r2 = 0.0;
    double b2 = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) {
        r2 += std::norm(ax[i] - b[i]);
        b2 += std::norm(b[i]);
    }
    return std::sqrt(r2) / std::max(std::sqrt(b2), 1.0e-300);
}

}  // namespace

// 初始化 PARDISO 控制参数；矩阵类型为复对称矩阵。
MklPardisoSolver::MklPardisoSolver() {
#ifdef BPFEM_USE_MKL
    initializeParameters();
#endif
}

// 释放 PARDISO 内部保存的符号分析和数值分解状态。
MklPardisoSolver::~MklPardisoSolver() {
#ifdef BPFEM_USE_MKL
    release();
#endif
}

#ifdef BPFEM_USE_MKL
// 设置 PARDISO 的零基 CSR、重排序和迭代改进参数。
void MklPardisoSolver::initializeParameters() {
    pt_.fill(nullptr);
    iparm_.fill(0);
    iparm_[0] = 1;
    iparm_[1] = 2;
    iparm_[7] = 2;
    iparm_[9] = 13;
    iparm_[10] = 1;
    iparm_[12] = 1;
    iparm_[17] = -1;
    iparm_[18] = -1;
    iparm_[26] = 1;
    iparm_[34] = 1;
}

// 判断新的矩阵是否可以复用当前符号分析结果。
bool MklPardisoSolver::hasSamePattern(
    MKL_INT n,
    const std::vector<int>& rowOffsets,
    const std::vector<int>& columnIndices) const {
    if (!analyzed_ || n != n_ || rowOffsets.size() != ia_.size()
        || columnIndices.size() != ja_.size()) {
        return false;
    }
    return std::equal(rowOffsets.begin(), rowOffsets.end(), ia_.begin())
        && std::equal(columnIndices.begin(), columnIndices.end(), ja_.begin());
}

// 执行 PARDISO phase=11，仅分析并缓存稀疏模式。
void MklPardisoSolver::analyzePattern(
    MKL_INT n,
    const std::vector<int>& rowOffsets,
    const std::vector<int>& columnIndices,
    std::vector<MKL_Complex16>& values) {
    n_ = n;
    ia_.assign(rowOffsets.begin(), rowOffsets.end());
    ja_.assign(columnIndices.begin(), columnIndices.end());

    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT mtype = 6;
    MKL_INT phase = 11;
    MKL_INT nrhs = 1;
    MKL_INT msglvl = 0;
    MKL_INT error = 0;
    MKL_INT idum = 0;
    MKL_Complex16 ddum{0.0, 0.0};

    const auto started = std::chrono::steady_clock::now();
    pardiso(pt_.data(), &maxfct, &mnum, &mtype, &phase, &n_, values.data(),
            ia_.data(), ja_.data(), &idum, &nrhs, iparm_.data(), &msglvl,
            &ddum, &ddum, &error);
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    statistics_.symbolicAnalysisSec += elapsed;
    if (error != 0) {
        release();
        throw std::runtime_error(
            "MKL PARDISO analysis failed with error " + std::to_string(error));
    }
    analyzed_ = true;
    ++statistics_.symbolicAnalysisCount;
}

// 释放 PARDISO 的全部内部内存，并清空当前分解标记。
void MklPardisoSolver::release() {
    if (analyzed_) {
        MKL_INT maxfct = 1;
        MKL_INT mnum = 1;
        MKL_INT mtype = 6;
        MKL_INT phase = -1;
        MKL_INT nrhs = 1;
        MKL_INT msglvl = 0;
        MKL_INT error = 0;
        MKL_INT idum = 0;
        MKL_Complex16 ddum{0.0, 0.0};
        pardiso(pt_.data(), &maxfct, &mnum, &mtype, &phase, &n_, &ddum,
                ia_.data(), ja_.data(), &idum, &nrhs, iparm_.data(), &msglvl,
                &ddum, &ddum, &error);
    }
    pt_.fill(nullptr);
    ia_.clear();
    ja_.clear();
    factorValues_.clear();
    n_ = 0;
    analyzed_ = false;
    factorized_ = false;
}
#endif

// 兼容旧的一次性求解接口：先分解矩阵，再对单个 RHS 回代。
SolveResult MklPardisoSolver::solve(
    const SparseMatrix& matrix,
    const std::vector<std::complex<double>>& rhs) {
#ifndef BPFEM_USE_MKL
    (void)matrix;
    (void)rhs;
    throw std::runtime_error(
        "MKL PARDISO solver was requested, but this binary was built without BPFEM_USE_MKL.");
#else
    factorize(matrix);
    SolveResult result = solveFactorized(rhs);
    result.residual = relativeResidual(matrix, result.field, rhs);
    return result;
#endif
}

// 缓存矩阵数值并执行 phase=22；相同稀疏模式不会重复 phase=11。
void MklPardisoSolver::factorize(const SparseMatrix& matrix) {
#ifndef BPFEM_USE_MKL
    (void)matrix;
    throw std::runtime_error(
        "MKL PARDISO factorization was requested, but this binary was built without BPFEM_USE_MKL.");
#else
    const MKL_INT n = static_cast<MKL_INT>(matrix.size());
    const auto& rowOffsets = matrix.rowOffsets();
    const auto& columnIndices = matrix.columnIndices();
    const auto& csrValues = matrix.values();

    const bool samePattern = hasSamePattern(n, rowOffsets, columnIndices);
    if (!samePattern && analyzed_) {
        release();
    }

    factorValues_.clear();
    factorValues_.reserve(csrValues.size());
    for (const auto& value : csrValues) {
        factorValues_.push_back(MKL_Complex16{value.real(), value.imag()});
    }

    if (!samePattern) {
        analyzePattern(n, rowOffsets, columnIndices, factorValues_);
    }

    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT mtype = 6;
    MKL_INT phase = 22;
    MKL_INT nrhs = 1;
    MKL_INT msglvl = 0;
    MKL_INT error = 0;
    MKL_INT idum = 0;
    MKL_Complex16 ddum{0.0, 0.0};

    factorized_ = false;
    iparm_[17] = -1;
    iparm_[18] = -1;
    const auto started = std::chrono::steady_clock::now();
    pardiso(pt_.data(), &maxfct, &mnum, &mtype, &phase, &n_,
            factorValues_.data(), ia_.data(), ja_.data(), &idum, &nrhs,
            iparm_.data(), &msglvl, &ddum, &ddum, &error);
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    statistics_.numericFactorizationSec += elapsed;
    if (error != 0) {
        release();
        throw std::runtime_error(
            "MKL PARDISO factorization failed with error " + std::to_string(error));
    }
    // phase=33 的迭代改进需要原始矩阵数值；phase=22 后重新装载输入值。
    for (std::size_t i = 0; i < csrValues.size(); ++i) {
        factorValues_[i] = MKL_Complex16{csrValues[i].real(), csrValues[i].imag()};
    }
    factorized_ = true;
    ++statistics_.numericFactorizationCount;
#endif
}

// 使用当前数值因子执行单 RHS phase=33；实现统一委托给批量入口。
SolveResult MklPardisoSolver::solveFactorized(
    const std::vector<std::complex<double>>& rhs) {
    auto results = solveFactorizedBatch({rhs});
    return std::move(results.front());
}

// 用一次 phase=33 回代多个列主序 RHS，减少 PARDISO 调用和线程调度开销。
std::vector<SolveResult> MklPardisoSolver::solveFactorizedBatch(
    const std::vector<std::vector<std::complex<double>>>& rhsList) {
#ifndef BPFEM_USE_MKL
    (void)rhsList;
    throw std::runtime_error(
        "MKL PARDISO batch solve was requested, but this binary was built without BPFEM_USE_MKL.");
#else
    if (!factorized_) {
        throw std::runtime_error(
            "MKL PARDISO solveFactorizedBatch called before a successful factorize");
    }
    if (rhsList.empty()) {
        return {};
    }
    for (const auto& rhs : rhsList) {
        if (rhs.size() != static_cast<std::size_t>(n_)) {
            throw std::runtime_error(
                "MKL PARDISO batch RHS size does not match the factorized matrix");
        }
    }

    const MKL_INT nrhs = static_cast<MKL_INT>(rhsList.size());
    const std::size_t valueCount = static_cast<std::size_t>(n_) * rhsList.size();
    rhsValues_.resize(valueCount);
    solutionValues_.resize(valueCount);
    auto& b = rhsValues_;
    auto& x = solutionValues_;
    for (std::size_t column = 0; column < rhsList.size(); ++column) {
        for (std::size_t row = 0; row < static_cast<std::size_t>(n_); ++row) {
            const auto value = rhsList[column][row];
            b[column * static_cast<std::size_t>(n_) + row] =
                MKL_Complex16{value.real(), value.imag()};
        }
    }

    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT mtype = 6;
    MKL_INT phase = 33;
    MKL_INT msglvl = 0;
    MKL_INT error = 0;
    MKL_INT idum = 0;
    MKL_INT mutableNrhs = nrhs;
    const MKL_INT savedRefinementSteps = iparm_[7];
    if (rhsList.size() > 1) {
        iparm_[7] = 0;
    }
    const auto started = std::chrono::steady_clock::now();
    pardiso(pt_.data(), &maxfct, &mnum, &mtype, &phase, &n_,
            factorValues_.data(), ia_.data(), ja_.data(), &idum, &mutableNrhs,
            iparm_.data(), &msglvl, b.data(), x.data(), &error);
    statistics_.factorizedRhsSolveSec += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    iparm_[7] = savedRefinementSteps;
    if (error != 0) {
        release();
        throw std::runtime_error(
            "MKL PARDISO batch solve failed with error " + std::to_string(error));
    }
    statistics_.factorizedRhsSolveCount += rhsList.size();
    ++statistics_.factorizedSolveCallCount;
    statistics_.batchRhsMax = std::max<std::uint64_t>(
        statistics_.batchRhsMax, static_cast<std::uint64_t>(rhsList.size()));

    std::vector<SolveResult> results(rhsList.size());
    for (std::size_t column = 0; column < rhsList.size(); ++column) {
        auto& field = results[column].field;
        field.resize(static_cast<std::size_t>(n_));
        results[column].iterations = 1;
        for (std::size_t row = 0; row < field.size(); ++row) {
            const auto value = x[column * static_cast<std::size_t>(n_) + row];
            field[row] = std::complex<double>(value.real, value.imag);
        }
    }
    return results;
#endif
}

// 显式清除当前符号分析和数值因子；累计统计保持不变。
void MklPardisoSolver::clearFactorization() {
#ifdef BPFEM_USE_MKL
    release();
#endif
}

}  // namespace fem
