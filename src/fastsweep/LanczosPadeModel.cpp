#include "bpfem/fastsweep/LanczosPadeModel.hpp"

#include "bpfem/fastsweep/LinearAlgebra.hpp"

#ifdef BPFEM_USE_MKL
#include <mkl_lapacke.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace fem::fastsweep {
namespace {

using Complex = LanczosPadeModel::Complex;
using Vector = LanczosPadeModel::Vector;

// 从目标向量中减去一个带复系数的基向量。
void subtractScaled(Vector& target, const Vector& basis, Complex scale) {
    for (std::size_t i = 0; i < target.size(); ++i) {
        target[i] -= scale * basis[i];
    }
}

// 计算向量差的二范数，供三项递推残差诊断使用。
double differenceNorm(const Vector& left, const Vector& right) {
    double squared = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        squared += std::norm(left[i] - right[i]);
    }
    return std::sqrt(squared);
}

// 计算候选向量对已有双正交基的最大相对投影。
double rightProjectionError(const Vector& candidate,
                            const std::vector<Vector>& leftBasis) {
    const double scale = std::max(norm2(candidate), 1.0e-300);
    double error = 0.0;
    for (const auto& left : leftBasis) {
        error = std::max(error, std::abs(bilinear(left, candidate)) / scale);
    }
    return error;
}

// 计算左候选向量对已有右基的最大相对投影。
double leftProjectionError(const Vector& candidate,
                           const std::vector<Vector>& rightBasis) {
    const double scale = std::max(norm2(candidate), 1.0e-300);
    double error = 0.0;
    for (const auto& right : rightBasis) {
        error = std::max(error, std::abs(bilinear(candidate, right)) / scale);
    }
    return error;
}

// 对右候选向量执行两遍选择性双正交化。
void reorthogonalizeRight(Vector& candidate,
                          const std::vector<Vector>& rightBasis,
                          const std::vector<Vector>& leftBasis) {
    for (int pass = 0; pass < 2; ++pass) {
        for (std::size_t i = 0; i < rightBasis.size(); ++i) {
            subtractScaled(candidate, rightBasis[i], bilinear(leftBasis[i], candidate));
        }
    }
}

// 对左候选向量执行两遍选择性双正交化。
void reorthogonalizeLeft(Vector& candidate,
                         const std::vector<Vector>& rightBasis,
                         const std::vector<Vector>& leftBasis) {
    for (int pass = 0; pass < 2; ++pass) {
        for (std::size_t i = 0; i < rightBasis.size(); ++i) {
            subtractScaled(candidate, leftBasis[i], bilinear(candidate, rightBasis[i]));
        }
    }
}

// 将算子像正交投影到当前双正交基的补空间。
void projectComplement(Vector& right,
                       Vector& left,
                       const std::vector<Vector>& rightBasis,
                       const std::vector<Vector>& leftBasis) {
    reorthogonalizeRight(right, rightBasis, leftBasis);
    reorthogonalizeLeft(left, rightBasis, leftBasis);
}

}  // namespace

// 用分离的左右算子构造标准双边 Lanczos 三项递推。
void LanczosPadeModel::build(const Vector& rightStart,
                             const Vector& leftOutput,
                             const Operator& applyRight,
                             const Operator& applyTranspose,
                             int order,
                             double breakdownTolerance,
                             bool retainReconstructionBasis,
                             bool retainLeftBasisForDiagnostics) {
    buildImpl(rightStart, leftOutput, applyRight, applyTranspose, {}, order,
              breakdownTolerance, retainReconstructionBasis,
              retainLeftBasisForDiagnostics);
}

// 用成对算子构造模型，使右、左移动可共享一次批量稀疏回代。
void LanczosPadeModel::buildPaired(const Vector& rightStart,
                                   const Vector& leftOutput,
                                   const Operator& applyRight,
                                   const PairedOperator& applyPair,
                                   int order,
                                   double breakdownTolerance,
                                   bool retainReconstructionBasis,
                                   bool retainLeftBasisForDiagnostics) {
    if (!applyPair) {
        throw std::runtime_error("LanczosPadeModel: paired operator is empty");
    }
    buildImpl(rightStart, leftOutput, applyRight, {}, applyPair, order,
              breakdownTolerance, retainReconstructionBasis,
              retainLeftBasisForDiagnostics);
}

// 生成 alpha/beta/gamma 三项系数，并在真正耦合击穿时尝试一次 2x2 look-ahead。
void LanczosPadeModel::buildImpl(const Vector& rightStart,
                                 const Vector& leftOutput,
                                 const Operator& applyRight,
                                 const Operator& applyTranspose,
                                 const PairedOperator& applyPair,
                                 int order,
                                 double breakdownTolerance,
                                 bool retainReconstructionBasis,
                                 bool retainLeftBasisForDiagnostics) {
    ready_ = false;
    dimension_ = 0;
    breakdownDetected_ = false;
    happyBreakdown_ = false;
    poleResidueReady_ = false;
    lastSolveSucceeded_ = true;
    biorthogonalityError_ = 0.0;
    tridiagonalLeakage_ = 0.0;
    finalCouplingMagnitude_ = 0.0;
    recurrenceResidual_ = 0.0;
    momentMatchingError_ = 0.0;
    poleResidueReconstructionError_ = 0.0;
    poleDecompositionSec_ = 0.0;
    lookAheadCount_ = 0;
    selectiveReorthogonalizationCount_ = 0;
    rightBasis_.clear();
    leftBasis_.clear();
    reducedOperator_.clear();
    reducedRightStart_.clear();
    reducedOutput_.clear();
    poleEigenvalues_.clear();
    poleWeights_.clear();

    if (rightStart.empty() || rightStart.size() != leftOutput.size()) {
        throw std::runtime_error("LanczosPadeModel: incompatible start vectors");
    }
    if (order < 1 || !applyRight) {
        throw std::runtime_error("LanczosPadeModel: invalid order or right operator");
    }
    if (!applyPair && !applyTranspose) {
        throw std::runtime_error("LanczosPadeModel: transpose operator is empty");
    }

    originalRightStart_ = rightStart;
    originalLeftOutput_ = leftOutput;
    const double rightNorm = norm2(rightStart);
    const double leftNorm = norm2(leftOutput);
    if (rightNorm <= std::numeric_limits<double>::min()
        || leftNorm <= std::numeric_limits<double>::min()) {
        throw std::runtime_error("LanczosPadeModel: start vector is numerically zero");
    }

    Vector firstRight = rightStart;
    for (auto& value : firstRight) {
        value /= rightNorm;
    }
    const Complex initialOverlap = bilinear(leftOutput, firstRight);
    if (std::abs(initialOverlap) <= breakdownTolerance * leftNorm) {
        throw std::runtime_error(
            "LanczosPadeModel: left/right start vectors have zero transfer moment");
    }
    Vector firstLeft = leftOutput;
    for (auto& value : firstLeft) {
        value /= initialOverlap;
    }

    std::vector<Vector> rightBasis;
    std::vector<Vector> leftBasis;
    rightBasis.reserve(static_cast<std::size_t>(order));
    leftBasis.reserve(static_cast<std::size_t>(order));
    rightBasis.push_back(std::move(firstRight));
    leftBasis.push_back(std::move(firstLeft));
    std::vector<Complex> alpha;
    std::vector<Complex> beta;
    std::vector<Complex> gamma;
    alpha.reserve(static_cast<std::size_t>(order));
    beta.reserve(static_cast<std::size_t>(std::max(0, order - 1)));
    gamma.reserve(static_cast<std::size_t>(std::max(0, order - 1)));

    // 让选择性再双正交化在双正交漂移进入结果精度前触发；该阈值
    // 仍允许绝大多数稳定步保持纯三项递推。
    const double reorthogonalizationTolerance = std::max(
        breakdownTolerance, 100.0 * std::numeric_limits<double>::epsilon());
    for (int j = 0; j < order; ++j) {
        Vector rightImage;
        Vector leftImage;
        if (j + 1 < order) {
            if (applyPair) {
                auto images = applyPair(rightBasis[static_cast<std::size_t>(j)],
                                        leftBasis[static_cast<std::size_t>(j)]);
                rightImage = std::move(images.first);
                leftImage = std::move(images.second);
            } else {
                rightImage = applyRight(rightBasis[static_cast<std::size_t>(j)]);
                leftImage = applyTranspose(leftBasis[static_cast<std::size_t>(j)]);
            }
        } else {
            rightImage = applyRight(rightBasis[static_cast<std::size_t>(j)]);
        }
        if (rightImage.size() != rightStart.size()
            || (j + 1 < order && leftImage.size() != rightStart.size())) {
            throw std::runtime_error("LanczosPadeModel: operator returned the wrong size");
        }

        const Complex diagonal = bilinear(
            leftBasis[static_cast<std::size_t>(j)], rightImage);
        alpha.push_back(diagonal);
        if (j + 1 == order) {
            break;
        }

        Vector rightCandidate = rightImage;
        Vector leftCandidate = leftImage;
        subtractScaled(rightCandidate, rightBasis[static_cast<std::size_t>(j)], diagonal);
        subtractScaled(leftCandidate, leftBasis[static_cast<std::size_t>(j)], diagonal);
        if (j > 0) {
            subtractScaled(rightCandidate, rightBasis[static_cast<std::size_t>(j - 1)],
                           gamma[static_cast<std::size_t>(j - 1)]);
            subtractScaled(leftCandidate, leftBasis[static_cast<std::size_t>(j - 1)],
                           beta[static_cast<std::size_t>(j - 1)]);
        }
        const Vector rawRightCandidate = rightCandidate;

        if (rightProjectionError(rightCandidate, leftBasis) > reorthogonalizationTolerance
            || leftProjectionError(leftCandidate, rightBasis) > reorthogonalizationTolerance) {
            projectComplement(rightCandidate, leftCandidate, rightBasis, leftBasis);
            ++selectiveReorthogonalizationCount_;
        }

        double rightCandidateNorm = norm2(rightCandidate);
        double leftCandidateNorm = norm2(leftCandidate);
        Complex coupling = bilinear(leftCandidate, rightCandidate);
        finalCouplingMagnitude_ = std::abs(coupling);
        double couplingScale = std::max(
            rightCandidateNorm * leftCandidateNorm, std::numeric_limits<double>::min());
        if (rightCandidateNorm <= breakdownTolerance * std::max(norm2(rightImage), 1.0)
            || leftCandidateNorm <= breakdownTolerance * std::max(norm2(leftImage), 1.0)) {
            happyBreakdown_ = true;
            breakdownDetected_ = true;
            break;
        }

        if (std::abs(coupling) <= breakdownTolerance * couplingScale) {
            // 最小 look-ahead：生成下一 Krylov 层，并从三个交叉配对中选择
            // 数值耦合最强的一对。若仍失败则显式终止，不制造虚假基向量。
            std::pair<Vector, Vector> lookImages = applyPair
                ? applyPair(rightCandidate, leftCandidate)
                : std::make_pair(applyRight(rightCandidate),
                                 applyTranspose(leftCandidate));
            Vector lookRight = std::move(lookImages.first);
            Vector lookLeft = std::move(lookImages.second);
            projectComplement(lookRight, lookLeft, rightBasis, leftBasis);
            ++lookAheadCount_;

            struct CandidatePair {
                Vector* right;
                Vector* left;
                double score;
            };
            auto score = [](const Vector& right, const Vector& left) {
                return std::abs(bilinear(left, right))
                    / std::max(norm2(right) * norm2(left), 1.0e-300);
            };
            CandidatePair candidates[] = {
                {&lookRight, &leftCandidate, score(lookRight, leftCandidate)},
                {&rightCandidate, &lookLeft, score(rightCandidate, lookLeft)},
                {&lookRight, &lookLeft, score(lookRight, lookLeft)}
            };
            const auto best = std::max_element(
                std::begin(candidates), std::end(candidates),
                [](const CandidatePair& a, const CandidatePair& b) {
                    return a.score < b.score;
                });
            if (best == std::end(candidates)
                || best->score <= breakdownTolerance) {
                breakdownDetected_ = true;
                break;
            }
            rightCandidate = *best->right;
            leftCandidate = *best->left;
            rightCandidateNorm = norm2(rightCandidate);
            leftCandidateNorm = norm2(leftCandidate);
            coupling = bilinear(leftCandidate, rightCandidate);
            couplingScale = std::max(
                rightCandidateNorm * leftCandidateNorm, 1.0e-300);
            if (std::abs(coupling) <= breakdownTolerance * couplingScale) {
                breakdownDetected_ = true;
                break;
            }
        }

        const Complex lower = std::sqrt(coupling);
        const Complex upper = coupling / lower;
        Vector nextRight = rightCandidate;
        Vector nextLeft = leftCandidate;
        for (auto& value : nextRight) {
            value /= lower;
        }
        for (auto& value : nextLeft) {
            value /= upper;
        }
        Vector represented = nextRight;
        for (auto& value : represented) {
            value *= lower;
        }
        recurrenceResidual_ = std::max(
            recurrenceResidual_,
            differenceNorm(rawRightCandidate, represented)
                / std::max(norm2(rightImage), 1.0e-300));
        beta.push_back(lower);
        gamma.push_back(upper);
        rightBasis.push_back(std::move(nextRight));
        leftBasis.push_back(std::move(nextLeft));
    }

    dimension_ = static_cast<int>(rightBasis.size());
    reducedOperator_.assign(static_cast<std::size_t>(dimension_ * dimension_),
                            Complex(0.0, 0.0));
    for (int i = 0; i < dimension_; ++i) {
        reducedOperator_[static_cast<std::size_t>(i * dimension_ + i)] =
            alpha[static_cast<std::size_t>(i)];
        if (i + 1 < dimension_) {
            reducedOperator_[static_cast<std::size_t>((i + 1) * dimension_ + i)] =
                beta[static_cast<std::size_t>(i)];
            reducedOperator_[static_cast<std::size_t>(i * dimension_ + i + 1)] =
                gamma[static_cast<std::size_t>(i)];
        }
    }

    reducedRightStart_.resize(static_cast<std::size_t>(dimension_));
    reducedOutput_.resize(static_cast<std::size_t>(dimension_));
    for (int i = 0; i < dimension_; ++i) {
        reducedRightStart_[static_cast<std::size_t>(i)] =
            bilinear(leftBasis[static_cast<std::size_t>(i)], rightStart);
        reducedOutput_[static_cast<std::size_t>(i)] =
            bilinear(leftOutput, rightBasis[static_cast<std::size_t>(i)]);
    }

    for (int i = 0; i < dimension_; ++i) {
        for (int j = 0; j < dimension_; ++j) {
            const Complex expected = i == j ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
            biorthogonalityError_ = std::max(
                biorthogonalityError_,
                std::abs(bilinear(leftBasis[static_cast<std::size_t>(i)],
                                  rightBasis[static_cast<std::size_t>(j)]) - expected));
        }
    }
    Complex projectedMoment(0.0, 0.0);
    for (int i = 0; i < dimension_; ++i) {
        projectedMoment += reducedOutput_[static_cast<std::size_t>(i)]
            * reducedRightStart_[static_cast<std::size_t>(i)];
    }
    const Complex exactMoment = bilinear(leftOutput, rightStart);
    momentMatchingError_ = std::abs(projectedMoment - exactMoment)
        / std::max(std::abs(exactMoment), 1.0e-300);

    if (retainLeftBasisForDiagnostics) {
        // 单元测试需要同时观察左右基，因此在小系统中保留两份。
        rightBasis_ = rightBasis;
        leftBasis_ = leftBasis;
    } else if (retainReconstructionBasis) {
        // 大网格只保留场重构所需右基，移动可避免 q 个全尺寸向量的峰值复制。
        rightBasis_ = std::move(rightBasis);
    }
    ready_ = true;
    buildPoleResidues();
}

// 用 Thomas 递推求解三对角 Padé 系统，病态时回退到稠密主元法。
LanczosPadeModel::Vector LanczosPadeModel::solveReduced(
    Complex parameter,
    int activeDimension) const {
    if (!ready_ || dimension_ == 0) {
        throw std::runtime_error("LanczosPadeModel::solveReduced called before build");
    }
    const int q = std::max(1, std::min(activeDimension, dimension_));
    Vector diagonal(static_cast<std::size_t>(q));
    Vector upper(static_cast<std::size_t>(std::max(0, q - 1)));
    Vector lower(static_cast<std::size_t>(std::max(0, q - 1)));
    Vector rhs(reducedRightStart_.begin(),
               reducedRightStart_.begin() + static_cast<std::ptrdiff_t>(q));
    for (int i = 0; i < q; ++i) {
        diagonal[static_cast<std::size_t>(i)] = Complex(1.0, 0.0)
            - parameter * reducedOperator_[static_cast<std::size_t>(i * dimension_ + i)];
        if (i + 1 < q) {
            upper[static_cast<std::size_t>(i)] =
                -parameter * reducedOperator_[static_cast<std::size_t>(i * dimension_ + i + 1)];
            lower[static_cast<std::size_t>(i)] =
                -parameter * reducedOperator_[static_cast<std::size_t>((i + 1) * dimension_ + i)];
        }
    }

    bool stable = true;
    constexpr double pivotTolerance = 1.0e-14;
    for (int i = 1; i < q; ++i) {
        const Complex pivot = diagonal[static_cast<std::size_t>(i - 1)];
        if (std::abs(pivot) <= pivotTolerance) {
            stable = false;
            break;
        }
        const Complex factor = lower[static_cast<std::size_t>(i - 1)] / pivot;
        diagonal[static_cast<std::size_t>(i)] -=
            factor * upper[static_cast<std::size_t>(i - 1)];
        rhs[static_cast<std::size_t>(i)] -= factor * rhs[static_cast<std::size_t>(i - 1)];
    }
    if (stable && std::abs(diagonal.back()) > pivotTolerance) {
        rhs.back() /= diagonal.back();
        for (int i = q - 2; i >= 0; --i) {
            if (std::abs(diagonal[static_cast<std::size_t>(i)]) <= pivotTolerance) {
                stable = false;
                break;
            }
            rhs[static_cast<std::size_t>(i)] =
                (rhs[static_cast<std::size_t>(i)]
                 - upper[static_cast<std::size_t>(i)] * rhs[static_cast<std::size_t>(i + 1)])
                / diagonal[static_cast<std::size_t>(i)];
        }
    } else {
        stable = false;
    }
    if (stable) {
        lastSolveSucceeded_ = true;
        return rhs;
    }

    Vector matrix(static_cast<std::size_t>(q * q), Complex(0.0, 0.0));
    for (int i = 0; i < q; ++i) {
        matrix[static_cast<std::size_t>(i * q + i)] = Complex(1.0, 0.0)
            - parameter * reducedOperator_[static_cast<std::size_t>(i * dimension_ + i)];
        if (i + 1 < q) {
            matrix[static_cast<std::size_t>(i * q + i + 1)] =
                -parameter * reducedOperator_[static_cast<std::size_t>(i * dimension_ + i + 1)];
            matrix[static_cast<std::size_t>((i + 1) * q + i)] =
                -parameter * reducedOperator_[static_cast<std::size_t>((i + 1) * dimension_ + i)];
        }
    }
    rhs.assign(reducedRightStart_.begin(),
               reducedRightStart_.begin() + static_cast<std::ptrdiff_t>(q));
    lastSolveSucceeded_ = denseSolve(matrix, rhs, q);
    if (!lastSolveSucceeded_) {
        throw std::runtime_error("LanczosPadeModel: reduced Padé system is singular");
    }
    return rhs;
}

// 用指定截断阶数评估三对角模型。
LanczosPadeModel::Complex LanczosPadeModel::evaluateAtOrder(
    Complex parameter,
    int activeDimension) const {
    const int q = std::max(1, std::min(activeDimension, dimension_));
    const auto coordinates = solveReduced(parameter, q);
    Complex value(0.0, 0.0);
    for (int i = 0; i < q; ++i) {
        value += reducedOutput_[static_cast<std::size_t>(i)]
            * coordinates[static_cast<std::size_t>(i)];
    }
    return value;
}

// 优先按极点-留数求和评估，特征分解不可用时使用三对角求解。
LanczosPadeModel::Complex LanczosPadeModel::evaluate(Complex parameter) const {
    if (poleResidueReady_) {
        Complex value(0.0, 0.0);
        for (std::size_t i = 0; i < poleEigenvalues_.size(); ++i) {
            value += poleWeights_[i]
                / (Complex(1.0, 0.0) - parameter * poleEigenvalues_[i]);
        }
        return value;
    }
    return evaluateAtOrder(parameter, dimension_);
}

// 以相邻 Padé 阶数差估计频带局部误差。
double LanczosPadeModel::relativeTruncationError(
    Complex parameter,
    int orderDrop) const {
    if (dimension_ <= 1) {
        return 0.0;
    }
    const int reducedOrder = std::max(1, dimension_ - std::max(1, orderDrop));
    const Complex full = evaluate(parameter);
    const Complex truncated = evaluateAtOrder(parameter, reducedOrder);
    return std::abs(full - truncated) / std::max(std::abs(full), 1.0e-12);
}

// 按右基线性组合低维坐标，恢复完整增广状态。
LanczosPadeModel::Vector LanczosPadeModel::reconstruct(Complex parameter) const {
    if (rightBasis_.empty()) {
        throw std::runtime_error("LanczosPadeModel: reconstruction basis was not retained");
    }
    const auto coordinates = solveReduced(parameter, dimension_);
    Vector state(rightBasis_.front().size(), Complex(0.0, 0.0));
    for (int j = 0; j < dimension_; ++j) {
        const Complex coefficient = coordinates[static_cast<std::size_t>(j)];
        const auto& basis = rightBasis_[static_cast<std::size_t>(j)];
        for (std::size_t i = 0; i < state.size(); ++i) {
            state[i] += coefficient * basis[i];
        }
    }
    return state;
}

// 标记带内且留数异常小的极点；仅用于诊断，不做筛除。
int LanczosPadeModel::countSuspectedSpuriousPoles(
    double parameterMin,
    double parameterMax) const {
    if (!poleResidueReady_) {
        return 0;
    }
    const double maxWeight = std::accumulate(
        poleWeights_.begin(), poleWeights_.end(), 0.0,
        [](double value, Complex weight) { return std::max(value, std::abs(weight)); });
    int count = 0;
    for (std::size_t i = 0; i < poleEigenvalues_.size(); ++i) {
        if (std::abs(poleEigenvalues_[i]) <= 1.0e-300) {
            continue;
        }
        const Complex pole = Complex(1.0, 0.0) / poleEigenvalues_[i];
        const double imaginaryTolerance = 1.0e-7 * std::max(1.0, std::abs(pole.real()));
        if (pole.real() >= parameterMin && pole.real() <= parameterMax
            && std::abs(pole.imag()) <= imaginaryTolerance
            && std::abs(poleWeights_[i]) <= 1.0e-8 * std::max(maxWeight, 1.0e-300)) {
            ++count;
        }
    }
    return count;
}

// 显式生成全阶和降阶矩序列，验证 Padé 的 2q 阶矩匹配性质。
double LanczosPadeModel::validateMomentMatching(
    const Operator& applyRight,
    int momentCount) const {
    if (!ready_ || !applyRight || momentCount < 1) {
        throw std::runtime_error("LanczosPadeModel: invalid moment validation request");
    }
    Vector fullMoment = originalRightStart_;
    Vector reducedMoment = reducedRightStart_;
    double maxError = 0.0;
    for (int order = 0; order < momentCount; ++order) {
        const Complex exact = bilinear(originalLeftOutput_, fullMoment);
        Complex reduced(0.0, 0.0);
        for (int i = 0; i < dimension_; ++i) {
            reduced += reducedOutput_[static_cast<std::size_t>(i)]
                * reducedMoment[static_cast<std::size_t>(i)];
        }
        maxError = std::max(
            maxError, std::abs(exact - reduced) / std::max(std::abs(exact), 1.0e-300));
        if (order + 1 < momentCount) {
            fullMoment = applyRight(fullMoment);
            Vector next(static_cast<std::size_t>(dimension_), Complex(0.0, 0.0));
            for (int i = 0; i < dimension_; ++i) {
                for (int j = 0; j < dimension_; ++j) {
                    next[static_cast<std::size_t>(i)] +=
                        reducedOperator_[static_cast<std::size_t>(i * dimension_ + j)]
                        * reducedMoment[static_cast<std::size_t>(j)];
                }
            }
            reducedMoment = std::move(next);
        }
    }
    return maxError;
}

// 对小三对角矩阵做特征分解，得到 h(t)=sum_i weight_i/(1-t*lambda_i)。
void LanczosPadeModel::buildPoleResidues() {
#ifdef BPFEM_USE_MKL
    const auto decompositionStarted = std::chrono::steady_clock::now();
    if (dimension_ == 0) {
        return;
    }
    const lapack_int n = static_cast<lapack_int>(dimension_);
    std::vector<MKL_Complex16> matrix(static_cast<std::size_t>(n * n));
    for (std::size_t i = 0; i < matrix.size(); ++i) {
        matrix[i] = MKL_Complex16{reducedOperator_[i].real(), reducedOperator_[i].imag()};
    }
    std::vector<MKL_Complex16> eigenvalues(static_cast<std::size_t>(n));
    std::vector<MKL_Complex16> rightEigenvectors(static_cast<std::size_t>(n * n));
    const lapack_int info = LAPACKE_zgeev(
        LAPACK_ROW_MAJOR, 'N', 'V', n, matrix.data(), n, eigenvalues.data(),
        nullptr, n, rightEigenvectors.data(), n);
    if (info != 0) {
        return;
    }

    std::vector<Complex> eigenvectorMatrix(static_cast<std::size_t>(n * n));
    poleEigenvalues_.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < dimension_; ++i) {
        poleEigenvalues_[static_cast<std::size_t>(i)] = Complex(
            eigenvalues[static_cast<std::size_t>(i)].real,
            eigenvalues[static_cast<std::size_t>(i)].imag);
        for (int j = 0; j < dimension_; ++j) {
            const auto value = rightEigenvectors[static_cast<std::size_t>(i * dimension_ + j)];
            eigenvectorMatrix[static_cast<std::size_t>(i * dimension_ + j)] =
                Complex(value.real, value.imag);
        }
    }
    Vector modalInput = reducedRightStart_;
    auto inverseSolve = eigenvectorMatrix;
    if (!denseSolve(inverseSolve, modalInput, dimension_)) {
        poleEigenvalues_.clear();
        return;
    }
    poleWeights_.assign(static_cast<std::size_t>(dimension_), Complex(0.0, 0.0));
    for (int mode = 0; mode < dimension_; ++mode) {
        Complex modalOutput(0.0, 0.0);
        for (int row = 0; row < dimension_; ++row) {
            modalOutput += reducedOutput_[static_cast<std::size_t>(row)]
                * eigenvectorMatrix[static_cast<std::size_t>(row * dimension_ + mode)];
        }
        poleWeights_[static_cast<std::size_t>(mode)] =
            modalOutput * modalInput[static_cast<std::size_t>(mode)];
    }

    const Complex probes[] = {
        Complex(-0.75, 0.0), Complex(-0.25, 0.0), Complex(0.0, 0.0),
        Complex(0.25, 0.0), Complex(0.75, 0.0)
    };
    for (Complex parameter : probes) {
        const Complex reference = evaluateAtOrder(parameter, dimension_);
        Complex poleValue(0.0, 0.0);
        for (std::size_t i = 0; i < poleEigenvalues_.size(); ++i) {
            poleValue += poleWeights_[i]
                / (Complex(1.0, 0.0) - parameter * poleEigenvalues_[i]);
        }
        poleResidueReconstructionError_ = std::max(
            poleResidueReconstructionError_,
            std::abs(reference - poleValue) / std::max(std::abs(reference), 1.0e-300));
    }
    poleResidueReady_ = poleResidueReconstructionError_ <= 1.0e-10;
    if (!poleResidueReady_) {
        poleEigenvalues_.clear();
        poleWeights_.clear();
    }
    poleDecompositionSec_ = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - decompositionStarted).count();
#endif
}

}  // namespace fem::fastsweep