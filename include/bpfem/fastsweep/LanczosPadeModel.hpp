#pragma once

#include <complex>
#include <functional>
#include <utility>
#include <vector>

namespace fem::fastsweep {

// 标准非对称双边 Lanczos-Padé 单输入单输出模型。
// 左右过程使用非共轭双线性积，并保持 W^T V=I。
class LanczosPadeModel {
public:
    using Complex = std::complex<double>;
    using Vector = std::vector<Complex>;
    using Operator = std::function<Vector(const Vector&)>;
    using PairedOperator = std::function<std::pair<Vector, Vector>(
        const Vector&, const Vector&)>;

    // 用独立的 G 与 G^T 算子构造 q 阶 [q-1/q] Padé 模型。
    void build(const Vector& rightStart,
               const Vector& leftOutput,
               const Operator& applyRight,
               const Operator& applyTranspose,
               int order,
               double breakdownTolerance = 1.0e-12,
               bool retainReconstructionBasis = true,
               bool retainLeftBasisForDiagnostics = false);

    // 用成对算子构造模型；直接求解器可借此批量回代右、左两个 RHS。
    void buildPaired(const Vector& rightStart,
                     const Vector& leftOutput,
                     const Operator& applyRight,
                     const PairedOperator& applyPair,
                     int order,
                     double breakdownTolerance = 1.0e-12,
                     bool retainReconstructionBasis = true,
                     bool retainLeftBasisForDiagnostics = false);

    // 计算 c^T(I-tG)^(-1)r0；优先使用已验证的极点-留数表达。
    Complex evaluate(Complex parameter) const;

    // 用 q 阶与 q-orderDrop 阶模型差异估计局部截断误差。
    double relativeTruncationError(Complex parameter, int orderDrop = 2) const;

    // 将 Padé 坐标按右 Lanczos 基恢复到完整增广状态。
    Vector reconstruct(Complex parameter) const;

    // 统计参数区间内小留数的疑似伪极点，不改变任何极点或求解结果。
    int countSuspectedSpuriousPoles(double parameterMin, double parameterMax) const;

    // 用显式全阶算子比较前 momentCount 个标量矩，供单元测试使用。
    double validateMomentMatching(const Operator& applyRight, int momentCount) const;

    int dimension() const { return dimension_; }
    bool ready() const { return ready_; }
    bool breakdownDetected() const { return breakdownDetected_; }
    bool happyBreakdown() const { return happyBreakdown_; }
    bool lastSolveSucceeded() const { return lastSolveSucceeded_; }
    bool poleResidueReady() const { return poleResidueReady_; }
    double biorthogonalityError() const { return biorthogonalityError_; }
    double tridiagonalLeakage() const { return tridiagonalLeakage_; }
    double finalCouplingMagnitude() const { return finalCouplingMagnitude_; }
    double recurrenceResidual() const { return recurrenceResidual_; }
    double momentMatchingError() const { return momentMatchingError_; }
    double poleResidueReconstructionError() const { return poleResidueReconstructionError_; }
    double poleDecompositionSec() const { return poleDecompositionSec_; }
    int lookAheadCount() const { return lookAheadCount_; }
    int selectiveReorthogonalizationCount() const {
        return selectiveReorthogonalizationCount_;
    }
    const std::vector<Vector>& rightBasis() const { return rightBasis_; }
    const std::vector<Vector>& leftBasis() const { return leftBasis_; }
    const std::vector<Complex>& reducedOperator() const { return reducedOperator_; }

private:
    // 共享实现负责标准三项递推；paired 为空时逐个调用左右算子。
    void buildImpl(const Vector& rightStart,
                   const Vector& leftOutput,
                   const Operator& applyRight,
                   const Operator& applyTranspose,
                   const PairedOperator& applyPair,
                   int order,
                   double breakdownTolerance,
                   bool retainReconstructionBasis,
                   bool retainLeftBasisForDiagnostics);

    // 求解截断后的三对角 (I-tT)y=W^T r0。
    Vector solveReduced(Complex parameter, int activeDimension) const;

    // 用指定阶数评估同一 Lanczos 序列。
    Complex evaluateAtOrder(Complex parameter, int activeDimension) const;

    // 在 MKL 构建中对 T 做特征分解并生成极点-留数权重。
    void buildPoleResidues();

    int dimension_ = 0;
    bool ready_ = false;
    bool breakdownDetected_ = false;
    bool happyBreakdown_ = false;
    bool poleResidueReady_ = false;
    mutable bool lastSolveSucceeded_ = true;
    double biorthogonalityError_ = 0.0;
    double tridiagonalLeakage_ = 0.0;
    double finalCouplingMagnitude_ = 0.0;
    double recurrenceResidual_ = 0.0;
    double momentMatchingError_ = 0.0;
    double poleResidueReconstructionError_ = 0.0;
    double poleDecompositionSec_ = 0.0;
    int lookAheadCount_ = 0;
    int selectiveReorthogonalizationCount_ = 0;
    Vector originalRightStart_;
    Vector originalLeftOutput_;
    std::vector<Vector> rightBasis_;
    std::vector<Vector> leftBasis_;
    std::vector<Complex> reducedOperator_;
    std::vector<Complex> reducedRightStart_;
    std::vector<Complex> reducedOutput_;
    std::vector<Complex> poleEigenvalues_;
    std::vector<Complex> poleWeights_;
};

}  // namespace fem::fastsweep