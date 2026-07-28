# WCAWE 代码实现规范 Skill

## 1. 文档地位

本 skill 规定 Slone、Lee 与 Lee（2003）WCAWE 在 BP-FEM 中的 C++ 实现边界、输入输出、索引、异常语义和验收条件。

数学依据：

- [WCAWE 数学推导 LaTeX](wcawe-theory-cn.tex)
- [WCAWE 数学推导 PDF](wcawe-theory-cn.pdf)
- [WCAWE 理论入口](wcawe-skill.md)

若本 skill 与论文公式冲突，以论文式 (7)–(9) 和数学推导 PDF 为准。实现不得以现有类名或旧结果作为偏离论文的理由。

## 2. 实现目标与非目标

### 2.1 必须实现

- 直接处理参数 τ 下的多项式系统 A(τ)x(τ)=b(τ)。
- 逐阶生成 WCAWE 候选向量；第 n 阶必须读取已生成的 V_{n-1} 和 U_{n-1}。
- 实现论文 P_U1/P_U2 连续主子块逆乘积作用。
- 使用 Hermitian MGS 更新论文上三角矩阵 U。
- 一次分解 A0，每个成功候选仅执行一次全阶 RHS 回代。
- 用生成的 V 构造现有 GalerkinReducedModel。
- 输出足以证明式 (7)、式 (8)、正交性和矩匹配成立的诊断。

### 2.2 本轮非目标

- 不实现多展开点 WCAWE、block RHS WCAWE 或自适应展开点。
- 不实现有损材料路径。
- 不显式求 U 子块逆。
- 不改变 GAWE/MGAWE/AWE 的数学行为。
- 不要求把在线 reduced model 强制写成 Padé 分子分母。

## 3. 强制数学约定

统一使用归一化展开参数 τ：

$$
A(\tau)=\sum_{r=0}^{d_A}\tau^r A_r,
\qquad
b(\tau)=\sum_{r=0}^{d_b}\tau^r b_r.
$$

当前无损端口适配器使用：

$$
\tau=\frac{k_0^2-\lambda_0}{\lambda_{scale}}.
$$

所有 A_r、b_r、展开点和有限差分验证必须使用同一个 τ。禁止把关于 ω、s=jω、k0² 的导数混在一次递推中。

核心关系：

$$
V_n=\widetilde V_n U_n^{-1},
\qquad
\widetilde V_n=V_nU_n.
$$

第 n 阶候选必须按数学推导中的论文式 (7) 生成。U 是下一阶递推输入，不是全部矩生成完成后的 QR 诊断量。

## 4. 目标文件与职责

新增：

- include/bpfem/fastsweep/PolynomialWcaweRecurrence.hpp
- src/fastsweep/PolynomialWcaweRecurrence.cpp
- include/bpfem/fastsweep/UpperTriangularBlockAction.hpp
- src/fastsweep/UpperTriangularBlockAction.cpp

修改：

- WellConditionedBasisBuilder.hpp/.cpp
- WcaweSweep.hpp/.cpp
- FastSweepDiagnostics.hpp/.cpp
- tests/fastsweep_tests.cpp
- CMakeLists.txt

不得把论文递推全部塞入 WcaweSweep.cpp。Sweep 层只负责工程适配、求解器会话、ROM 投影、输出和日志。

## 5. 通用多项式核心接口

建议接口如下；实现可微调命名，但不能改变数据语义。

~~~cpp
namespace fem::fastsweep {

using Complex = std::complex<double>;
using Vector = std::vector<Complex>;
using LinearOperator = std::function<Vector(const Vector&)>;
using RhsCoefficientFunction = std::function<Vector(std::size_t order)>;
using SolveFunction = std::function<Vector(const Vector& rhs)>;

struct PolynomialWcaweModel {
    std::size_t fullDimension = 0;
    std::size_t matrixDegree = 0;
    std::size_t rhsDegree = 0;

    // applyA0(x) = A0*x，用于递推残差；不负责求解。
    LinearOperator applyA0;

    // applyHigherOrder[r-1](x) = Ar*x，r=1...matrixDegree。
    std::vector<LinearOperator> applyHigherOrder;

    // order > rhsDegree 时必须返回 fullDimension 个零。
    RhsCoefficientFunction rhsCoefficientAt;
};

struct WcaweBuildOptions {
    int requestedOrder = 12;
    double breakdownTolerance = 1.0e-12;
    double recurrenceTolerance = 1.0e-10;
    double orthogonalityTolerance = 1.0e-10;
    int reorthogonalizationPasses = 2;
};

enum class WcaweTerminationReason {
    ReachedRequestedOrder,
    HappyBreakdown,
    NearSingularUpperTriangularFactor,
    SolverFailure,
    NonFiniteValue
};

struct WcaweStepDiagnostics {
    int order = 0;
    double candidateNormBeforeMgs = 0.0;
    double diagonalAbs = 0.0;
    double recurrenceResidual = 0.0;
    double basisRelationResidual = 0.0;
    double orthogonalityError = 0.0;
    int reorthogonalizationPasses = 0;
    std::uint64_t triangularSolveCount = 0;
};

struct WcaweBuildResult {
    std::vector<Vector> basis;
    std::vector<Complex> upperTriangularU; // row-major, stride=requestedOrder
    int requestedOrder = 0;
    int achievedOrder = 0;
    WcaweTerminationReason terminationReason{};
    std::vector<WcaweStepDiagnostics> steps;

    double maxRecurrenceResidual = 0.0;
    double maxBasisRelationResidual = 0.0;
    double maxOrthogonalityError = 0.0;
    double minUpperDiagonalAbs = 0.0;
    std::uint64_t triangularSolveCount = 0;
};

class PolynomialWcaweRecurrence {
public:
    static WcaweBuildResult build(
        const PolynomialWcaweModel& model,
        const SolveFunction& solveAtExpansion,
        const WcaweBuildOptions& options = {});
};

} // namespace fem::fastsweep
~~~

### 5.1 输入契约

PolynomialWcaweModel：

- fullDimension 必须大于零。
- matrixDegree 至少为 1；P1 适配器取 1，通用数学单测必须覆盖 2。
- applyHigherOrder.size() 必须等于 matrixDegree。
- 所有 operator 输入输出长度必须等于 fullDimension。
- rhsCoefficientAt(k) 必须返回 fullDimension 个元素；k>rhsDegree 时返回零向量。
- 不允许回调返回 NaN/Inf。

solveAtExpansion：

- 语义固定为求解 A0*x=rhs。
- 不得在每次调用内重新装配 A0。
- WcaweSweep 必须用 FactorizedSolveSession 包装，使 PARDISO 只分解一次。
- 返回 SolveResult 失败时，适配 lambda 必须抛出带阶数上下文的异常。

WcaweBuildOptions：

- requestedOrder >= 1。
- reorthogonalizationPasses 只能为 1 或 2，默认 2。
- tolerance 必须有限且大于零。

### 5.2 输出契约

WcaweBuildResult：

- basis.size() 等于 achievedOrder。
- 每个 basis 列长度等于 fullDimension。
- upperTriangularU 固定分配 requestedOrder²，row-major；未生成区域保持零。
- U 的有效 leading achievedOrder×achievedOrder 区域必须上三角。
- 成功达到目标时 terminationReason=ReachedRequestedOrder。
- breakdown 后保留已经成功生成的基，但不得跳过失败列继续递推。
- achievedOrder>=1 时结果可以交给 GalerkinReducedModel；零维必须抛错。

## 6. 一基索引与 C++ 存储

论文 n、m、t 使用一基索引。U 使用 row-major 固定 stride=q：

~~~cpp
U[(rowOneBased - 1) * q + (colOneBased - 1)]
~~~

对论文

$$
P_{U_w}(n,m)=\prod_{t=w}^{m}
U_{[t:n-m+t-1,\;t:n-m+t-1]}^{-1},
$$

令 r=n-m。每个子块尺寸均为 r×r；一基起点 t 转成零基 t-1。

计算 P_{U_w}(n,m)e_r 时，矩阵乘积写作 B_w^{-1}...B_m^{-1}，因此作用到向量的三角求解顺序必须反向：

~~~cpp
y = e_r;
for (t = m; t >= w; --t) {
    // B_t = U[t : n-m+t-1, t : n-m+t-1]，论文一基闭区间。
    // 解 B_t * z = y，然后 y = z。
    y = solveUpperTriangularSubblock(U, q, t - 1, r, y);
}
return y;
~~~

必须用随机上三角小矩阵把该结果与显式高精度参考比较。只检查 U=I 不足以发现乘积顺序错误。

## 7. WellConditionedBasisBuilder 规范

旧接口 append(moment) 的语义必须改为 appendCandidate(candidate)，因为输入必须是已经包含论文校正项的 WCAWE 候选，而不是传统 AWE 矩。

建议结果：

~~~cpp
struct WcaweAppendResult {
    bool accepted = false;
    int column = 0;
    double preNorm = 0.0;
    double diagonalAbs = 0.0;
    double basisRelationResidual = 0.0;
    double orthogonalityError = 0.0;
    int passes = 0;
};
~~~

强制行为：

- MGS 内积使用 hdot，即 v_j^H*candidate。
- 两遍 MGS 时，第二遍系数累加到同一个 U(j,n)，不能覆盖第一遍。
- U(n,n)=MGS 后范数。
- 候选副本只在当前步骤保留，用于检查 candidate≈V_n*U(:,n)。
- 若 U(n,n)<=breakdownTolerance*max(preNorm,1)，终止整条单点递推。
- 不允许将该列计为 deflation 后继续 n+1，因为后续 P_U 需要非奇异 U。
- 不允许列交换或 pivot；论文递推依赖固定阶次。

## 8. PolynomialWcaweRecurrence 离线流程

1. 校验 model/options。
2. 固定分配 q×q 的 U；basis.reserve(q)。
3. n=1：求解 A0*vtilde1=b0，appendCandidate。
4. 对 n=2...q：
   - 对每个 RHS 阶 m 计算 e1^T P_U1(n,m)e_{n-m}。
   - 计算 -A1*v_{n-1}。
   - 对 m=2...min(dA,n-1) 计算 -A_m*V_{n-m}*P_U2(n,m)e_{n-m}。
   - 在一个全阶 rhs scratch 中累加所有项。
   - 调用 solveAtExpansion(rhs) 一次得到候选。
   - 在覆盖候选前保留当前候选副本或即时重构残差。
   - appendCandidate 并记录式 (7)、式 (8) 和正交性诊断。
5. 达到 q 或发生 breakdown 后返回结果。

每个成功阶只允许一次 A0 回代。P_U 三角回代是 q 维小系统，不得调用 PARDISO。

## 9. 当前 P1 工程适配

WcaweSweep::buildOffline 必须改为：

1. buildAffineSystem() 和 buildPortVectors()。
2. buildLosslessLinearization() 得到 matrixAtExpansion、rhs0、rhs1、lambda0、lambdaScale。
3. 构造 PolynomialWcaweModel：
   - applyA0(x)=matrixAtExpansion.multiply(x)。
   - matrixDegree=1。
   - applyHigherOrder[0](x)=applyLosslessFirstOrderMatrix(...,x)。
   - rhsDegree=1；rhsCoefficientAt(0)=rhs0，(1)=rhs1，其余为零。
4. 创建 FactorizedSolveSession(solver,matrixAtExpansion,config)。
5. solveAtExpansion 调用 session.solve(rhs)，并验证 converged、field 长度和有限值。
6. 调用 PolynomialWcaweRecurrence::build()。
7. 将 result.basis move 给 GalerkinReducedModel::build()。
8. 保存 U 和 step diagnostics；删除 buildAweMoments() 主路径。

当前 GalerkinReducedModel 可以继续投影 K/M/portVectors 并在线使用精确端口频率律。它属于“论文基 + 更精确 reduced 频率律”的工程扩展，但必须满足：在展开点处及 P1 导数层面与离线 A0/A1、b0/b1 一致。

## 10. 所有权、速度与内存

- WCAWE 核心只持有 basis、U、一个 candidate、一个 rhs 和少量 operator scratch。
- 禁止同时保存完整传统 AWE moments 和 WCAWE basis。
- basis 移交 ROM 时使用 std::move。
- A0 只允许一次 numeric factorization。
- q 阶 P1 WCAWE 应为 q 个 factorized RHS、q 次依赖串行回代调用；不同阶不能批量，因为后阶依赖前阶 U。
- full-size scratch 应复用容量，避免每个来源项重复分配 N 个 Complex。
- 目标内存复杂度：N*q 个基元素 + O(N) scratch + O(q²) reduced 数据。

## 11. 异常与终止语义

直接抛 std::invalid_argument：

- 空模型、维度不一致、次数与 operator 数量不一致、非法阈值。

直接抛 std::runtime_error：

- A0 求解失败、回调返回错误长度、出现 NaN/Inf、零阶候选失败。

返回部分有效 ROM 并记录 terminationReason：

- n>1 时 U(n,n) 近零或 happy breakdown。

工程层不允许静默回退 direct。以下情况必须终止 WCAWE、写明原因并建议显式改用 direct：

- 有损材料。
- 展开频带跨越端口截止点，局部端口导数不可靠。
- achievedOrder 低于配置的最小可用阶数。
- 递推残差或矩匹配超过硬阈值。

日志必须明确区分“正常达到目标阶数”“happy breakdown”“近奇异 U”和“数值失败”。允许返回部分有效 ROM 的 breakdown 必须写入 achievedOrder 和 terminationReason，不能继续下一阶。

## 12. diagnostics.json 新字段

新增并实际写入：

- wcawe_recurrence_residual_max
- wcawe_basis_relation_residual_max
- wcawe_orthogonality_error
- wcawe_min_u_diagonal
- wcawe_u_diagonal_ratio
- wcawe_triangular_solve_count
- wcawe_breakdown_order
- wcawe_termination_reason
- wcawe_moment_matching_error
- wcawe_reorthogonalization_count
- wcawe_breakdown_detected
- wcawe_termination_reason
- wcawe_identity_limit_error（仅单测/诊断构建）

basis_condition.csv 改为：

- order
- awe_moment_condition_proxy（仅诊断小样或流式代理，不得要求保存全部 AWE 矩）
- wcawe_basis_orthogonality_error
- u_diagonal_abs
- u_diagonal_ratio
- recurrence_residual
- basis_relation_residual

旧字段 X≈VR 不能继续作为论文主验收项。

## 13. 单元测试与硬阈值

### 13.1 UpperTriangularBlockAction

- q=2...8 的随机复上三角 U，对所有合法 n,m,w 比较显式参考。
- 相对误差 <=1e-13。
- 覆盖非单位对角、强非正规上三角和接近但未达到奇异的情况。
- 非法下标、零对角和越界必须抛异常。

### 13.2 U=I 退化

- 使用 dA=2、db=2 的小型复多项式系统。
- WCAWE 强制 U=I 时逐列等于 PolynomialMomentRecurrence。
- 每列相对误差 <=1e-12。

### 13.3 论文式 (7) 与式 (8)

- 每阶 recurrenceResidual <=1e-11。
- 每阶 basisRelationResidual <=1e-12。
- max ||V_n^H V_n-I|| <=1e-11。
- U 有效区域下三角元素严格为零。

### 13.4 矩匹配

- 用解析可逆的小型 dA=2 系统计算 full transfer Taylor 系数。
- reduced model 至少匹配前 q 个系统矩。
- 相对误差 <=1e-10。
- 必须加入一个论文朴素式 (6) 失败的反例，证明测试能区分真正 WCAWE。

### 13.5 Breakdown

- happy breakdown：返回部分 ROM，不抛异常。
- 首列零向量：抛 runtime_error。
- U 对角近零：停止且 achievedOrder 不包含失败列。
- breakdown 后不得继续访问下一阶 P_U。

### 13.6 求解器生命周期

P1、q=12：

- symbolic analysis count=1（模式首次出现时）。
- numeric factorization count=1。
- factorized RHS vector count=12。
- factorized solve call count=12。
- 优化前后 A0 回代解相对差异 <=1e-12。

## 14. 工程集成验收

### 14.1 五点烟测

命令继续使用：

~~~powershell
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 0 --max-sweep-points 5 --sweep wcawe --wcawe-order 12 --linear-solver direct --no-write-all-fields --out result
~~~

result/result_WCAWE 必须包含：

- s_parameters.csv
- run.log
- run.json
- timing.json
- diagnostics.json
- basis_condition.csv
- field_last.vtu

### 14.2 展开点精确性

95 GHz 展开点：

- WCAWE 与 direct 的复数 S11、S21 相对误差 <=1e-10。
- reduced/full 残差 <=1e-10。
- P1 A1*x 和 b1 的有限差分相对误差 <=1e-7。

### 14.3 当前 30 万网格、101 点基准

固定：

- current.ngmesh，basis-order=0。
- 90–100 GHz，101 点。
- HFSS 参考 `IOStructure_S_parameters.csv`。
- PARDISO Release，wcawe-order=12。
- 只比较 S11/S21，主指标为幅值相对 L2；深零点另报绝对误差。

最低验收：

- HFSS S11 幅值相对 L2 <=1.0%。
- HFSS S21 幅值相对 L2 <=1.2%。
- 不得比相同基和端口定义下的当前 GAWE q=12 相对 L2 恶化超过 max(0.01 个百分点, GAWE误差的 0.1%)。
- 峰值内存不超过当前 GAWE/WCAWE 基线的 105%。
- q=12 必须为 1 次数值分解、12 个 RHS、12 次回代调用。

### 14.4 高阶病态压力测试

WCAWE 的优势不能只用 q=12 简单算例判断。至少测试 q={12,20,30,40}：

- WCAWE 正交误差保持 <=1e-10，或在 U 近奇异时明确停止。
- 不允许输出 NaN/Inf 或把失败列计入 ROM。
- 与传统 AWE/GAWE 比较矩匹配停滞阶数。
- 在达到同一 direct/HFSS 精度时，记录有效阶数、离线时间和内存。
- 只有在“更晚停滞、更低高阶误差或更少有效列达到同精度”至少一项成立时，才宣称 WCAWE 在该算例有优势。

## 15. 完成定义

代码只有同时满足以下条件才可将当前实现标记为论文 WCAWE：

1. WcaweSweep 不再调用 generateLosslessMoments() 生成完整传统矩主路径。
2. 存在独立 PolynomialWcaweRecurrence 和 UpperTriangularBlockAction。
3. 第 n 阶候选确实使用第 n-1 阶 U。
4. 式 (7)、式 (8)、正交性、矩匹配测试全部通过。
5. Debug/Release 构建及全部 CTest 通过。
6. 5 点烟测和 101 点工程基准通过。
7. diagnostics.json 能证明递推而不只是证明 QR 后基正交。
8. 文档、类注释和日志不再把“传统 AWE 矩事后 MGS”描述为 WCAWE。
