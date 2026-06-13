# WCAWE Skill：良条件 AWE

## 适用目标

WCAWE（Well-Conditioned Asymptotic Waveform Evaluation）用于解决传统 AWE 矩向量病态的问题。
它直接处理 FEM 产生的多项式矩阵方程，不需要通过线性化引入额外未知量，也避免了传统 AWE 高阶矩
条件数快速恶化。

WCAWE 的核心不是外部频段偏好策略，而是用非奇异上三角系数矩阵把 AWE 矩向量组合成良条件基。
当该系数矩阵取单位阵时，WCAWE 退化为传统 AWE；当该系数矩阵取 modified Gram-Schmidt 的
正交化系数时，WCAWE 得到稳定的 Arnoldi-like 基。

## 输入

- `polynomialModel`：多项式矩阵方程 `A(t)x(t)=b(t)`，含 `A_0...A_d` 和 `b_0...b_d`。
- `expansionPoint`：展开波数/频率，初版使用单点；多点版本可复用 MGAWE 的展开点管理。
- `order`：目标 WCAWE 基阶数。
- `aweRecurrence`：传统 AWE 系统矩递推接口，用于定义矩匹配目标。
- `orthogonalizer`：modified Gram-Schmidt + reorthogonalization。
- `triangularConfig`：上三角系数矩阵 `R` 的构造、存储和奇异性阈值。
- `projectionConfig`：用良条件基构造 Galerkin reduced model 的投影约定。

## 输出

- `wellConditionedBasis V`：良条件正交/近正交基。
- `triangularCoefficients R`：连接传统 AWE 矩向量和 WCAWE 基的非奇异上三角矩阵。
- `correctionTerms`：WCAWE 递推中的校正项与每阶来源记录。
- `reducedModel`：投影后的 Galerkin 降阶模型。
- `conditionCurves`：传统 AWE 矩基与 WCAWE 基的条件数随阶数变化曲线。
- `diagnostics`：正交性、矩匹配保持性、direct 抽检误差和无额外未知量检查结果。

## 算法步骤

1. 用 `PolynomialMomentRecurrence` 定义传统 AWE 矩序列 `{x_0, x_1, ...}` 的目标递推。
2. 初始化 WCAWE 候选向量，使其张成与对应 AWE 矩相同的子空间。
3. 对每一阶候选向量执行 modified Gram-Schmidt，得到新的基向量和上三角系数。
4. 根据已生成的上三角系数矩阵 `R` 构造 WCAWE 校正项，使新基仍保持 AWE 的矩匹配性质。
5. 检查 `R` 的对角元和条件数；若接近奇异，停止增加阶数并报告饱和。
6. 用良条件基 `V` 构造 Galerkin reduced model。
7. 在线扫频时求解小系统并输出 S 参数；需要场时使用 `x_hat=V x_r`。
8. 输出 AWE 原始矩基与 WCAWE 基的条件数曲线，证明良条件化确实发生。

## 实现建议

- `WellConditionedBasisBuilder` 应同时保存 `V` 和 `R`，因为可靠性诊断需要二者。
- WCAWE 的递推必须显式调用传统 AWE 矩目标，不能只做普通 Krylov 正交化后宣称等价。
- 初版不需要 Padé；WCAWE 更适合作为 Galerkin reduced model 的稳定基生成器。
- 多点 WCAWE 可以后续作为 MGAWE 的局部向量生成器，不在第一版文档实现范围内。

## 风险点

- 若只正交化 AWE 矩向量而不加入 WCAWE 校正项，会破坏论文中的矩匹配性质。
- `R` 接近奇异时继续增阶会制造虚假精度。
- WCAWE 比传统 AWE 构造稍贵，但应通过减少病态和减少无效高阶向量来换取整体可靠性。
- 当前项目的端口 `beta(omega)` 非多项式，需要先在展开点局部多项式化或通过仿射导数接口处理。

## 可靠性验证

- **条件数对比**：输出传统 AWE 矩基和 WCAWE 基的 `cond(V_n)` 曲线，WCAWE 应明显更稳定。
- **MGS 正交性**：检查 `||V^T V-I||` 或当前 bilinear 内积下的等价指标。
- **矩匹配保持性**：小系统上验证 WCAWE reduced model 与传统 AWE 匹配同阶矩。
- **Arnoldi 对照**：在线性化的小型多项式系统上，WCAWE 结果应接近 Arnoldi，但不增加全空间未知量。
- **无额外未知量检查**：诊断中记录 WCAWE 基向量长度仍为原 FEM DOF，不使用线性化扩维系统。
- **宽带 direct FEM 对比**：在 BP filter 抽检点上比较 direct，目标同 MGAWE：`max delta <= 0.05 dB`。
