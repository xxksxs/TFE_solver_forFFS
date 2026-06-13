# AWE / GAWE / MGAWE / WCAWE 快速扫频算法族

本文档组记录 AWE、GAWE、MGAWE、WCAWE 四类 FEM 快速扫频算法的工程 skill。它们的共同目标是：
在少量展开点上提取矩信息，构造低阶传递模型，从而避免对每个频率点都进行一次完整 FEM 求解。

本轮文档定义算法、输入输出、模块边界和可靠性验证；当前工程已逐步接入 `--sweep awe|gawe|mgawe|wcawe`。

## 论文对齐后的定位

三篇论文给出的分工很清楚：

1. **AWE**：在单个展开波数/频率附近对 FEM 解向量做 Taylor 矩递推，再把感兴趣的标量响应
   转成 Padé 有理函数。输出通常是输入阻抗、S 参数、远场或某个场分量。
2. **GAWE**：在单个展开点生成 AWE 矩向量，正交化后构造 Galerkin 降阶空间。它不走 Padé 标量外推，而是在线求解小型 reduced FEM 系统。
3. **MGAWE**：把多个展开点产生的 AWE/GAWE 向量一起正交化，形成统一 Galerkin 降阶空间。
   重点是宽带精度和残差正交，不是多个局部 ROM 的简单拼接。
4. **WCAWE**：全称是 **Well-Conditioned Asymptotic Waveform Evaluation**。它通过非奇异上三角
   系数矩阵把传统 AWE 矩向量组合为良条件基；当系数矩阵取单位阵时退化为 AWE，当系数来自
   modified Gram-Schmidt 时得到稳定的 Arnoldi-like 过程。

因此，WCAWE 的 `W` 表示良条件化，不应设计为外部频率/端口偏好策略。

## 推荐模块边界

后续实现时建议新增独立 `fastsweep` 模块，`src/sweep/` 只保留 CLI 策略适配层：

```text
include/bpfem/fastsweep/
  AffineModelView.hpp
  PolynomialMomentRecurrence.hpp
  PadeApproximant.hpp
  GalerkinReducedModel.hpp
  WellConditionedBasisBuilder.hpp
  ReliabilityDiagnostics.hpp
  AweSweep.hpp
  GaweSweep.hpp
  MgaweSweep.hpp
  WcaweSweep.hpp

src/fastsweep/
  PolynomialMomentRecurrence.cpp
  PadeApproximant.cpp
  GalerkinReducedModel.cpp
  WellConditionedBasisBuilder.cpp
  ReliabilityDiagnostics.cpp
  AweSweep.cpp
  GaweSweep.cpp
  MgaweSweep.cpp
  WcaweSweep.cpp
```

## 共享模块职责

- `AffineModelView`：只读包装 `FEMAssembler::buildAffineSystem()`。提供 `A(omega)`、`b(omega)`、
  输出投影、端口传播常数和可选频率导数。
- `PolynomialMomentRecurrence`：处理 `A(k)=A0+kA1+k^2A2+...` 或项目仿射形式的系统矩递推。
- `PadeApproximant`：由标量/端口传递函数矩生成 Padé 分子、分母、极点和留数。
- `GalerkinReducedModel`：保存全局基 `V`，投影后的矩阵/RHS/输出，并执行在线小系统求解。
- `WellConditionedBasisBuilder`：保存上三角系数矩阵、MGS 系数和 WCAWE 校正项，生成良条件基。
- `ReliabilityDiagnostics`：统一记录矩匹配误差、残差、条件数、极点位置、正交性和 direct 抽检误差。

## Skill 列表

- [AWE skill](awe-skill.md)：单展开点系统矩递推 + Padé 极点留数近似。
- GAWE：单展开点系统矩递推 + 正交化 Galerkin ROM；工程上可复用 MGAWE 的单点路径。
- [MGAWE skill](mgawe-skill.md)：多个展开点同时构造统一 Galerkin 降阶空间。
- [WCAWE skill](wcawe-skill.md)：良条件 AWE，通过上三角正交化系数稳定矩基。
- [统一校验方法](validation.md)：三类算法进入实现前后的数值可靠性和工程验收。

## 与现有 ALPS/Krylov 文档的关系

当前 `--sweep alps` 已实现单展开点、块 shift-and-invert Krylov、复对称 Galerkin ROM。它和
GAWE/MGAWE/WCAWE 的工程基础高度重合：都需要仿射系统、展开点因式分解、基生成、正交化、投影和
S 参数还原。

建议演进顺序：

1. 抽出 `AffineModelView`、`GalerkinReducedModel`、`ReliabilityDiagnostics`。
2. 把现有 ALPS 中的单点 Galerkin ROM 逻辑迁到 fast-sweep 基建。
3. 实现 MGAWE 的多展开点统一 Galerkin 空间。
4. 实现 WCAWE 的良条件基构造，作为单点/多点矩基的稳定版本。
5. 最后补显式 Padé AWE，因为它对 Hankel 条件数和 spurious pole 处理最敏感。
