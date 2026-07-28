# AWE / GAWE / MGAWE / WCAWE 统一校验方法

任何快速扫频算法进入实现前后，都必须同时通过数学单测、论文一致性检查、项目集成测试、HFSS 基准对比和失效诊断。
本文件定义校验口径；其中 HFSS 文件是 BP filter 的外部物理基准，direct FEM 结果是工程回归和定位问题用的内部基准。

WCAWE 的论文递推、接口和硬阈值以 [WCAWE 代码实现规范](wcawe-implementation-skill.md) 为准；本文件的通用条款不得覆盖该规范。

## 0. BP Filter 外部基准

BP filter 的首要校验文件为项目根目录下的 `wg_bp_filter_S_parameters.csv`，该文件来自 HFSS：

- 频率列：`Freq [GHz]`。
- 反射幅度列：`mag(S(1,1)) []`。
- 传输幅度列：`mag(S(2,1)) []`。

所有 sweep 输出的 `s_parameters.csv` 必须先换算出：

- `|S11| = sqrt(S11_real^2 + S11_imag^2)`。
- `|S21| = sqrt(S21_real^2 + S21_imag^2)`。
- `Sij_dB = 20 log10(|Sij|)`。

验收优先级：

1. 与 HFSS CSV 在相同 GHz 频点对齐比较，这是最终验收口径。
2. 与 direct FEM 输出比较，这是快速定位实现误差的内部口径。
3. 若 HFSS 与 direct FEM 出现差异，必须先记录差异来源，不得只用 direct FEM 结果替代 HFSS 验收。

推荐输出 `hfss_comparison.csv`，至少包含：

- `freq_GHz`。
- `hfss_abs_s11`、`solver_abs_s11`、`abs_err_s11`、`db_err_s11`。
- `hfss_abs_s21`、`solver_abs_s21`、`abs_err_s21`、`db_err_s21`。

BP filter 的默认工程目标：

- direct FEM 对 HFSS：全频带 `max |abs_err| <= 2e-2`，且通带附近 `max |db_err| <= 0.2 dB`。
- AWE 单展开点：展开点必须与 HFSS/direct 同阶吻合；窄带内按 3.1 验收。
- GAWE/MGAWE 或自适应多点 AWE：旧 BP filter 的 40--43 GHz 全带按 3.2 验收。
- WCAWE：当前 IOStructure 的 90--100 GHz、`IOStructure_S_parameters.csv` 按 3.3 验收。

## 1. 数学单测

### 1.1 AWE 矩递推

构造小型多项式矩阵系统：

`A(t) = A0 + t A1 + t^2 A2`，
`b(t) = b0 + t b1`。

要求：

- `PolynomialMomentRecurrence` 生成的 `{x0, x1, ...}` 满足系数比较得到的递推式。
- `x_j` 与有限差分得到的 `d^j x/dt^j / j!` 一致。
- 二阶矩阵项和 RHS 导数项必须进入递推。

### 1.2 Padé 矩匹配

构造已知有理函数：

`H(t) = (1 + 0.2 t) / (1 - 0.3 t + 0.05 t^2)`。

要求：

- AWE 生成的 `[L/M]` Padé 与 Taylor 矩匹配到 `L+M` 阶。
- 误差阈值：`max |m_ref - m_fit| <= 1e-10`。
- 输出 Padé 分母根，并检查极点是否落入目标频带。

### 1.3 Galerkin 投影一致性

对随机复对称小矩阵测试：

- `K_r = V^T K V`、`M_r = V^T M V`、`m_{p,r}=V^T m_p` 维度正确。
- full matrix 复对称时，reduced matrix 保持同一投影约定下的复对称。
- reduced residual 在 `span(V)` 内正交：`V^T(A V x_r - b) ≈ 0`。

### 1.4 WCAWE 论文递推

构造 dA=2、db=2 的小型复多项式系统，直接验证论文式 (7)–(9)：

- 强制 U=I 时逐列恢复传统 AWE 矩，相对误差不高于 1e-12。
- 候选与基满足 Vtilde_n≈V_n*U_n，相对 Frobenius 误差不高于 1e-12。
- 论文式 (7) 的逐阶递推残差不高于 1e-11。
- Hermitian MGS 正交误差不高于 1e-11。
- reduced model 至少匹配前 q 个系统矩，相对误差不高于 1e-10。
- 必须包含一个能区分真正 WCAWE 与论文错误式 (6) 的反例。
- U 接近奇异时必须终止单点递推；不得跳过失败列继续生成下一阶。
- 基向量长度保持为原 FEM DOF，不允许通过线性化扩维。

X≈VR 只能作为普通 QR 诊断，不能替代式 (7)、式 (8) 和矩匹配验收。

## 2. 工程烟测

通用旧 BP filter 烟测先生成 direct 内部基准，并用 `wg_bp_filter_S_parameters.csv` 做 HFSS 外部基准对比；WCAWE 当前工程改用第 3.3 节的 `IOStructure_S_parameters.csv`：

```powershell
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 5 --sweep direct --no-write-all-fields --out result
```

现有 ALPS 对照：

```powershell
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 5 --sweep alps --alps-order 12 --no-write-all-fields --out result
```

未来新增算法后的目标命令示例：

```powershell
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 5 --sweep awe --awe-order 8 --no-write-all-fields --out result
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 5 --sweep gawe --gawe-order 12 --no-write-all-fields --out result
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 5 --sweep mgawe --mgawe-points 3 --mgawe-order 20 --no-write-all-fields --out result
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 0 --max-sweep-points 5 --sweep wcawe --wcawe-order 12 --linear-solver direct --no-write-all-fields --out result
```

## 3. 精度验收

### 3.1 AWE 窄带

- 频带：围绕中心频率的小范围，例如 41.0--42.0 GHz。
- 展开点：频带中心。
- 阈值：相对 HFSS CSV 的 `max delta |S11|, |S21| <= 0.01 dB`；同时抽检 direct FEM，确认误差不是由 FEM/HFSS 基准差异引入。
- 若扫到 40--43 GHz 全带，允许失败，但必须通过极点位置、残差或 direct 抽检说明超出单点有效带宽。

### 3.2 GAWE / MGAWE 宽带

- 频带：40--43 GHz。
- GAWE 初始展开点：中心点；MGAWE 初始展开点：左端、中心、右端。
- 阈值：相对 HFSS CSV 的 `max delta |S11|, |S21| <= 0.05 dB`。
- 展开点处必须接近 direct；频带内最大误差应随 GAWE 阶数、MGAWE 展开点数量或局部阶数增加而下降。
- 必须记录每个展开点贡献的保留列数和 deflation 列数。

### 3.3 WCAWE 良条件宽带

论文一致性首先按 1.4 和 [代码实现规范](wcawe-implementation-skill.md) 验收。当前大型工程基准固定为：

- 零阶 current.ngmesh，290856 个四面体、364932 个边未知量。
- 90–100 GHz，101 点，中心展开点 95 GHz。
- PARDISO Release，wcawe-order=12。
- HFSS 外部参考为 `IOStructure_S_parameters.csv`，只比较 S11/S21。
- 主指标使用幅值相对 L2；深零点额外报告绝对幅值误差。

最低门槛：

- 95 GHz 展开点相对 direct 的复数 S11/S21 误差不高于 1e-10。
- HFSS S11 幅值相对 L2 不高于 1.0%。
- HFSS S21 幅值相对 L2 不高于 1.2%。
- q=12 时 1 次数值分解、12 个 RHS、12 次依赖串行回代。
- 峰值内存不超过当前 GAWE/WCAWE 基线的 105%。

另做 q={12,20,30,40} 高阶压力测试。只有在更晚停滞、更低高阶误差或以更少有效列达到同精度至少一项成立时，才宣称 WCAWE 相对 GAWE/AWE 有优势。

2026-07-28 实测：q=12、numerical port 的 HFSS S11/S21 幅值相对 L2 为 0.774%/0.990%，95 GHz 复数 direct 相对误差为 7.9e-13/9.2e-14，峰值内存 2900.7 MB，1 次数值分解和 12 次回代。q=20/30/40 均达到目标阶数且正交误差小于 5e-14，但宽带误差未单调改善，因此当前算例只证明高阶递推稳定，不声明相对 GAWE/MGAWE 的精度优势。

## 4. 残差和物理检查

每个算法都必须输出：

- 代数残差：
  `rho(omega) = ||A(omega) x_hat(omega) - b(omega)|| / ||b(omega)||`。
- S 参数 HFSS 基准误差。
- S 参数 direct 抽检误差。
- 无源/无损双端口的被动性检查：`|S11|^2 + |S21|^2 <= 1 + 1e-3`。
- 端口传播常数检查：若目标频率靠近截止，提示拆分频带或回退 direct。
- reduced system 条件数或 LDL 主元诊断。

## 5. 失效诊断

必须把以下情况作为显式失败或警告写入日志：

- AWE 的 Padé Hankel 条件数过高，例如 `cond > 1e12`。
- Padé 极点落在目标频带内或离采样点过近。
- GAWE/MGAWE 全局基正交性恶化或 deflation 过多。
- WCAWE 上三角系数矩阵 `R` 接近奇异。
- WCAWE 未能改善传统 AWE 矩基条件数。
- ROM 小系统奇异或条件数异常。
- residual 高但 S 参数表面上正常。
- lossless 以外的材料进入尚未支持的快速扫频路径。

## 6. 回归记录

每次算法变更至少保存：

- `s_parameters.csv`：S 参数结果。
- `hfss_comparison.csv`：与 `wg_bp_filter_S_parameters.csv` 对齐后的外部基准误差。
- `diagnostics.json`：算法参数、展开点、阶数、ROM 维度、retained/deflated 列数、Padé pivot ratio、WCAWE 式 (7) 递推残差、式 (8) 关系残差、U 对角、正交性、矩匹配、无源性和 reduced solve 状态。
- `timing.json`：offline、basis、projection、online 的耗时。
- `basis_condition.csv`：WCAWE 基条件数诊断；必须包含 AWE 矩基条件代理、WCAWE 正交基条件代理、U 对角元、递推残差、基关系残差和正交性误差。
- `benchmark_summary.csv`：`scripts/compare_with_hfss.py --batch-root` 生成的 direct/ALPS/AWE/GAWE/MGAWE/WCAWE 批量 HFSS 对比汇总。
- `direct_comparison.csv`：direct 抽检频点上的幅度误差、相位误差和残差。

这些文件让 AWE、GAWE、MGAWE、WCAWE 即使将来抽离成独立库，也能继续复用本工程的 BP filter 基准。
