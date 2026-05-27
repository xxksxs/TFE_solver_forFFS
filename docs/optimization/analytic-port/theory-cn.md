# 超限元端口建模 (Transfinite Element, TFE)

## 1. 背景与动机

数值端口模（Numerical Port Mode, NPM）在端口面网格上求解二维 H(curl) 本征值问题：

```
∇_t × (μ_r⁻¹ ∇_t × E_t) - k_c² E_t = 0     在 Γ_p
n × E_t = 0                                 在 ∂Γ_p
```

这一路径由 `PortModeSolver` 实现，对**任意**端口截面均适用，但要做一次小规模 LAPACK 广义本征求解 + 后续的 mass 投影。在 BP filter 基准上：每端口约 1100 dofs，对 LAPACKE_dsygv 一次开销 ~1 s，共两端口 ~2 s。

**TFE 的核心思想**：当端口截面是规则几何（矩形 / 圆 / 同轴）时，端口模可写成解析形式，直接代入投影即可，**无需求本征值**。这在大规模问题（端口面 dof 数 > 10⁴）或扫频 + h 自适应（端口拓扑频繁变化）场景下能省下显著时间。

本 MVP 仅实现**矩形波导 + TE10 主模**，因为 BP filter 例题正是这种几何，物理可与 NPM 直接对照验证。

## 2. 数学模型（矩形 TE10）

设矩形端口截面在端口本地坐标 (u, v) 下满足

```
0 ≤ u ≤ a,    0 ≤ v ≤ b,    a ≥ b
```

其中 `a` 是宽度（长边），`b` 是高度（短边）。TE10 模解析形式：

```
e_t(u, v) = √(2 / (a·b)) · sin(π u / a) · v̂
```

其中 `v̂` 是沿短边的单位向量。该归一化保证

```
∫_Γ |e_t|² dS = 1
```

与 `PortModeSolver` 的 L2 模式形状归一化一致。

截止波数为

```
k_c² = (π / a)²
```

传播常数沿用工程约定 `β = √(k₀² - k_c²)`。

## 3. 端口几何探测

每个端口面的 (`origin`, `axisU`, `axisV`, `normal`, `a`, `b`) 由 `TransfinitePortBuilder::detectRectangular(faceId)` 自动估计：

1. 收集 `faceId` 上所有 surface triangles。
2. 累积面元法向（外积加权）→ `normal`。
3. 顶点投影到端口平面，做平面内 PCA：协方差矩阵的两个特征向量给出主轴。
4. 沿主轴的 bounding box 给出 `a`（长）、`b`（短）。
5. `axisU` 取长边方向，`axisV` 取短边方向，使 TE10 的 E 沿 `v̂`。

这一探测对**严格矩形**的端口面是无误差的；对略偏的端口（曲面或几何变形）会产生 PCA 拟合残差，进而影响 TFE 端口模的形状。

## 4. 离散与 FE 投影

NPM 路径里端口模通过本征求解直接得到 FE 空间的 dof：`dofs^T M dofs = 1`，并把 `couplingWeights = M·dofs` 喂给 `FEMAssembler` / `ResultExtractor`。TFE 路径要让外部模块感知不到差异，必须输出**相同形状的 PortMode 结构体**。

具体步骤（与 NPM 严格对应）：

1. **DOF 收集**：与 `PortModeSolver` 完全相同——枚举端口三角形的边 / 一阶面 DOF，按 outline 边（仅被一个端口三角形使用 → PEC 外轮廓）排除。
2. **装配端口面 mass M**：在 Wandzura 7 点求积下装配 `M_ij = ∫ N_i · N_j dS`。这一步是 SPD，规模 ~1000，Cholesky O(n³) 可控（每端口一次）。
3. **解析 RHS**：`b_i = ∫ N_i · e_t^analytic dS`，同一份求积。
4. **解 mass 系统**：`dofs = M⁻¹ b`（端口模在 FE 空间的最佳 L2 投影系数）。
5. **L2 重归一化**：

   ```
   α² = dofs^T M dofs   （Bessel 不等式：≤ 1，仅当 e_t 在 FE 空间内时取等）
   dofs ← dofs / α
   b    ← b / α
   ```

   这一步保证 FE 重构的端口模 `e_t^FE = Σ dofs_i N_i` 满足 `∫|e_t^FE|² dS = 1`。
6. **输出**：

   - `PortMode::couplingWeights = b`（即 `M·dofs`，与 NPM 接口完全相同）。
   - `PortMode::edgeDofs = dofs`（仅供 logging，下游不读其值）。
   - `PortMode::quadrature` 用 `e_t^FE`（不是 `e_t^analytic`），以保证 `poyntingPowerIntegral` 与 `applyWavePorts` / `ResultExtractor` 用同一份模式形状。

第 5 步的 L2 重归一**至关重要**：`b = ∫ N · e_t^analytic dS` 编码的是 `e_t^analytic` 在 FE 空间的投影，但 `dofs^T M dofs ≤ 1`（投影系数的 L2 在 FE 空间度量下小于 `e_t^analytic` 的 L2）。不归一会让 S 参数缩放偏差，曾在第一次实现中出现 `S21 = +0.13 dB` 的非物理结果。

## 5. 数据契约

```cpp
PortMode = {
    int    faceId;
    double cutoffWavenumberSquared;          // = (pi/a)^2
    vector<pair<int, double>> edgeDofs;       // dofs (L2 normalized in FE space)
    vector<pair<int, double>> couplingWeights;// = M * dofs
    vector<PortQuadraturePoint> quadrature;   // (e_t^FE, normal, area*w_qp)
};
```

`FEMAssembler::applyWavePorts`、`ResultExtractor::portProjection`、`poyntingPowerIntegral` 三处都不感知端口模式的来源——TFE 把数据塞进 `PortModeSolver` 的 cache 后，整个下游管线照常运行。

## 6. 与 NPM 的物理等价性

理论上：当端口截面**严格矩形**且网格在端口面上能精确表达解析 TE10 时，TFE 与 NPM 给出**相同**的 `couplingWeights`（仅相差一个全局相位/符号，由 mass-Cholesky 与本征求解的归一化约定决定）。

实际：

- 网格离散（边数有限）让 NPM 求得的本征模在 FE 空间内是 `e_t^analytic` 的 L2 最佳近似（带离散误差）。
- TFE 直接给出 `e_t^analytic` 的同种 L2 最佳近似。
- 两者在严格矩形几何下等价；非矩形几何下 NPM 跟踪实际几何，TFE 受 PCA 拟合限制。

实测对照（BP filter，5 个频点）：

| 频点 | NPM S21_dB | TFE S21_dB | 偏差 |
|---|---|---|---|
| 40.0 GHz | -10.910 | -10.900 | 0.010 |
| 40.75 GHz | -0.000160 | -0.000325 | < 0.001 |
| 41.5 GHz | -0.0692 | -0.0718 | 0.003 |
| 42.25 GHz | -12.636 | -12.653 | 0.017 |
| 43.0 GHz | -29.32 | -29.34 | 0.02 |

S11 在 deep nulls（如 40.75 GHz S11 ≈ -44 dB）上 TFE 与 NPM 偏差可达 3 dB，这是端口模"形状"在反射极深点的极端敏感性，TFE 的 PCA 矩形拟合会在这里放大几何误差。

## 7. 性能

- **NPM**：每端口约 1 s（LAPACK 广义本征 + mass 投影），两端口 ~2 s 一次。
- **TFE**：端口几何探测 < 1 ms，mass 装配 + Cholesky ~50 ms，总耗时与 NPM 相同量级。MVP 阶段未见明显加速；但当 NPM 求本征因为高阶 hierarchical 让端口面 dof 增长到 10⁴ 时，TFE 的 Cholesky 仍是 O(n³) 但常数更小（无本征求解），并且天然支持多模、多端口快速重建（仅改 `(a, b, mode index)`）。

主要价值不在 BP filter 这种小算例，而在以下场景：

- **大端口面**：端口面网格几何加密 / 高阶基函数让 NPM 端口面 dof 数到 5×10⁴ 以上，本征求解 LAPACK 单线程 O(n³) 开销超过 100 s。
- **多模端口**：NPM 取前 N 个本征模；TFE 直接生成 TE_mn / TM_mn 各阶解析模。
- **几何参数化扫描**：端口尺寸 (a, b) 改变只需重新评估解析公式；NPM 必须重新装配 + 求本征。

## 8. 适用范围与限制

| 场景 | TFE | NPM |
|------|-----|-----|
| 矩形波导主模 | ✓（MVP） | ✓ |
| 矩形波导高阶模 | 计划中 | ✓ |
| 圆波导 / 同轴线 | 计划中 | ✓ |
| 任意截面 | ✗ | ✓ |
| Floquet 周期端口 | ✗ | 计划中 |
| 端口几何明显非矩形 | 拟合误差大 | 受网格影响 |
| 工程容差 0.05 dB 通带验证 | ✓ | ✓ |
| S11 < -40 dB deep nulls | 偏差 ~ 1 ~ 3 dB | 准确 |

## 9. CLI 与代码入口

```
--port-method numerical|tfe       # 默认 numerical（保持与历史行为一致）
```

代码：

| 角色 | 文件 |
|------|------|
| 公共类型 | `include/bpfem/tfe/TransfinitePortBuilder.hpp` |
| 实现 | `src/tfe/TransfinitePortBuilder.cpp` |
| 注入接口 | `PortModeSolver::setPrecomputed(faceId, PortMode)` |
| 应用层分流 | `src/app/Application.cpp`（按 `Options::portMethod` 调用 TFE，再让 `PortModeSolver` 走 cache-first） |

## 10. 待办

| 任务 | 优先级 | 估算 |
|------|--------|------|
| 多模 TE_mn / TM_mn（矩形） | 高 | M |
| 圆波导 TE_nm / TM_nm（Bessel 函数） | 中 | M |
| 同轴线 TEM | 中 | S |
| 端口几何 user override（直接给 a, b, axes） | 中 | S |
| 与 ALPS 协同：在端口拟合参数（a, b）扫描中复用 ROM | 低 | L |
| 验证基准入 CI（与 NPM 自动对比） | 高 | S |

## 11. 与 ALPS 的关系

ALPS 模块（`fem::mor::AlpsSweep`）依赖 `FEMAssembler::buildAffineSystem` 给出 `(K, M, m_p, k_c²)`。TFE 替换的是 `m_p` 和 `k_c²` 的来源，**不影响** ALPS 的下游逻辑。`--port-method tfe --sweep alps` 同时启用是合法组合。
