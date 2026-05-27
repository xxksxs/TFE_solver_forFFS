# 物理原语

本文说明工程中最常见的电磁物理概念。修改端口、材料或频率相关代码前应先阅读本文。

## 频率

- **内部单位**：Hz。
- **日志显示**：可以转换为 GHz，仅用于显示。
- **主要来源**：AEDT 工程中的 driven-modal sweep。
- **主要消费者**：`FEMAssembler`、`PortModeSolver`、`ResultExtractor`。

频率会决定自由空间波数：

```text
k0 = 2 * pi * f / c0
```

## 波数与传播常数

- **自由空间波数 `k0`**：由频率和光速计算。
- **截止波数**：由端口本征模计算得到。
- **传播常数 `beta`**：通常由 `sqrt(k0^2 - kc^2)` 得到。

如果频率低于端口截止频率，传播常数会变为虚数或不应被当作正常传播模式处理。

## 材料

材料由 AEDT 文件解析得到，核心属性包括：

- **相对介电常数**：影响电场质量项。
- **相对磁导率**：影响磁场 curl-curl 项。
- **电导率**：引入损耗项。

材料数据应集中保存在共享类型中，避免在 FEM 装配中硬编码材料参数。

## 波端口

波端口是由 AEDT port 定义和 NGMesh faceId 共同确定的边界。

波端口包含：

- **faceId**：用于匹配网格表面三角形。
- **激励幅度**：通常来自端口功率或等效幅度设置。
- **激励相位**：以角度或弧度转换后参与复数源项。
- **是否激励**：非激励端口仍应作为吸收边界参与装配。

## 端口模式与功率归一化

端口数值本征模 `e_t` 在 `PortModeSolver` 中以 L2 范数 `∫|e_t|² dS = 1` 归一化。这只确定了模式的"形状"，与"P 瓦特入射"之间还差一个频率相关的幅度因子。

端口模式可以来自三条路径：

- **NPM (`--port-method numerical`，默认)**：`PortModeSolver::solve(faceId)` 在端口面三角网上求解二维 H(curl) 广义本征值问题 `K_port v = k_c^2 M_port v`（LAPACKE_dsygv，无 MKL 时回落 Jacobi），取最低非伪本征对作为端口主模。适用于任意截面，是默认路径。
- **APM (`--port-method analytic`)**：`bpfem::apm::AnalyticPortBuilder` 用解析 TE10（矩形）形式直接生成端口模，再投影到 FE 端口 dof + L2 归一。仅适用于规则矩形截面，作为闭式参考路径保留。
- **TFE (`--port-method tfe`)**：`bpfem::tfe::TransfiniteElementBuilder` 在端口面上做同样的二维 H(curl) 本征解，但保留 `--tfe-modes-per-port N` 个最低非伪本征对。每个模通过 `setMultiMode` 注入 `PortModeSolver::multiCache_`，`WavePortBC::apply` 对多模并入秩-1 端口算子。

三条路径输出同一个 `PortMode` 数据契约（见 `fem.md`）；下游模块完全等价对待。

波端口的功率归一化通过显式的坡印廷面积分给出：

```text
P_unit(f) = ½ ∫ Re{ E_t × H_t* } · n̂ dS
```

对 TE/TEM 主模，磁场由电场重构：

```text
H_t = (β / ω μ₀) · (n̂ × E_t),   β = √(k₀² − k_c²)
```

一瓦特入射对应的模幅度因子记为

```text
s(f) = 1 / √P_unit(f)
```

由 `powerNormalizationFactor(mode, f)` 提供。`PortMode::quadrature` 缓存了端口三角形上的求积点（位置场值 `E_t`、法向 `n̂`、权重 `area * w_qp`），供 `poyntingPowerIntegral` 数值累加。求积规则为 7 点 Wandzura（degree 5），`w_qp` 是按 `Σ w_qp = 1` 归一的重心坐标权重；零阶路径下任何 ≥ degree 2 的规则都够用，一阶路径下必须 ≥ degree 4。

约定：

- **截止端口** (β² ≤ 0)：`s(f) = 0`，不可激励，S 参数置零。
- **入射 RHS**：FEMAssembler 用 `√P_W · s(f)` 作为入射模幅度，使域内电场量纲与 `magnitudeW` 一致。
- **S 参数投影**：ResultExtractor 把每个端口的投影除以该端口的 `s_p(f)`，得到 √W 量纲的波幅度 `b`，再做 S 参数比值。两端口截面不同时 `s` 不同，公式仍成立。

## S 参数

S 参数由求解后的边自由度投影到端口模式得到。

- **`S11`**：输入端口反射。
- **`S21`**：输出端口传输。
- **dB 值**：由复幅度模值换算。

S 参数依赖端口模式归一化（L2 模式形状 + 坡印廷功率归一化）、边方向和投影定义。不要只修改 `ResultExtractor` 而不检查 `FEMAssembler` 中的端口装配逻辑，也不要只改幅度因子而不同步 `PortModeSolver::computeMode` 末尾的求积点缓存。
