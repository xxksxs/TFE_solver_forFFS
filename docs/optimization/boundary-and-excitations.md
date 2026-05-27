# 边界条件与激励

本文给出边界条件与激励源的扩展计划。当前求解器仅支持 PEC + 数值波端口主模激励，距离商业级有较大缺口。

## 1. 现状

- **PEC**：`FEMAssembler::collectPecConstrainedDofs` 把非端口边界全部视为零切向 Dirichlet。
- **波端口**：`PortModeSolver` 求 2D 端口本征模 → L2 形状归一 + 坡印廷功率归一 → `FEMAssembler::applyWavePorts` 装配吸收项与 RHS。
- **激励**：仅波端口主模，幅值 √P_W，相位由 AEDT。

## 2. 边界条件扩展矩阵

| 边界 | 物理含义 | 实现路径 |
|------|---------|----------|
| **PMC** | 完美磁壁 | H(curl) 自然边界，等价于不施加任何约束；用作对称面。 |
| **阻抗 / SIBC** | 表面阻抗 Z_s 或 Leontovich 边界 | 在端口/边界面上加面积分项 `(jωμ/Z_s) ∫ (n×E)·(n×W) dS`。 |
| **薄层 SIBC** | 多层金属/介质薄层等效 | 等效张量阻抗（双侧 Leontovich）或 *Karlsson-Kristensson* 透射条件。 |
| **ABC (1st/2nd-order Sommerfeld)** | 一阶/二阶吸收边界 | 在边界面加 `jk0 ∫ (n×E)·(n×W) dS`（一阶）或加二阶曲率项。 |
| **PML (Perfectly Matched Layer)** | 各向异性吸收层 | 在外围网格区使用 *complex coordinate stretching* 或 *Berenger split* 张量材料；`materials.md` 中扩张量材料后即可装配。 |
| **Floquet / 周期** | 单胞 + 主从映射 + 相移 | `EdgeTopology` 加 *master-slave constraint*，主自由度 = 从自由度 × e^{-jk·d}。 |
| **对称面 (PEC/PMC)** | 1/2 / 1/4 / 1/8 模型 | 利用 PEC/PMC 自然条件实现，自动检测对称面 faceId。 |
| **辐射面 (天线问题)** | 远场吸收 | PML 优先；ABC 作为低阶 fallback。 |
| **完美匹配端口 (Wave port match)** | 端口本征模匹配吸收 | 已实现，保持。 |

## 3. 激励扩展矩阵

| 激励类型 | 物理 | 实现路径 |
|---------|------|----------|
| **多模波端口** | 端口可激励 / 接收前 N 阶 TE/TEM 模 | `PortMode` 改为 `std::vector<PortMode>`；S 参数变为多模矩阵。 |
| **集总端口 (Lumped port)** | 跨两个金属面注入电压源/电流源 | 在端口截面（通常一段二维条带）施加 `E·t = V/L` 的均匀切向源 + 内阻；S 参数按 `Z_0` 归一。 |
| **入射平面波 (Plane wave)** | 散射 / RCS / 透射 | 整体散射场公式 (SF/TF) 或 *total-field/scattered-field* 分区；引入入射场体源 `(jωμ_r σ + ω²μ_0(ε - ε_inc))·E_inc`。 |
| **电压 / 电流源** | 馈电点 / 偶极子 | 在某一边或一组边上指定 `∫ E·t dl = V` 或 `J = J_0`。 |
| **电流环 / 偶极子探针** | 天线馈电 | 集中体电流源 RHS。 |
| **近场 / 远场等价源** | 二次散射 | Huygens 表面 + 等效电/磁流。 |
| **Floquet 端口 (周期天线)** | 单元辐射 / TX-RX | Floquet 模式 + Floquet 本征展开。 |

## 4. 数据契约扩展（`include/bpfem/core/Types.hpp`）

```cpp
enum class BoundaryKind {
    PEC, PMC, Impedance, ThinLayer,
    AbcFirstOrder, AbcSecondOrder, Pml,
    Periodic, Symmetry, Radiation
};

struct BoundaryDefinition {
    int faceId = -1;
    BoundaryKind kind = BoundaryKind::PEC;
    std::complex<double> impedance{};        // for Impedance / SIBC
    double thickness = 0.0;                  // for ThinLayer
    Vec3 periodicTranslation{};              // for Periodic
    std::complex<double> phaseShift{};       // for Periodic / Floquet
};

enum class ExcitationKind {
    WavePortMode, LumpedPort, PlaneWave,
    VoltageSource, CurrentSource, FieldImport
};

struct ExcitationDefinition {
    int id = -1;
    ExcitationKind kind = ExcitationKind::WavePortMode;
    int faceId = -1;
    int modeIndex = 0;
    std::complex<double> magnitude;          // sqrt(W) or V or A
    Vec3 incidentDirection{};                // plane wave
    Vec3 incidentPolarization{};             // plane wave
    double referenceImpedance = 50.0;        // lumped port
};
```

`ProjectDefinition` 同步增加 `std::vector<BoundaryDefinition> boundaries` 与 `std::vector<ExcitationDefinition> excitations`。

## 5. 装配影响

- `FEMAssembler` 重构为 `assembleVolume` + `assembleBoundary(BoundaryDefinition)` + `assembleExcitation(ExcitationDefinition)`，按多态分发。
- PEC 约束统一通过 `imposeZeroDirichlet` 走批量路径（已有），避免和阻抗边界耦合。
- ABC / PML / SIBC 的边界项是面积分质量项，按对称矩阵单写规则装配（见 `primitives/invariants.md`）。
- 平面波散射场公式必须引入 `E_inc(r, f)` 体源 RHS；要求 `Mesh` 有体单元的材料对 `(ε_r, μ_r)` 与背景对比。

## 6. 多激励 / 多端口扫频

- 当前一次扫频只算一个激励；商业级要求"一次装配，多个 RHS 同时求解"。
- `MklPardisoSolver` 应支持 `nrhs > 1`；`SolveResult` 升级为 `std::vector<SolveResult>` 或 `MultiRhsSolveResult`。
- S 参数提取改为对每个端口/激励一次投影，得到 N×N 散射矩阵。

## 7. 验收基准

- **PMC**：矩形腔 PEC↔PMC 频率对换。
- **阻抗 / SIBC**：解析铜壁波导衰减常数。
- **ABC**：自由空间偶极子辐射阻抗误差 < 5%。
- **PML**：PML 反射 < −60 dB（在工作频段中心）。
- **周期**：FSS 单元周期阵列与无限阵列解析解对照。
- **集总端口**：50 Ω 同轴线 S11 < −50 dB。
- **平面波**：金属球 RCS 与 Mie 解析在 0.5 dB 内重合。

## 8. 任务清单

| 任务 | 优先级 | 估算 |
|------|--------|------|
| PMC / 对称面 | 高 | S |
| ABC（1 阶） | 高 | M |
| PML（UPML 张量） | 高 | L |
| 多模波端口 + N×N 散射矩阵 | 高 | M |
| 集总端口 + Z₀ 重归一 | 高 | M |
| 阻抗 / SIBC 边界 | 中 | M |
| 平面波 (SF/TF) 激励 + RCS | 中 | L |
| 周期 / Floquet 主从 | 中 | L |
| 薄层 SIBC | 低 | L |
| 多 RHS 一次求解 | 高 | M |
