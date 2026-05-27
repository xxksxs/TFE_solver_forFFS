# 商业级求解器优化路线（顶层入口）

本文是把当前 BP-FEM 频域求解器演进到商业级 FEM 电磁求解器（对标 HFSS / CST MWS / FEKO / COMSOL RF）所需的**完整优化路线图**的顶层入口。详细内容已按主题拆分到 `optimization/` 目录，本文件仅保留导航与摘要。

## 适用对象

- **项目负责人**：用于制定 2 ~ 4 期迭代计划与资源预算。
- **算法工程师**：用于按主题选取最近一期要落地的数值/求解器特性。
- **架构 / 平台工程师**：用于规划并行化、I/O、API、UI 边界。
- **测试 / 验证工程师**：用于规划基准库与回归框架。

## 推荐入口

- **总览与差距分析**：`optimization/overview.md`
- **数值方法升级**：`optimization/numerical-methods.md`
- **边界条件与激励**：`optimization/boundary-and-excitations.md`
- **材料模型扩展**：`optimization/materials.md`
- **线性求解器与扫频**：`optimization/solvers-and-sweeps.md`
- **并行化与性能**：`optimization/parallelism-and-performance.md`
- **几何与网格**：`optimization/geometry-and-mesh.md`
- **后处理与结果**：`optimization/post-processing.md`
- **多物理与瞬态**：`optimization/multiphysics-and-transient.md`
- **可用性 / API / UI**：`optimization/usability.md`
- **验证、测试、CI**：`optimization/verification-and-ci.md`
- **分期路线图**：`optimization/roadmap-phases.md`

## 当前能力与差距摘要

当前求解器已具备：AEDT/NGMesh 解析、零阶 + 一阶层次 Nedelec H(curl) 棱元、curl-curl 体装配、PEC 切向约束、数值波端口 + 坡印廷功率归一化、对称上三角 CSR + MKL PARDISO 直接求解、BiCGSTAB fallback、S11/S21 提取、VTU 场输出。

距离商业级求解器最关键的差距集中在十二个方向：

1. **基函数阶数受限**：仅支持零阶 + 一阶层次，缺二阶/三阶/可变阶 hp。
2. **边界条件单一**：只有 PEC + 波端口，缺 PMC、SIBC、ABC、PML、Floquet/周期、辐射边界、对称面、集总端口。
3. **激励单一**：只有波端口模式激励，缺集总端口、电压源、入射平面波、电流源、近场/远场等价源、线天线馈电。
4. **材料模型单一**：标量、频率独立，缺各向异性张量、Debye/Drude/Lorentz 多极色散、铁氧体/旋磁、薄层/表面粗糙、温度相关、非线性。
5. **网格自适应缺失**：无误差估计驱动的 h / p / hp 自适应，无目标量误差控制，无非协调网格。
6. **几何能力缺失**：无 CAD（STEP/IGES/Parasolid/ACIS）导入、无几何核（OCCT/Parasolid 集成）、无网格生成器（依赖外部 NGMesh）、无曲边/曲面单元。
7. **求解器手段单薄**：直接 PARDISO + 弱 BiCGSTAB，缺 AMS/HX/多重网格/AMG 等谱等价预处理、缺 DDM/FETI-DP 区域分解、缺特征模/谐振腔本征求解、缺 AWE/Krylov MOR 快速扫频。
8. **并行化几乎为零**：单线程装配、单节点求解，无 OpenMP/TBB、无 MPI、无 GPU。
9. **后处理薄弱**：只输出 S11/S21 + VTU，缺 Touchstone / Y / Z / 群时延 / 远场方向图 / 增益 / 轴比 / RCS / 端口去嵌 / 端口重归一化 / Q 值 / 谐振参数 / 表面电流。
10. **多物理与时域缺失**：无瞬态时域 (TD-FEM / DGTD)、无热-电耦合、无电路 (SPICE/Touchstone) 联仿、无优化与参数扫描框架。
11. **可用性差**：仅 CLI、无脚本 API（Python/Lua）、无 GUI/工程文件、无项目管理、无作业调度接口、无许可机制。
12. **验证体系不完整**：无规范化基准库（矩形腔、介质腔、矩形/圆波导、同轴线、贴片天线、SRR/EBG、Mie 散射、…）、无 CI、无收敛性回归。

**演进哲学**：先把"可信"做透（基函数阶数、边界、扫频、误差估计、并行装配），再把"通用"做开（CAD/网格/材料/天线后处理/瞬态），最后做"产品化"（GUI/API/分布式/优化/许可）。

## 分期路线图（详见 `optimization/roadmap-phases.md`）

| 阶段 | 周期 | 目标 | 关键交付 |
|------|------|------|---------|
| Phase 1 | 1 ~ 2 季 | 数值核心可信 | 二阶/三阶 hierarchical 基、PML/ABC、误差估计 + h-adapt、AMS 预处理、Touchstone/远场/Y-Z 输出、基准库 |
| Phase 2 | 2 ~ 3 季 | 工程通用性 | 集总端口 + 平面波 + 电压源、各向异性 + Debye/Drude/Lorentz、AWE/Krylov MOR 快速扫频、CAD 导入 + 内置 Netgen/OCCT 网格、OpenMP 装配 + 多线程 PARDISO |
| Phase 3 | 3 ~ 4 季 | 大规模并行 | DDM/FETI-DP、MPI、GPU 装配/求解、特征模求解器、热-电耦合、瞬态 DGTD 原型 |
| Phase 4 | 4 ~ 6 季 | 产品化 | Python 脚本 API、工程文件格式、GUI/结果浏览、优化框架、许可与诊断、商业级 CI |

## 跨模块影响摘要

任何一项优化都需要在这条不变量链上做一致性更新，详见 `primitives/invariants.md` 与 `optimization/overview.md`：

```
faceId 匹配 → 边方向 → 局部基函数族 → 端口模式 → 装配 → 约束 → 求解器接口 → 后处理投影 → 输出
```

若新增基函数阶数或边界条件，必须同步检查 `EdgeTopology / PortModeSolver / FEMAssembler / ResultExtractor / OutputWriter` 五处。

## 维护约定

- 在 `optimization/` 中新增主题文档时，应同步更新本文件的"推荐入口"与 `optimization/README.md`。
- 每项优化落地后，应在 `optimization/roadmap-phases.md` 标注完成状态，并在 `validation/` 中补充对应的回归基准。
- 新增数据契约前先更新 `include/bpfem/core/Types.hpp` 与 `primitives/invariants.md`。
