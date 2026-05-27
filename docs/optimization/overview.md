# 总览与差距分析

本文给出当前 BP-FEM 频域求解器与商业级 3D 电磁 FEM 求解器（HFSS、CST MWS Frequency-domain、FEKO、COMSOL RF、EMPro 等）之间的差距，并提出优化策略。

## 当前能力（事实清单）

通过通读 `src/`、`include/bpfem/`、`docs/primitives/` 与 `docs/architecture/` 可以确认：

- **几何**：无 CAD 输入，依赖外部 NGMesh 离线生成的网格。
- **网格**：仅四面体一阶单元，背景体被剔除，无曲面/曲边、无网格生成器、无网格自适应。
- **基函数**：H(curl) Nedelec/Whitney 棱元；零阶（每四面体 6 dof）+ 一阶层次（每四面体 20 dof）。
- **物理**：标量、频率独立的 ε_r、μ_r、σ；curl-curl 体装配；时谐 e^{+iωt} 频域。
- **边界**：PEC（切向场零 Dirichlet）+ 数值波端口（端口本征模 + 坡印廷归一化）。
- **激励**：仅波端口主模激励，幅值 √P_W、相位由 AEDT 给出。
- **求解器**：`SparseMatrix` 复对称上三角 CSR；`MklPardisoSolver` 直接法（pt 状态复用 symbolic）；`BiCGStabSolver` 迭代 fallback（无预处理）。
- **扫频**：均匀线性扫频，逐频点重装-求解；MKL 路径在稀疏图不变时复用 symbolic。
- **后处理**：S11 / S21（CSV）、`field_last.vtu` 与可选逐频点 VTU；PEC 切向幅值诊断字段。
- **并行**：单线程；MKL PARDISO 内部并行依赖 oneAPI 运行时配置。
- **可用性**：仅 CLI；无脚本 API、无 GUI、无工程文件格式。
- **验证**：以 `wg_bp_filter` 单一波导带通滤波器为示例；无规范化基准库与回归 CI。

## 商业级求解器的最小特性集（基线对比）

下表是商业级 3D 电磁 FEM 频域求解器普遍具备的能力。"现状"列只反映本工程：

| 维度 | 商业级最小要求 | 本工程现状 |
|------|---------------|------------|
| 基函数阶数 | 至少 2 阶向量元，普遍支持 hp 自适应 | 0 阶 + 1 阶层次，无自适应 |
| 网格自适应 | 残差 / 目标量驱动的 h / p / hp 自适应 | 无 |
| 几何输入 | STEP / IGES / Parasolid / ACIS + 内置布尔 | 无，仅外部 ngmesh |
| 网格生成 | 内置或集成 (Netgen/Gmsh/OCCT) | 无 |
| 边界条件 | PEC / PMC / 阻抗 / SIBC / ABC / PML / Floquet 周期 / 对称面 / 辐射 | 仅 PEC + 数值波端口 |
| 激励 | 波端口（多模）/ 集总端口 / 电压源 / 入射平面波 / 电流源 / 偶极子 | 仅波端口主模 |
| 材料 | 各向异性张量 / Debye-Drude-Lorentz 多极色散 / 铁氧体 / 薄层 / 表面阻抗 / 温度相关 / 非线性 | 标量、频率独立、各向同性、无损扩展 |
| 直接求解 | 多前端 / 多波前 / 64-bit indexing / out-of-core | MKL PARDISO 单节点 |
| 迭代求解 | AMS / HX 预处理、AMG、ILU、GMRES / FGMRES / BiCGStab | BiCGSTAB 无预处理 |
| 区域分解 | DDM (FETI-DP / IETI / OS-DDM) | 无 |
| 扫频 | 离散 / 插值 / AWE / Krylov 多端口 MOR / 自适应频带 | 仅均匀离散 |
| 本征求解 | 谐振腔 / 周期结构 / 特征模 (CMA) | 无 |
| 并行 | OpenMP + MPI + GPU；混合 | 单线程，仅 PARDISO 内部并行 |
| 后处理 | S/Y/Z/Touchstone、远场方向图、增益、轴比、RCS、近场、Q、群时延、表面电流、传输线参数 | S11/S21 + VTU |
| 多物理 | 时域 (DGTD/TD-FEM) / 热-电 / 电路联仿 | 无 |
| 优化 | DOE / 参数扫描 / 梯度 + 伴随 / 拓扑优化 | 无 |
| API/UI | Python/Lua/COM API + GUI + 工程文件 | 仅 CLI |
| 验证 | 数百个内部基准 + 收敛性 + 长期 CI | 单一示例 |

## 差距优先级矩阵

横轴：用户感知（"做这件事用户立刻发现"），纵轴：实现成本。优先做"高感知 / 中等成本"的项目。

```
高感知
   │
   │  集总端口 / 平面波       Python API
   │  PML / ABC               GUI
   │  Touchstone + 远场       优化框架
   │  二阶/三阶 hierarchical   分布式 DDM
   │  AWE 快速扫频            GPU
   │  hp 自适应               热电耦合
   │  AMS / HX 预处理         瞬态 DGTD
   │  各向异性 / 色散
   │  内置网格 (Netgen/OCCT)
   │  OpenMP 装配
低感知
   └─────────────────────────────────────────
        低成本                   高成本
```

## 演进哲学

1. **先把"可信"做透**：Phase 1 的目标是任何一个用户跑出来的 S 参数能与 HFSS / CST 在 0.1 dB 以内重合，并自动诊断收敛性。这意味着先补**高阶基函数 + PML/ABC + 误差估计 + AMS 预处理 + Touchstone/远场/Y-Z 输出 + 基准库**。
2. **再把"通用"做开**：Phase 2 解决"够不到"的工程问题——集总端口 / 平面波 / 电压源、各向异性 / 色散、AWE 快速扫频、内置 Netgen 网格 + CAD 导入、OpenMP 装配。
3. **最后做"规模与产品化"**：Phase 3 / 4 加 DDM、MPI、GPU、瞬态、热-电耦合、Python API、GUI、优化器、CI。

## 不变量传播

任何优化必须沿如下不变量链一次性更新：

```
faceId 匹配
  → 边/面方向（按全局点号小到大）
    → 局部 DOF 族（0/1 阶 → 高阶）
      → 端口模式（L2 形状 + 坡印廷功率归一化）
        → 体装配 + 边界装配（对称上三角，单写）
          → 求解器 SolveResult 契约
            → 端口投影 + 后处理（用相同的 M_port * e_mode）
              → CSV / VTU / Touchstone / 远场输出
```

新增基函数阶数或边界类型时，必须把这条链上的 5 个模块（`EdgeTopology / PortModeSolver / FEMAssembler / ResultExtractor / OutputWriter`）一次性更新，并同步 `primitives/invariants.md` 与 `validation/`。

## 与本目录其他文档的关系

- 数值核心：`numerical-methods.md`、`boundary-and-excitations.md`、`materials.md`
- 求解层：`solvers-and-sweeps.md`、`parallelism-and-performance.md`
- 工程支撑：`geometry-and-mesh.md`、`usability.md`
- 输出与验证：`post-processing.md`、`verification-and-ci.md`
- 拓展：`multiphysics-and-transient.md`
- 计划：`roadmap-phases.md`
