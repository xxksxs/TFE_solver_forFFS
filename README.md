# BP-FEM 频域求解器

本目录包含一个基于 Visual Studio/CMake 的 C++17 频域 FEM 求解器原型，用于 `wg_bp_filter.aedt` 与 `current.ngmesh` 波导带通滤波器示例。

如需按推荐路径阅读工程和拆分开发任务，请从 `docs/overview.md` 开始；若需要更细的专题入口，请进入 `docs/overview/README.md`。

## 功能概览

- **AEDT 工程读取**
  - 材料定义
  - 设计变量
  - driven-modal 频率扫频
  - 波端口 faceId 和源激励
- **NGMesh 网格读取**
  - body 信息
  - 点坐标
  - 表面 facet
  - 四面体体网格单元
- **频域 FEM 求解**
  - 零阶 Nedelec/Whitney H(curl) 棱元
  - 可选一阶层次棱元，每个四面体 20 个局部基函数
  - 复数稀疏矩阵装配
  - Intel oneMKL PARDISO 直接求解器
  - BiCGSTAB fallback
- **结果输出**
  - `results/s_parameters.csv`
  - `results/field_last.vtu`
  - 可选逐频点 VTU 场文件

## 使用 Visual Studio 或 CMake 构建

可直接用 Visual Studio 的 CMake 集成打开本目录，也可以运行：

```powershell
cmake -S . -B build
cmake --build build --config Release
```

如果 Intel oneMKL 可用，CMake 默认启用 `BPFEM_USE_MKL` 并链接 PARDISO 稀疏直接求解器。若要强制使用迭代 fallback：

```powershell
cmake -S . -B build -DBPFEM_USE_MKL=OFF
```

## 模块布局

求解器组织为可复用静态库和轻量命令行入口。Phase 1 接口化（详见
`docs/optimization/strategy-interfaces/plan.md`）后，运行时分发全部走工厂 +
抽象接口；新增策略只动新文件 + 工厂的一行注册。

```text
include/bpfem/
  app/        命令行应用门面
  core/       常量、数学函数、共享数据类型、日志和工具
  fem/        Nedelec 边拓扑和 curl-curl FEM 装配（含 NPM）
  bc/         IBoundaryCondition 抽象 + WavePortBC + Absorbing/Impedance/FiniteConductor 占位
  linalg/     SparseMatrix + ISparseSolver / IPreconditioner 抽象 + 后端实现
  sweep/      ISweepStrategy 抽象 + DirectSweep + AlpsSweep + AweSweep + GaweSweep + MgaweSweep + WcaweSweep
  factory/    SparseSolverFactory + SweepStrategyFactory（按 Options 路由）
  apm/        解析端口模 (Analytic Port Mode)
  tfe/        超限元端口模 (Transfinite Element)
  mor/        [DEPRECATED] AlpsSweep 旧路径，转发到 sweep/
  io/         AEDT 与 NGMesh 输入前端
  post/       S 参数提取和 VTU/CSV 输出
src/                       镜像 include 目录
src/main.cpp               轻量入口
```

- **`bp_fem_core`**：可复用静态求解器库。
- **`bp_fem_solver`**：链接 `bp_fem_core` 的 CLI 程序。

## 运行

```powershell
.\build\Release\bp_fem_solver.exe --aedt wg_bp_filter.aedt --mesh current.ngmesh --out results
```

默认使用零阶棱元。若要启用一阶层次棱元：

```powershell
.\build\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 1 --out results_basis1_first
```

快速检查可运行：

```powershell
.\build\Release\bp_fem_solver.exe --max-sweep-points 3 --max-iterations 100
```

如需为每个频点写出 VTU 场文件：

```powershell
.\build\Release\bp_fem_solver.exe --write-all-fields
```

如需关闭逐频点 VTU（例如长扫频时减少磁盘占用）：

```powershell
.\build\Release\bp_fem_solver.exe --no-write-all-fields
```

VTU 抽样密度由 `--field-output-order 1|2|3` 控制：1 = 每四面体 4 顶点（与历史行为兼容）；2 = 4 顶点 + 6 边中点（一阶基函数默认）；3 = 再加 4 面心（用于核对一阶 face bubble）。该选项与 `--basis-order` 解耦，仅影响 VTU 抽样不影响 S 参数 CSV。

## 快速扫频（ALPS）

如需把扫频时间从 N 次直接求解降到一次离线 + N 次廉价投影：

```powershell
.\build\Release\bp_fem_solver.exe --basis-order 1 --sweep alps --alps-krylov-order 30 --no-write-all-fields
```

该路径在 BP filter 基准上 101 频点扫频从 ~22 min 降到 ~104 s（加速 ~12×），与 `--sweep direct` 在所有频点的 S 参数偏差 < 1e-10。当前为 MVP 版本（单展开点 + 复对称 Galerkin + 仅 lossless 材料）；详见 `docs/optimization/alps-sweep/plan.md`。

## 快速扫频（AWE / GAWE / MGAWE / WCAWE）

单点 Padé AWE：

```powershell
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 101 --sweep awe --awe-order 8 --no-write-all-fields --out results_awe
```

单点 Galerkin AWE（GAWE）：

```powershell
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 101 --sweep gawe --gawe-order 12 --no-write-all-fields --out results_gawe
```

多点 Galerkin AWE（MGAWE）：

```powershell
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 101 --sweep mgawe --mgawe-points 3 --mgawe-order 4 --no-write-all-fields --out results_mgawe
```

良条件 AWE（WCAWE）：

```powershell
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 101 --sweep wcawe --wcawe-order 12 --no-write-all-fields --out results_wcawe
```

- `--gawe-order`：单展开点生成并正交化的 AWE 矩向量数量，用于构造一个 Galerkin ROM。
- `--gawe-expansion`：GAWE 展开频率；不指定时默认取扫频区间中心。
- `--mgawe-points`：在扫频区间内自动均匀选择展开点，例如 3 点对应左端、中心、右端。
- `--mgawe-order`：每个展开点生成的局部 AWE 矩向量数量。
- MGAWE 会把所有局部矩向量正交化成一个统一 Galerkin ROM，不是多个局部模型拼接。
- `--wcawe-order`：单展开点 WCAWE 的目标良条件基维度。
- WCAWE 会输出 `basis_condition.csv`，记录传统 AWE 矩基与 WCAWE 正交基的条件曲线。
- 每次运行都会输出 `run.log`、`run.json` 和 `timing.json`，可用于和 direct / ALPS / AWE 做时间与峰值内存对比。

## 端口建模选择（NPM / APM / TFE）

CLI `--port-method` 控制端口模建立方式：

```powershell
.\build\Release\bp_fem_solver.exe --basis-order 1                                              # NPM (默认), 二维 H(curl) 数值本征主模
.\build\Release\bp_fem_solver.exe --basis-order 1 --port-method analytic                       # APM, 单 TE10 解析闭式
.\build\Release\bp_fem_solver.exe --basis-order 1 --port-method tfe --tfe-modes-per-port 5     # TFE, 数值多模
```

- `numerical`（默认）：在端口面三角网上求解 2D H(curl) 广义本征 `K_port v = k_c² M_port v`，取最低非伪本征对作为主模。任意截面（矩形 / 脊形 / 圆形 / 偏心）都可用。模块 `bpfem::fem::PortModeSolver`。
- `analytic`：闭式 TE10 + FE L2 投影；仅适用于 PCA 可拟合为矩形的端口。作为闭式参考路径保留。模块 `bpfem::apm`。详见 `docs/optimization/analytic-port/`。
- `tfe`：与 numerical 同一二维本征解，但保留 N 个最低物理本征对，每个模并入秩-1 端口算子。模块 `bpfem::tfe`。详见 `docs/optimization/transfinite-element/`。
- `--tfe-modes-per-port 1` 时 TFE 与 numerical 字节级等价。
- 与 `--sweep alps` 兼容（按 (face, mode) 展开为虚拟端口）。

## 求解器后端与预条件（CLI）

Phase 1 接口化之后线性求解器和预条件可以在运行时切换，不需要重新编译：

```powershell
# 默认（推荐）：MKL 构建走 PARDISO，否则 BiCGSTAB（与 main 分支字节级一致）
.\build\Release\bp_fem_solver.exe

# GMRES + ILU(0) （插件式占位，当前 H(curl) 频域系统上不收敛）
.\build\Release\bp_fem_solver.exe --linear-solver gmres --precon ilu0 --tolerance 1e-7 --max-iterations 1500 --gmres-restart 50

# BiCGSTAB + Jacobi（同样不收敛，框架占位）
.\build\Release\bp_fem_solver.exe --linear-solver bicgstab --precon jacobi
```

- `--linear-solver auto|direct|bicgstab|gmres`：默认 `auto`，按构建期 `BPFEM_USE_MKL`
  路由。`direct` 在无 MKL 构建上会 throw，`bicgstab` / `gmres` 始终可用。
- `--precon none|jacobi|ilu0`：默认 `none`。仅对迭代法生效；`direct` 时静默忽略。
- `--gmres-restart 30`：GMRES 重启长度，默认 30。

**重要**：`bicgstab` / `gmres` + 任何当前实现的预条件（`jacobi` / `ilu0`）在频域
H(curl) indefinite 系统（`A = K - k₀² M + jβ m mᵀ`）上**不收敛**。BP filter 实测：

| 配置 | 1500 iter 后残差 |
|---|---|
| GMRES(50) + ILU(0) | 0.98 |
| BiCGSTAB + Jacobi (4000 iter) | 0.91 |

这是 H(curl) 频域 Maxwell 的领域共识：indefinite 谱跨原点，常规预条件不能消除
原点附近的奇点。真正能用的预条件是 AMS (Hiptmair-Xu) / Schwarz-DDM /
shifted Laplacian + multigrid，工程成本 ~2000+ 行且需离散梯度算子等额外结构，
当前未实现。**生产请走 PARDISO**；GMRES / ILU(0) / Jacobi 框架作为后续 AMS 等
H(curl) 专用预条件接入的 plug-in 占位保留。

加新求解器后端 / 新预条件 / 新边界条件 / 新扫频策略的标准步骤见
`docs/optimization/strategy-interfaces/example.md`。

## 文档入口

- **工程阅读总览**：`docs/overview.md`
- **项目树说明**：`docs/project-tree.md`
- **工程总览专题**：`docs/overview/README.md`
- **架构说明**：`docs/architecture/README.md`
- **原语文档索引**：`docs/primitives/README.md`
- **测试说明**：`docs/testing/README.md`
- **验证说明**：`docs/validation/README.md`
- **任务拆分指南**：`docs/task-splitting/README.md`
- **开发流程**：`docs/development-workflow/README.md`

## 备注

数值内核默认使用零阶四面体 Whitney/Nedelec 棱元，并支持可选一阶层次棱元、全局 H(curl) 自由度、curl-curl 体装配、PEC 切向边界约束、数值波端口模态边投影、S 参数 CSV 输出和重构矢量场 VTU 输出。
