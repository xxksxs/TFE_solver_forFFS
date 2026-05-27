# 模块布局

本工程采用"可复用静态库 + 轻量命令行入口"的结构。Phase 1 接口化 (见
`docs/optimization/strategy-interfaces/plan.md`) 之后，*策略选择* 全部走
工厂 + 抽象接口，新增策略只动新文件 + 工厂的一行注册，不再触碰 Application
和核心装配。

## 顶层目标

- **`bp_fem_core`**：单一静态库，包含解析、FEM 装配、边界条件、线性代数、
  扫频策略、降阶模型、后处理。
- **`bp_fem_solver`**：链接 `bp_fem_core` 的 CLI 程序。

## 目录结构

```text
include/bpfem/
  app/        命令行应用门面 (Options + runApplication)
  core/       常量、数学、共享类型、日志、工具
  io/         AEDT 与 NGMesh 输入前端
  fem/        Nedelec 边拓扑、curl-curl 装配、PortModeSolver (NPM)
  bc/         IBoundaryCondition 抽象 + WavePortBC + Absorbing/Impedance/FiniteConductor 占位
  linalg/     稀疏矩阵 + ISparseSolver / IPreconditioner 抽象 + 后端实现
  sweep/      ISweepStrategy 抽象 + DirectSweep + AlpsSweep
  factory/    SparseSolverFactory + SweepStrategyFactory（按 Options 路由）
  apm/        解析端口模 (Analytic Port Mode)
  tfe/        超限元端口模 (Transfinite Element)
  mor/        [DEPRECATED] AlpsSweep 旧路径，转发到 sweep/
  post/       S 参数提取 + VTU/CSV 输出

src/                       镜像 include 目录
src/main.cpp               轻量入口，调用 fem::runApplication
```

## 模块职责

- **`app`**：CLI 参数解析；端到端 setup + 调用工厂 + 调用 sweep。**不**做策略
  分发——策略选择全部委托给 `factory/*Factory`。
- **`core`**：共享常量、共享类型、日志、工具函数。无外部依赖。
- **`io`**：把 AEDT 和 NGMesh 文件转换为内部数据结构。**不**做 FEM 计算。
- **`fem`**：电磁离散、`PortModeSolver` (NPM 数值端口模) 和体积装配。
  `FEMAssembler::setBoundaryConditions(...)` 接收 `IBoundaryCondition` 列表，
  `assemble` / `sparsityPattern` / `buildAffineSystem` 通过列表分发，**不**硬编码
  wave port 概念。
- **`bc`**：边界条件抽象层。`WavePortBC` 是当前唯一的具体实现，封装现 NPM/APM/TFE
  的多模波端口装配。`AbsorbingBC` / `ImpedanceBC` / `FiniteConductorBC` 是已注册
  到构建系统的占位骨架，等待具体实现。
- **`linalg`**：稀疏矩阵存储 (`SparseMatrix`)、抽象求解器 (`ISparseSolver`) 和抽象
  预条件 (`IPreconditioner`)、具体后端 (`PardisoBackend` / `BiCGStabBackend` /
  `PreconJacobi`)。CMake 宏 `BPFEM_USE_MKL` 仅在 `factory/SparseSolverFactory.cpp`
  内部决定 Auto 路径选择，其它代码不再直接 `#ifdef`。
- **`sweep`**：扫频策略抽象。`DirectSweep` (per-freq assemble + solve)；
  `AlpsSweep` (block Krylov MOR)。`SweepContext` 携带一个 `onFieldSolved`
  回调让 Application 提供 VTU 写出钩子，避免 Sweep 与 OutputWriter 直接耦合。
- **`factory`**：根据 `Options` 实例化具体策略。**唯一**包含 `BPFEM_USE_MKL`
  分支的求解器路由代码就在这里。
- **`apm` / `tfe`**：端口模式建立路径。把得到的 `PortMode` / `MultiPortMode`
  注入 `PortModeSolver` 缓存；`WavePortBC::apply` 在装配时读取。
- **`mor`**：deprecated 别名层；保留一个 release 周期。新代码应直接 include
  `bpfem/sweep/AlpsSweep.hpp`。
- **`post`**：S 参数提取 (`ResultExtractor`) 和场输出 (`OutputWriter`)。

## 入口约束

`src/main.cpp` 应保持精简，只调用 `fem::runApplication(argc, argv)`。
`Application::runApplication` 主体本身已瘦身：策略分发占 < 10 行，剩余是
io / port-mode setup / 结果导出。
