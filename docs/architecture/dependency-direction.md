# 依赖方向

本文说明模块间允许的依赖关系，目标是避免循环依赖和职责混乱。Phase 1
接口化之后，依赖图从"Application 直调具体类"变为"Application 调工厂 → 工厂
返回接口"。

## 当前依赖结构

```text
app
 |-- io
 |-- fem
 |-- bc
 |-- linalg          (仅 ISparseSolver 接口；具体后端通过 factory 拿到)
 |-- sweep           (仅 ISweepStrategy 接口；具体策略通过 factory 拿到)
 |-- factory
 |-- apm
 |-- tfe
 |-- post

factory
 |-- bc              (Application 注册 BCs 时；factory 本身不实例化 BCs，
 |                   只实例化 solver / sweep；BCs 由 Application 自己 push_back)
 |-- linalg
 |-- sweep
 |-- app             (只为读 Options 字段；不调用 runApplication)

bc
 |-- core
 |-- fem             (AssemblyContext 引用 PortModeSolver / EdgeTopology)
 |-- linalg          (SparseMatrix + SparsePatternBuilder)

sweep
 |-- core
 |-- fem             (AlpsSweep 直接调 assembler.buildAffineSystem)
 |-- linalg          (ISparseSolver；ALPS offline 通过接口注入求解器)
 |-- post            (DirectSweep 调 ResultExtractor)

linalg
 |-- core            (Types / Constants)
 (仅此；不依赖 fem / app)

fem
 |-- core
 |-- linalg
 |-- bc              (FEMAssembler 持有 IBoundaryCondition 列表)

mor (deprecated)
 |-- sweep           (转发别名)

apm
 |-- core
 |-- fem

tfe
 |-- core
 |-- fem
 |-- apm             (复用 RectPortGeometry 和几何探测器)

post
 |-- core
 |-- fem

io
 |-- core            (仅此)
```

## 设计原则

- **`core` 是最底层**：共享类型、常量、日志、工具，被几乎所有其它模块使用。
- **`io` 只填充数据结构**：不做 FEM 计算，不依赖 `fem`。
- **`linalg` 不知 `fem` / `app`**：求解器 / 预条件不知道 H(curl)，不知道 CLI。
  它们看到的就是 `SparseMatrix * x = b`。
- **`bc` 不知 `app`**：边界条件不依赖 CLI；只依赖 fem / linalg 的对外接口。
- **`sweep` 是 fem + linalg + post 的编排层**：策略类不知道 CLI，但需要装配
  和求解能力。
- **`factory` 是 Options 到接口的路由**：是唯一允许 include `app/Application.hpp`
  来读 Options 字段的非 app 模块。它可以 include `linalg/*Backend.hpp`
  和 `sweep/*Sweep.hpp` 等具体类型来 `make_unique`。
- **`app` 装配整个流水线**：唯一允许同时 include `factory/*` + `bc/*` +
  `post/*` + `io/*` 的模块。不直接依赖具体 `linalg` / `sweep` 实现。

## 加新模块时的判断

1. 它属于输入、装配、求解、扫频、还是输出？
2. 它能不能写成接口 + 实现的形式？如果能，把接口放在自己的目录，让
   factory 实例化具体实现。
3. 不让底层模块 include 上层模块的头。例如 `linalg` 不应该 include `fem`。
4. 出现循环时，把共享类型下沉到 `core`；把策略选择上推到 `factory`。

## 具体策略加新条目时的最小改动集

- **新增稀疏求解后端**（例如 GMRES、SuperLU、AMD/HIP 路径）：
  1. `include/bpfem/linalg/<X>Backend.hpp` + `src/linalg/<X>Backend.cpp` 实现
     `ISparseSolver`。
  2. `include/bpfem/app/Application.hpp` 的 `enum LinearSolverKind` 加一项。
  3. `src/factory/SparseSolverFactory.cpp` 的 `switch` 加一个 case。
  4. `src/app/Application.cpp::parseOptions` 的 `--linear-solver` 字符串映射
     加一项。
  - **不**触碰 `BiCGStabSolver.cpp` / `MklPardisoSolver.cpp` / `Application` 主流程 / `FEMAssembler`。

- **新增扫频策略**（例如 AAA 有理插值、自适应频点）：
  1. `include/bpfem/sweep/<X>Sweep.hpp` + `src/sweep/<X>Sweep.cpp` 实现
     `ISweepStrategy`。
  2. `enum SweepStrategy` + `SweepStrategyFactory` switch + `--sweep` 字符串映射。
  - **不**触碰 `DirectSweep.cpp` / `AlpsSweep.cpp` / `Application` 主流程。

- **新增边界条件**（例如真正的 ABC、ImpedanceBC 实现）：
  1. `include/bpfem/bc/<X>BC.hpp` + `src/bc/<X>BC.cpp` 实现 `IBoundaryCondition`。
  2. 在 `Application::runApplication` 注册阶段 `bcs.push_back(make_shared<X>())`。
  - 注册策略 (面 ID 来自哪里、什么参数) 是 Application 的职责，不由 factory 处理；
    边界条件本身没有"几个备选项二选一"的需求，每个 BC 都可以独立 push 到列表里。
  - **不**触碰 `FEMAssembler.cpp` / 其它 BCs。

- **新增预条件**（例如 SSOR、ILU(0)、AMS）：
  1. `include/bpfem/linalg/Precon<X>.hpp` + `src/linalg/Precon<X>.cpp` 实现
     `IPreconditioner`。
  2. `enum PreconditionerKind` + `factory/SparseSolverFactory.cpp::makePreconditioner`
     switch + `--precon` 字符串映射。
  - **不**触碰 `BiCGStabSolver.cpp` 或其它预条件。
