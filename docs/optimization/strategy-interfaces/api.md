# Phase 1 接口契约

本文是 `IBoundaryCondition` / `ISparseSolver` / `IPreconditioner` /
`ISweepStrategy` 四个接口的契约说明：参数语义、调用顺序、前后置条件、
错误模型。和 `plan.md` 的关系：plan.md 是设计意图（why），本文是契约
（what & when）。

## `IBoundaryCondition`

定义于 `include/bpfem/bc/IBoundaryCondition.hpp`。

### 调用顺序

```text
注册期    Application 构造 BC 列表 (push_back into vector<shared_ptr<IBoundaryCondition>>)
        Application 调 FEMAssembler::setBoundaryConditions(std::move(bcs))

第一次 assemble() / sparsityPattern() / buildAffineSystem()：
  for each bc in bcs:
    bc.declareSparsity(builder, ctx)         <-- 装 SparsePattern
  ... volume curl-curl 装配 ...
  for each bc in bcs:
    bc.apply(matrix, rhs, ctx)               <-- 加 BC 贡献
  matrix.imposeZeroDirichlet(constrainedDofs, rhs)

ALPS 路径还会调：
  for each bc in bcs:
    auto pieces = bc.affineContributions(ctx)   <-- 收集 rank-1 端口贡献
```

### `apply(matrix, rhs, ctx)`

- **前置条件**：`matrix` 已经是带稀疏 pattern 的对象，体积装配已加入；
  `rhs` 已经分配为 `topology.edgeCount()` 大小、初值 0；`ctx.constrainedDofs`
  已 finalize。
- **职责**：把这一条边界对系统矩阵和 RHS 的贡献加上。**只能 add，不能
  remove**——zero-Dirichlet 投影由 `FEMAssembler` 在最后一步统一处理。
- **必须**：跳过 `ctx.constrainedDofs` 中的 `(row, col)` 对（任意一端在该
  集合都跳）。否则 `imposeZeroDirichlet` 会清掉，但代码失去显式约束安全网。
- **可以**：每次调用做不同的工作（频率相关），不需要幂等。
- **错误**：通过 `throw std::runtime_error` 报告无法装配的情况；返回值无意义。

### `declareSparsity(builder, ctx)`

- **前置条件**：`builder` 是 `SparsePatternBuilder` 的实例，包装当前 pattern
  的 `columnsByRow` 集合；`ctx.constrainedDofs` 已 finalize。
- **职责**：声明 `apply` 将要写入的所有 `(row, col)` 对。
- **必要性**：跳过此步会导致 `SparseMatrix::add` 在 `apply` 阶段触发"未声明
  位置"路径，该路径会破坏 PARDISO 的 symbolic-factorization 缓存复用，
  扫频整体变慢 5–10×。
- **顺序无关**：`builder.add(row, col)` 内部去重 + upper-tri 规范化；可以
  重复声明同一 pair，可以以任意顺序声明。
- **可以**为空实现：占位边界（如 `AbsorbingBC` 当前）的 `declareSparsity`
  返回不做事是合法的，只要其 `apply` 也不写新 `(row, col)`（或干脆 throw）。

### `affineContributions(ctx) -> vector<AffinePortContribution>`

- **前置条件**：与 `apply` 相同。
- **职责**：把 BC 装入的每个秩-1 算子拆成 `(coupling, k_c^2, faceId,
  projectPortIndex, isExcitationMode)`。仅当 BC 能写成
  `+j β m m^T` 形式时返回非空；ABC / SIBC 等不能拆分的返回空。
- **顺序约束**：返回顺序定义 ALPS 中的"虚拟端口索引"。多 BC 时所有 BC 的
  返回拼接成一个大列表。当前 `WavePortBC` 的内部顺序是 `(projectPortIndex
  ascending, modeIndex ascending)`；保持这个顺序对 ALPS 内 `dominantVirtualByProject`
  的正确性是关键。
- **每个项目端口至多一个 `isExcitationMode = true`**。多端口激励后续支持。

### 默认实现

基类 `affineContributions` 默认返回空 vector，所以非端口型 BC 不需要重写。

## `ISparseSolver`

定义于 `include/bpfem/linalg/ISparseSolver.hpp`。

### 调用顺序

```text
factory::makeSparseSolver(options) -> unique_ptr<ISparseSolver>
solver.rememberPatternForReuse(true)    [可选；默认 true 行为]
loop:
  solver.solve(A, b, cfg)
```

### `solve(A, b, cfg) -> SolveResult`

- **前置条件**：`A.size() == b.size()`；`A` 已 compress（被 `ISparseSolver`
  内部触发 `A.values()` 即 compress）。
- **必须返回**：
  - `field` 大小等于 `b.size()`。
  - `iterations`：迭代解为实际迭代次数；直接解返回 1 作为哨兵。
  - `residual`：相对残差 `||A x - b|| / max(||b||, 1e-30)`。
- **直接后端**忽略 `cfg.maxIterations` / `cfg.tolerance`。
- **迭代后端**遵守 `cfg`；`maxIterations` 用尽时返回当前残差，不抛错（让
  Application / Sweep 决定是否接受）。
- **错误**：因子分解失败、奇异矩阵等抛 `runtime_error`。

### `rememberPatternForReuse(bool)`

- 直接后端实现 reuse；迭代后端 noop。
- 默认 true，即扫频默认沿用 PARDISO symbolic 缓存。

### `name() -> const char*`

仅用于日志。

## `IPreconditioner`

定义于 `include/bpfem/linalg/IPreconditioner.hpp`。

### 调用顺序

```text
solver.setPreconditioner(make_shared<PreconJacobi>())   [可选]
loop:
  precon.setup(A)             <-- 在 BiCGSTAB::solve 内部调一次
  while iter:
    precon.apply(x, y)        <-- 每次需要 M^{-1} 时调
```

### `setup(A)`

- **前置条件**：`A` 必须在后续 apply() 调用期间保持有效。
- **职责**：从 `A` 提取计算 `M^{-1}` 所需的数据（对角线、ILU 因子等）。
- **频率**：当前实现每次 `BiCGStabSolver::solve` 调一次 setup，下一次 solve
  会再次 setup。这是为了简化生命周期管理；如果未来引入扫频间复用，需要
  在 IPreconditioner 上加 `clearCache()` / `dimension()`。

### `apply(x, y)`

- **前置条件**：setup() 已被调用过。
- **必须**：`y` 大小 = `x` 大小（实现自己 resize 即可）。
- **允许**：in-place (`&x == &y`)。
- **数学**：`y = M^{-1} x`。当前 Jacobi 实现的对角接近 0 的项把 `M^{-1}` 那
  位置存为 0（degenerate row 透传 x），避免 NaN。

## `ISweepStrategy`

定义于 `include/bpfem/sweep/ISweepStrategy.hpp`。

### 调用顺序

```text
factory::makeSweepStrategy(opts, project, assembler, portModeSolver)
  -> unique_ptr<ISweepStrategy>
SweepContext ctx { project, assembler, portModeSolver, extractor, solver, log };
ctx.linearMaxIterations = options.maxIterations;
ctx.linearTolerance     = options.tolerance;
ctx.onFieldSolved       = ... [optional, see below];
result = sweep.run(frequencies, ctx);
```

### `run(frequencies, ctx) -> SweepResult`

- **必须返回**：`points` 是 `SParameterPoint` 列表，按频率升序排序。
- **`lastEdgeDofs`**：直接式策略填最后频点的解；ALPS 留空（因为没有保留
  全空间基）。
- **`lastFrequencyHz`**：与 `lastEdgeDofs` 对应；空时填 0。
- **错误**：奇异系统、收敛失败抛 `runtime_error`，由 Application 顶层 catch。

### `SweepContext::onFieldSolved`

- 直接式策略每个频点解完调一次；ALPS 不调。
- Application 用它写 per-freq VTU。
- 默认空 closure；策略实现必须接受空的 sink（用 `if (ctx.onFieldSolved)` 守
  门）。
- 不要在 sink 内做长时间工作 — 它在求解循环内部，会拖慢扫频。

### `name() -> const char*`

仅用于日志（"direct" / "alps"）。

## `SparsePatternBuilder`

辅助类，定义于 `include/bpfem/bc/IBoundaryCondition.hpp`。

- `add(row, col)`：声明一个非零位置。`row > col` 自动 swap 到 upper-tri；
  index 越界或在 constrainedDofs 中静默跳过。
- 不暴露删除接口。pattern 阶段是只增不减。

## 错误模型一览

| 接口方法 | 错误处理 |
|---|---|
| `IBoundaryCondition::apply` | throw runtime_error |
| `IBoundaryCondition::declareSparsity` | throw runtime_error |
| `IBoundaryCondition::affineContributions` | 返回空 vector，不抛 |
| `ISparseSolver::solve` | 直接：分解失败 throw；迭代：超 max-iter 不抛 |
| `IPreconditioner::setup` | throw runtime_error |
| `IPreconditioner::apply` | 不抛（内部静默 fall back to identity） |
| `ISweepStrategy::run` | throw runtime_error |
| `SparseSolverFactory::makeSparseSolver` | 选 direct 但无 MKL 不抛构造，solve 时才抛 |
| `SweepStrategyFactory::makeSweepStrategy` | 不抛 |

## 默认行为不变性

下表是 main 分支字节级兼容的命令组合。Phase 1 完成后所有这些组合的 CSV
和 VTU 输出与 Phase 1 之前的 main 分支应等价（max\|Δ\| < 1e-12）。

| 命令 | 选用的 BC / 求解器 / 扫频 |
|---|---|
| `bp_fem_solver` (无任何选项) | WavePortBC + PARDISO + DirectSweep |
| `bp_fem_solver --port-method tfe` | + TFE-1 数值多模端口 |
| `bp_fem_solver --sweep alps` | + AlpsSweep |
| `bp_fem_solver --linear-solver bicgstab` | + BiCGStabBackend (no precon) |
| `bp_fem_solver --linear-solver bicgstab --precon jacobi` | + PreconJacobi |

`--precon` 默认 `none`，`--linear-solver` 默认 `auto`。任何未传新选项的
现有脚本应当与 Phase 1 之前完全一致。
