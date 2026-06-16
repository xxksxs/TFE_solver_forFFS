# 策略接口化 + 工厂分发（架构演进 Phase 1）

> 状态：草案 / PR 计划。等同于 issue 描述 + task split + design review。
> 目标 commit 颗粒度：每个 Step 一个独立 PR，互不阻塞。

## 0. Why

当前代码目录上已经按 `core / io / fem / linalg / mor / post / apm / tfe` 分组，但**运行时分发**全部写死在 `Application::runApplication`：

- `if (options.portMethod == PortMethod::Analytic) {...} else {...}`
- `#ifdef BPFEM_USE_MKL` 选 PARDISO，否则选 BiCGSTAB（编译期，不可运行时切换）
- `if (options.sweepStrategy == Alps) {...} else {direct loop 写死在 Application}`
- 后处理循环耦合在求解循环里（每个频点立刻 `extractor.extract(...)`、立刻 `OutputWriter::writeVTU(...)`）

**症状**：

| 想加什么 | 现在要改的文件 |
|---|---|
| 阻抗边界 (Z_s) | `FEMAssembler.cpp` 加分支 + `Types.hpp` 加配置 + `AEDTParser` 解析 + `Application` 分发，4 个文件 |
| ILU(0) 预条件 BiCGSTAB | `BiCGStabSolver.cpp` 改签名 + `Application.cpp` 改调用，至少 2 个文件 |
| AAA 有理插值扫频 | 新建文件 + `Application.cpp` 加 enum 分支，2 个文件 |
| 后处理重新计算 group delay | 必须重跑全部频点（解向量没存盘） |

我们不打算解决最后一项（那是 Phase 2，独立后处理 exe），但把前 3 项的扩展成本从"改 N 个文件 + 阅读 Application.cpp 上下文"压到"加 1 个 .cpp + 注册到工厂"。

## 1. Non-goals

- **不**把 `bp_fem_core` 拆成多个静态库 / DLL。单一静态库 + 抽象接口已经足够。多 lib 引入 ABI / `Types.hpp` 共享 / DLL 边界 / `/MD` vs `/MT` 一致性等额外负担，等真正出现第二个 binary 消费者再做。
- **不**改 CLI 用户接口。所有现有 `--port-method analytic|tfe` / `--sweep direct|alps` / `--alps-krylov-order N` 行为字节级保留。
- **不**引入新外部依赖。只在 C++ 标准库范围内做接口化（`std::unique_ptr` / `std::function`）。

## 2. 目标架构

```
include/bpfem/
├── core/                      不变
├── io/                        不变
├── fem/                       EdgeTopology, FEMAssembler 保留；FEMAssembler 内部不再
│                              直接 know 'wave port'：调用 IBoundaryCondition::apply()
│   └── PortModeSolver.hpp     不变
├── bc/                        新建 ── 边界条件抽象层
│   ├── IBoundaryCondition.hpp ── apply(matrix, rhs, freq, ctx)
│   ├── WavePortBC.hpp         ── 现 applyWavePorts 抽出（含 NPM/APM/TFE 多模分支）
│   ├── PecBC.hpp              ── PEC zero-Dirichlet（现 imposeZeroDirichlet 的封装）
│   ├── AbsorbingBC.hpp        ── 1st-order ABC（占位空实现 + 单测）
│   ├── ImpedanceBC.hpp        ── Z_s 边界（占位空实现 + 单测）
│   └── FiniteConductorBC.hpp  ── Leontovich SIBC（占位空实现 + 单测）
├── linalg/                    SparseMatrix 不变
│   ├── ISparseSolver.hpp      ── solve(A, b) → SolveResult
│   ├── PardisoBackend.hpp     ── 现 MklPardisoSolver 实现 ISparseSolver
│   ├── BiCGStabBackend.hpp    ── 现 BiCGStabSolver 实现 ISparseSolver
│   ├── IPreconditioner.hpp    ── apply(M^{-1}, x) → y（占位 + Jacobi 实现）
│   └── PreconJacobi.hpp       ── 第一个具体实现，用于 BiCGStabBackend
├── sweep/                     新建 ── 扫频策略抽象层
│   ├── ISweepStrategy.hpp     ── run(freqs, ctx) → vector<SParameterPoint>
│   ├── DirectSweep.hpp        ── 现 Application 内 direct loop 抽出
│   ├── AlpsSweep.hpp          ── sweep::AlpsSweep + ISweepStrategy
│   └── AdaptiveSweep.hpp      ── 占位（频点自适应、AAA 插值 等留空）
├── factory/                   新建
│   ├── BoundaryConditionFactory.hpp
│   ├── SparseSolverFactory.hpp
│   └── SweepStrategyFactory.hpp
├── mor/                       逐步弃用 ── 等 Step 4 完成后移除（保 alpsdir/AlpsSweep 别名一周）
└── post/                      不变（Phase 1 不动后处理）

src/                           镜像 include 目录
src/bc/                        新增
src/sweep/                     新增（含 DirectSweep.cpp 与从 mor/ 平移过来的 AlpsSweep.cpp）
src/factory/                   新增
src/app/Application.cpp        瘦身，从分发器变成"读 Options → 调工厂 → 运行 sweep"
```

## 3. 关键接口草案

### 3.1 `IBoundaryCondition`

```cpp
// include/bpfem/bc/IBoundaryCondition.hpp
namespace fem::bc {

struct AssemblyContext {
    double frequencyHz;
    const Mesh& mesh;
    const ProjectDefinition& project;
    const EdgeTopology& topology;
    const PortModeSolver& portModeSolver;  // 仅 WavePortBC 用
    const std::unordered_set<int>& constrainedDofs;
};

class IBoundaryCondition {
public:
    virtual ~IBoundaryCondition() = default;

    // 装配阶段：把这条边界对系统矩阵和 RHS 的贡献加上。
    virtual void apply(SparseMatrix& matrix,
                       std::vector<std::complex<double>>& rhs,
                       const AssemblyContext& ctx) const = 0;

    // 稀疏 pattern 预声明：边界产生的非零位置（哪些 (row, col) 对会被 apply 写入）。
    // 装配前调用一次，让 FEMAssembler 把 pattern 合入全局 SparsityPattern。
    virtual void declareSparsity(SparsePatternBuilder& builder,
                                 const AssemblyContext& ctx) const = 0;

    // 仅 lossless ALPS 路径需要：返回此边界对应的 affine 数据
    // (rank-1 端口对：m_p、k_{c,p}^2、是否激励)。其他边界返回空。
    // 这是为了让 ISweepStrategy::Alps 不必硬编码 "wave port" 概念。
    virtual std::vector<AffinePortContribution> affineContributions(
        const AssemblyContext& ctx) const { return {}; }
};

}  // namespace fem::bc
```

`AffinePortContribution` 字段直接复用现有 `FEMAssembler::AffineSystem` 内的并行数组语义：

```cpp
struct AffinePortContribution {
    std::vector<std::pair<int, double>> coupling;  // m_p
    double cutoffSquared;                          // k_{c,p}^2
    int faceId;
    int projectPortIndex;
    bool isExcitationMode;
};
```

### 3.2 `ISparseSolver`

```cpp
// include/bpfem/linalg/ISparseSolver.hpp
namespace fem::linalg {

struct SolverConfig {
    int maxIterations = 400;
    double tolerance = 1.0e-7;
    // 直接求解器忽略以上字段；迭代求解器使用。
};

class ISparseSolver {
public:
    virtual ~ISparseSolver() = default;
    virtual SolveResult solve(const SparseMatrix& A,
                              const std::vector<std::complex<double>>& b,
                              const SolverConfig& cfg = {}) = 0;

    // 现 Pardiso 已经支持的 symbolic 复用：当稀疏 pattern 不变时，
    // 后续 solve 跳过 reorder/analysis。直接求解器实现；迭代器返回 noop。
    virtual void rememberPatternForReuse(bool enable) { (void)enable; }
};

}  // namespace fem::linalg
```

### 3.3 `ISweepStrategy`

```cpp
// include/bpfem/sweep/ISweepStrategy.hpp
namespace fem::sweep {

struct SweepContext {
    const ProjectDefinition& project;
    const FEMAssembler& assembler;
    const PortModeSolver& portModeSolver;
    linalg::ISparseSolver& solver;
    Logger& log;
};

struct SweepResult {
    std::vector<SParameterPoint> points;
    std::vector<std::complex<double>> lastEdgeDofs;  // 最后频点解；direct 路径填，ALPS 路径空
    double lastFrequencyHz = 0.0;
};

class ISweepStrategy {
public:
    virtual ~ISweepStrategy() = default;
    virtual SweepResult run(const std::vector<double>& frequencies,
                            const SweepContext& ctx) = 0;
};

}  // namespace fem::sweep
```

### 3.4 工厂

```cpp
// include/bpfem/factory/SparseSolverFactory.hpp
namespace fem::factory {
std::unique_ptr<linalg::ISparseSolver> makeSparseSolver(const Options& opts);
}
```

工厂内部按 `Options` 字段路由：

```cpp
// src/factory/SparseSolverFactory.cpp
std::unique_ptr<linalg::ISparseSolver> makeSparseSolver(const Options& opts) {
#ifdef BPFEM_USE_MKL
    if (opts.linearSolver == LinearSolverKind::Direct) {
        return std::make_unique<linalg::PardisoBackend>();
    }
#endif
    return std::make_unique<linalg::BiCGStabBackend>();
}
```

`LinearSolverKind { Direct, BiCGStab }` 是 `Options` 新增字段，CLI 加 `--linear-solver direct|bicgstab` 但默认仍为"有 MKL 用 PARDISO，无 MKL 用 BiCGSTAB"，行为兼容。

## 4. PR 拆分（每条 = 一个 PR）

### Step 1 — 引入 `ISparseSolver` 抽象（最小 PR，零行为变更）

**目标**：`MklPardisoSolver` / `BiCGStabSolver` 都实现 `ISparseSolver`，`Application.cpp` 用接口指针调用。

**文件改动**：

```
include/bpfem/linalg/ISparseSolver.hpp        新增
include/bpfem/linalg/PardisoBackend.hpp       新增（包装现 MklPardisoSolver，实现 ISparseSolver）
include/bpfem/linalg/BiCGStabBackend.hpp      新增（包装现 BiCGStabSolver）
include/bpfem/factory/SparseSolverFactory.hpp 新增
src/factory/SparseSolverFactory.cpp           新增
src/app/Application.cpp                       solver 局部变量改成 unique_ptr<ISparseSolver>
include/bpfem/app/Application.hpp             Options 加 LinearSolverKind（默认 Direct）
CMakeLists.txt                                src/factory/SparseSolverFactory.cpp 加入
```

**不改动**：
- 现有 `MklPardisoSolver` / `BiCGStabSolver` 类签名保持，让 `PardisoBackend` 用组合而不是继承（避免动旧代码）
- AlpsSweep 离线阶段通过 `ISparseSolver` 注入求解器，不再内部直接 new `MklPardisoSolver`

**验证**：
```powershell
.\build_mkl\Release\bp_fem_solver.exe --max-sweep-points 3 --basis-order 1 --no-write-all-fields
```
S 参数与 main 分支字节级一致（diff CSV）。

**估时**：0.5 day

---

### Step 2 — 引入 `IBoundaryCondition` 抽象，搬迁 wave port

**目标**：`FEMAssembler::applyWavePorts` 从成员函数变成"遍历注册的 `IBoundaryCondition` 列表，逐条调 `.apply`"。第一条具体边界 = `WavePortBC`，行为与现有完全一致（包括 NPM/APM/TFE 多模分支）。

**文件改动**：

```
include/bpfem/bc/IBoundaryCondition.hpp       新增
include/bpfem/bc/WavePortBC.hpp               新增（封装现 applyWavePorts）
src/bc/WavePortBC.cpp                         新增
include/bpfem/bc/PecBC.hpp                    新增（封装现 imposeZeroDirichlet 调用）
src/bc/PecBC.cpp                              新增
include/bpfem/fem/FEMAssembler.hpp            构造函数加 std::vector<unique_ptr<IBoundaryCondition>>
                                              applyWavePorts 标 [[deprecated]] 但保留实现
src/fem/FEMAssembler.cpp                      assemble() 改成遍历 bcs_
src/app/Application.cpp                       构造时把 PEC + WavePortBC 注册进 assembler
                                              （根据 project.ports / collectPecConstrainedDofs 推断）
CMakeLists.txt                                src/bc/*.cpp 加入
```

**关键细节**：
- `WavePortBC` 内部不要再 own 一个 `PortModeSolver`，从 `AssemblyContext` 取，避免双缓存。
- `affineContributions()` 实现也搬到 `WavePortBC`，让 `FEMAssembler::buildAffineSystem` 不再写死 "wave port"。

**验证**：现有 BP filter 全套对照（5 频点 ALPS-direct 等价 1e-11、TFE-1 vs APM 1e-3）必须不退化。

**估时**：1 day（核心是把 `applyWavePorts` 内部逻辑无脑迁过去 + `sparsityPattern` 里的端口 column 声明也走 declareSparsity）

---

### Step 3 — `IPreconditioner` + Jacobi / ILU(0) 占位（已落地，结论为负）

**目标**：让 `BiCGStabBackend` / `GmresBackend` 接受可选预条件，提供 `PreconJacobi`（对角）和 `PreconIlu0`（零填充 ILU）两个具体实现。

**文件改动**：

```
include/bpfem/linalg/IPreconditioner.hpp      已落地
include/bpfem/linalg/PreconJacobi.hpp         已落地
include/bpfem/linalg/PreconIlu0.hpp           已落地
include/bpfem/linalg/GmresSolver.hpp          已落地（restarted GMRES + Givens + 右预条件）
include/bpfem/linalg/GmresBackend.hpp         已落地（ISparseSolver 包装）
src/linalg/PreconJacobi.cpp                   已落地
src/linalg/PreconIlu0.cpp                     已落地
src/linalg/GmresSolver.cpp                    已落地
src/linalg/GmresBackend.cpp                   已落地
src/linalg/BiCGStabSolver.cpp                 加 setPreconditioner / 右预条件 BiCGSTAB
src/factory/SparseSolverFactory.cpp           解析 --precon none|jacobi|ilu0；--linear-solver gmres
include/bpfem/app/Application.hpp             Options 加 PreconditionerKind / gmresRestart
```

**验证（实测，BP filter @ 41.5 GHz, 325k DOF）**：

| 配置 | iter | residual | S 参数有效 |
|---|---|---|---|
| BiCGSTAB + Jacobi | 4000 | 0.91 | 否 |
| GMRES(50) + ILU(0) | 1500 | 0.98 | 否 |
| PARDISO 直接解 | 1 | ~1e-15 | ✓ |

**结论**：Jacobi / ILU(0) / 对角归一化（数学上等价于左 Jacobi）**对当前装配的频域 H(curl) indefinite 系统 `A = K - k₀² M + jβ m mᵀ` 都不收敛**。这是 H(curl) 频域 Maxwell 的领域共识：indefinite 谱跨原点，常规预条件不能消除原点附近的奇点。

**真正能跑通的预条件**：AMS (Hiptmair-Xu) / Schwarz-DDM / shifted Laplacian + multigrid，工程成本 ~2000+ 行且需要离散梯度算子等额外结构，当前未实现。

**框架价值仍然成立**：`IPreconditioner` 抽象 + `GmresSolver` + `PreconIlu0` / `PreconJacobi` 作为后续 AMS 等 H(curl) 专用预条件的 plug-in 接口保留。下游加 AMS 只需新建 `PreconAms.hpp/cpp` + factory 一行注册，不动现有任何代码。

**生产推荐**：默认走 PARDISO（`--linear-solver auto` 或 `--linear-solver direct`）。

**估时**：0.5 day（框架）+ 已花约 0.5 day（GMRES + ILU(0) 实现 + 验证 + 文档）= 1 day 实际

---

### Step 4 — `ISweepStrategy` + DirectSweep + AlpsSweep 迁移

**目标**：把 `Application::runApplication` 内 196–251 行的 direct loop 抽成 `DirectSweep`；`sweep::AlpsSweep` 实现 `ISweepStrategy`。

**文件改动**：

```
include/bpfem/sweep/ISweepStrategy.hpp        新增
include/bpfem/sweep/DirectSweep.hpp           新增（接管 Application 内 direct loop）
src/sweep/DirectSweep.cpp                     新增
include/bpfem/sweep/AlpsSweep.hpp             新增（重命名自 mor/AlpsSweep.hpp）
src/sweep/AlpsSweep.cpp                       从 src/mor/ 移过来，实现 ISweepStrategy
include/bpfem/factory/SweepStrategyFactory.hpp 新增
src/factory/SweepStrategyFactory.cpp          新增
src/app/Application.cpp                       Application 缩到 ~50 行：
                                                  - 解 Options
                                                  - 装配 io 和 fem
                                                  - 注册 BCs
                                                  - 工厂建 solver / sweep
                                                  - 调 sweep.run() → 写 CSV / VTU
include/bpfem/mor/AlpsSweep.hpp               留薄壳别名，标 [[deprecated]]，下个 release 删
CMakeLists.txt                                src/sweep/*.cpp 加入；ALPS 位于 src/sweep/AlpsSweep.cpp
```

**验证**：BP filter 5 频点 direct vs 101 频点 ALPS 同时跑，CSV 与 main 分支字节一致。

**估时**：1 day（搬代码 + Application 重写约 60 行，无逻辑变更）

---

### Step 5 — 文档与 examples

**目标**：把架构文档同步到 `docs/optimization/strategy-interfaces/`，加一个最小 plugin 示范。

**文件改动**：

```
docs/optimization/strategy-interfaces/plan.md       本文档（已存在）
docs/optimization/strategy-interfaces/api.md        IBoundaryCondition / ISparseSolver / ISweepStrategy 的字段意义
docs/optimization/strategy-interfaces/example.md    用 ImpedanceBC 占位实现演示如何加新边界
docs/architecture/strategy-interfaces.md            Application 层的指针图（替换现有"Application 直接调 Pardiso"那张图）
docs/task-splitting/boundary-condition-tasks.md     新增（拆 BC 类的开发任务模板）
docs/task-splitting/linear-solver-tasks.md          新增
docs/task-splitting/sweep-strategy-tasks.md         新增
README.md                                            "模块布局" 段同步
```

**估时**：0.5 day

---

## 5. 总成本

| Step | 工作量 | 阻塞下一步？ |
|---|---|---|
| 1. ISparseSolver | 0.5 d | 是（4 依赖） |
| 2. IBoundaryCondition | 1 d | 否 |
| 3. IPreconditioner | 0.5 d | 否 |
| 4. ISweepStrategy | 1 d | 是 |
| 5. Docs | 0.5 d | 否 |
| **总计** | **3.5 d** | |

Step 2 与 Step 3 互不阻塞，可并行；Step 4 必须 Step 1 之后做。

## 6. 兼容性矩阵

| 现有调用 | 改动后 | 用户感知 |
|---|---|---|
| `bp_fem_solver --port-method analytic` | 不变 | 无 |
| `bp_fem_solver --port-method tfe --tfe-modes-per-port 5` | 不变 | 无 |
| `bp_fem_solver --sweep direct` | 走 DirectSweep | 无 |
| `bp_fem_solver --sweep alps --alps-krylov-order 30` | 走 AlpsSweep | 无 |
| 新增 `--linear-solver bicgstab` | 强制 BiCGSTAB（即便有 MKL）| 新功能 |
| 新增 `--precon jacobi` | 启用预条件 | 新功能 |
| 老用户 default | PARDISO 直接、direct sweep、analytic 端口 | 与 main 字节一致 |

## 7. 风险与回退

- **风险**：Step 2 把 `applyWavePorts` 拆出去时，可能破坏 `sparsityPattern` 里端口列的声明顺序，导致 PARDISO symbolic 缓存命中失败、扫频整体变慢 5–10×。
  **回退**：保留 `IBoundaryCondition::declareSparsity` 接口，让 WavePortBC 在装配前先声明完所有端口列；和现版 `FEMAssembler::sparsityPattern` 完全等价。CI 加一项"50 freq 扫频时间不超过 baseline ×1.05"。

- **风险**：Step 4 把 ALPS 从 `mor/` 迁到 `sweep/` 是逻辑等价改名，但下游 Python / 第三方代码若直接 include `bpfem/mor/AlpsSweep.hpp` 会断。
  **回退**：保留 `bpfem/mor/AlpsSweep.hpp` 一个 release 周期，内容就是 `#include "bpfem/sweep/AlpsSweep.hpp"`；deprecation warning 提示迁移。

- **风险**：抽象层引入虚函数开销。
  **评估**：所有虚调用都在"每频点 1 次"或"每装配 1 次"尺度，相对一次 PARDISO 因式分解（O(n^{1.5}) 操作 = 几亿次 flop）开销是 1e-9，可忽略。已在 plan 中明示，避免后续"虚函数开销变性能问题"的迷信回滚。

## 8. 验收标准（Acceptance）

完成 Step 1–5 后：

- [ ] `Application::runApplication` 主体 ≤ 80 行，分发逻辑完全在 factory / sweep / bc 三类工厂内
- [ ] 加一个新的边界类型（即便是 1 行打印 placeholder）只需 *新建 2 个文件*（`Foo.hpp`、`Foo.cpp`）+ *改 1 个文件*（factory 注册），不需要动 `Application.cpp` / `FEMAssembler.cpp`
- [ ] 加一个新的迭代求解器 / 预条件类似上一条
- [ ] BP filter 5 频点 direct + 101 频点 ALPS 测试与 main 分支 S 参数 max\|Δ\| < 1e-10 dB
- [ ] 文档 5 件齐
- [ ] CI 跑一遍 baseline，无新增 warning（`/W4 /permissive-`）

## 9. Phase 2 预告（不在本 PR 范围）

完成 Phase 1 后，下一阶段的两件事：

1. **解向量序列化** + 独立 `bp_fem_post.exe`：用 HDF5 存 `(freq, edgeDofs, sparams, port_modes)`，post-processor 读取后做 group delay / passivity / 远场 / 多 VTU 输出，不再重跑求解器。
2. **可选**：把 `bp_fem_core` 拆成多个静态库（仅当 Phase 2 完成后真的有第二个 binary 消费者时）。
