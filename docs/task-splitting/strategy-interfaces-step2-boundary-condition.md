# Step 2：引入 `IBoundaryCondition` 抽象，搬迁 wave port

> Phase 1 / Step 2 of `docs/optimization/strategy-interfaces/plan.md`。可与 Step 3 并行。Step 4 依赖此 PR。

## 目标

把 `FEMAssembler::applyWavePorts` / `imposeZeroDirichlet` 从硬编码成员函数变成"装配阶段遍历注册的 `IBoundaryCondition` 列表"。第一次发布两个具体实现：`WavePortBC`（封装现 NPM/APM/TFE 多模逻辑）+ `PecBC`（封装现 PEC 边）。新增空骨架 `AbsorbingBC` / `ImpedanceBC` / `FiniteConductorBC`，仅含构造函数和单测占位，不实现 apply。

新增任意边界类型只需：1) 新建 hpp/cpp 实现 `IBoundaryCondition`；2) 在 Application 的"build BCs"段调一次 `register`。**不**碰 `FEMAssembler.cpp` / `Application.cpp` 主体逻辑。

## 优先阅读文件

- `../optimization/strategy-interfaces/plan.md`（必读）
- `../primitives/fem.md`
- `../primitives/invariants.md`
- `../../src/fem/FEMAssembler.cpp`（重点：`applyWavePorts` 全部，约第 270–330 行；`sparsityPattern` 中端口列声明段）
- `../../src/post/ResultExtractor.cpp`（投影必须用同一个 `couplingWeights`）
- `../../include/bpfem/fem/PortModeSolver.hpp`（`MultiPortMode` / `setMultiMode` 接口）

## 可能修改文件

```
include/bpfem/bc/IBoundaryCondition.hpp       新增
include/bpfem/bc/WavePortBC.hpp               新增
src/bc/WavePortBC.cpp                         新增（迁出现 applyWavePorts 全部逻辑）
include/bpfem/bc/PecBC.hpp                    新增
src/bc/PecBC.cpp                              新增（封装 collectPecConstrainedDofs + imposeZeroDirichlet）
include/bpfem/bc/AbsorbingBC.hpp              新增（空 apply）
include/bpfem/bc/ImpedanceBC.hpp              新增（空 apply）
include/bpfem/bc/FiniteConductorBC.hpp        新增（空 apply）
src/bc/AbsorbingBC.cpp                        新增（throw "not implemented"）
src/bc/ImpedanceBC.cpp                        新增（throw "not implemented"）
src/bc/FiniteConductorBC.cpp                  新增（throw "not implemented"）

include/bpfem/fem/FEMAssembler.hpp            构造函数 + setBoundaryConditions(...)
src/fem/FEMAssembler.cpp                      assemble() / sparsityPattern() 改成走 bcs_
                                              applyWavePorts 标 [[deprecated]] 但保留兼容路径

src/app/Application.cpp                       构造时把 PEC + WavePortBC 注册进 assembler
CMakeLists.txt                                src/bc/*.cpp 加入
```

## 接口签名（与 plan.md §3.1 一致）

```cpp
namespace fem::bc {

struct AssemblyContext {
    double frequencyHz;
    const Mesh& mesh;
    const ProjectDefinition& project;
    const EdgeTopology& topology;
    const PortModeSolver& portModeSolver;
    const std::unordered_set<int>& constrainedDofs;
};

struct AffinePortContribution {
    std::vector<std::pair<int, double>> coupling;  // m_p
    double cutoffSquared;
    int faceId;
    int projectPortIndex;
    bool isExcitationMode;
};

class IBoundaryCondition {
public:
    virtual ~IBoundaryCondition() = default;
    virtual void apply(SparseMatrix& matrix,
                       std::vector<std::complex<double>>& rhs,
                       const AssemblyContext& ctx) const = 0;
    virtual void declareSparsity(SparsePatternBuilder& builder,
                                 const AssemblyContext& ctx) const = 0;
    virtual std::vector<AffinePortContribution> affineContributions(
        const AssemblyContext& ctx) const { return {}; }
};

}
```

`SparsePatternBuilder` 是新建的轻量 helper，包装 `std::vector<std::unordered_set<int>>& columnsByRow`，让 BC 不需要直接接触 `FEMAssembler` 的私有数据结构。

## 验证命令

```powershell
cmake --build build_mkl --config Release

# 5 频点 direct - 确认与 main 字节一致
.\build_mkl\Release\bp_fem_solver.exe --max-sweep-points 5 --basis-order 1 --no-write-all-fields --out results_step2_direct

# 101 频点 ALPS（TFE multi-mode 走 affine 路径）
.\build_mkl\Release\bp_fem_solver.exe --max-sweep-points 101 --basis-order 1 --port-method tfe --tfe-modes-per-port 5 --sweep alps --alps-krylov-order 30 --no-write-all-fields --out results_step2_alps

# APM 单频点对照
.\build_mkl\Release\bp_fem_solver.exe --max-sweep-points 1 --basis-order 1 --port-method analytic --no-write-all-fields --out results_step2_apm
```

三个 CSV 与 main 分支对应输出 max\|Δ\|S\| < 1e-12（机器精度）。

## 不变量 / 风险提示

- **`PortMode` 数据契约不变**：APM、TFE、NPM 三条路径仍通过 `PortMode { faceId, cutoffWavenumberSquared, edgeDofs, couplingWeights, quadrature }` 与下游交付。`couplingWeights = M_port · v` 的定义在 `PortModeSolver` 内做完，BC 只 *读取*。
- **稀疏 pattern 必须一致**：装配前 `FEMAssembler::sparsityPattern` 调用每个 BC 的 `declareSparsity()`，让端口列、约束行的非零位置在 pattern 阶段就已经声明完。否则首次 `assemble` 调用会 hit `SparseMatrix::add` 的"未声明 (row, col)"路径，PARDISO symbolic 缓存就废了，扫频会慢 5–10×。
- **affine 多模展开**：`buildAffineSystem` 现走 `WavePortBC::affineContributions`，按 (face, mode) 展开虚拟端口。返回顺序必须严格匹配现 `buildAffineSystem` 顺序（projectIdx 升序、模 idx 升序），否则 `AlpsSweep` 的 `dominantVirtualByProject[]` 索引会错位。
- **CLI 兼容性**：本 PR **不**新增任何 CLI 选项。`AbsorbingBC` / `ImpedanceBC` / `FiniteConductorBC` 只是占位骨架，实际启用逻辑在后续单独 PR 加。
- **`PecBC` 与 `cachedConstrainedDofs`**：现 `FEMAssembler::cachedConstrainedDofs` 是私有 mutable 缓存。迁到 `PecBC` 后，`PecBC` 内部持有这个缓存（首次 `apply` 时构造），`AssemblyContext` 字段 `constrainedDofs` 由 `PecBC::declareSparsity` 在 pattern 阶段提供给 `FEMAssembler`。

## 估时

1 day。核心是把 `applyWavePorts` 内部 NPM/APM/TFE 多模分支 + `sparsityPattern` 内端口列声明无脑迁过去；不要趁机重构其他东西。

## 输出物

- 1 个 PR
- 三个对照 CSV 与 main 分支字节一致（或 max\|Δ\| < 1e-12）
- `Application::runApplication` 中 wave-port 装配段从约 30 行降到 ~5 行（`assembler.setBoundaryConditions(...)`）
- 三个空 BC 骨架（编译通过，不被注册）
