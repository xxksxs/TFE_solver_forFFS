# Step 1：引入 `ISparseSolver` 抽象

> Phase 1 / Step 1 of `docs/optimization/strategy-interfaces/plan.md`。最小 PR、零行为变更。后续步骤依赖此 PR。

## 目标

把"Application.cpp 内 `#ifdef BPFEM_USE_MKL` 编译期选 PARDISO 否则 BiCGSTAB"换成"运行时通过 `ISparseSolver*` 多态分发"。结果：

- 增加 `--linear-solver direct|bicgstab` CLI 选项
- 默认行为字节级与 main 分支一致（有 MKL 用 PARDISO，无 MKL 用 BiCGSTAB）
- 加新求解器后端只需新建 2 个文件（Backend.hpp / Backend.cpp）+ 在 factory 注册

## 优先阅读文件

- `../optimization/strategy-interfaces/plan.md`（必读，完整上下文）
- `../../include/bpfem/linalg/MklPardisoSolver.hpp`
- `../../include/bpfem/linalg/BiCGStabSolver.hpp`
- `../../src/app/Application.cpp`（看现 `#ifdef BPFEM_USE_MKL` 分支位置，约第 220–260 行）
- `../../src/mor/AlpsSweep.cpp`（注意：内部直接 `new MklPardisoSolver`，本步骤**不**改它，留到 Step 4）

## 可能修改文件

```
include/bpfem/linalg/ISparseSolver.hpp        新增
include/bpfem/linalg/PardisoBackend.hpp       新增（包装 MklPardisoSolver）
include/bpfem/linalg/BiCGStabBackend.hpp      新增（包装 BiCGStabSolver）
src/linalg/PardisoBackend.cpp                 新增
src/linalg/BiCGStabBackend.cpp                新增
include/bpfem/factory/SparseSolverFactory.hpp 新增
src/factory/SparseSolverFactory.cpp           新增
include/bpfem/app/Application.hpp             Options 加 LinearSolverKind 枚举（默认 Direct）
src/app/Application.cpp                       solver 局部变量改成 unique_ptr<ISparseSolver>
                                              parseOptions 加 --linear-solver
CMakeLists.txt                                src/linalg/PardisoBackend.cpp
                                              src/linalg/BiCGStabBackend.cpp
                                              src/factory/SparseSolverFactory.cpp
```

## 接口签名（与 plan.md §3.2 一致）

```cpp
namespace fem::linalg {

struct SolverConfig {
    int maxIterations = 400;
    double tolerance = 1.0e-7;
};

class ISparseSolver {
public:
    virtual ~ISparseSolver() = default;
    virtual SolveResult solve(const SparseMatrix& A,
                              const std::vector<std::complex<double>>& b,
                              const SolverConfig& cfg = {}) = 0;
    virtual void rememberPatternForReuse(bool /*enable*/) {}
};

}
```

`PardisoBackend` 用组合包装现 `MklPardisoSolver`（不要继承，避免动旧类签名）。`BiCGStabBackend` 同理。

## 验证命令

```powershell
cmake --build build_mkl --config Release
.\build_mkl\Release\bp_fem_solver.exe --max-sweep-points 3 --basis-order 1 --no-write-all-fields --out results_step1
.\build_mkl\Release\bp_fem_solver.exe --max-sweep-points 3 --basis-order 1 --no-write-all-fields --out results_step1_bicg --linear-solver bicgstab
```

第一条 CSV 应与 main 分支 `--max-sweep-points 3` 输出字节一致。第二条与 main 分支 `BPFEM_USE_MKL=OFF` 编译输出字节一致。

## 不变量 / 风险提示

- **`MklPardisoSolver` 的 symbolic 缓存**：现有代码中扫频每个频点稀疏 pattern 不变（同一 `assembler.assemble` 返回同一 `SparsityPattern`），PARDISO 内部跳过 reorder。`PardisoBackend` 必须**持有同一个 `MklPardisoSolver` 实例跨 solve 调用**，不能每次 solve 都 new 一个。`rememberPatternForReuse(true)` 的语义就是它。
- **Options 二进制兼容**：在 `Options` 末尾追加新字段，不要插中间，避免下游构造代码 misalign。
- **CLI 兼容**：不传 `--linear-solver` 时行为与 main 完全一致。
- **mor::AlpsSweep 不动**：它内部 `new MklPardisoSolver()` 留到 Step 4 一起重构。本 PR 别碰。
- **#ifdef 不能进 ISparseSolver.hpp 公共头**：`SparseSolverFactory.cpp` 是唯一允许出现 `#ifdef BPFEM_USE_MKL` 的地方。

## 估时

0.5 day。

## 输出物

- 1 个 PR
- main 分支 CSV diff 全 0
- `--linear-solver bicgstab --tolerance 1e-9 --max-iterations 800` 在 BP filter 单频点上能收敛（max iter 内 residual < tol）
