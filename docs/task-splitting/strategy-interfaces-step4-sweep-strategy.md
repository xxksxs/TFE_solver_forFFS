# Step 4：`ISweepStrategy` 抽象 + DirectSweep / AlpsSweep 迁移

> Phase 1 / Step 4 of `docs/optimization/strategy-interfaces/plan.md`。**依赖 Step 1**。可与 Step 5 并行。

## 目标

把 `Application::runApplication` 内 ~196–251 行的 direct loop 抽成 `sweep::DirectSweep`；`sweep::AlpsSweep` 实现 `ISweepStrategy`。`Application::runApplication` 主体缩到 ~50–80 行。旧 `bpfem/mor/AlpsSweep.hpp` 口径已过期。

新增任意扫频策略只需 1 个 hpp + 1 个 cpp + 在 `SweepStrategyFactory` 注册。

## 优先阅读文件

- `../optimization/strategy-interfaces/plan.md`（必读）
- `../../src/app/Application.cpp`（重点：196–280 行 direct loop、ALPS 调用段）
- `../../src/sweep/AlpsSweep.cpp`（ALPS-Galerkin MVP 实现）
- `../../include/bpfem/mor/AlpsSweep.hpp`
- `../optimization/alps-sweep/plan.md`
- Step 1 的成果（`ISparseSolver`）

## 可能修改文件

```
include/bpfem/sweep/ISweepStrategy.hpp        新增
include/bpfem/sweep/DirectSweep.hpp           新增
src/sweep/DirectSweep.cpp                     新增（接管 Application 内 direct loop）
include/bpfem/sweep/AlpsSweep.hpp             新增（重命名自 mor/AlpsSweep.hpp，加实现 ISweepStrategy）
src/sweep/AlpsSweep.cpp                       从 src/mor/ 移过来，buildOffline + evaluate 改成 run() 一次跑完
include/bpfem/factory/SweepStrategyFactory.hpp 新增
src/factory/SweepStrategyFactory.cpp          新增
src/app/Application.cpp                       重写 runApplication，从分发 → 调工厂 + 调 sweep.run()
include/bpfem/mor/AlpsSweep.hpp               留薄壳：#include "bpfem/sweep/AlpsSweep.hpp" + [[deprecated]]
src/sweep/AlpsSweep.cpp                       实现 buildOffline / evaluate / reconstructField
CMakeLists.txt                                src/sweep/*.cpp 加入构建；
                                              src/factory/SweepStrategyFactory.cpp 加
```

## 接口签名（与 plan.md §3.3 一致）

```cpp
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

}
```

`DirectSweep` 实现：每个频点 `ctx.assembler.assemble(...)` + `ctx.solver.solve(...)` + extractor.extract + 写 VTU；与现 Application 里那段逻辑等价。

`AlpsSweep` 实现：`run()` 内部先做现 `buildOffline()`，然后循环每个 freq 调 `evaluate()`。`buildOffline()` / `evaluate()` 保留为 public 给单测，但顶层入口是 `run()`。

## 验证命令

```powershell
cmake --build build_mkl --config Release

# direct 5 频点
.\build_mkl\Release\bp_fem_solver.exe --max-sweep-points 5 --basis-order 1 --no-write-all-fields --out result

# ALPS 101 频点 (TFE-5 multi-mode)
.\build_mkl\Release\bp_fem_solver.exe --max-sweep-points 101 --basis-order 1 --port-method tfe --tfe-modes-per-port 5 --sweep alps --alps-krylov-order 30 --no-write-all-fields --out result

# direct 单频点 (TFE-5)
.\build_mkl\Release\bp_fem_solver.exe --max-sweep-points 1 --basis-order 1 --port-method tfe --tfe-modes-per-port 5 --sweep direct --no-write-all-fields --out result
```

三个 CSV 与 main 分支对应输出 max\|Δ\|S\| < 1e-12（机器精度）。

`Application::runApplication` 行数缩到 ≤ 80 行，且不再有任何"if Sweep == Alps"分支。

## 不变量 / 风险提示

- **VTU 写出策略不变**：DirectSweep 内部仍按 `options.writeAllFields` 决定逐频点写 VTU；ALPS 路径仍不写 VTU（ROM 不持有解向量，这是已记录的 MVP 限制）。
- **`SweepResult::lastEdgeDofs`**：DirectSweep 填，AlpsSweep 留空。Application 拿到结果后判空决定是否写 `field_last.vtu`。
- **`buildFrequencies`** 仍在 Application 层调用，把 `std::vector<double>` 传给 `sweep.run()`。这是为了让"频点列表"这件事独立于策略 —— 未来 AdaptiveSweep 可以无视输入列表自己决定频点（接口允许返回不同长度的 `points`）。
- **离线求解器注入**：AlpsSweep::buildOffline 接收 `ISparseSolver&` 和 `SolverConfig`。PARDISO 后端仍可复用 pattern 缓存，但 ALPS 不再硬编码 `MklPardisoSolver`。
- **mor 别名期**：保留 `bpfem/mor/AlpsSweep.hpp` 一个 release 周期：
  ```cpp
  // include/bpfem/mor/AlpsSweep.hpp
  #pragma once
  #include "bpfem/sweep/AlpsSweep.hpp"
  namespace fem::mor {
      using AlpsSweep [[deprecated("use fem::sweep::AlpsSweep")]] = fem::sweep::AlpsSweep;
      using AlpsOptions [[deprecated("use fem::sweep::AlpsOptions")]] = fem::sweep::AlpsOptions;
  }
  ```

## 估时

1 day。

## 输出物

- 1 个 PR
- 三个对照 CSV 与 main 字节一致
- `Application::runApplication` ≤ 80 行（实际行数附在 PR 描述）
- Step 5 的文档可以紧接更新
