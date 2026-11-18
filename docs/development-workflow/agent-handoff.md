# Agent 快速接手指南

## 先做三件事

1. 在包含 `CMakeLists.txt`、`src/` 和 `docs/` 的 BP-FEM 仓库根目录运行 `git status --short --branch`、`git diff --stat`。核对当前分支和未提交改动，确认本次任务涉及哪些文件，保留其他任务的改动。
2. 阅读[修改规范](change-policy.md)和[变更记录](../changes/README.md)。从下表选读任务相关的模块和权威页。
3. 用当前源码构建，并跑相关测试和最小算例。已有二进制可能对应旧工作树。记录构建配置、后端、版本和实际结果。

| 任务 | 源码起点 | 接手时先读 |
|---|---|---|
| CLI、输入、输出 | `src/app/Application.cpp`、`src/io/`、`src/post/` | [工程参考](../engineering-reference.md) |
| FEM、端口和 TFE | `src/fem/FEMAssembler.cpp`、`src/fem/PortModeSolver.cpp`、`src/bc/WavePortBC.cpp`、`src/tfe/` | [TFE 专题](../optimization/transfinite-element/README.md)及[端口验证](../validation/port-validation.md) |
| AWE、MGAWE、WCAWE | `src/sweep/`、`src/fastsweep/`、`tests/fastsweep_tests.cpp` | [公式复现](../optimization/awe-family/paper-reproduction.md) |
| 线性求解和预条件 | `src/factory/SparseSolverFactory.cpp`、`src/linalg/`、`tests/sparse_solver_tests.cpp` | [求解器选择](../architecture/solver-selection.md) |
| 构建和验证 | `CMakeLists.txt`、`tests/` | [测试入口](../testing/README.md) |

## 工作树快照 · 2026-11-18

以下是接手线索。每次开工以实际的 `git status`、源码和变更记录为准。

- AWE 逐场 Padé、MGAWE 2001/自适应模式与高阶 WCAWE 已在工作树接入。WCAWE 的离线和在线共用一套多项式。做矩验证时使用 `--wcawe-check-moments`。做性能计时时关闭该验证并单独记其开销。
- 当前 CLI 的 `tfe` 路径是在端口网格上求多个数值本征模并装配对应波端口项。教科书中的模态自由度替换和 Schur 推导在 TFE 专题中作为理论与设计参考。
- 工作树含 MUMPS、Schwarz、MKL SpMV、PARDISO out-of-core 及稀疏求解器测试改动。先核对当前编译条件、依赖和测试结果。早期 BiCGSTAB/GMRES 的测量仅代表对应旧配置。
- 项目有并行的未提交改动。提交、切换分支或集成之前，先确认这些改动的归属。历史测量不能代替当前构建验证。

## 最小复跑

在仓库根目录用 Visual Studio 多配置构建。PARDISO 需要找到 oneMKL，缺 DLL 时先核对实际运行的可执行文件路径和该构建的运行库配置。

```powershell
cmake -S . -B build_agent_pardiso -DBPFEM_USE_MKL=ON
cmake --build build_agent_pardiso --config Release
ctest --test-dir build_agent_pardiso -C Release --output-on-failure
.\build_agent_pardiso\Release\bp_fem_solver.exe --aedt IOStructure.aedt --mesh current.ngmesh --basis-order 0 --port-method numerical --sweep direct --linear-solver direct --max-sweep-points 1 --no-write-all-fields --out result_agent_smoke
```

检查 `result_agent_smoke/result_DIRECT/` 下的 `s_parameters.csv`、`run.json` 和 `timing.json`。单点烟测使用工程设置频率 95 GHz。宽带基准为同一 AEDT 和网格的 90–100 GHz、101 点。定量精度比较使用同网格 Direct 的复数 S 参数。仓库的 `IOStructure_S_parameters.csv` 提供 HFSS 幅值参考。对比应为每次运行指定独立 `--out`，同一策略会重建其结果子目录。复现算法命令见[公式复现](../optimization/awe-family/paper-reproduction.md)。

## 离开前

把修改点、验证命令与实测结果、失败或未核实项写到 [docs/changes](../changes/README.md)。需交给下一人继续做的事项，用[交接模板](handoff-template.md)注明当前状态。
