# 工程与自动化参考

本文承接根目录 README 中不适合首次阅读的工程细节，供维护者、测试脚本和 AI 编码助手查询。项目用途、动机和主要能力请先阅读 [README](../README.md)。

## 输入、求解与输出数据流

```text
AEDT project + NGMesh tetrahedral mesh
        |
        v
AEDTParser + NGMeshParser + PortFaceResolver
        |
        v
ProjectDefinition + Mesh + EdgeTopology
        |
        v
PortModeSolver / AnalyticPortBuilder / TfeSystem
        |
        v
FEMAssembler -> sparse K/M/port operators + RHS
        |
        +--> DirectSweep -> full-order solves
        |
        +--> AWE / GAWE / MGAWE / WCAWE / ALPS -> reduced model
        |
        v
ResultExtractor + OutputWriter
        |
        v
result/result_<SWEEP>/
```

HFSS CSV 只作为外部参考，不参与有限元矩阵装配。网格和自由度来自 `current.ngmesh`。

## 构建配置

项目要求 CMake 3.20、C++17 和支持 C++17 的编译器。Windows 推荐 Visual Studio 2022。

### oneMKL/PARDISO

```powershell
cmake -S . -B build_pardiso -DBPFEM_USE_MKL=ON
cmake --build build_pardiso --config Release
```

`BPFEM_USE_MKL` 默认开启。CMake 找不到 MKL 时会输出警告并退回 BiCGSTAB；需要 PARDISO 时应从 oneAPI 环境配置，或显式设置 `MKL_DIR`。

### 无 MKL

```powershell
cmake -S . -B build -DBPFEM_USE_MKL=OFF
cmake --build build --config Debug
```

### VTU 兼容选项

`BPFEM_OUTPUT_PARAVIEW=ON` 输出 ParaView 风格的 VTK Lagrange tetrahedron。默认关闭时，一阶/二阶场输出采用更兼容 VisIt 的单元类型。该选项不影响 S 参数。

## CLI 基本约定

```powershell
bp_fem_solver.exe `
  --aedt <project.aedt> `
  --mesh <mesh.ngmesh> `
  --basis-order 0|1 `
  --port-method numerical|analytic|numerical-multimode|tfe-block `
  --linear-solver auto|pardiso|bicgstab|gmres|mumps `
  --sweep direct|awe|gawe|mgawe|wcawe|alps `
  --out result
```

主要参数：

| 参数 | 含义 |
|---|---|
| `--aedt` | AEDT 工程文件路径 |
| `--mesh` | NGMesh 网格路径 |
| `--basis-order` | `0` 为 Whitney/Nédélec 零阶棱元，`1` 为一阶层次棱元 |
| `--max-sweep-points` | 从工程扫频中最多取多少个频点 |
| `--port-method` | 数值单模、解析 TE10、全棱边数值多模或原生 TFE 分块。历史 `tfe`、`transfinite` 映射为 `numerical-multimode`。见[TFE 专题](optimization/transfinite-element/README.md)。 |
| `--port-modes` | 全棱边多模和分块 TFE 的每端口模式数量，兼容 `--tfe-modes-per-port`。 |
| `--linear-solver` | 自动路由、PARDISO/MUMPS 直接后端或 GMRES/BiCGSTAB 迭代后端。历史 `direct` 映射为 `pardiso`。 |
| `--precon` | `none`、`jacobi`、`ilu0`、`schwarz`，或 oneMKL 构建下专供 GMRES 使用的 `pardiso-single`。见[求解器选择](architecture/solver-selection.md)。 |
| `--sweep` | `direct` 表示逐频全阶求解，线性后端由 `--linear-solver` 决定。其余策略使用降阶扫频。 |
| `--write-all-fields` | 为每个频点写 VTU |
| `--no-write-all-fields` | 只保留最后频点场，长扫频推荐 |
| `--field-output-order` | 控制 VTU 采样阶数，不改变 S 参数 |
| `--out` | 统一结果根目录，通常为 `result` |

参数解析和默认值以 `src/app/Application.cpp` 为准。新增参数时必须同步 CLI 帮助、README/本文和 smoke test。

## 扫频策略

当前 AWE 逐场 Padé（`--awe-field-pade`）、MGAWE 固定点与自适应（`--mgawe-expansion`、`--mgawe-adaptive`、`--mgawe-2001`）以及 WCAWE 独立矩检查（`--wcawe-check-moments`）的参数和限制见[公式复现](optimization/awe-family/paper-reproduction.md)。

### Direct

```powershell
.\build_pardiso\Release\bp_fem_solver.exe `
  --aedt IOStructure.aedt --mesh current.ngmesh `
  --sweep direct --linear-solver pardiso --out result
```

### AWE

```powershell
.\build_pardiso\Release\bp_fem_solver.exe `
  --aedt IOStructure.aedt --mesh current.ngmesh `
  --sweep awe --awe-order 8 --out result
```

### GAWE

```powershell
.\build_pardiso\Release\bp_fem_solver.exe `
  --aedt IOStructure.aedt --mesh current.ngmesh `
  --sweep gawe --gawe-order 12 --out result
```

### MGAWE

```powershell
.\build_pardiso\Release\bp_fem_solver.exe `
  --aedt IOStructure.aedt --mesh current.ngmesh `
  --sweep mgawe --mgawe-points 3 --mgawe-order 4 --out result
```

`--mgawe-points` 控制展开点数量，`--mgawe-order` 是每个展开点的局部阶数。所有局部矩向量进入同一个全局 Galerkin 空间，不是多个局部模型拼接。

### WCAWE

```powershell
.\build_pardiso\Release\bp_fem_solver.exe `
  --aedt IOStructure.aedt --mesh current.ngmesh `
  --sweep wcawe --wcawe-order 12 --out result
```

`--wcawe-order` 是目标良条件基维数。实现必须遵守 `docs/optimization/awe-family/wcawe-implementation-skill.md`：候选向量在求解前加入上三角校正，禁止退化为“传统 AWE 矩事后 QR”。

### ALPS

```powershell
.\build_pardiso\Release\bp_fem_solver.exe `
  --aedt IOStructure.aedt --mesh current.ngmesh `
  --sweep alps --alps-order 6 --out result
```

`--alps-order` 表示每个标量 Padé 模型的阶数，不是端口色散矩阵的阶数。旧参数 `--alps-krylov-order` 仅作为弃用别名保留。

## 端口方法

- `numerical`：在端口表面三角网格上求解二维 H(curl) 广义本征问题，适用于任意可解析截面，是默认和工程校验路径。
- `analytic`：矩形波导 TE10 闭式模式，只适用于可可靠识别为矩形的端口。
- `tfe`：保留多个数值端口本征模并组装多模端口算子；`--tfe-modes-per-port` 控制每个物理端口保留的模式数。

端口面识别、方向、激励端口和 S 参数索引必须联合验证，不能只根据 HFSS 曲线形状判断正确性。

## 线性求解器约定

- `auto`：MKL 构建选择 PARDISO，否则使用 BiCGSTAB。
- `direct`：强制 PARDISO；无 MKL 构建应明确报错。
- `bicgstab`、`gmres`：迭代后端，保留 Jacobi/ILU0 扩展接口。

当前工作树还提供可选 MUMPS 后端、PARDISO out-of-core、Schwarz 预条件和 MKL SpMV 入口。构建依赖和组合限制见[求解器选择](architecture/solver-selection.md)。常规预条件器对频域 H(curl) 不定系统的效果需按具体算例验证。生产基准和 fast-sweep 离线阶段优先使用 PARDISO。支持分解复用的后端通过 `IFactorizedSparseSolver` 暴露 `factorize`、`solveFactorized` 和批量 RHS 能力。

## 结果目录与文件

`--out result` 会映射到固定算法子目录：

```text
result/
  result_DIRECT/
  result_ALPS/
  result_AWE/
  result_GAWE/
  result_MGAWE/
  result_WCAWE/
```

每次运行只清空当前算法对应的子目录。不要按频点数或时间戳创建新的结果目录。

通用输出：

- `s_parameters.csv`
- `field_last.vtu`
- `run.log`
- `run.json`
- `timing.json`
- `diagnostics.json`

WCAWE 额外输出 `basis_condition.csv`。逐频点场只有在 `--write-all-fields` 时写出，长扫频默认关闭，避免磁盘膨胀。

## HFSS 比较

- `wg_bp_filter_S_parameters.csv` 对应 `wg_bp_filter.aedt`，频带为 40–43 GHz。
- `IOStructure_S_parameters.csv` 对应 `IOStructure.aedt`，频带为 90–100 GHz。

旧 BP filter 的批量比较入口：

```powershell
python scripts/compare_with_hfss.py `
  wg_bp_filter_S_parameters.csv `
  result/result_BENCHMARK `
  --batch-root result
```

IOStructure 专用对比逻辑位于 `scripts/compare_ios_fastsweeps.mjs`。主要指标使用 S11/S21 幅值相对 L2 误差；深零点应同时报告绝对幅值误差，避免逐点相对误差失真。

## 模块边界

```text
include/bpfem/ + src/
  app/         CLI 与端到端应用流程
  core/        共享类型、日志、计时与运行报告
  io/          AEDT、NGMesh 和端口面解析
  fem/         H(curl) 拓扑、端口模和 FEM 装配
  bc/          边界条件
  linalg/      稀疏矩阵、PARDISO、迭代后端和预条件器
  fastsweep/   矩递推、WCAWE、Lanczos-Padé 和统一 ROM
  sweep/       Direct/AWE/GAWE/MGAWE/WCAWE/ALPS 策略
  factory/     求解器和扫频策略运行时路由
  post/        S 参数、CSV 和 VTU 输出
  apm/         解析端口模
  tfe/         多模超限元端口
```

公共数据结构优先放在 `include/bpfem/core/Types.hpp`。新后端应实现既有接口并通过 factory 注册，不要在 `Application.cpp` 或 sweep 中硬编码具体求解器。

## 测试与提交前检查

```powershell
cmake --build build --config Debug
cmake --build build_pardiso --config Release
ctest --test-dir build -C Debug --output-on-failure
ctest --test-dir build_pardiso -C Release --output-on-failure
```

修改 fast-sweep 数学过程后，至少检查：

- 展开点处与 Direct 的复数 S11/S21 一致性。
- 矩递推或 Lanczos 递推残差。
- 基正交性、矩匹配和 breakdown 语义。
- PARDISO 数值分解次数与 RHS 回代次数。
- HFSS 相对 L2 误差、总时间和峰值内存。

详细测试入口见 `docs/testing/README.md`，外部校验见 `docs/validation/README.md`，AWE-family 的强制条件见 `docs/optimization/awe-family/validation.md`。

## 文档阅读顺序

1. `README.md`：项目用途和快速开始。
2. `docs/overview.md`：工程数据流和扩展入口。
3. `docs/architecture/README.md`：模块依赖和风险点。
4. `docs/primitives/README.md`：物理量、网格、FEM 与线性代数约定。
5. `docs/optimization/awe-family/README.md` 或 `docs/optimization/alps-sweep/README.md`：算法专题。
6. `docs/testing/README.md` 与 `docs/validation/README.md`：验证和回归。

面向任务拆分和自动化修改的约定保留在 `docs/task-splitting/` 与 `docs/development-workflow/`，不应重新堆回根 README。
