# 求解器选择

本文说明应用层如何选择线性求解器、预条件和扫频策略。Phase 1 接口化后，
应用层通过工厂分发。构建宏也用于条件编译求解器后端。

## 选择逻辑

### 线性求解器（`--linear-solver`）

- **`auto`**（默认）：`BPFEM_USE_MKL` 已定义时返回 `PardisoBackend`，否则返回
  `BiCGStabBackend`。这条路径与 Phase 1 之前的编译期选择字节级一致。
- **`direct`**：始终返回 `PardisoBackend`。若构建时 `BPFEM_USE_MKL=OFF`，
  在 `solve()` 调用时抛出 runtime_error，提示重新构建或回退到 bicgstab。
- **`bicgstab`**：始终返回 `BiCGStabBackend`，即使有 MKL。
- **`gmres`**：复数 GMRES，默认重启长度 30，用 `--gmres-restart` 调整。
- **`mumps`**：可选复数直接后端。需使用 `-DBPFEM_USE_MUMPS=ON` 配置，支持 `--mumps-mode incore|ooc|blr`。
- PARDISO 可用 `--pardiso-mode incore|ooc`。迭代后端支持 `--iterative-spmv native|mkl`，MKL 路径要求 oneMKL 构建。

实现位于 `src/factory/SparseSolverFactory.cpp::makeSparseSolver(options)`。
此处负责运行时后端选择和选项组合校验。

### 预条件（`--precon`）

- **`none`**（默认）：无预条件，BiCGSTAB / GMRES 走原始未修改的迭代序列。
- **`jacobi`**：对向量应用复对角逆，极小或缺失的对角项使用保护处理。
- **`ilu0`**：零填充 ILU 因式分解。

- **`schwarz`**：重叠子域预条件。局部后端由 `--schwarz-local-solver pardiso|mumps` 选择，要求相应构建支持。局部 MUMPS 仅支持 `incore` 或 `blr`。

以下是旧版 BP filter 试验记录。适用范围为当时的配置，详细参数待核实。

| 配置 | iter | residual |
|---|---|---|
| BiCGSTAB + Jacobi | 4000 | 0.91 |
| GMRES(50) + ILU(0) | 1500 | 0.98 |

频域 H(curl) 不定系统的收敛取决于离散、预条件器和频点，应报告残差、迭代上限及未收敛点。当前新增后端的性能与发布状态待验证。

预条件只对 `bicgstab` / `gmres` 生效。选中 PARDISO 或 MUMPS 时，工厂忽略预条件并输出警告。

**当前推荐**：默认 PARDISO；GMRES + ILU(0) 框架作为后续 AMS 等更复杂预条件的 plug-in 接口保留，**不**作为现阶段 BP filter 求解路径推荐。

### 扫频策略（`--sweep`）

- **`direct`**（默认）：`sweep::DirectSweep`。每个频点重装配 + 直接求解 + 后处理。
  唯一支持 per-frequency VTU 输出的路径。
- **`alps`**：单点双边 Lanczos–Padé，见 [ALPS 专题](../optimization/alps-sweep/README.md)。
- **`awe`、`gawe`、`mgawe`、`wcawe`**：详见[工程参考](../engineering-reference.md)和[公式复现](../optimization/awe-family/paper-reproduction.md)。

实现位于 `src/factory/SweepStrategyFactory.cpp::makeSweepStrategy(...)`。

## CMake 构建标志

`BPFEM_USE_MKL` 控制是否尝试启用 Intel oneMKL PARDISO：

- **`BPFEM_USE_MKL=ON`（默认）且找到 MKL**：定义编译宏 `BPFEM_USE_MKL` 并链接
  `MKL::MKL`。`PardisoBackend` 走真实 PARDISO；`auto` 默认选 PARDISO。
- **`BPFEM_USE_MKL=OFF`**：不定义宏，不链接 MKL。`auto` 选 BiCGSTAB。
  显式 `--linear-solver direct` 在 solve() 时抛错。
- **`BPFEM_USE_MKL=ON` 但未找到 MKL**：CMake 打印 warning 并相当于 `OFF`。

## `SolveResult` 契约

所有 `ISparseSolver` 实现返回相同语义的 `SolveResult`：

- `field`：复数边自由度向量。
- `iterations`：迭代次数（直接求解器返回 1 作为哨兵）。
- `residual`：相对残差 `||A x - b|| / ||b||`。

修改求解器后端时必须保持此契约稳定。

## 修改求解器 / 预条件 / 扫频时的要求

- **保持 `--linear-solver auto` / `--precon none` / `--sweep direct` 默认行为
  与 main 分支字节级一致**，否则下游脚本会回归。
- **新增条目走"接口实现 + 工厂枚举 + CLI 字符串映射"三步法**，不要直接在
  `Application.cpp` 里加 `if/else`。详见
  `docs/architecture/dependency-direction.md` 第 "具体策略加新条目时的最小改动集"
  段。
- **保持 MKL 和非 MKL 两种构建均可编译**。`PardisoBackend` 的非 MKL fallback
  必须在 solve() 时抛 runtime_error 而不是构造时；构造侧要做到无副作用，
  让 `factory` 始终能 `make_unique`。
- **CMake 修改后至少验证一次 `-DBPFEM_USE_MKL=OFF` 构建**。
- **求解结果差异先看残差**：BiCGSTAB 默认 tol = 1e-7 与 PARDISO 直接解差到
  1e-7 级别属于正常；只有当残差 < tol 而 S 参数仍大幅偏离参考时才追根因。

## 扩展示例

加新条目的可执行示例见
`docs/optimization/strategy-interfaces/example.md`（用 Impedance BC + GMRES
backend 的虚构组合演示完整流程）。
