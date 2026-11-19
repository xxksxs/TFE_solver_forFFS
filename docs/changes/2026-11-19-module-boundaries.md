# 端口、扫频与线性求解职责核对

## 审阅范围

以 WCAWE 残差校正版 `e0d1823` 为基线，沿输入、端口、本征模及体装配、边界条件、线性后端、预条件、扫频、降阶数值组件和输出报告核对调用关系。重点检查是否将使用某个底层组件等同于采用该组件的求解策略。当前分层见[模块布局](../architecture/module-layout.md)和[求解器选择](../architecture/solver-selection.md)。

## 调整

- **线性与扫频选择**。`--linear-solver pardiso` 明确选中 PARDISO，`direct` 是历史兼容别名。MUMPS 同属直接后端。`--sweep direct` 表示逐频求解全阶系统，可使用直接或迭代后端。`auto` 在工厂统一解析。
- **预条件器与直接因子**。单精度 PARDISO 预条件器独立为 `PreconPardisoSingle`，实现 `IPreconditioner`，共享底层 `MklPardisoSolver`。`PardisoBackend.hpp` 保留旧包含方式的兼容性。Schwarz 的 MUMPS 子域求解仍由预条件器工厂组装。
- **端口表示**。`numerical-multimode` 通过 `PortModeSolver` 和 `WavePortBC` 使用全棱边坐标。`tfe-block` 由 `TfeSystem` 执行模态迹替换。`tfe` 和 `transfinite` 保留为全棱边多模的别名，`--port-modes` 和旧 `--tfe-modes-per-port` 等价。
- **残差校正与降阶**。通用回代校正移到 `linalg/FactorizedSolveSession.hpp` 的 `solveWithRefinement`。WCAWE 传递 `1e-10` 验收门槛，并负责构基、递推与降阶投影。常规 `solve`、`solveBatch` 路径保持各自契约。
- **报告分层**。CLI 参数转换移到 `app/RunReportOptions.cpp`，`core/RunReport.cpp` 负责 JSON 写出。`run.json` 和 `timing.json` 增加实际执行的 `port_formulation`、`linear_solver_family`、`active_preconditioner`、`active_iterative_spmv` 和 `active_schwarz_local_solver`。保留旧字段和输出目录名。直接后端收到迭代 SpMV 请求时给出提示，实际配置字段与警告一致。

## 验证口径

在相互独立的构建目录执行 Release 和 CTest，覆盖 MKL、MUMPS、两者同时启用和两者同时关闭四种组合。CLI 使用 `tests/fixtures/tfe-waveguide.aedt` 和 `tfe-waveguide.ngmesh`，检查旧别名与新名称、自动后端选择、两种使用 MUMPS 的方式、混合精度 GMRES 和分块 TFE。IOStructure 使用 90–100 GHz 41 点、95 GHz 展开点、16 阶 BLR WCAWE，复核回代计数、递推残差和 S 参数。此改动的测试目标为分类、接线和数值回归。

四种 Release 构建共通过 16 项 CTest，12 个 CLI 对照达到预期，其中包含改前基准和一个预期拒绝的 `tfe-block`＋WCAWE 组合。旧 `tfe`、`transfinite` 与 `numerical-multimode` 的 S 参数 CSV 字节一致。全局 MUMPS 记录 `direct`，GMRES＋Schwarz 记录 `iterative` 并单列本地 MUMPS。BLR WCAWE 在全带构足 16 阶，一次分解、32 次回代，最大递推残差 `4.58e-12`，复数 S 参数各分量与改前结果一致。修改页的相对链接及 `core`、`linalg` 的项目内部 include 方向均已核查。
