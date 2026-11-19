# 模块布局

`bp_fem_solver` 是命令行入口，链接承载各模块实现的 `bp_fem_core` 静态库。下面的三项选择由不同模块处理。

| 选择 | 负责的问题 | 源码入口 |
| --- | --- | --- |
| `--port-method` | 端口的离散表示和边界装配 | `PortModeSolver`、`AnalyticPortBuilder`、`WavePortBC`、`TfeSystem` |
| `--sweep` | 逐频求解全阶系统，或用展开点构造降阶模型 | `SweepStrategyFactory`、`ISweepStrategy` |
| `--linear-solver` / `--precon` | 全阶稀疏方程的求解后端及外层迭代的预条件 | `SparseSolverFactory`、`ISparseSolver`、`IPreconditioner` |

`--sweep direct` 对各频点保留全阶问题，可搭配 PARDISO、MUMPS、GMRES 和 BiCGSTAB。WCAWE 通过同一个线性接口构造离线矩，在线求解降阶矩阵。原生 `tfe-block` 替换端口迹自由度，支持全阶扫频和分段二阶 MGAWE。历史 `tfe` 入口对应 `numerical-multimode` 全棱边路径。[求解器选择](solver-selection.md) · [TFE 专题](../optimization/transfinite-element/README.md)

## 源码职责

| 目录 | 职责 |
| --- | --- |
| `app/` | 解析 CLI、串起读取到导出的流程、记录请求参数和实际配置。 |
| `factory/` | 将配置解析为稀疏后端和扫频策略。 |
| `core/` | 常量、共享数据、日志与通用报告序列化。CLI 选项到报告的转换在 `app/RunReportOptions.cpp`。 |
| `io/` | AEDT 与 NGMesh 输入及端口面映射。 |
| `fem/` | 棱元拓扑、体矩阵和数值波导本征模。`FEMAssembler` 通过 `IBoundaryCondition` 添加边界项。 |
| `bc/` | 棱边坐标下的边界项。`WavePortBC` 支持单模及多模秩一更新。 |
| `apm/` | 矩形端口几何识别和解析 TE10 模。 |
| `tfe/` | `TfeSystem` 在模态迹坐标装配并恢复全棱边场。旧 `TransfiniteElementBuilder` 名称保留作 API 兼容。 |
| `linalg/` | 稀疏矩阵、直接和迭代后端、预条件、因子复用及真实残差校正。PARDISO 与 MUMPS 可作为全局直接后端，也可当 Schwarz 的局部直接后端。 |
| `fastsweep/` | 多项式、矩递推、正交化、Padé 和 Galerkin 等扫频数值组件。 |
| `sweep/` | 组合同一端口模型与指定线性后端，管理展开点、构基和在线求值。通过 `SweepContext` 回调传递场与统计量。 |
| `post/` | S 参数提取及 VTU、CSV 输出。 |
| `mor/` | 旧 ALPS 头文件的兼容转发。 |

新增配置或改动接线时，同时核对实际后端、已启用的预条件与端口表示。详见[依赖方向](dependency-direction.md)和[修改规范](../development-workflow/change-policy.md)。
