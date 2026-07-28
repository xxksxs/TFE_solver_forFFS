# BP-FEM 微波频域有限元求解器

BP-FEM 是一个面向微波器件宽带分析的 C++17 频域有限元求解器。它可以读取 AEDT 工程定义和 NGMesh 四面体网格，建立三维 Maxwell 有限元方程，计算端口 S 参数与电磁场，并用多种模型降阶算法加速扫频。

项目当前主要服务于波导滤波器和端口器件研究，同时也是一个用于比较 AWE、GAWE、MGAWE、WCAWE 与 ALPS 等快速扫频方法的实验平台。

## 为什么做这个项目

对一个包含数十万自由度的三维电磁模型，传统扫频需要在每个频点重新装配或求解大型稀疏线性系统。频点增加后，数值分解时间和峰值内存很快成为主要成本。

这个项目希望回答三个问题：

1. 能否用开放、可检查的有限元流程复现商业软件的 S 参数结果？
2. 能否在保持可控误差的前提下，用模型降阶减少宽带扫频中的大型求解次数？
3. AWE-family 与 Lanczos-Padé 方法在精度、时间、内存和数值稳定性上分别适合什么场景？

因此，项目不仅输出最终曲线，也记录矩阵分解、回代、降阶维数、正交性、残差、运行时间和峰值内存，方便对算法做可重复的比较。

## 主要功能

- **三维频域 FEM**：支持零阶 Whitney/Nédélec 棱元，并提供可选的一阶层次棱元。
- **工程与网格输入**：读取 AEDT 中的材料、变量、扫频和端口定义，以及 NGMesh 四面体网格。
- **端口建模**：支持数值端口模、矩形波导解析模和多模 TFE 端口。
- **线性求解**：支持 Intel oneMKL PARDISO，并保留 BiCGSTAB、GMRES 与预条件器扩展接口。
- **扫频方法**：提供 Direct、AWE、GAWE、MGAWE、WCAWE 和 ALPS。
- **结果输出**：生成 S 参数 CSV、最后频点场 VTU、运行日志、计时和数值诊断文件。
- **HFSS 校验**：仓库包含两套 HFSS S 参数参考数据，可用于比较 S11 和 S21。

| 方法 | 主要特点 | 更适合的场景 |
|---|---|---|
| Direct | 每个频点求解完整 FEM 系统 | 基准结果、少量频点 |
| AWE | 单点矩展开与 Padé 近似 | 展开点附近的窄带问题 |
| GAWE | 单点矩空间正交化并进行 Galerkin 投影 | 比传统 AWE 更稳定的单点降阶 |
| MGAWE | 多个展开点组成统一 Galerkin 空间 | 共振较多、频带较宽的问题 |
| WCAWE | 在矩递推过程中构造良条件基 | 高阶单点展开、传统 AWE 易病态时 |
| ALPS | 双边 Lanczos-Padé 降阶 | 追求较低在线扫频成本的实验路径 |

## WCAWE 的优势

WCAWE 的全称是 **Well-Conditioned Asymptotic Waveform Evaluation**。它解决的是传统 AWE 在阶数升高后，矩向量迅速趋于线性相关、条件数恶化的问题。

本项目中的 WCAWE 按 Slone 等人在 2003 年提出的递推形式实现。新基向量不是在完整 AWE 矩生成后再统一做 QR，而是在每一阶求解前，就利用已有正交基和上三角矩阵构造校正项。这样做有几个实际好处：

- 基在生成过程中保持良条件，高阶递推更不容易因矩向量病态而失效。
- 不需要扩展全阶未知量，内存主体仍约为 `N × q` 的降阶基。
- 单展开点只需对展开矩阵做一次数值分解，`q` 阶模型随后执行 `q` 个 RHS 回代。
- 相比需要多个展开点的 MGAWE，单点 WCAWE 通常需要更少的全阶矩阵分解。
- 会输出条件数、递推残差、基关系残差和正交误差，数值稳定性可以被直接检查。

WCAWE 的优势主要是**高阶数值稳定性和单点离线成本**，并不意味着它在所有宽带问题上都比 MGAWE 更精确。当前 90–100 GHz、约 36.5 万自由度的算例中，`q=12` 的 WCAWE 对 HFSS 的 S11/S21 幅值相对 L2 误差约为 `0.774%/0.990%`；继续提高到 `q=20/30/40` 时递推仍然稳定，但宽带误差没有单调下降。对于跨越多个共振的宽频带，MGAWE 的多展开点覆盖通常更有优势。

WCAWE 的数学推导和代码约定见 [WCAWE 文档](docs/optimization/awe-family/wcawe-skill.md)。

## 快速开始

### 1. 构建

推荐使用 Visual Studio 2022、CMake 和 Intel oneMKL：

```powershell
cmake -S . -B build_pardiso -DBPFEM_USE_MKL=ON
cmake --build build_pardiso --config Release
```

没有 oneMKL 时可以构建迭代求解版本：

```powershell
cmake -S . -B build -DBPFEM_USE_MKL=OFF
cmake --build build --config Release
```

### 2. 运行一次 WCAWE 扫频

下面的命令使用仓库内的 IOStructure 工程、当前零阶网格和 PARDISO，对 90–100 GHz 执行 101 点扫频：

```powershell
.\build_pardiso\Release\bp_fem_solver.exe `
  --aedt IOStructure.aedt `
  --mesh current.ngmesh `
  --basis-order 0 `
  --port-method numerical `
  --linear-solver direct `
  --sweep wcawe `
  --wcawe-order 12 `
  --max-sweep-points 101 `
  --no-write-all-fields `
  --out result
```

结果写入 `result/result_WCAWE/`。再次运行 WCAWE 时只会清空并重建这个子目录，不会不断创建带频点数的新目录。

常用结果文件包括：

- `s_parameters.csv`：复数 S 参数及幅值。
- `field_last.vtu`：最后频点的重构场，可用 ParaView 或 VisIt 查看。
- `run.log`、`run.json`：运行过程与环境信息。
- `timing.json`：分解、回代、离线建模、在线扫频和峰值内存。
- `diagnostics.json`：残差、正交性、矩匹配和算法终止原因。
- `basis_condition.csv`：WCAWE 各阶基条件与递推诊断。

## 示例与参考数据

- `wg_bp_filter.aedt` 与 `wg_bp_filter_S_parameters.csv`：40–43 GHz 波导带通滤波器示例。
- `IOStructure.aedt`、`current.ngmesh` 与 `IOStructure_S_parameters.csv`：90–100 GHz、约 36.5 万未知量的当前大网格基准。

HFSS 数据用于外部校验，不是求解器的网格来源。求解器直接读取 `current.ngmesh` 建立自己的有限元拓扑和矩阵。

## 文档

- [工程与命令参考](docs/engineering-reference.md)：完整 CLI、输出约定、模块布局和开发者信息。
- [项目阅读入口](docs/overview.md)：运行流程、扩展点和文档地图。
- [AWE-family](docs/optimization/awe-family/README.md)：AWE、GAWE、MGAWE、WCAWE 理论与验证。
- [ALPS](docs/optimization/alps-sweep/README.md)：Lanczos-Padé 路径与当前实现状态。
- [测试](docs/testing/README.md)：构建、单元测试和烟测。
- [验证](docs/validation/README.md)：端口、残差、S 参数和场结果检查。

## 当前定位

BP-FEM 目前是研究与算法验证性质的求解器，不是商业电磁软件的完整替代品。快速扫频路径目前重点覆盖无损材料、单端口激励和 P1 端口色散模型；涉及有损材料、更复杂端口色散或生产级鲁棒性时，应优先使用 Direct 路径进行交叉验证。
