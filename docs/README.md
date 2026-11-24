# BP-FEM 项目文档

本目录是项目用法、算法、开发规范和变更记录的统一入口。新建正式项目文档时放在这里，并更新对应专题索引。

BP-FEM 是读取 AEDT 工程与 NGMesh 网格的 C++17 三维频域 Maxwell 有限元求解器。它支持零阶和一阶层次棱元，输出 S 参数、场和数值诊断。项目用于微波器件验证和快速扫频研究。

主要算例有 IOStructure（90–100 GHz，当前主仓库网格）和波导滤波器（40–43 GHz，使用其配套网格）。HFSS 导出的 CSV 是对照数据。关于算例与参考文件，见[工程参考](engineering-reference.md)。

新 Agent 先读[快速接手指南](development-workflow/agent-handoff.md)，核对当前工作树并复跑最小算例。

## 快速开始

在仓库根目录构建。使用 PARDISO 时需要已配置的 oneMKL。找不到 MKL 时，默认构建可以回退到迭代后端。构建与运行的完整参数见[工程参考](engineering-reference.md)。

```powershell
cmake -S . -B build_pardiso -DBPFEM_USE_MKL=ON
cmake --build build_pardiso --config Release
ctest --test-dir build_pardiso -C Release --output-on-failure
.\build_pardiso\Release\bp_fem_solver.exe --aedt IOStructure.aedt --mesh current.ngmesh --sweep direct --linear-solver direct --max-sweep-points 1 --out result_smoke
```

该示例使用 Visual Studio 多配置生成器。每次运行会重建指定输出根目录下该策略的结果子目录。做对照实验时使用独立的 `--out`。

## 查阅入口

| 需要查什么 | 维护入口 |
|---|---|
| 构建、CLI、输入和输出 | [工程参考](engineering-reference.md) |
| 调用链、模块、依赖、线性求解器 | [架构](architecture/README.md) · [求解器选择](architecture/solver-selection.md) · [文档地图](overview/document-map.md) |
| 物理、装配、线性代数和不变量 | [原语](primitives/README.md) |
| AWE、MGAWE、WCAWE 现行公式接线和命令 | [公式复现](optimization/awe-family/paper-reproduction.md) |
| ALPS 现行实现 | [ALPS 专题](optimization/alps-sweep/README.md) |
| 端口理论与规划 | [解析端口](optimization/analytic-port/README.md) · [TFE 专题](optimization/transfinite-element/README.md) |
| 构建、烟测及数值验收 | [测试](testing/README.md) · [验证](validation/README.md) |
| 修改代码或文档 | [修改规范](development-workflow/change-policy.md) · [开发流程](development-workflow/README.md) · [变更记录](changes/README.md) |
| 历史方案、未完成计划和迁移候选 | [优化路线](optimization/README.md) · [历史档案](history/README.md) · [迁移记录](development-workflow/document-migration-plan.md) · [参考资料](references/README.md) |

阅读技术页时区分**现行参考**、**理论或设计**、**历史基准**与**待核实**。现行行为使用当前源码和可复现验证确认。旧配置与理论推导仍可作为研究资料。
