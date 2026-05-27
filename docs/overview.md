# 工程阅读总览

本文是阅读和扩展 BP-FEM 频域求解器工程的主入口。详细总览内容已拆分到 `overview/` 目录，但本文件仍保留核心导航。

## 适用对象

- **新开发者**：用于快速了解模块顺序和主数据流。
- **AI 编码模型**：用于在修改代码前选择最小必要上下文。
- **代码审查者**：用于从输入文件追踪到 S 参数和 VTU 场输出。

## 推荐入口

- **总览索引**：`overview/README.md`
- **推荐阅读顺序**：`overview/reading-order.md`
- **文档地图**：`overview/document-map.md`
- **运行时数据流**：`overview/runtime-flow.md`
- **主要扩展点**：`overview/extension-points.md`
- **安全规则**：`overview/safety-rules.md`

## 专题文档入口

- **架构说明**：`architecture/README.md`
- **原语文档**：`primitives/README.md`
- **测试说明**：`testing/README.md`
- **验证说明**：`validation/README.md`
- **任务拆分**：`task-splitting/README.md`
- **开发流程**：`development-workflow/README.md`
- **商业级优化路线**：`optimization-roadmap.md`（专题目录 `optimization/README.md`）

## 运行时数据流摘要

```text
wg_bp_filter.aedt + current.ngmesh
        |
        v
AEDTParser + NGMeshParser
        |
        v
ProjectDefinition + Mesh
        |
        v
EdgeTopology + PortModeSolver
        |
        v
FEMAssembler -> SparseMatrix + RHS
        |
        v
MklPardisoSolver 或 BiCGStabSolver
        |
        v
ResultExtractor + OutputWriter
        |
        v
s_parameters.csv + field_last.vtu
```

## 快速安全规则

- **保持数据契约**：共享结构变化应先反映到 `include/bpfem/core/Types.hpp`。
- **不要孤立修改端口**：端口模式、边界装配和 S 参数投影必须保持一致。
- **求解器修改后必须验证**：至少跑一个单频点，检查残差和 S 参数。
- **保持 `main.cpp` 精简**：应用流程应放在 `src/app/Application.cpp`。
