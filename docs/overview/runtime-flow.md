# 运行时数据流

本文说明求解器从输入文件到输出结果的主数据路径。

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

## 输入阶段

- `AEDTParser` 读取材料、扫频、端口和激励。
- `NGMeshParser` 读取点、面三角形、facet 和四面体。

## 离散阶段

- `EdgeTopology` 构建全局边自由度。
- `PortModeSolver` 计算端口表面的数值端口模式。
- `FEMAssembler` 装配矩阵和 RHS。

## 求解阶段

- MKL 可用时使用 `MklPardisoSolver`。
- 否则使用 `BiCGStabSolver` fallback。

## 输出阶段

- `ResultExtractor` 提取 S 参数。
- `OutputWriter` 写出 CSV 和 VTU 文件。
