# 主要扩展点

本文说明常见功能修改应从哪里入手。

## 新增输入字段

优先修改：

- `../../include/bpfem/core/Types.hpp`
- `../../src/io/AEDTParser.cpp`
- `../../src/io/NGMeshParser.cpp`

共享结构变化必须同步更新所有消费者。

## 修改网格拓扑

优先修改：

- `../../src/fem/EdgeTopology.cpp`
- `../../include/bpfem/fem/EdgeTopology.hpp`

修改后检查 FEM 装配、端口模式和输出场重构。

## 修改端口物理

优先修改：

- `../../src/fem/PortModeSolver.cpp`
- `../../src/fem/FEMAssembler.cpp`
- `../../src/post/ResultExtractor.cpp`

端口模式、边界装配和 S 参数提取必须保持一致。

## 修改求解器

优先修改：

- `../../src/linalg/SparseMatrix.cpp`
- `../../src/linalg/MklPardisoSolver.cpp`
- `../../src/linalg/BiCGStabSolver.cpp`
- `../../src/app/Application.cpp`

保持 `SolveResult` 契约稳定。

## 修改输出

优先修改：

- `../../src/post/ResultExtractor.cpp`
- `../../src/post/OutputWriter.cpp`

修改输出字段时同步更新测试和文档。

修改 VTU 节点排布、新增高阶单元类型支持时，必须保持 `OutputWriter::buildLagrangeTetTable` 与 `OutputWriter::reconstructLagrangePointField` 的节点顺序一致：边内点按 VTK 边序、面内点按 VTK 面序、体内点按重心坐标 lex 顺序。`--field-output-order` 与 `Cells/HigherOrderDegrees` 必须同步更新；只改其中一处会让 ParaView 渲染出错。

VTU 中只输出 `E_real / E_imag / H_real / H_imag` 四个矢量场。`writeVTU` 因此需要传入 `ProjectDefinition&`（取材料 μ_r）和当前 `frequencyHz`（构造 `H = j ∇×E / (ω μ_r μ_0)`）。新增字段时务必把派生量留给 ParaView Calculator，不要把模值/dB/切向分量等扩展回 VTU——这是上一轮"清理场图输出"的明确不变量。
