# 线性代数原语

本文说明稀疏矩阵、右端项和求解器契约。修改 `src/linalg`、`CMakeLists.txt` 或求解器选择逻辑前应先阅读本文。

## SparseMatrix

`SparseMatrix` 是装配阶段使用的复数稀疏矩阵。

当前存储方式：

- 装配阶段可使用预构建 CSR 稀疏图，按固定 `(row, col) -> values[index]` 位置累加矩阵上三角非零元。
- 未提供预构建图时，仍可按行使用临时哈希表累加并在 `compress()` 时生成 CSR。
- 由于 FEM 主矩阵按双线性型对称装配，`add(r, c, v)` 会把下三角项折叠到上三角。
- 压缩后导出 0-based CSR：`rowOffsets`、`columnIndices` 和 `values`。
- 支持 `add`、`multiply`、单个 `imposeDirichlet` 和批量零约束 `imposeZeroDirichlet`。

该格式只存储非零元，并利用对称性减少存储量。`multiply` 会从上三角 CSR 恢复完整矩阵乘法；MKL PARDISO 路径直接读取 0-based CSR 和复对称矩阵类型。

`FEMAssembler` 会缓存固定稀疏图。只要约束 DOF 集合不变，扫频时每个频点复用同一套 `rowOffsets` 和
`columnIndices`，只重新填充 `values`。

## RHS 向量

RHS 是与矩阵同尺寸的复数向量。

主要来源：

- 波端口激励。
- 非零 Dirichlet 约束对 RHS 的修正。
- PEC 零 Dirichlet 约束会通过 mask 跳过装配并把受约束 RHS 项保持为 0。

RHS 的长度必须等于 `EdgeTopology::edgeCount()`。这里的 `edgeCount()` 表示全局 H(curl) 自由度数量；启用一阶层次棱元后，它会大于几何边数量。

## SolveResult

所有求解器都应返回同一个结构：

- **field**：求解得到的复数边自由度。
- **iterations**：迭代次数或直接求解器的等效标记。
- **residual**：相对残差。

下游模块不应依赖具体求解器类型。

## MKL PARDISO 求解器

`MklPardisoSolver` 使用 Intel oneMKL PARDISO 直接求解复数稀疏线性系统。

特点：

- 适合当前较难收敛的频域 Maxwell 矩阵。
- 使用 0-based CSR，与 `SparseMatrix` 压缩结果一致。
- 使用复对称矩阵类型，因此只需要上三角存储。
- 同一 CSR 稀疏结构下复用 symbolic analysis；扫频时通常只重复数值分解和求解。
- 通过 CMake 选项 `BPFEM_USE_MKL` 控制。

## BiCGSTAB fallback

`BiCGStabSolver` 是无 MKL 环境下的迭代 fallback。

注意：

- 高残差时结果不可信。
- 增加迭代次数不一定能解决病态矩阵问题。
- 修改 fallback 时也必须保证 MKL 路径仍能编译。

## 求解器修改清单

- 保持 `SolveResult` 契约不变。
- 保持 MKL 和非 MKL 两种构建可用。
- 求解后计算并报告残差。
- 先用小系统或单频点验证，再做全频扫。
