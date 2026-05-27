# 架构风险点

本文列出架构层面的高风险耦合点。修改这些区域时应同时阅读相关专题文档和验证文档。

## 端口模式耦合

`PortModeSolver`、`bc::WavePortBC`、`FEMAssembler` 和 `ResultExtractor` 必须在
以下方面一致：

- faceId 匹配。
- 边方向：全局边按点号小到大，体单元和端口面的 `localNodes` 存储时也应按
  该方向排序。
- L2 模态归一化（模式形状）。
- 坡印廷功率归一化（依赖频率与 β，由 `powerNormalizationFactor` 给出）。
- 激励和投影定义。

`WavePortBC::apply` 和 `ResultExtractor::portProjection` 共用 `PortModeSolver`
缓存的 `couplingWeights = M_port v`，因此修改任何一处都要回归两条路径。

不要只修改其中一个模块而不验证其他模块。

## Dirichlet 约束副作用

PEC 零 Dirichlet 约束通过 DOF mask 批量施加。`FEMAssembler` 应先收集约束 DOF（`cachedConstrainedDofs()` lazy 缓存），体装配和端口装配跳过这些 DOF，
再调用 `SparseMatrix::imposeZeroDirichlet` 设置对角项。矩阵内部只存储对称上三角 CSR，因此修改 PEC 或约束逻辑时应检查：

- 被约束的边集合是否正确。
- 端口边是否被排除。
- RHS 是否在约束 DOF 上保持为 0。
- 上三角半存储中的对应列项是否被一次性清除。
- 修改端口集合 / 表面三角形 / faceId 解析后必须把 `constrainedDofsCached_` 与 `elementCacheBuilt_` 重置，否则后续频点会复用过期约束。

## 对称装配双计数

`SparseMatrix::add` 会把 `(row, col)` 规范到上三角位置。对称局部矩阵不能用完整双循环直接写入，
否则 `(i, j)` 和 `(j, i)` 会落到同一个上三角条目并被重复累加。体装配和端口吸收项都必须只写
局部上三角。

## 单位体系

频率、坐标、波数和材料参数必须保持单位一致。单位错误会直接造成截止频率、传播常数和场尺度错误。

## 场重构

VTU 输出来自边自由度重构，不是节点标量未知量。修改输出时必须理解边元场的重构方式。

## 求解器路径

MKL 直接求解和 BiCGSTAB fallback 两条路径都应保持可构建。Phase 1 之后唯一
允许出现 `BPFEM_USE_MKL` 求解器分支的文件是
`src/factory/SparseSolverFactory.cpp`；其它代码通过 `ISparseSolver*` 接口
访问。新增求解器后端不要在其它文件加 `#ifdef BPFEM_USE_MKL`。

## 边界条件分发

`FEMAssembler` 的 `assemble` / `sparsityPattern` / `buildAffineSystem` 都通过
注册的 `IBoundaryCondition` 列表分发。不要把新边界条件硬编码进 FEMAssembler，
应该新建 `bc/` 下的类。详见 `../optimization/strategy-interfaces/example.md`。

## 扫频与求解器解耦

`ISweepStrategy` 实现可以使用 `SweepContext::solver`（DirectSweep 这样做）或
内部建一个直接求解器（AlpsSweep 这样做，因为 ALPS 的 offline 因式分解需要
独立生命周期）。两条路径不要混用——sweep 内部的求解器不要把状态泄漏到
ctx.solver。
