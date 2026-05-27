# 常见失败特征

本文记录数值验证中常见失败现象和优先排查方向。

## 残差高且 S 参数异常

优先检查：

- 求解器收敛。
- 矩阵条件。
- Dirichlet 约束。
- RHS 激励。

## 全域零场

可能原因：

- RHS 激励缺失。
- 端口边被 PEC 约束。
- 端口模式无效。
- 端口 faceId 匹配失败。

## 传输接近数值零

可能原因：

- 端口投影失败。
- faceId 错误。
- 高残差。
- 模态归一化（L2 模式形状或坡印廷功率归一化）错误。
- 物理上确实处于阻带。

## 端口附近场不物理

可能原因：

- 模态归一化（L2 模式形状或坡印廷功率归一化）不一致。
- 边方向符号错误。
- 端口吸收项或激励项不一致。

## 端口本征求解 fatal `LAPACKE_dsygv info > N`

`LAPACKE_dsygv` 的 `info > N`（或 `info > 端口面 DOF 数`）表示端口面 mass 矩阵 Cholesky 失败，**M 不是正定的**。
对一阶 hierarchical 路径，这通常对应以下两类问题：

- 局部基函数族线性相关。最常见的是 `FaceFirst1` 写成 `lambda_a * N_bc - lambda_b * N_ac`，它恒等于 `-FaceFirst0 = -lambda_c * N_ab`。两份 face DOF 完全相反，让局部 mass 秩亏。
- 端口三角面求积阶不足。3 点 (degree 2) 规则对一阶 hierarchical mass（degree 4 被积函数）截断过大。

优先检查：

- `PortModeSolver::evaluateBasis` 与 `FEMAssembler::evaluateBasis` 的 face bubble 公式是否一致且线性独立。
- `PortModeSolver::computeMode` 中端口三角形求积是否 ≥ degree 4（当前为 7 点 Wandzura, degree 5）。

## 场图呈现块状或单元级异常振荡

可能原因：

- 对称上三角矩阵装配时把非对角项重复加入。
- 四面体重心坐标梯度符号与边方向约定不一致。
- 网格中存在质量很差的四面体，导致 `grad(lambda)` 被放大。
- VTU 点场平均前的单元顶点重构值已经异常。

优先检查 `FEMAssembler` 的局部矩阵循环范围、端口吸收项外积装配，以及 `tetraGradients` 是否满足
零阶 Whitney 边基函数的线积分约定。

## 构建成功但未使用 MKL

可能原因：

- CMake 没找到 oneMKL。
- `BPFEM_USE_MKL` 未开启。
- 构建目录缓存了旧配置。

优先检查 CMake 配置输出和 `BPFEM_USE_MKL` 状态。
