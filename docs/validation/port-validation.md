# 端口验证

端口验证用于确认波端口定义、端口模式和激励装配是否一致。

## faceId 一致性

确认每个 AEDT 波端口 faceId 都存在于 NGMesh 表面三角形中。

期望：

- 每个 `PortDefinition::faceId` 至少对应一个 `SurfaceTriangle`。
- 端口 face 上的边不会被 PEC 约束清零。

## 模态合理性

对每个端口模式检查：

- 边 DOF 映射非空。
- 端口截面外轮廓边应被视为 PEC 边界，不参与端口本征模自由度。
- 零阶端口三角单元应使用 3 个边基函数。
- `--basis-order 1` 时，端口三角单元应使用 6 个边基函数和 2 个面基函数。
- `--basis-order 1` 时，两个面基函数必须线性独立，参考 `../primitives/fem.md` 中 `FaceFirst0 / FaceFirst1` 的公式。
- `--basis-order 1` 时，端口模式应包含一阶边 DOF 和端口面 DOF，而不是只包含零阶边 DOF。
- 端口面 H(curl) 本征问题求积规则应至少为 degree 4（实现上为 7 点 Wandzura, degree 5）。
- 范数有限且非零。
- 截止波数与端口尺寸物理上相符。
- 频率高于截止时传播常数为实数。

## 激励合理性

对激励端口检查：

- 入射幅度使用配置的幅度和相位。
- 端口边上的 RHS 非零。
- RHS 和吸收项应使用端口质量矩阵投影权重 `M_port * e_mode`，不能直接使用本征向量裸 DOF。
- 一阶层次棱元下，一阶边 DOF 和端口面 DOF 的 RHS/吸收项也应参与端口边界装配。
- 所有端口都装配吸收边界项，而不仅仅是激励端口。

## 跨模块一致性

以下模块必须共享同一端口约定：

- `../../src/fem/PortModeSolver.cpp`
- `../../src/fem/FEMAssembler.cpp`
- `../../src/post/ResultExtractor.cpp`
