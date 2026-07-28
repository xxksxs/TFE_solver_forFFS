# 关键不变量

本文列出跨模块必须保持一致的约束。任何涉及解析、网格、FEM、求解器或后处理的修改都应检查这些不变量。

## faceId 必须匹配

AEDT 中的端口 faceId 必须能在 NGMesh 表面三角形中找到。

对于 AEDT 中采用 Objects(...) 定义的 sheet-object 波端口，必须满足：

- NGMesh 保留对应 object body 的完整六向包围盒；
- sheet 包围盒与一个且仅一个保留网格面的包围盒重合；
- PortFaceResolver 在建立 EdgeTopology 前已经解析出有效 faceId。

如果匹配失败：

- 端口模式可能为空。
- RHS 激励可能为空。
- 端口可能被误当作 PEC 边界。
- S 参数会失去物理意义。

## 边方向必须一致

以下模块必须共享同一套边方向约定：

- `EdgeTopology`
- `PortModeSolver`
- `FEMAssembler`
- `ResultExtractor`
- `OutputWriter`

任何符号错误都可能造成端口投影抵消、场方向错误或 S 参数异常。

当前约定是全局边方向按全局点号小到大。`EdgeTopology` 和 `PortModeSolver` 生成边引用时，应把
`localNodes` 也调整为这个方向；`FEMAssembler`、`ResultExtractor` 和 `OutputWriter` 不应再各自猜测
或重新定义局部边方向。

## 端口不能被 PEC 约束

波端口边界应参与端口吸收和激励装配，不能再被 PEC 零切向场约束。

检查点：

- 端口 faceId 是否被识别。
- 端口边是否从 PEC 集合中排除。
- PEC 约束是否在端口装配后正确应用。

## 矩阵和向量尺寸必须一致

以下对象尺寸必须等于全局边数量：

- SparseMatrix 行数。
- RHS 向量长度。
- 求解结果 field 长度。

如果新增自由度类型，应重新设计数据结构，而不是混入当前边 DOF 向量。

## 对称矩阵项只能装配一次

`SparseMatrix` 以对称上三角形式存储矩阵，并在矩阵-向量乘法和 MKL PARDISO 路径中按对称矩阵解释。

因此：

- 体单元局部矩阵只装配局部上三角。
- 波端口吸收项的模态外积只装配一次非对角项。
- 不要同时加入 `(i, j)` 和 `(j, i)`，否则非对角耦合会被加倍。

该错误可能在残差仍然很低时产生错误场图，因为求解器精确求解的是已经装错的矩阵。

## 局部基函数族必须线性独立

每个单元（四面体或端口三角形）内的局部基函数集合必须**线性独立**。这看似平凡，但一阶 hierarchical
H(curl) 元最常见的实现错误就是把两个 face bubble 写成代数上线性相关的形式，例如
`lambda_a * N_bc - lambda_b * N_ac`，它恒等于 `-lambda_c * N_ab`。

后果链：

- 局部 mass 矩阵秩亏。
- 端口面广义本征问题 `K x = λ M x` 中 `M` 失去正定性 → `LAPACKE_dsygv` 报 `info > N`。
- 体内局部矩阵秩亏，但 PEC 约束 + 直接求解可能掩盖错误，导致看似有结果但物理失真。

正确做法：

- `FaceFirst0 = lambda_c * N_ab`、`FaceFirst1 = lambda_a * N_bc`（或任意 `{lambda_x * N_yz}` 三选二）。
- `PortModeSolver` 与 `FEMAssembler` 中的 `evaluateBasis` 必须共享同一份基函数定义。

## 求积阶必须匹配基函数阶

数值积分阶不足会让局部 mass 在数值上偏离正定，对一阶 hierarchical 元尤其明显。约束如下：

- 端口三角面（2D）求积 ≥ degree 4（实现上使用 7 点 Wandzura, degree 5）。
- 四面体（3D）求积 ≥ degree 4（实现上使用 11 点 Keast）。

修改基函数阶或装配项时，必须同步检查 `PortModeSolver` 与 `FEMAssembler` 各自的求积规则。

## 元素缓存与频率因子分离必须保持一致

`FEMAssembler` 的体装配把每个四面体的局部矩阵分解为：

- 几何项 `K_e[a,b] = ∫_T (∇×N_a)·(∇×N_b) dV`、`M_e[a,b] = ∫_T N_a·N_b dV`，**只在第一次 `assemble()` 时按 11 点 Keast 求积一次**。
- 频率/材料因子 `invMu_e`、`εc_e(ω) = ε_r,e − jσ_e/(ωε_0)`，每频点重新计算。

这条不变量的含义：

- 任何修改 `evaluateBasis` 的 PR 必须假定 `K_e, M_e` 会被缓存——即基函数行为在第一次装配后不能再变。
- 若基函数依赖材料（如各向异性张量、Debye/Drude 多极），单一标量 `invMu / εc` 不足以表达 `A_local`，必须扩展 `ElementMatrices` 为多分量缓存（每个张量分量 / 极各一份），或暂时关闭缓存走原路径。
- `EdgeTopology` 改变（自适应 h-refine / 升 p）必须使 `elementCacheBuilt_ = false`，否则缓存与新拓扑不匹配。
- PEC 约束 dof 集合（`cachedConstrainedDofs_`）也按 lazy 缓存；端口集合或表面三角形改变时同样要 invalidate。

## 求解器输出语义必须一致

MKL PARDISO 和 BiCGSTAB fallback 都必须返回同一语义的 `SolveResult`。

下游模块只能依赖：

- `field`
- `iterations`
- `residual`

不能依赖求解器内部数据结构。

## 单位必须一致

频率、坐标、波数和材料参数必须使用一致单位体系。

常见风险：

- NGMesh 头部单位与实际坐标值解释不一致。
- 频率显示用 GHz，但内部计算必须用 Hz。
- 截止频率和传播常数对长度单位非常敏感。

## 端口装配和端口提取必须一致

波端口相关的三个环节不能单独修改：

1. `PortModeSolver` 计算端口模式（NPM 数值本征求解，或 `bpfem::apm::AnalyticPortBuilder` 通过 `setPrecomputed` 注入解析单模，或 `bpfem::tfe::TransfiniteElementBuilder` 通过 `setMultiMode` 注入解析多模）。
2. `FEMAssembler` 使用端口模式装配边界项和 RHS。多模端口（`PortModeSolver::multiMode(faceId) != nullptr`）会循环每个模贡献一个秩-1 端口算子。
3. `ResultExtractor` 使用端口模式提取 S 参数，仍读取 `PortModeSolver::solve(faceId)` 返回的主模 `couplingWeights`（多模缓存里 `excitationModeIndex` 对应的那个）。

如果只改其中一个环节，S 参数可能看似有数值但无物理意义。

端口模式的边自由度和弱式投影权重不是同一个对象。`PortModeSolver` 求得本征场 DOF 后，应通过
端口质量矩阵得到投影权重 `M_port * e_mode`。`FEMAssembler` 和 `ResultExtractor` 必须使用同一组
投影权重，否则端口附近场分布和 S 参数都会失真。

无论是 NPM、APM、TFE 三条路径，都必须输出**相同的 `PortMode` 数据契约**：

- `couplingWeights[i] = ∫ N_i · e_t dS = M·dofs`（FE 空间下的投影）。
- `quadrature.modeFieldValue` 用 FE 重构的 `e_t^FE`，**不是**解析 `e_t^analytic`。
- `e_t^FE` 在 FE 空间下 L2 归一：`dofs^T M dofs = 1`。

APM 与 TFE 路径必须装配端口面 mass 矩阵并解出 `dofs = M⁻¹ b`，再做 L2 重归一；不可直接把解析投影 `b = ∫ N · e_t^analytic dS` 当作 `couplingWeights` 而不归一，否则 `quadrature` 与 `couplingWeights` 的归一化口径会失配，引发非物理 |S21| > 1。

TFE 多模路径还要遵守 *AffineSystem 多模展开规则*：每个 (faceId, modeIndex) 在 `affine_.portCoupling`、`portFaceIds`、`portCutoffSquared`、`projectPortIndex`、`isExcitationMode` 五个并行数组中各占一项；ALPS 用虚拟端口索引装配 ROM、用 `dominantVirtualByProject` 还原物理端口主模做 S 参数提取。修改任一处必须同步更新另外四处。

## 高残差结果不可信

如果求解残差很高：

- 场图不应被当作物理结果。
- S 参数不应作为验证依据。
- 应优先检查求解器、矩阵条件数、边界条件和 RHS。
