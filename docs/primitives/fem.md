# FEM 原语

本文说明工程中的 H(curl) 棱元概念。修改 `EdgeTopology`、`PortModeSolver` 或 `FEMAssembler` 前应先阅读本文。

## 全局边自由度

本工程默认使用零阶 Nedelec/Whitney 棱元。未知量不是节点标量，而是网格边上的切向电场自由度。命令行参数 `--basis-order 0|1` 可选择零阶或一阶层次棱元。

关键点：

- 零阶棱元每条唯一网格边对应一个全局自由度。
- 零阶每个四面体有 6 个局部基函数。
- 一阶层次棱元每个四面体有 20 个局部基函数，其中包含零阶的 6 个边基函数。
- 全局边方向统一为全局点号小到大。
- 局部边引用在存储时即调整为全局点号小到大，体装配和端口模式都使用同一方向。
- 后处理不能把边自由度当作节点标量使用。

## 局部边引用

局部边引用记录：

- 当前四面体内已经按全局点号小到大排序的两个局部节点编号。
- 对应全局边编号。
- 方向符号字段仍保留在数据结构中，但边引用应在生成时定向到全局方向，因此零阶和一阶边引用通常使用 `+1`。

体装配、端口模式投影和场重构必须共享这套存储约定。不要只在某一个模块中交换局部端点。

## Nedelec/Whitney 基函数

零阶和一阶层次棱元基函数都属于 H(curl) 空间，适合求解频域 Maxwell 方程。

它的作用：

- 保证切向场连续。
- 支持 curl-curl 弱式装配。
- 使 PEC 边界可通过切向边自由度约束实现。

四面体重心坐标梯度必须与边方向约定一致。对已经按全局点号小到大存储的局部边 `a -> b`，零阶 Whitney 基函数
`lambda_a * grad(lambda_b) - lambda_b * grad(lambda_a)` 沿该局部边的线积分应为 `+1`。如果
`tetraGradients` 的符号反了，体单元场重构会和端口面投影的方向语义不一致。

一阶边基函数使用 `(lambda_a - lambda_b) * N_ab`。交换端点时，`lambda` 差和零阶边基函数都会变号，
因此一阶边基函数本身对端点交换不变；仍建议存储为全局点号小到大，以便和零阶边、端口投影保持同一约定。

一阶面基函数（`LocalDofKind::FaceFirst0` 与 `FaceFirst1`）使用如下两个线性独立的 face bubble：

- `FaceFirst0`：`lambda_c * N_ab`
- `FaceFirst1`：`lambda_a * N_bc`

注意 `lambda_c * N_ab + lambda_a * N_bc + lambda_b * N_ca ≡ 0` 是恒等关系，因此从 `{lambda_c * N_ab,
lambda_a * N_bc, lambda_b * N_ca}` 中**任取两个即线性独立**，但任何形如 `lambda_a * N_bc - lambda_b *
N_ac` 的组合在代数上等于 `-lambda_c * N_ab`，与 `FaceFirst0` 完全相反，会让端口面 mass 矩阵秩亏，
触发 `LAPACKE_dsygv info > N` 一类的 Cholesky 失败。`PortModeSolver::evaluateBasis` 与
`FEMAssembler::evaluateBasis` 必须采用同一套面基函数定义。

## 体装配项

FEM 主矩阵通常包含：

- **curl-curl 项**：代表磁场能量或旋度项。
- **质量/材料项**：代表介电响应和损耗。

装配时需要同时考虑材料参数、单元体积、局部边方向、基函数阶数和复数频域系数。一阶层次基函数使用四面体数值积分装配局部矩阵。

数值积分阶必须匹配被积函数的多项式阶：

- 零阶 Whitney 基：被积函数 ≤ degree 2，体单元采用 11 点 Keast (degree 4) 已余量充足。
- 一阶 hierarchical 基：mass 项被积函数最高 degree 4，curl-curl 项最高 degree 2，体单元 11 点 Keast 仍然适用。
- 端口三角面（2D）：一阶 hierarchical mass 被积达到 degree 4，**必须**采用至少 degree 4 的求积规则；当前实现使用 7 点 Wandzura (degree 5)。曾经的 3 点 (degree 2) 规则会让端口面 mass 矩阵失去正定性，进而让 `LAPACKE_dsygv` 失败。

### 几何与频率因子分离（element cache）

每个四面体的局部矩阵在频域内可以分解为两个**与频率无关**的实矩阵：

```
K_e[a,b] = ∫_T (∇×N_a)·(∇×N_b) dV     （仅依赖几何 + 基函数）
M_e[a,b] = ∫_T  N_a·N_b dV             （仅依赖几何 + 基函数）
```

体装配在每个频率点退化为标量组合：

```
A_local(a,b)(ω) = invMu_e · K_e[a,b] − k0(ω)² · εc_e(ω) · M_e[a,b]
```

其中 `invMu_e = 1/μ_r,e`、`εc_e(ω) = ε_r,e − jσ_e/(ωε_0)` 是按元素缓存的材料。`FEMAssembler::ensureElementCache()` 在第一次 `assemble()` 时一次性预计算 `(K_e, M_e)` + 材料数据，按 upper-triangle packed 存放。这把内层求积的 `Nq × dof²` 操作转为 `dof²` 标量 FMA，扫频时单点装配代价跨频点摊销后约 6× ~ 8× 加速。

修改材料模型（例如引入色散）或基函数族时，必须扩展 `ElementMatrices`：色散材料的 `εc(ω)` 已经按 ω 在 `assemble()` 中计算，但若引入张量 ε / μ，则 `K_e, M_e` 不再能用单个标量缩放，需要拆分为多分量缓存（参考 `optimization/materials.md`）。

当前 `SparseMatrix` 只存储对称上三角 CSR。体单元局部矩阵和波端口边界项都是对称贡献，因此
装配时每个非对角矩阵项只能写入一次。不要用完整双循环把 `(i, j)` 和 `(j, i)` 都加入
`SparseMatrix::add`，否则上三角半存储会把非对角耦合项加倍，导致场图出现非物理的局部振荡。

## PEC 约束

PEC 边界通过把非端口边界上的切向电场自由度设为 0 实现。

注意：

- 端口边不能被 PEC 约束。
- PEC 约束的是切向电场；法向电场可由表面电荷支持，不要求完整 `|E|` 为 0。
- `FEMAssembler` 先收集 PEC 约束 DOF mask，体装配和端口装配跳过这些 DOF。
- `SparseMatrix::imposeZeroDirichlet` 在一次矩阵遍历后设置所有零 Dirichlet 对角项，避免逐约束边扫描全矩阵。

## 波端口边界项

波端口边界项由数值端口模式提供。

端口模式计算规则：

- 从三维端口 face 上抽取表面三角网格。
- 端口截面外轮廓边对应波导金属壁，应在二维端口本征问题中施加切向零场约束。
- 在端口二维网格上建立并求解端口本征问题。
- 把端口二维模式按端口三角单元 H(curl) 基函数投影回三维全局 DOF。
- 端口边界项和 RHS 使用质量矩阵投影权重 `M_port * e_mode`，而不是直接使用本征向量裸 DOF。

端口三角单元局部基函数：

- 零阶端口单元使用 3 个边基函数。
- 一阶端口单元使用 6 个边基函数和 2 个面基函数。两个面基函数必须线性独立（见上文 *Nedelec/Whitney 基函数* 一节中关于 `FaceFirst0 / FaceFirst1` 的公式）。
- 一阶端口投影会同时写入三维零阶边 DOF、一阶边 DOF 和端口面 DOF，避免高阶自由度在端口激励和 S 参数提取中被遗漏。
- 端口面 H(curl) 本征问题的求积规则必须 ≥ degree 4（实现上为 7 点 Wandzura, degree 5）。

涉及模块：

- `PortModeSolver`：计算端口面本征模式和边 DOF。
- `FEMAssembler`：使用端口模式装配吸收项和激励 RHS。
- `ResultExtractor`：使用同一端口模式做 S 参数投影。

端口相关修改必须跨这三个模块一起检查。
