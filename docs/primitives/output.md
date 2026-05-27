# 输出原语

本文说明后处理阶段的核心输出概念。修改 `ResultExtractor` 或 `OutputWriter` 前应先阅读本文。

## 边场向量

求解器输出的 `field` 是复数边自由度向量，不是节点场。

含义：

- 每个值对应一条全局边的切向积分自由度。
- 输出场可视化需要从边自由度重构到点或单元。
- 后处理时必须使用 `EdgeTopology` 提供的边和单元关系。

## 端口投影

S 参数通过把求解场投影到端口模式得到。

涉及一致性：

- 使用与装配相同的端口 faceId。
- 使用与 `PortModeSolver` 相同的边方向约定。
- 零阶端口三角单元使用 3 个边基函数。
- 一阶端口三角单元使用 6 个边基函数和 2 个面基函数。
- 一阶层次棱元下，端口模式应同时包含零阶边 DOF、一阶边 DOF 和端口面 DOF。
- L2 模态归一化保证 `e_t` 的形状（`∫|e_t|² dS = 1`）；坡印廷功率归一化保证不同频率/不同截面端口下入射与出射在 √W 量纲一致。
- `ResultExtractor` 应对每个端口用其自己的 `powerNormalizationFactor` 折算回 √W 后再做 S 参数比值。
- S 参数投影应使用与端口 RHS 相同的质量矩阵投影权重 `M_port * e_mode`。不要用裸模态 DOF 做欧氏点积投影。

## CSV 输出

`OutputWriter` 写出 `s_parameters.csv`。

典型列包括：

- `freq_Hz`
- `S11_real`
- `S11_imag`
- `S11_dB`
- `S21_real`
- `S21_imag`
- `S21_dB`

CSV 输出应保证每个求解频点对应一行。

## VTU 输出

VTU 文件用于可视化电场分布。

常见文件：

- `field_last.vtu`：最后一个频点的场。
- `field_<frequency>Hz.vtu`：开启全频点输出时的逐频点场文件。

### 单元类型与高阶节点

为了忠实表现一阶 hierarchical 棱元在边和面上的非平凡分布，VTU 输出在线性 `VTK_TETRA`（cell type 10）+ 顶点平均之外提供两条不同的高阶路径，由 CMake 选项 `BPFEM_OUTPUT_PARAVIEW` 切换：

#### 默认（VisIt-friendly，`BPFEM_OUTPUT_PARAVIEW=OFF`）

- **`--field-output-order 1`**：每个四面体 4 顶点，cell type **10 (VTK_TETRA)**。零阶基函数下默认。
- **`--field-output-order 2`**：每个四面体 10 节点（4 顶点 + 6 边中点），cell type **24 (VTK_QUADRATIC_TETRA)**。一阶基函数下默认。
- `--field-output-order 3` 在该路径下会运行时报错，需要切换到 ParaView 路径。
- VisIt ≥ 2.x 与 ParaView ≥ 4.x 都能直接渲染。

#### 可选（ParaView 高保真，`BPFEM_OUTPUT_PARAVIEW=ON`）

- **`--field-output-order 1`**：4 顶点，cell type **10**。
- **`--field-output-order 2`**：10 节点，cell type **71 (VTK_LAGRANGE_TETRAHEDRON)** + `HigherOrderDegrees=(2,2,2)`。可看见 EdgeFirst bubble 在边方向上的分布。
- **`--field-output-order 3`**：20 节点（再加 4 面心 + 0 体内），cell type **71** + `HigherOrderDegrees=(3,3,3)`。可看见 FaceFirst bubble 的面贡献。
- 需要 ParaView ≥ 5.5 / VisIt ≥ 3.0 的 Lagrange 元渲染支持；老版 VisIt 不识别 cell type 71 / 会报错或拒绝加载。

#### 通用约束（两条路径共有）

- **节点共享（默认）**：H(curl) Nédélec 棱元在 cell 边界上切向连续；同一条几何边或同一个三角面在邻接 cell 中给出的重构值在切向上严格相等、在法向上仅 O(ε) 量级差异。`OutputWriter` 默认按几何实体共享：
  - 角点 → 共享于网格顶点
  - 边中点（order=2 / order=3）→ 共享于几何边
  - 面中心点（order=3）→ 共享于几何三角面
  - 体内点（order=3）→ cell-local（未共享，几何位置不重合）
  邻接 cell 给出的多个值用算术平均合并。在 BP filter 例题（48,710 tet，basis-1 + order=2 输出）下，从 *cell-local* 的 452,080 节点降到 *共享* 的 76,166 节点（缩减 5.9×）；VTU 体积从 ~149 MB 降到 ~28 MB。
- VTU 体积仍随 p_out 大致线性增长：order=2 ≈ order=1 的 ~7×，order=3 又再 ~3× （相对 order=2）。
- 如需缩小 VTU 文件用于 long-running 扫频时的逐频点输出，搭配 `--no-write-all-fields` 只保留 `field_last.vtu`。

### 可视化字段

VTU 中只输出**两个矢量场，各分实部和虚部**，共 4 个 `DataArray`：

- **`E_real`**：电场实部，单位 V/m。
- **`E_imag`**：电场虚部，单位 V/m。
- **`H_real`**：磁场实部，单位 A/m。由 Faraday 定律按 `H = j ∇×E / (ω μ_r μ_0)` 在每个 Lagrange 节点重构得到（时谐约定 `e^{+jωt}`）。
- **`H_imag`**：磁场虚部，单位 A/m。

注意：

- VTU 中的场来自 H(curl) 棱元重构，不是直接节点未知量。
- 复频域解需要同时看实部和虚部；只看 `E_real`/`H_real` 任一份会随相位旋转出现"局部为零"。完整瞬时场是 `Re{(E_real + j·E_imag) · e^{+jωt}}`。
- 模值、幅值、对数刻度、PEC 切向分量等派生量不再写入 VTU。可在 ParaView 中用 *Calculator* / *Python Calculator* / *Threshold* 滤镜按需计算，例如：
  - `mag(E_real)` → 实部矢量长度
  - `sqrt(E_real.E_real + E_imag.E_imag)` → 复模 |E|
  - `abs(E_real_n) + abs(E_imag_n)`（其中 `E_n` 为指定法向投影）→ 切向幅值，用于 PEC 检查
- H 场依赖材料（`μ_r`），跨材料界面会出现物理跳变。VTU 中每个 cell 独立持有 Lagrange 节点，这种跳变会被忠实绘出，**不是可视化 bug**。
- 如果场看起来不传播，应先检查求解残差和端口激励，而不是只修改可视化。

### 何时升 `--field-output-order`

- 默认按 `basis-order + 1` 自动选择，零阶 → p_out=1，一阶 → p_out=2。
- 仅当需要在场图中验证一阶 face bubble（FaceFirst0 / FaceFirst1）的面贡献是否合理时才显式设置 `--field-output-order 3`。
- p_out 与求解器基函数阶数无关：`--field-output-order` 改变的只是 VTU 抽样密度，不会改变 S 参数 CSV。

## 输出修改清单

- 修改 S 参数前，检查 `PortModeSolver` 和 `FEMAssembler` 的端口约定。
- 修改 VTU 前，确认边元重构公式和单位。
- 添加新输出列时，同步更新测试和文档。
- 输出异常时先确认求解残差是否足够低。
