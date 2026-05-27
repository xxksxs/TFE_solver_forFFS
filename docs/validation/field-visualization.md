# 场可视化验证

VTU 场文件用于定性检查电场和磁场分布，但必须结合求解残差和端口激励判断。

## 检查内容

打开 VTU 后检查：

- `E_real`、`E_imag` 在驱动端口外不应恒为零。
- `H_real`、`H_imag` 应与 `E` 在主模区域呈现 `H ⊥ E ⊥ k` 的几何关系。
- 场型应符合波导传播直觉。
- 强不连续应能由几何或边界条件解释。
- 端口面应呈现激励/吸收行为，而不是 PEC 零场行为。
- 检查 PEC 墙面是否干净时，在 ParaView 中用 Calculator 计算 `E_real - (E_real . n) * n` 的模长（n 为表面法向），不再以单独的 DataArray 提供。

## 注意事项

- VTU 输出来自边自由度重构，不是节点标量未知量。
- 一阶 hierarchical 路径下默认 `--field-output-order 2`：每个四面体输出 10 个 Lagrange 节点（VTK type 71），ParaView/VisIt 按 p=2 多项式渲染。EdgeFirst bubble 在边中点处的非平凡分布可以直接看见。
- 当需要核对 FaceFirst bubble 的面贡献时显式 `--field-output-order 3`，每个四面体输出 20 节点；VTU 体积约比 p_out=1 大 5×，仅按需启用。
- 同一基函数阶下不同 `--field-output-order` 输出的 S 参数 CSV 完全相同；该选项**只**改变 VTU 抽样密度。
- 不同单元的 Lagrange 节点**不共享**，因此在 ParaView 中默认会看到 cell 边界处的小幅跳变。这是 H(curl) 棱元的物理特性（切向连续，法向允许跳变），不是可视化 bug。
- 复频域场必须同时看实部与虚部；任一时刻的瞬时场为 `Re{(real + j·imag) · e^{+jωt}}`，单看实部会被相位旋转误导。
- H 场依赖材料 `μ_r`，跨材料界面跳变是物理结果。
- ParaView ≥ 5.5 / VisIt ≥ 3.0 才完整支持 VTK Lagrange Tet (cell type 71)；旧版本可能把高阶节点当作"游离点"绘制成散点。
- 如果场看起来不传播，应先检查残差、端口激励和 PEC 约束。
- 可视化异常不一定说明输出代码错误，也可能是上游求解失败。

## 相关文档

- `../primitives/output.md`
- `../primitives/invariants.md`
- `s-parameters.md`
