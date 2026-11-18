# AWE、MGAWE 与 WCAWE 公式复现

在支持的无损 FEM 工程上，三种方法使用局部坐标 t=(k−k₀)/k₀。
对于 MGAWE 和 WCAWE，这也是论文变量 (s−s₀)/s₀（s=j2πf）。
缩放后的第 r 阶系数是原坐标系数乘以相应尺度的 r 次幂。
GAWE 和 ALPS 保留原有展开变量和在线模型。

| 论文 | 公式 | C++ 实现 |
|---|---|---|
| Reddy、Deshpande，NASA/CR-97-206261 | (21)–(22) Taylor 系数和场矩 | PolynomialPortMomentBuilder::buildTaylorSeries、PolynomialMomentRecurrence::generatePolynomialMoments |
| 同上 | (24)–(26) Padé 线性方程与求值 | PadeApproximant::build；AweSweep 的 --awe-field-pade 逐场分量求值后再提取 S 参数 |
| Slone 等，2001 MGAWE，第 III、V 节；2002 MGAWE，式 (9)、(12)–(16) | 局部矩、统一空间及复对称投影 | MgaweSweep::build、GalerkinReducedModel::build 和 solveReduced |
| 同上 | (20)–(24) 局部空间耗尽；(19)、(25)–(27) 残差驱动选点 | --mgawe-adaptive，MgaweSweep.cpp 中的 localSpaceExhausted 和 nextExpansion |
| Slone、Lee、Lee，2003 WCAWE | (7)–(9) 校正递推和 MGS | PolynomialWcaweRecurrence、UpperTriangularBlockAction、WellConditionedBasisBuilder |
| 同上 | (2)、(4)–(5)、Algorithm 1 | WcaweSweep 将同一组高阶矩阵及 RHS 系数用于离线构基和在线 Galerkin 求解 |

参考原文：[AWE](https://ntrs.nasa.gov/api/citations/19980006278/downloads/19980006278.pdf) · [MGAWE 2002](https://bpb-us-w2.wpmucdn.com/u.osu.edu/dist/3/47088/files/2017/06/slon02a-1m1cnmt.pdf) · [WCAWE 2003](https://bpb-us-w2.wpmucdn.com/u.osu.edu/dist/3/47088/files/2017/06/slon03b-2med5w4.pdf)。
[MGAWE 2001 作者公开全文](https://www.researchgate.net/publication/3011305_Multipoint_Galerkin_asymptotic_waveform_evaluation_for_model_order_reduction_of_frequency_domain_FEM_electromagnetic_radiation_problems)。

2001 年原文从约化解的 Taylor 残差取基。设残差定义为 b−Ax_rom。
若前 n−1 阶已经匹配，则第 n 阶满足 A₀⁻¹rₙ=xₙ−x_rom,n。
后一项已在当前全局子空间内，所以对原始 xₙ 和校正残差候选做同样的全局 MGS，会得到相同的新方向。
这给出当前局部原始矩加全局正交化与原文残差构基之间的等价关系，前提是低阶矩匹配及约化 A₀ 可逆。
2001 版本的在线阶段使用最后一个展开点的多项式，2002 版本使用最近的展开点。

## 数学及参数约定

- Taylor 模型保留 λ=k² 的多项式系数，并以标量 Taylor 级数近似非多项式端口频率律。
  在每个传播支上，矩阵系数来自 sqrt(λ−kc²)，激励来自 λ^(1/4)(λ−kc²)^(1/4)。
  低于截止时保留 Direct 的既有导纳 jk₀ 和零入射激励。
- 展开点恰在截止处会报错。MGAWE 在与目标频点相同的截止分支上选最近的已有展开多项式。
  AWE 逐场路径和单点 WCAWE 在跨截止的频带上需要分段运行。
- AWE 默认保留端口投影 Padé 用法。--awe-field-pade 使用 [q−1/q] 构造场分量 Padé，
  在出现奇异子问题时报告自由度索引。恒零场分量直接求值为零。
- MGAWE 固定模式使用均匀展开点。--mgawe-2001 按 2001 年论文使用最后一个展开点的多项式。
  可重复指定 --mgawe-expansion <Hz> <order>，逐点给出不同阶数。
  --mgawe-adaptive 将 --mgawe-points 和 --mgawe-order
  当作计算上限，局部容差为 1e−6，扫频残差容差为 1e−2。日志与 diagnostics.json
  记录停止原因、实际局部阶数和最大全阶多项式相对残差。触及计算上限时，结果仍可供误差分析。
- WCAWE 以 q−1 阶系数生成 q 个矩所需的基，模型次数至少为 1。
  --wcawe-check-moments 用同一 A₀ 分解额外回代独立的传统矩，并检验其到基的相对投影残差。
  此选项会增加回代次数和内存。full_moment_projection_error 是该检查值。
  旧 wcawe_moment_matching_error 字段仍表示内部递推及基关系误差。
- diagnostics.json 的 max_passivity_error 计算已输出端口功率超出 1 的正值，
  因此适用于只报告部分端口的单入射归一化情形。完整能量守恒需要所有传播端口的数据。

## 命令

以下示例以项目根目录为当前目录。先构建 oneMKL/PARDISO 的 Release 可执行文件。

```powershell
$solver = '.\build_pardiso\Release\bp_fem_solver.exe'
$common = @('--aedt','IOStructure.aedt','--mesh','current.ngmesh','--basis-order','0',
            '--port-method','numerical','--linear-solver','direct','--max-sweep-points','101',
            '--no-write-all-fields','--out','result_paper')
& $solver @common --sweep direct
& $solver @common --sweep awe --awe-field-pade --awe-order 8
& $solver @common --sweep mgawe --mgawe-adaptive --mgawe-points 8 --mgawe-order 8
& $solver @common --sweep wcawe --wcawe-order 12 --wcawe-check-moments
```
