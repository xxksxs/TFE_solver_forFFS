# AWE、MGAWE 与 WCAWE 公式复现

在支持的无损 FEM 工程上，三种方法使用局部坐标 t=(k−k₀)/k₀。
对于 MGAWE 和 WCAWE，这也是论文变量 (s−s₀)/s₀（s=j2πf）。
缩放后的第 r 阶系数是原坐标系数乘以相应尺度的 r 次幂。
GAWE 和 ALPS 保留原有展开变量和在线模型。

| 论文 | 公式 | C++ 实现 |
|---|---|---|
| Reddy、Deshpande，NASA/CR-97-206261 | (21)–(22) Taylor 系数和场矩 | PolynomialPortMomentBuilder::buildTaylorSeries、PolynomialMomentRecurrence::generatePolynomialMoments |
| 同上 | (24)–(26) Padé 线性方程与求值 | PadeApproximant::build；AweSweep 的 --awe-field-pade 逐场分量求值后再提取 S 参数 |
| Slone 等，2001 MGAWE，第 III、V 节；2002 MGAWE，式 (9)、(12)–(16) | 局部矩、统一空间和在线投影。2001 使用共轭转置，2002 使用普通转置 | MgaweSweep::build、GalerkinReducedModel::build 和 solveReduced |
| 同上 | (20)–(24) 局部空间耗尽；(19)、(25)–(27) 残差驱动选点 | --mgawe-adaptive，MgaweSweep.cpp 中的 localSpaceExhausted 和 nextExpansion |
| Slone、Lee、Lee，2003 WCAWE | (7)–(9) 校正递推和 MGS | PolynomialWcaweRecurrence、UpperTriangularBlockAction、WellConditionedBasisBuilder |
| 同上 | (2)、(4)–(5)、Algorithm 1 | WcaweSweep::projectPolynomial 预投影每阶矩阵和右端项，solveReduced 用 Horner 求各频率的降阶系统 |

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
  首点按式 (18) 取扫频序列中间索引，偶数频点时取两个中点中靠后的一个。后续选点按式 (19)、(27)
  对未收敛子带的端点索引加权并向上取整，频带边界取实际首末点。每次扩基后重新检查已接受频点，
  因为程序用最终统一 ROM 输出整段扫频。
- WCAWE 以最高 q−1 阶系数生成 q 个矩所需的基。在线多项式次数至少为 2，以精确保留体积分中的 k² 项。
  构基回代使用用户要求的容差和 1e−10 中较严的一项。回代失败或递推残差超标时，本次运行报错。
  --wcawe-check-moments 用同一 A₀ 分解额外回代独立的传统矩，并检验其到基的相对投影残差。
  此选项会增加回代次数和内存。full_moment_projection_error 是该检查值。
  旧 wcawe_moment_matching_error 字段仍表示内部递推及基关系误差。
- diagnostics.json 的 max_passivity_error 计算已输出端口功率超出 1 的正值，
  因此适用于只报告部分端口的单入射归一化情形。完整能量守恒需要所有传播端口的数据。

## 工程模式的真实算子残差

- `--sweep mgawe --mgawe-adaptive-true-residual` 启用工程原型。局部基来自端口 Taylor 矩，在线 ROM 直接投影求值无损 FEM 的真实端口频率律。真实残差按 $\|b(f)-A(f)Vz(f)\|_\infty/\|b(f)\|_\infty$ 计算。Taylor 残差使用同一降阶场和最近的同截止分支多项式，各自按其右端项归一化。残差采样只做稀疏矩阵向量乘法和端口秩一作用。
- `--mgawe-true-residual-tolerance <值>` 接受有限正数，默认 `1e-3`。`--mgawe-points` 限制全阶展开点数，`--mgawe-order` 限制每点矩数。首点取靠近中点的可用网格点，之后选择真实相对残差最高且尚未展开的可用频点。最终按 `stopping_reason`、`unmet_frequencies_hz` 和 `residual_samples` 检查达标情况。达到预算或找不到新候选时仍输出当前近似，进程退出码为 0，`true_residual_acceptance` 为 `unmet`。阈值判断的是方程残差，S 参数精度通过同网格 Direct 单独验收。
- 任一保留端口模满足 $|k^2-k_c^2|\le10^{-8}\max(k^2,k_c^2)$ 时，该频点标为 `near_cutoff`，从 Taylor 展开候选中排除。工程模式在线使用与 Direct 相同的分支和截止约定。缺少同分支多项式或频点处在此保护区时，Taylor 残差为 `null`。
- 真右端项为零时记录 `zero_rhs`，相对残差为 `null`，并用零场和零值占位 S 参数输出，避免把无入射波的频点计入验收。全频带右端项为零或所有非零右端项频点都落在截止保护区时，分别给出对应错误，进程退出码为 1。仍需采用现有无损波端口和至少双端口设置。
- `--sweep wcawe --wcawe-residual-diagnostics` 为论文 WCAWE 的同一个降阶场增加逐点真实和 Taylor 残差。WCAWE 继续使用单展开点、现有递推和在线多项式，跨截止需分段。展开点右端项为零时在稀疏分解前给出错误。诊断开启后会增加残差评估的时间和内存。

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

同网格评估降阶误差时，将 Direct、MGAWE、WCAWE 输出写到同一个结果根目录，各策略生成自己的子目录。运行以下命令得到复数 S11/S21 的相对误差曲线、逐频点 CSV 和 JSON 汇总。逐点指标为 `|S_ROM - S_Direct| / |S_Direct|`。参考值为零时只报告绝对误差。

```powershell
node scripts/compare_ios_fastsweeps.mjs --direct --results-root result_paper
```

见 [MGAWE 与 WCAWE 审查及验证记录](../../changes/2026-11-20-mgawe-wcawe-audit.md)。
