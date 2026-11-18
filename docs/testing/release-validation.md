# 版本发布验收规范

在隔离分支合并待发布改动，再从该分支重新构建。归档 Git 提交、编译器、依赖版本、构建选项、算例及网格、线程数和每轮命令。每组运行使用独立 `--out`，保留 `run.json`、`timing.json`、`s_parameters.csv` 及适用时的模态 CSV。验收结果记在 [docs/changes](../changes/README.md)。日常改动使用[回归清单](regression-checklist.md)。

## 验收矩阵

| 项目 | 方法与门槛 |
| --- | --- |
| Release 构建 | 分别验证应用启用 MKL、启用 MUMPS、同时启用及同时关闭四种组合。每种执行 CTest。若 ZMUMPS 使用 MKL BLAS，应用关闭 MKL 仅表示禁用应用自身 PARDISO 与 MKL SpMV。 |
| 基本物理与 CLI | 用 20 mm × 12 mm × 20 mm 的单端激励均匀波导，在 9–10 GHz 五点跑基阶 0/1。确认 TE10 截止与解析值约 7.495 GHz 相符，并检查 `abs(S11) ≤ 0.02`、`abs(S21) ≥ 0.99` 及主模功率守恒。完整相位、谱和模式验证由 `tfe_tests` 覆盖。 |
| 直接后端与迭代校正 | 同一网格、频率和端口模型，以双精度 PARDISO 为基准。要求实际全阶迭代解的真实相对残差不超过设定容差。每个复数 S11/S21 的绝对差应满足 `abs(S − S_ref) ≤ 1e-5 + 1e-3 × abs(S_ref)`。小振幅通道同时报告绝对差。 |
| 降阶精度 | IOStructure 使用同端口模型的 Direct 比较复数 S11 和 S21，分别计算全带相对 L₂。已选发布配置两路均应 ≤1%。另外报告全带各通道最大 dB 偏差及发生频率。 |
| 深谷 | 独立 Direct 在 93.960–94.030 GHz 加密到 5 MHz。比较原生 TFE 与 MGAWE 的采样谷频和谷深。采样谷频差应 ≤5 MHz，谷深差应 ≤0.1 dB。启用精确区间时，再验证区间内请求频点直接求解并达到真实残差要求。 |
| 失败保护 | 强制迭代一步且容差 `1e-12`。应退出 1，`run.json` 为 `fatal`，保留 `iteration_limit` 与真实残差，且该回合无 S 参数 CSV。 |
| 工程成本 | 同一算例上统一线程和输出选项，比较总耗时和峰值工作集。涉及运行时波动的代表路径至少三轮，报告中位数。外存除精度外报告节省的内存及附加耗时。仅单次测量的配置注明次数。 |

复数相对 L₂ 定义为 `sqrt(Σ|S − S_ref|² / Σ|S_ref|²)`。dB 差用幅值的 20 倍常用对数。HFSS CSV 的幅值对照在共同频点进行。对照结果要注明网格、端口编号、参考面及归一化，确实未核对的条件写入待验证项。

## 可复跑命令

从仓库根运行。构建命令与 ZMUMPS 依赖见[构建检查](build-checks.md)。无 MKL 的波导对照可用 `--linear-solver gmres --precon ilu0 --gmres-restart 60 --max-iterations 1500`。串行 ZMUMPS 构建可替换为 `--linear-solver mumps --mumps-mode ooc`。夹具 `tests/fixtures/tfe-waveguide.ngmesh` 的 `coords` 以米存储，按 AEDT 导出的 NGMesh 约定处理，`user_unit_name` 保存为模型显示单位。

```powershell
$solver = '.\build_combined\Release\bp_fem_solver.exe'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$wave = @('--aedt', 'tests/fixtures/tfe-waveguide.aedt',
          '--mesh', 'tests/fixtures/tfe-waveguide.ngmesh', '--port-method', 'tfe-block',
          '--tfe-modes-per-port', '2', '--no-write-all-fields')
& $solver @wave --basis-order 1 --sweep direct --linear-solver direct --out result/wave-ref
& $solver @wave --basis-order 1 --sweep direct --linear-solver gmres --precon pardiso-single --tolerance 1e-7 --out result/wave-mixed
& $solver @wave --basis-order 0 --sweep direct --linear-solver gmres --max-iterations 1 --tolerance 1e-12 --out result/wave-failure
```

```powershell
$env:MKL_NUM_THREADS = '4'
$env:OMP_NUM_THREADS = '4'
$io = @('--aedt', 'IOStructure.aedt', '--mesh', 'current.ngmesh', '--basis-order', '0',
        '--port-method', 'tfe-block', '--tfe-modes-per-port', '3', '--no-write-all-fields')
$rom = @('--sweep', 'mgawe', '--tfe-subband-width', '2e9', '--mgawe-points', '1')
& $solver @io --max-sweep-points 81 --sweep direct --linear-solver direct --out result/io-ref
& $solver @io --max-sweep-points 81 @rom --tfe-exact-band 93.96e9 94.03e9 --linear-solver direct --out result/io-rom
& $solver @io --max-sweep-points 81 @rom --tfe-exact-band 93.96e9 94.03e9 --linear-solver gmres --precon pardiso-single --iterative-spmv mkl --max-iterations 120 --tolerance 1e-7 --out result/io-rom-mixed
```

在工程副本中将扫频设为 93.960–94.030 GHz、等距 15 点。命令使用原 `current.ngmesh` 和该副本，分别跑完整 Direct 和局部 MGAWE，读取 S11 最小的采样点、谷深和最大全阶残差。外存与混合精度滤波器命令见[混合精度回归](../changes/2026-11-18-iterative-mixed-precision.md)。数值单模端口的既有快扫频入口见[AWE 系列复现](../optimization/awe-family/paper-reproduction.md)。

41 点数值端口兼容性使用 `--aedt IOStructure.aedt --mesh current.ngmesh --basis-order 0 --port-method numerical --max-sweep-points 41 --no-write-all-fields`。先用 `--sweep direct --linear-solver direct` 建立基准，再按表中各行设置，分别指定独立 `--out`。这些阶数对应本次回归网格。

| 方法 | 选项 |
| --- | --- |
| AWE | `--sweep awe --awe-field-pade --awe-order 17 --linear-solver direct` |
| MGAWE | `--sweep mgawe --mgawe-order 6 --mgawe-points 3 --linear-solver direct` |
| WCAWE | `--sweep wcawe --wcawe-order 16 --linear-solver direct` |
| WCAWE 混合精度 | `--sweep wcawe --wcawe-order 16 --linear-solver gmres --precon pardiso-single --iterative-spmv mkl --tolerance 1e-7 --max-iterations 120` |

## 固化

具备门槛记录后，提交代码、规范与变更记录，并在固定提交上建立附注标签。只推送验收过的分支和对应标签。构建中间产物、日志、临时网格和 VTU 留在本地或归档系统。Git 工作流见[仓库管理约定](../development-workflow/git-workflow.md)。
