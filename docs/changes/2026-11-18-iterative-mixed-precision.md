# 混合精度 GMRES 与报告容差

## 目标与范围

恢复三模 `wg_bp_filter` 上的工程迭代求解。旧配置在 41.5 GHz 使用 GMRES(30)、4 个 Schwarz 子域、一层重叠、局部 MUMPS 和 MKL SpMV。独立复测 120 步后的真实相对残差为 0.470939，退出码为 1。验收沿用 `1e-7` 的真实相对残差门槛，并和双精度 PARDISO 对照复数 S 参数。

## 实现与数值约定

- `--linear-solver gmres --precon pardiso-single` 新增显式的混合精度路径，要求 oneMKL。预条件器用 `mtype=6` 和 `iparm[27]=1` 准备全局单精度复数对称因子，PARDISO 内部精化关闭。分析和分解沿用现有后端生命周期。同一个展开点的多个右端项可共享因子，换矩阵时重建。
- 外层用双精度 flexible GMRES。每一步保留实际用于 Arnoldi 的 `z_j = M⁻¹v_j`，用 `Zy` 更新解。原始矩阵和 RHS 的物理模型、精度及收敛判据保持一致。原生和 MKL SpMV 路径都支持此选项。原有 Schwarz 和默认求解器选择沿用现有配置。 [oneMKL PARDISO 数组精度](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-2/pardiso-data-type.html)
- `run.json` 中的线性求解、MUMPS BLR 和降阶丢弃容差使用可往返解析的十进制表示。字段仍为字符串，`1e-7` 和 `1e-10` 等小容差能完整保留。单测对线性容差与非整幂 BLR 容差做回读检查。小型复数系统测试覆盖混合精度迭代、同矩阵多 RHS 和换矩阵后重建。
- 主要代码位于 `src/linalg/MklPardisoSolver.cpp`、`PardisoBackend.cpp`、`GmresSolver.cpp`、`src/factory/SparseSolverFactory.cpp`、`src/app/Application.cpp`、`src/core/RunReport.cpp` 和对应头文件。CLI 说明见 [求解器选择](../architecture/solver-selection.md)。

## 复现命令

使用 Windows x64 Release，MSVC 19.51、oneMKL 2025.3 ILP64、MUMPS 5.9.2。MUMPS 为串行 MPI 接口构建。用启用 MKL 和 MUMPS 的可执行文件运行。以仓库根目录为当前目录。

```powershell
$exe = 'C:/path/to/bp_fem_solver.exe'
$case = Resolve-Path '../../example/bp_filter'
$env:MKL_NUM_THREADS = '4'
$env:OMP_NUM_THREADS = '4'
$common = @('--aedt', "$case/wg_bp_filter.aedt", '--mesh', "$case/current.ngmesh",
            '--basis-order', '0', '--port-method', 'tfe', '--tfe-modes-per-port', '3',
            '--sweep', 'direct', '--no-write-all-fields', '--tolerance', '1e-7')
$iter = @('--linear-solver', 'gmres', '--precon', 'pardiso-single', '--iterative-spmv', 'mkl',
          '--gmres-restart', '30', '--max-iterations', '120')
& $exe @common @iter --max-sweep-points 1 --out ./result/single-point
& $exe @common @iter --max-sweep-points 21 --out ./result/mixed-21
& $exe @common --max-sweep-points 21 --linear-solver direct --out ./result/pardiso-21
& $exe @common --max-sweep-points 21 --linear-solver mumps --mumps-mode blr --mumps-blr-tolerance 1e-2 --out ./result/blr-21
& $exe @common --max-sweep-points 21 --linear-solver mumps --mumps-mode ooc --out ./result/ooc-21
& $exe @common --max-sweep-points 1 --linear-solver gmres --precon schwarz --schwarz-local-solver mumps --schwarz-subdomains 4 --schwarz-overlap 1 --max-iterations 120 --iterative-spmv mkl --out ./result/original-schwarz
```

读取各运行目录下 `result_DIRECT/run.json` 的真实相对残差、迭代数、耗时和峰值工作集。按频率对齐 `s_parameters.csv`，对复数值计算最大绝对差和整带相对 L₂ 误差。

## 验证结果

41.5 GHz 的混合精度路径以 2 步达到 1.04e−9，退出码 0。原 Schwarz 配置复测得到 0.470939，退出码 1，且无 S 参数 CSV。保留失败门槛和诊断。

多频点为 40–43 GHz 等距 21 点，共 62,306 个未知量。所有表中路径均完成 21 点。PARDISO 与混合精度 GMRES 各运行 3 次，表中显示耗时与峰值工作集的中位数，残差和 S 参数差取三次最大值。MUMPS 各测一次，耗时包含导入、装配、线性求解及输出。

| 路径 | 外层迭代步数 | 最大真实相对残差 | 峰值 MiB | 总耗时 s | 复数 S11 / S21 最大绝对差 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 双精度 PARDISO 直接 | — | 1.16e−13 | 302.6 | 7.97 | 基线 |
| 双精度 flexible GMRES + 单精度 PARDISO | 2–3 | 9.05e−8 | 270.3 | 7.07 | 3.17e−8 / 2.85e−8 |
| MUMPS BLR 直接，压缩容差 1e−2 | — | 4.22e−12 | 306.6 | 77.43 | 1.21e−12 / 1.70e−12 |
| MUMPS 外存直接 | — | 4.91e−12 | 204.4 | 64.80 | 1.24e−12 / 1.71e−12 |

混合精度路径的复数整带相对 L₂ 误差在三次中最大为 S11 1.05e−8、S21 9.68e−9。MUMPS BLR 为 S11 4.42e−13、S21 5.41e−13，MUMPS 外存为 S11 4.50e−13、S21 5.37e−13。所有实际比较通过 `|S − S_ref| ≤ 1e-5 + 1e-3 × |S_ref|`。

MKL + MUMPS 构建及无 MKL/MUMPS 构建各通过 CTest 3/3。命令行输出已验证 `tolerance` 为 `1e-07`，默认 `mumpsBlrTolerance` 为 `1e-10`，指定 `1e-2` 时记为 `0.01`。MUMPS 单次耗时受环境波动影响，另一次相同配置的耗时为 BLR 34.9 秒、外存 42.1 秒。混合精度与 PARDISO 的三次原始耗时分别为 6.90 / 11.22 / 7.07 秒和 7.21 / 21.09 / 7.97 秒；峰值分别为 270.0 / 272.2 / 270.3 MiB 和 287.0 / 302.6 / 304.8 MiB。

## 限制与后续验收

混合精度路径依赖全局单精度分解，具备实际双精度 Krylov 迭代和误差修正。MUMPS BLR 直接与外存直接分别使用压缩因子和磁盘因子。初期实验也把全局 MUMPS BLR 用作 GMRES 预条件器，在该工程单点只需 1 步，且未观察到明确的峰值优势，因此没有在交付补丁中新增该 CLI 路径。

原 Schwarz 仍未达到验收门槛，此次由新增混合精度配置提供经验证的迭代解。21 点峰值中位数相对 PARDISO 降低 10.7%。IOStructure 及高阶降阶展开点的混合精度表现需要独立验证。
