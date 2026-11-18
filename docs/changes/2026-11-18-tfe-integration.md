# 原生 TFE 与 MGAWE 集成

## 目标与实现

在桌面当前工作树的隔离快照上集成端口迹替换、TE/TM 模、分块 Direct、二阶 Chebyshev MGAWE 及模态 S 参数。原生方程使用教材 Zhu–Cangellaris §6.3.2 式 (6.139)、IVEC 式 (3)–(6) 和专利实施方式步骤 1–8。经 PEC 消元后使用 $[x_E;\alpha]$ 及 $Q=\operatorname{diag}(I,B)$，体项是 $Q^TDQ$，出射波为 $\alpha-a$。默认的数值端口及旧 `tfe` 入口保持原有模型。

CMake 延用现有 PARDISO、MUMPS 和 Fortran 运行库部署目标，向既有测试增加 `tfe_tests`。应用经工厂选择原生策略。Direct、代理 Direct 和 MGAWE 离线矩求解共用 `FactorizedSolveSession` 的真实残差与失败状态。失败时保存相应频点及残差，终止计算。MGAWE 在线原生残差衡量近似质量，单独报告；代理残差写入频点日志。入口见 [TFE 专题](../optimization/transfinite-element/README.md)。

## 构建与 CLI

本机 MSVC 19.51 Release，Ninja，各构建目录独立。CTest 注册 `fastsweep_tests`、`io_parser_tests`、`sparse_solver_tests` 和 `tfe_tests`。

| 构建 | 验证 |
|---|---|
| `BPFEM_USE_MKL=OFF, BPFEM_USE_MUMPS=OFF` | 四项 CTest 通过；GMRES + ILU(0) 原生 TFE Direct/MGAWE 均跑通 5 点均匀波导。 |
| `BPFEM_USE_MKL=ON, BPFEM_USE_MUMPS=OFF` | 四项 CTest 通过；PARDISO 原生 TFE Direct/MGAWE 均跑通同一波导及 IOStructure。 |
| `BPFEM_USE_MKL=ON, BPFEM_USE_MUMPS=ON` | 使用本机串行 ZMUMPS 与 MKL sequential 时四项 CTest 通过。对本机另一套 threaded ZMUMPS/BLAS 组合，既有 `sparse_solver_tests` 未完成。 |
| `BPFEM_USE_MKL=OFF, BPFEM_USE_MUMPS=ON` | 四项 CTest 通过；MUMPS 及 GMRES + Schwarz（本地 MUMPS）均跑通原生 Direct/MGAWE。该构建的 ZMUMPS 自身链接 MKL sequential BLAS，程序的 PARDISO/MKL SpMV 功能关闭。 |

构建命令为 `cmake -S . -B <build> ...`、`cmake --build <build>`、`ctest --test-dir <build> --output-on-failure`。MUMPS 构建还须按实际 ZMUMPS 链接 ABI 提供 `BPFEM_MUMPS_INCLUDE_DIR`、`BPFEM_MUMPS_LIBRARIES` 与 `BPFEM_MUMPS_RUNTIME_DLLS`。均匀波导采用 20 mm × 12 mm 截面、20 mm 长、9–10 GHz 五点、每端口两模。CLI 选 `--port-method tfe-block --sweep direct|mgawe` 和相应后端。强制 GMRES 只迭代一次且容差 `1e-12` 时，Direct 返回 1，`run.json` 为 `fatal` 并保留 `iteration_limit` 及绝对/相对残差，该次未产生 S 参数 CSV。

## IOStructure 全带和局部深谷

共同设置为仓库 `IOStructure.aedt` 与配套 `current.ngmesh`，0 阶基、每端口 3 模、PARDISO、4 线程、仅末频场输出。全带 90–100 GHz 等距 41 点。误差指标为同网格精确 TFE Direct 下的复数相对 L2。耗时含进程内导入、装配、求解和输出，峰值为工作集。每组为单次实测。

| 全带配置 | S11 L2 | S21 L2 | 94 GHz S11 | 总用时 / 峰值 |
|---|---:|---:|---:|---:|
| 精确 TFE Direct（独立验收基准） | — | — | −46.0084 dB | 130.60 s / 2217.79 MiB |
| 二阶 MGAWE，4 GHz 子带，每带 1 展开点（独立验收） | 0.6186% | 0.1231% | −44.7852 dB | 49.01 s / 2287.40 MiB |
| 二阶 MGAWE，2 GHz 子带，每带 1 展开点 | 0.07469% | 0.01304% | −46.0748 dB | 69.36 s / 2282.66 MiB |
| 二阶 MGAWE，2 GHz 子带加深谷区间 Direct | 0.07457% | 0.01304% | −46.0084 dB | 59.15 s / 2285.64 MiB |

4 GHz 配置的代理 Direct 在 94 GHz 为 −44.7655 dB，说明该处约 1.22 dB 的谷深偏差主要来自二阶插值。收窄到 2 GHz 后，94 GHz 偏差降至 0.0664 dB。命令的扫频参数为 `--tfe-modes-per-port 3 --max-sweep-points 41 --no-write-all-fields --sweep mgawe --tfe-subband-width 2e9 --mgawe-points 1`。

另外复制 AEDT 并将频带改为 93.960–94.030 GHz，间隔 5 MHz，共 15 点。精确 Direct 在 93.985 GHz 采得 −53.5278 dB，窄区间 MGAWE 在同一采样点为 −53.5278 dB。谷深差 $2.36\times10^{-6}$ dB，采样谷频差 0 MHz，谷频的分辨率为 5 MHz。Direct 为 55.98 s / 2211.41 MiB，局部单展开点 MGAWE 为 44.26 s / 2278.27 MiB。这组局部模型独立拟合，精确 Direct 密集频点用于定位全带 250 MHz 采样未捕获的谷底。

在上述 2 GHz 命令后追加 `--tfe-exact-band 93.96e9 94.03e9`，41 点网格上只有 94 GHz 走精确 TFE，五个子带展开点和一个精确点共六次分解。94 GHz 的原生相对残差为 $1.06\times10^{-14}$，与独立 Direct 的谷深差低于 $3\times10^{-11}$ dB。该选项适用于原生 MGAWE，频点仍由 AEDT 定义，频率以 Hz 输入。局部加密时用相同 mesh 和局部 AEDT 频率区间。各次计时有机器负载波动，表中两种 2 GHz 的耗时差不能归因于额外一次精确解。

## 范围与待验证

上述 1% 结果限 IOStructure 41 点、3 模和所列配置。整个 10 GHz 频带的独立密集网格、连续频率谷值、HFSS 参考面/重归一化和 IVEC pillbox 尚未验证。截面支持均匀、无损且单一外缘的 PEC 波导。精确截止处的教材归一化有奇异性，CLI 给出诊断。对象型 Radiation 及跨网格简并模配准仍需验证。
