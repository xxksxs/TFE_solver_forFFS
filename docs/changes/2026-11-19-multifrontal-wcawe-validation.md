# 多波前＋WCAWE 首轮工程验证

在 `2026.11.18-tfe-solver` 的 `c8b4ed7` 版本上，通过 `--linear-solver mumps` 选择复数对称多波前直接求解。WCAWE 的 `FactorizedSolveSession` 在固定展开点分解一次，后续矩逐个调用 `solveFactorized`。MUMPS 和 PARDISO 共用 WCAWE 的矩递推、正交化及在线投影。四项 CTest 通过。

## 复现

建立启用 oneMKL 和 MUMPS 的 Release 构建。在仓库根运行，从同一可执行文件依次测试。在本机使用 MSVC 19.51、oneMKL 2025.3、串行 ZMUMPS 5.9.2。PARDISO 使用四线程。启用相同环境变量的各例仍受实际库编译方式约束。

```powershell
$exe = '.\build_combined\Release\bp_fem_solver.exe'
$env:MKL_NUM_THREADS = '4'
$env:OMP_NUM_THREADS = '4'
$common = @('--aedt','IOStructure.aedt','--mesh','current.ngmesh',
    '--basis-order','0','--port-method','numerical','--max-sweep-points','41',
    '--sweep','wcawe','--wcawe-order','16','--wcawe-expansion','95e9',
    '--tolerance','1e-7','--no-write-all-fields')
& $exe @common --linear-solver direct --out result/wcawe-pardiso
& $exe @common --linear-solver mumps --mumps-mode incore --out result/wcawe-mumps
& $exe @common --linear-solver mumps --mumps-mode ooc --out result/wcawe-mumps-ooc
& $exe @common --linear-solver direct --pardiso-mode ooc --out result/wcawe-pardiso-ooc
& $exe @common --linear-solver mumps --mumps-mode blr --mumps-blr-tolerance 1e-14 --out result/wcawe-mumps-blr
```

IOStructure 采用 90–100 GHz 的 41 点网格。检查 `result_WCAWE/diagnostics.json` 的有效维数、终止原因和递推残差，再对齐 `s_parameters.csv` 的频率，比较复数 S11/S21。用 `timing.json` 核对分解/回代计数、进程总耗时和峰值工作集。正式计时关闭 `--wcawe-check-moments`，每组单次实测。

带通滤波器采用 40–43 GHz、21 点、41.5 GHz 展开点，同样使用零阶数值单模端口。MUMPS 外存与 PARDISO 均完成 16 阶构基，两路 S 参数最大复数绝对差低于 `1.1e-12`。额外执行 `--wcawe-check-moments` 后，仍为一次分解，回代 32 次，独立参考矩的最大投影误差为 `4.84e-14`。

三模端口配置的该滤波器频带跨越高次模截止。单展开点 WCAWE 在该频带用两种后端均触发 `No polynomial expansion on the requested port cutoff branch`。数值单模对照避免了这项条件，三模全带需要分段。

## IOStructure 实测

| 后端和模式 | 有效维数 | 数值分解 / 回代 | 端到端秒数 | 峰值工作集（MiB） |
| --- | ---: | ---: | ---: | ---: |
| 双精度 PARDISO 内存内 | 16 | 1 / 16 | 41.66 | 2639.6 |
| MUMPS 内存内 | 16 | 1 / 16 | 48.65 | 2460.1 |
| MUMPS 外存 | 16 | 1 / 16 | 70.88 | 1543.5 |
| 双精度 PARDISO 外存 | 16 | 1 / 16 | 69.42 | 2707.7 |
| MUMPS BLR，压缩容差 `1e-14`，保护版 | 16 | 1 / 16 | 52.27 | 2471.0 |

MUMPS 内存内及外存与 PARDISO 内存内的复数 S11/S21 最大绝对差均低于 `9e-13`。相对全阶 Direct 的两路复数误差分别约为 0.000194% 和 0.05644%。MUMPS 外存峰值工作集降低 41.5%，16 次回代共 27.62 秒，内存内模式为 3.96 秒。

BLR 压缩容差 `1e-8` 时，首个 RHS 真实相对残差为 `2.60e-7`，超过 `1e-7`，退出 1。压缩容差 `1e-12` 时，第二阶未通过 `1e-10` 的递推残差门槛，原发布版输出一维结果并退出 0。该条为失败复现，不参与工况性能比较。

## 递推失败保护

在隔离分支 `experiment/multifrontal-wcawe` 中，`WcaweSweep` 将 `ToleranceExceeded` 按离线数值失败上报，在写 S 参数前终止。保留递推残差 `1e-10` 的门槛。用 `1e-12` 压缩容差复测时，日志报告 95 GHz 展开点第二个 RHS 递推残差超限，进程退出 1、`run.json` 状态为 `fatal`，结果目录没有 S 参数 CSV。

保护版的带通滤波器 MUMPS 外存复测构足 16 维基，独立参考矩最大投影误差仍为 `4.84e-14`，一次分解、32 次回代。4 项 CTest 均通过。IOStructure 的 BLR 试验在压缩容差 `1e-14` 下完成 16 阶，最大递推相对残差 `3.09e-11`，相对 PARDISO 的复数 S11/S21 最大绝对差低于 `5.1e-13`。两档压缩容差通过与失败的情况限于本算例和本机库构建。
