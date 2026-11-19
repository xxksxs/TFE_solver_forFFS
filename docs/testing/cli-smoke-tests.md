# CLI Smoke Test

本文说明命令行程序的最小运行检查。

> 算例配对：主仓库现行 `current.ngmesh` 对应 `IOStructure.aedt`。下文滤波器命令中的网格路径需换为与滤波器工程配套的网格。旧路径、默认选项和历史输出待复核。运行对照时分别指定 `--out`。

## 帮助信息

```powershell
.\build\Release\bp_fem_solver.exe --help
```

期望结果：

- 程序成功退出。
- 用法文本包含 `--aedt`、`--mesh`、`--out`、`--max-sweep-points`、`--max-iterations`、`--tolerance`、`--basis-order`、`--field-output-order`、`--write-all-fields`、`--sweep`、`--alps-order`、弃用别名 `--alps-krylov-order`、`--port-method`（`numerical|analytic|tfe`）与 `--tfe-modes-per-port`。

## 单频点求解

原生分块 TFE 的运行命令和两种扫频的验证见 [TFE 专题](../optimization/transfinite-element/README.md) 与 [集成变更记录](../changes/2026-11-18-tfe-integration.md)。

```powershell
.\build\Release\bp_fem_solver.exe --aedt wg_bp_filter.aedt --mesh current.ngmesh --out result --max-sweep-points 1
```

期望结果：

- 工程和网格被正确导入。
- 日志显示完成一个频点求解。
- 日志显示 `Nedelec basis order: 0`。
- 日志显示 `Local basis functions per tetrahedron: 6`。
- 写出 `result/result_DIRECT/s_parameters.csv`。
- 开启场输出时写出 `result/result_DIRECT/field_last.vtu`。

## 一阶层次棱元烟测

```powershell
.\build\Release\bp_fem_solver.exe --aedt wg_bp_filter.aedt --mesh current.ngmesh --out result --max-sweep-points 1 --basis-order 1
```

期望结果：

- 程序成功退出。
- 日志显示 `Nedelec basis order: 1`。
- 日志显示 `Field output order: 2`（一阶基函数下默认 p_out=2）。
- 日志显示 `Local basis functions per tetrahedron: 20`。
- 日志显示端口模式 DOF 数量大于零阶路径，以确认一阶端口边 DOF 和端口面 DOF 参与投影。
- 写出 `result/result_DIRECT/s_parameters.csv`。
- 写出 `result/result_DIRECT/field_last.vtu`，VTU 体积比零阶路径明显增大（每四面体 10 节点）。
- 在通带频点上检查 `|S11|² + |S21|² ≤ 1 + ε`（无源情形下能量守恒，ε ≈ 1e-3 数值误差量级）。

## 高阶可视化烟测

```powershell
.\build\Release\bp_fem_solver.exe --aedt wg_bp_filter.aedt --mesh current.ngmesh --out result --max-sweep-points 1 --basis-order 1 --field-output-order 3
```

期望结果：

- 程序成功退出。
- 日志显示 `Field output order: 3`。
- 写出的 VTU 体积约为 p_out=2 的 2×（每四面体 20 节点：4 顶点 + 12 边内点 + 4 面心）。
- 写出的 `s_parameters.csv` 与 p_out=2 的相同（`--field-output-order` 不影响求解结果）。
- 在 ParaView (≥ 5.5) 打开后，cell 类型应被识别为 *VTK_LAGRANGE_TETRAHEDRON*，按 p=3 多项式渲染。
- VTU 中 `PointData` 应仅包含 `E_real / E_imag / H_real / H_imag` 四个矢量数组，不再有模值或 PEC 切向标量。

## ALPS 快速扫频烟测

```powershell
.\build\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 5 --sweep alps --alps-order 12 --no-write-all-fields --out result
```

期望结果：

- 程序成功退出。
- 日志显示 `ALPS expansion frequency: <band-center> GHz`。
- 日志显示 `ALPS ROM dimension:` 数值（默认 q=30、Np=2 时通常为 60，可能因去冗（deflation）略小）。
- 离线阶段（`buildOffline`）耗时占主导，在线阶段 `evaluate` 极快（单点亚秒）。
- 写出 `result/result_ALPS/s_parameters.csv`、`result/result_ALPS/diagnostics.json` 和 `result/result_ALPS/field_last.vtu`；ALPS 只重构最后一个频点的场，不逐频点写 `field_*.vtu`。
- 同等参数下 `--sweep alps` 与 `--sweep direct` 在每个频点的 |S| 偏差应远小于 0.01 dB（实际通常 < 1e-10）。

如果材料含 σ > 0 的导率，会立即抛错并提示改用 `--sweep direct`。这是 MVP 阶段的 lossless 限定，详见 `../optimization/alps-sweep/plan.md` §7。

## 全棱边数值多模烟测

```powershell
.\build\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 5 --port-method numerical-multimode --port-modes 5 --no-write-all-fields --out result
```

期望结果：

- 程序成功退出。
- 日志显示 `Port-mode method: numerical-multimode (full edge coordinates, 5 modes/port)`。
- 日志记录端口的矩形尺寸及 `retained numerical modes`，报告的 `port_formulation` 为 `numerical-multimode`。
- `--port-modes 1` 时与 `--port-method numerical` 的 S 参数一致。两者共用二维 H(curl) 本征求解。
- 用同一设置的 `--port-method tfe --tfe-modes-per-port 5` 复跑，核对 S 参数与当前名称一致。

非矩形截面使用端口面三角网求本征模。矩形检测用于尺寸报告。

## 求解后端与预条件报告

- 用 `--sweep direct --linear-solver mumps` 检查 `linear_solver_family=direct` 和 `linear_solver_backend=MUMPS`。
- 用 `--sweep direct --linear-solver gmres --precon schwarz --schwarz-local-solver mumps` 检查 `linear_solver_family=iterative`、`active_preconditioner=schwarz`、`active_schwarz_local_solver=mumps`。
- `--linear-solver pardiso` 与兼容别名 `direct` 使用同一种线性后端。与 `--sweep direct` 组合时可分别核对这两层的选择。

## NPM（数值端口模，默认）烟测

```powershell
.\build\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 5 --no-write-all-fields --out result
```

期望结果：

- 程序成功退出。等价命令 `--port-method numerical`。
- 日志显示 `Port-mode method: numerical (2D H(curl) port eigensolve, dominant mode)`。
- 日志显示每个端口 `numerical mode: k_c^2=… (k_c=… rad/m)`，矩形端口下 $k_c\approx\pi/a$。
- 写出的 `s_parameters.csv` 与 `--port-method tfe --tfe-modes-per-port 1` 字节级等价。
- 与 `--port-method analytic` 在 TE10 主导频段偏差 ~ 1e-4 dB（端口网格对解析形状的几何拟合误差，物理上 NPM 更准）。
- 在通带频点上 `|S11|² + |S21|² ≤ 1 + ε`，ε ≈ 1e-3。

如果端口本征求解 fatal `LAPACKE_dsygv info > N`，参考 `../validation/failure-signatures.md` 中
"端口本征求解 fatal `LAPACKE_dsygv info > N`" 一节排查面基函数线性独立性与端口求积阶。

## APM（解析端口模，闭式 TE10 参考）烟测

```powershell
.\build\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 5 --port-method analytic --no-write-all-fields --out result
```

期望结果：

- 程序成功退出。
- 日志显示 `Port-mode method: analytic (TE10 closed-form, rectangular only)`。
- 日志显示每个端口的 `rectangular dimensions: a=… mm, b=… mm`。
- 写出的 `s_parameters.csv` 与默认 `--port-method numerical` 在通带 / 阻带平坦区域偏差通常 < 0.01 dB（端口网格几何拟合误差量级）。
- 在通带频点上 `|S11|² + |S21|² ≤ 1 + ε`，ε ~ 1e-2。
- 非矩形端口面会立即抛错：`face N is not a recognizable rectangular port`；遇此情况请改用默认 `--port-method numerical` 或 `--port-method tfe`。

## 失败时优先检查

- 输入文件路径是否存在。
- 构建配置是否与命令路径一致。
- MKL 是否被正确找到或是否已 fallback。
- 求解残差是否过高。
