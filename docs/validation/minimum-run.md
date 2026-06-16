# 最小验证运行

本文说明确认数值流程可运行的最小命令。

## 单频点运行

在扫频中心附近运行一个频点：

```powershell
.\build\Release\bp_fem_solver.exe --aedt wg_bp_filter.aedt --mesh current.ngmesh --out result --max-sweep-points 1
```

## 检查项

- 求解残差。
- 默认日志应显示 `Nedelec basis order: 0` 和 `Local basis functions per tetrahedron: 6`。
- `result/result_DIRECT/s_parameters.csv`。
- `result/result_DIRECT/field_last.vtu`。

## 一阶层次棱元检查

修改高阶基函数、拓扑或装配后，应额外运行：

```powershell
.\build\Release\bp_fem_solver.exe --aedt wg_bp_filter.aedt --mesh current.ngmesh --out result --max-sweep-points 1 --basis-order 1
```

期望日志显示 `Nedelec basis order: 1` 和 `Local basis functions per tetrahedron: 20`。

## 何时使用

以下修改后至少运行一次单频点验证：

- FEM 装配。
- 端口模式。
- 求解器。
- S 参数提取。
- 场重构或 VTU 输出。
