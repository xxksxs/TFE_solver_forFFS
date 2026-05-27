# 解析端口模 (Analytic Port Mode, APM) 专题

本目录维护**解析端口模**方法的资料。模块名 `bpfem::apm`、CLI `--port-method analytic`。

> 这部分内容曾经放在 `../transfinite-element/` 目录、命名为 `bpfem::tfe`，但实际上不是真正的超限元方法。重命名后这里专放 APM，真正 TFE 在 `../transfinite-element/`。

## 文档列表

- **`theory-cn.md`**：原理、与 NPM 对比、实测数据、限制、待办。

## MVP 状态

CLI 入口（默认）：

```powershell
.\build\Release\bp_fem_solver.exe --basis-order 1 --port-method analytic
```

仅支持 TE10 主模，做"解析模 + FE L2 投影 + 单秩-1 端口算子"。等价于"换了端口模来源的 NPM"。

## 与 TFE 的关系

APM = TFE 在 `--tfe-modes-per-port 1` 时的退化情形。两条路径的数值结果在所有频点上等价到 1e-13 量级（同样的端口面 dof 集合、同样的解析 TE10 公式、同样的 mass-projection 流程）。

如需多模端口建模，使用 `--port-method tfe`，详见 `../transfinite-element/`。
