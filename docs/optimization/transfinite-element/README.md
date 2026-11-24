# 超限元（Transfinite Element, TFE）专题

本目录包含 TFE 理论、历史实现方案和工程入口。现行 CLI `--port-method tfe` 使用端口网格的数值 H(curl) 本征模，并把多个模式的波端口项加到全局系统。完整模态自由度替换与 Schur 推导见理论稿。

> 注：本目录早期版本（理论文档与 README）描述的实际是"解析端口模"（APM）。该内容已迁出到 `../analytic-port/` 子目录，并对应模块 `bpfem::apm`、CLI `--port-method analytic`。

## 文档列表

- **`tfe-theory-cn.tex`** / **`tfe-theory-cn.pdf`**：理论文档（中文）。包含连续问题、解析模张量积基、增广系统矩阵、Schur 消元、与现有 NPM/APM 的关系、与 ALPS 协同。第 6.1 节"Schur 消元"开头有一个 remark 说明早期推导的因子约定问题与现行实现的对应关系。
- **`tfe-plan.md`**：开发计划与落地记录。包含已合入改动、CLI、性能、不变量、待办、限制。
- **`README.md`**：当前文件。

## MVP 状态

第一迭代已合入主干。CLI 入口：

```powershell
.\build\Release\bp_fem_solver.exe --basis-order 1 --port-method tfe --tfe-modes-per-port 5
```

每个端口面同时纳入 `--tfe-modes-per-port` 个数值本征模（按截止波数排序，最低保留模式 index 0），每个模贡献一个秩-1 端口算子并入全局 H(curl) 矩阵。激励仅施加在主模上；其它模式作为吸收。

## 与 APM 的关系

| 维度 | APM (`--port-method analytic`) | TFE (`--port-method tfe`) |
|---|---|---|
| 模式数 | 仅 TE10 | 数值多模 |
| 端口算子 | 单秩-1 | $\sum_{(m,n)} \mathrm{rank-1}$ |
| 与 ALPS 兼容 | ✓ | ✓（自动展开为虚拟端口） |
| 物理来源 | 解析模 + FE L2 投影 | 端口网格本征解 + 多模波端口项 |
| 实现复杂度 | 端口面 mass + Cholesky | 数值本征解 + 多模 RHS/LHS 分支 |

详见 `tfe-plan.md` 与 `tfe-theory-cn.tex`。
