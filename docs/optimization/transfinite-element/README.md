# 超限元（Transfinite Element, TFE）专题

本目录包含 TFE 理论、历史实现方案和工程入口。分块路径使用模态系数作为端口未知量。
## 原生分块装配

`--port-method tfe-block --sweep direct` 使用教材 §6.3.2 的端口迹替换。
有效端口迹由模态坐标表示，全阶未知量由内部棱元系数和总模态幅值组成。
`--sweep mgawe` 在同一系统上按截止划分子带，对矩阵作二阶 Chebyshev
插值，在多点生成局部矩并作复对称 Galerkin 降阶。

- 保留模数默认使用 AEDT 的 `NumModes`，可用 `--port-modes N` 覆盖。`--tfe-modes-per-port` 为兼容别名。
- 默认使用工程中的有效入射通道并按独立列求解。`--tfe-source PORT:MODE`
  指定单个单位入射，模态编号从 1 开始。`--tfe-all-sources` 求所有传播列。
- `s_parameters_modal.csv` 按频率、输入端口/模态和输出端口/模态列出结果。
  兼容 CSV 与 VTU 使用首个入射列。兼容 S21 取首个其他物理端口的主模。
- 支持均匀无损、单一外缘的 PEC 波导截面。精确截止返回诊断。
  AEDT 的显式面 Radiation 边界用一阶 ABC，阻抗取无损背景介质。
  截止点邻域用精确残差检查插值质量。
- 用同一批频点运行 Direct、MGAWE 和带 `--tfe-proxy-direct` 的 MGAWE 命令，
  可单独测量二阶代理插值误差和降阶误差。代理 Direct 写入
  `result_TFE_PROXY_DIRECT`。可用 `--tfe-subband-width HZ` 限制每个二阶
  插值区间的宽度。设置后按实际扫频采样点划分子带。
- `--tfe-exact-band LOW_HZ HIGH_HZ` 可重复指定闭区间，区间内已有的扫频
  频点使用精确 TFE Direct，区间外仍使用同一二阶代理。该选项不会加入额外频点。
  精确点记录各自的线性求解状态，并计入数值分解次数。

后端、IOStructure 全带与深谷工程结果见[联合验收](../../changes/2026-11-18-solver-release-validation.md)，独立分支的历史验证见[TFE 集成记录](../../changes/2026-11-18-tfe-integration.md)。

以下段落记录全棱边数值多模的实现。`tfe` 和 `transfinite` 是这条路径的历史别名。


> 注：本目录早期版本（理论文档与 README）描述的实际是"解析端口模"（APM）。该内容已迁出到 `../analytic-port/` 子目录，并对应模块 `bpfem::apm`、CLI `--port-method analytic`。

## 文档列表

- **`tfe-theory-cn.tex`** / **`tfe-theory-cn.pdf`**：理论文档（中文）。包含连续问题、解析模张量积基、增广系统矩阵、Schur 消元、与现有 NPM/APM 的关系、与 ALPS 协同。第 6.1 节"Schur 消元"开头有一个 remark 说明早期推导的因子约定问题与现行实现的对应关系。
- **`tfe-plan.md`**：开发计划与落地记录。包含已合入改动、CLI、性能、不变量、待办、限制。
- **`README.md`**：当前文件。

## 全棱边数值多模

CLI 入口如下。

```powershell
.\build\Release\bp_fem_solver.exe --basis-order 1 --port-method numerical-multimode --port-modes 5
```

每个端口面保留至多 `--port-modes` 个数值本征模，按截止波数排序。`PortModeSolver` 求本征模，`WavePortBC` 把各模的秩一端口项加到全棱边矩阵。激励施加在主模上。

## 与 APM 的关系

| 维度 | APM (`--port-method analytic`) | 全棱边数值多模 (`--port-method numerical-multimode`) |
|---|---|---|
| 模式数 | 仅 TE10 | 数值多模 |
| 端口算子 | 单秩-1 | $\sum_{(m,n)} \mathrm{rank-1}$ |
| 与 ALPS 兼容 | ✓ | ✓（自动展开为虚拟端口） |
| 物理来源 | 解析模 + FE L2 投影 | 端口网格本征解 + 多模波端口项 |
| 实现复杂度 | 端口面 mass + Cholesky | 数值本征解 + 多模 RHS/LHS 分支 |

详见 `tfe-plan.md` 与 `tfe-theory-cn.tex`。
