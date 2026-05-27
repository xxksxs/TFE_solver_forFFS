# ALPS 快速扫频专题目录

本目录维护"自适应 Lanczos-Padé 快速扫频 (ALPS, Adaptive Lanczos-Padé Sweep)"的理论资料。ALPS 是计算电磁学（CEM）领域历史最悠久的 FEM 频域快速扫频方法之一，由 Zhao 与 Lee 等人在 1999 ~ 2002 年提出，至今仍是商业级求解器（HFSS、CST、FEKO）"频率插值扫频/快速扫频"功能的主要数学骨架之一。

## 文档列表

- **`theory-cn.tex`** / **`theory-cn.pdf`**：基于 LaTeX 的中文理论文档。使用 `ctexart` + `xelatex theory-cn.tex` 编译。包含：原理概述、状态空间形式化、非对称 Lanczos 过程、Padé via Lanczos 定理、单点/多点自适应算法、对 Maxwell-FEM 的具体落地、误差指示器、与 Krylov MOR 的关系、收敛性与稳定性。
- **`plan.md`**：开发计划与 MVP 落地记录（CLI、模块布局、性能、不变量、待办）。
- **`README.md`**：当前文件。

## MVP 状态

第一迭代已合入主干（2026R1 时刻）。CLI 入口：

```powershell
.\build\Release\bp_fem_solver.exe --basis-order 1 --sweep alps --alps-krylov-order 30
```

在 BP filter 基准上 101 频点扫频从 ~22 min 降到 ~104 s（加速 ~12×），与 `--sweep direct` 在所有频点的 S 参数偏差 < 1e-10。详细说明见 `plan.md`。

## 与现有文档的关系

- 顶层入口：`../../optimization-roadmap.md`，`../solvers-and-sweeps.md` 中已经登记 *AWE / Krylov MOR* 等快速扫频候选；ALPS 现在作为补充候选与之并列。
- 平行专题：`../krylov-mor-sweep/` 为 Krylov MOR 的理论与开发计划。两者的数学关系见本目录 `theory-cn.tex` §10。
- 不变量：`../../primitives/invariants.md`，`../../architecture/risk-points.md`。

## 评审与放行

MVP（单展开点 + 复对称 Galerkin + lossless 限定）已合入主干。后续多展开点、自适应、有损材料、极点-留数解析、ROM 序列化等改动应通过 PR review 后扩展，并在 `plan.md` §5 待办表中标注完成状态。
