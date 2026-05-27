# 快速扫频模块（Krylov MOR）专题目录

本目录维护"频域 FEM 快速扫频"模块的理论与开发资料。新增此目录是因为 Krylov MOR 是一个独立、跨求解器/装配/后处理的功能，单独成文便于评审与版本演进。

## 文档列表

- **`theory.tex`** / **`theory.pdf`**：基于 LaTeX 的英文理论模型文档。`pdflatex theory.tex` 直接编译。
- **`theory-cn.tex`** / **`theory-cn.pdf`**：中文翻译版（与英文版一一对应）。使用 `ctexart` + `xelatex theory-cn.tex` 编译。任选一份作为权威文档；中英两份内容必须保持同步，新增/修订一份后必须同步另一份。
- **`plan.md`**：开发计划。包含目标、非目标、模块布局、CLI 与配置、API 草案、风险、验收标准、分阶段交付、文档同步清单。
- **`README.md`**：当前文件。
- **`reviews/`**（按需创建）：评审反馈与决议记录。

## 评审与放行

- 第一步：阅读 `theory.tex` 与 `plan.md`，在 `reviews/round-1.md` 中记录修改意见。
- 第二步：作者按反馈修订两文。直到双方在评审记录中签字后，才允许进入实现阶段。
- 第三步：实现阶段开始后，本目录补充 `progress.md` 跟踪每个里程碑的完成情况，并把代码改动指向 `plan.md` 中具体任务。

## 与现有文档的关系

- 顶层入口：`../../optimization-roadmap.md`（在 *Phase 2* 与 *solvers-and-sweeps.md* 中已经登记 Krylov MOR）。
- 不变量约束：`../../primitives/invariants.md`（端口归一化、对称装配、单位、稀疏图复用）。
- 求解器选择：`../../architecture/solver-selection.md`（MOR 不替代直接/迭代求解器，是一种"扫频策略"）。
- 验证基准：基准库一旦建立后，新增 *Cavity h-convergence + Sweep MOR* 案例到 `../../validation/`。
