# 任务拆分文档索引

本目录用于把工作拆分给不同开发者或 AI 模型，避免每个执行者都必须阅读整个工程。

## 通用任务说明模板

分配任务时建议提供：

- **目标**：本次需要改变什么。
- **优先阅读文件**：完成任务所需的最小上下文。
- **可能修改文件**：预期修改范围。
- **验证命令**：最小有效检查。
- **风险提示**：不能破坏的不变量。

## 文档列表

- **`parser-tasks.md`**：AEDT 与 NGMesh 解析任务。
- **`fem-assembly-tasks.md`**：体装配、材料处理、PEC 和波端口边界任务。
- **`port-mode-tasks.md`**：端口本征模、归一化和投影任务。
- **`solver-tasks.md`**：稀疏矩阵、直接求解器和 fallback 求解器任务。
- **`post-processing-tasks.md`**：S 参数提取和 VTU/CSV 输出任务。
- **`documentation-tasks.md`**：文档维护任务。
- **`model-routing.md`**：推荐模型分工。

## 架构演进 Phase 1：策略接口化（已落地）

接口化的 5 步全部合入主干，参考主计划 `../optimization/strategy-interfaces/plan.md`。
当前可作为后续 PR 的范式：每个 step 是一个独立 PR。

- **`strategy-interfaces-step1-sparse-solver.md`**：`ISparseSolver` 抽象 + Pardiso/BiCGSTAB Backend + 工厂。✓
- **`strategy-interfaces-step2-boundary-condition.md`**：`IBoundaryCondition` 抽象 + WavePort 迁移 + 三个空 BC 骨架。✓
- **`strategy-interfaces-step3-preconditioner.md`**：`IPreconditioner` 抽象 + Jacobi 实现。✓
- **`strategy-interfaces-step4-sweep-strategy.md`**：`ISweepStrategy` 抽象 + DirectSweep / AlpsSweep 迁移。✓
- **`strategy-interfaces-step5-docs.md`**：架构文档、API 契约、plugin 示例。✓

## 加新策略时使用的开发任务模板

按"接口实现 + 工厂枚举 + CLI 字符串映射"三步法，参考下面对应的步骤文档：

- 加新边界条件 → 参考 step2 + `../optimization/strategy-interfaces/example.md`
- 加新线性求解器后端 → 参考 step1 模板
- 加新预条件 → 参考 step3 模板
- 加新扫频策略 → 参考 step4 模板

接口契约（参数 / 调用顺序 / 错误模型）见
`../optimization/strategy-interfaces/api.md`。
