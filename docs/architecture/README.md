# 架构文档索引

本目录将工程架构说明拆分为多个专题文档。Phase 1 接口化（见
`../optimization/strategy-interfaces/plan.md`）后，运行时分发全部走工厂 +
抽象接口；架构文档已同步反映这一点。

## 阅读顺序

1. **模块布局**：`module-layout.md`，理解静态库 / 可执行 / 源码目录。
2. **依赖方向**：`dependency-direction.md`，含"加新策略最小改动集"。
3. **求解器选择**：`solver-selection.md`，含 CLI 路由表与 CMake 构建标志。
4. **风险点**：`risk-points.md`，跨模块耦合不变量。

## 文档列表

- **`module-layout.md`**：顶层目标、模块职责、`bc` / `sweep` / `factory` 的位置。
- **`dependency-direction.md`**：依赖方向、设计原则、加新条目的最小改动清单。
- **`solver-selection.md`**：CLI 路由 (`--linear-solver`、`--precon`、`--sweep`)、
  CMake `BPFEM_USE_MKL` 行为、`SolveResult` 契约。
- **`risk-points.md`**：端口模式 / Dirichlet / 单位 / 场重构 / 求解器路径不变量。

## 与其它文档的关系

- 设计意图与 PR 拆分：`../optimization/strategy-interfaces/plan.md`
- 接口契约（参数、调用顺序、错误模型）：`../optimization/strategy-interfaces/api.md`
- 加新策略示例：`../optimization/strategy-interfaces/example.md`
- 任务模板：`../task-splitting/strategy-interfaces-step{1..5}-*.md`
