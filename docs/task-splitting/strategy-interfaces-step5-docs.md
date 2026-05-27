# Step 5：架构文档与 plugin 示例

> Phase 1 / Step 5 of `docs/optimization/strategy-interfaces/plan.md`。可在 Step 1 完成后开始（用占位接口；Step 2/4 完成后再补充示例代码）。

## 目标

把 Phase 1 的成果固化到文档：

- API 字段意义（IBoundaryCondition / ISparseSolver / ISweepStrategy）
- 加新策略的 step-by-step 指南
- 把 README 与 architecture/* 的"Application 直接调 PARDISO"过期描述更新

## 优先阅读文件

- `../optimization/strategy-interfaces/plan.md`
- 完成后的 Step 1–4 PR
- `../README.md`
- `../architecture.md`
- `../task-splitting/README.md`

## 可能修改 / 新增文件

```
docs/optimization/strategy-interfaces/api.md             新增（接口字段意义、调用顺序）
docs/optimization/strategy-interfaces/example.md         新增（如何加 ImpedanceBC：完整可跑代码）
docs/architecture/strategy-interfaces.md                 新增（架构图，替换 Application 直调 PARDISO 的旧图）
docs/architecture.md                                     "求解器装配链路" 段同步
docs/task-splitting/boundary-condition-tasks.md          新增（开发新 BC 的任务模板，沿用 task-splitting 风格）
docs/task-splitting/linear-solver-tasks.md               新增
docs/task-splitting/sweep-strategy-tasks.md              新增
docs/task-splitting/README.md                            新增三条索引
README.md                                                "模块布局" 段同步、加 --linear-solver / --precon CLI 段
```

## 验证命令

```powershell
# 文档构建：本仓库无 docs CI，靠人工 review。这一项验证清单：
# 1. 所有新增 .md 文件 ≤ 200 行（避免单一文件臃肿）
# 2. example.md 给出的 ImpedanceBC plugin 真能编译过：
cmake --build build_mkl --config Release --target bp_fem_core
# 3. 验证 README CLI 段与 src/app/Application.cpp parseOptions 实际选项一致（grep 对比）
```

## 不变量 / 风险提示

- **不要重复 plan.md**：plan.md 是设计意图，api.md 是契约（参数语义、前后置条件、错误模型）。两者目标不同，避免内容拷贝。
- **example.md 的 plugin 示例必须能编译**：写完后实际放进 `src/bc/ImpedanceBC.cpp` 做最小实现（哪怕只是 PEC + warning），运行 BP filter 不退化。
- **架构图**：用 ASCII 画即可，不必引入 mermaid。本仓库 docs 没有渲染管线。

## 估时

0.5 day。

## 输出物

- 1 个 PR
- 7 个新文档 / 4 个更新
- 一份 example.md 内的 plugin 真能编译并跑通 BP filter
