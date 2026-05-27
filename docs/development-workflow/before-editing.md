# 编辑前准备

开始修改代码或文档前，应先确定任务范围和最小必要上下文。

## 推荐步骤

1. 阅读 `../overview.md`，确认工程主数据流。
2. 在 `../task-splitting.md` 中确认任务类别。
3. 只读取完成任务需要的最小文件集合。
4. 阅读 `../primitives/README.md` 和 `../primitives/invariants.md` 中的关键不变量。
5. 涉及数值结果时，同时阅读 `../validation.md`。

## 上下文选择原则

- **解析器任务**：优先读取 `include/bpfem/core/Types.hpp` 和 `src/io/`。
- **FEM 任务**：优先读取 `src/fem/` 和 `../primitives/fem.md`。
- **求解器任务**：优先读取 `src/linalg/` 和 `../primitives/linear-algebra.md`。
- **后处理任务**：优先读取 `src/post/` 和 `../primitives/output.md`。
- **文档任务**：优先读取当前文档索引和相关专题文档。

## 风险判断

如果任务涉及端口、边方向、单位体系或求解器路径，应先写出影响范围，再开始修改。
