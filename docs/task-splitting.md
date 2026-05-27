# 任务拆分指南

本文件保留为任务拆分文档的兼容入口。详细内容已拆分到 `task-splitting/` 目录。

## 推荐入口

- **索引**：`task-splitting/README.md`
- **解析器任务**：`task-splitting/parser-tasks.md`
- **FEM 装配任务**：`task-splitting/fem-assembly-tasks.md`
- **端口模式任务**：`task-splitting/port-mode-tasks.md`
- **线性求解任务**：`task-splitting/solver-tasks.md`
- **后处理任务**：`task-splitting/post-processing-tasks.md`
- **文档任务**：`task-splitting/documentation-tasks.md`
- **模型分工**：`task-splitting/model-routing.md`

## 通用任务说明模板

分配任务时建议提供：

- **目标**：本次需要改变什么。
- **优先阅读文件**：完成任务所需的最小上下文。
- **可能修改文件**：预期修改范围。
- **验证命令**：最小有效检查。
- **风险提示**：不能破坏的不变量。
