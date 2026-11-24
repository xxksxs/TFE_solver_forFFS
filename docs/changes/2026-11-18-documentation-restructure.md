# 文档体系重构全批次

- 目标与范围：为 `docs/` 建立统一入口、修改规范与变更留痕。同步源码已经证实的端口、扫频、求解器、构建和测试说明。
- 权威页：[文档入口](../README.md)、[工程参考](../engineering-reference.md)、[求解器选择](../architecture/solver-selection.md)、[修改规范](../development-workflow/change-policy.md)。
- 接口与数值约定：本次是文档变更。盘点后工作树出现的 MUMPS、Schwarz 和新增稀疏求解器测试按当前代码注明，性能和发布状态待另行确认。
- 验证：相对路径检查、补丁空白检查和工作树差异检查。测试和性能数据引用已有记录，本批未重新构建或运行求解器。
- 迁移和合并：根 README 正文归入文档入口和历史页，版本日志及根目录推导资料集中归档。旧实验口径和已有 PDF 删除保留历史状态，详见[迁移记录](../development-workflow/document-migration-plan.md)。
- 兼容性：根 README 保留入口，专题兼容索引保留。已迁移路径见迁移记录。
- 版本管理：为正式的 [测试文档目录](../testing/README.md) 添加 Git 忽略例外。CTest 生成目录仍忽略。
