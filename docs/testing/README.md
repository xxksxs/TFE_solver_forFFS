# 测试文档索引

当前 CMake 注册 `fastsweep_tests`、`io_parser_tests` 和 `sparse_solver_tests`。使用 `ctest --test-dir <构建目录> -C Release --output-on-failure` 执行 Release 测试，多配置生成器需要 `-C`。测试覆盖以当前源码为准。

io_parser_tests 覆盖 AEDT/NGMesh 输入和端口面映射。除 CTest 默认的小型夹具外，
还可对真实工程执行只解析、不装配、不求解的轻量检查：

    .\build\Debug\io_parser_tests.exe ..\IOStructure.aedt .\current.ngmesh

本目录说明构建检查、CLI smoke test、输出文件检查和回归检查。本工程目前以构建检查和 smoke test 为主，还没有完整单元测试框架。

## 阅读顺序

1. **构建检查**：阅读 `build-checks.md`。
2. **CLI smoke test**：阅读 `cli-smoke-tests.md`。
3. **输出文件检查**：阅读 `output-file-checks.md`。
4. **未来单元测试**：阅读 `future-unit-tests.md`。
5. **回归清单**：阅读 `regression-checklist.md`。

## 文档列表

- **`build-checks.md`**：默认构建、fallback 构建和 MKL 构建。
- **`cli-smoke-tests.md`**：帮助信息和单频点运行检查。
- **`output-file-checks.md`**：CSV 和 VTU 输出检查。
- **`future-unit-tests.md`**：建议补充的单元测试方向。
- **`regression-checklist.md`**：数值修改前后的最小回归清单。
