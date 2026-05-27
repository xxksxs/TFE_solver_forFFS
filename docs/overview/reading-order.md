# 推荐阅读顺序

本文给出新开发者或 AI 模型进入工程时的推荐阅读路径。

## 基础路径

1. **根目录说明**：先读 `../../README.md`，了解构建、运行和模块布局。
2. **项目树说明**：读 `../project-tree.md`，识别源码、文档、输入资产和生成目录。
3. **应用入口**：读 `../../src/main.cpp`，再读 `../../src/app/Application.cpp`。
4. **共享类型**：修改解析器、FEM、求解器或后处理前，先读 `../../include/bpfem/core/Types.hpp`。
5. **输入前端**：读 `../../src/io/AEDTParser.cpp` 和 `../../src/io/NGMeshParser.cpp`。
6. **FEM 核心**：读 `../../src/fem/EdgeTopology.cpp`、`../../src/fem/PortModeSolver.cpp` 和 `../../src/fem/FEMAssembler.cpp`。
7. **线性求解**：读 `../../src/linalg/SparseMatrix.cpp`、`../../src/linalg/MklPardisoSolver.cpp` 和 `../../src/linalg/BiCGStabSolver.cpp`。
8. **后处理**：读 `../../src/post/ResultExtractor.cpp` 和 `../../src/post/OutputWriter.cpp`。

## 按任务跳读

- **只改文档**：优先读 `document-map.md`、`../project-tree.md` 和相关专题目录。
- **只改测试说明**：优先读 `../testing/README.md`。
- **只改验证流程**：优先读 `../validation/README.md`。
- **只改原语说明**：优先读 `../primitives/README.md`。
