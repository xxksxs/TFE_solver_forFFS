# 项目树说明

本文档说明 `bp_filter` 示例工程的目录结构、关键文件职责，以及哪些目录属于源码、文档、输入资产或生成产物。

## 项目树

```text
bp_filter/
├── CMakeLists.txt
├── README.md
├── wg_bp_filter.aedt
├── current.ngmesh
├── include/
│   └── bpfem/
│       ├── app/
│       │   └── Application.hpp
│       ├── core/
│       │   ├── Constants.hpp
│       │   ├── Logger.hpp
│       │   ├── Math.hpp
│       │   ├── Types.hpp
│       │   └── Utilities.hpp
│       ├── fem/
│       │   ├── EdgeTopology.hpp
│       │   ├── FEMAssembler.hpp
│       │   └── PortModeSolver.hpp
│       ├── io/
│       │   ├── AEDTParser.hpp
│       │   └── NGMeshParser.hpp
│       ├── linalg/
│       │   ├── BiCGStabSolver.hpp
│       │   ├── MklPardisoSolver.hpp
│       │   └── SparseMatrix.hpp
│       └── post/
│           ├── OutputWriter.hpp
│           └── ResultExtractor.hpp
├── src/
│   ├── app/
│   │   └── Application.cpp
│   ├── core/
│   │   ├── Logger.cpp
│   │   ├── Math.cpp
│   │   ├── Types.cpp
│   │   └── Utilities.cpp
│   ├── fem/
│   │   ├── EdgeTopology.cpp
│   │   ├── FEMAssembler.cpp
│   │   └── PortModeSolver.cpp
│   ├── io/
│   │   ├── AEDTParser.cpp
│   │   └── NGMeshParser.cpp
│   ├── linalg/
│   │   ├── BiCGStabSolver.cpp
│   │   ├── MklPardisoSolver.cpp
│   │   └── SparseMatrix.cpp
│   ├── post/
│   │   ├── OutputWriter.cpp
│   │   └── ResultExtractor.cpp
│   └── main.cpp
├── docs/
│   ├── overview.md
│   ├── architecture.md
│   ├── primitives.md
│   ├── testing.md
│   ├── validation.md
│   ├── task-splitting.md
│   ├── development-workflow.md
│   ├── project-tree.md
│   ├── overview/
│   ├── architecture/
│   ├── primitives/
│   ├── testing/
│   ├── validation/
│   ├── task-splitting/
│   └── development-workflow/
├── build*/
└── results*/
```

`build*/` 和 `results*/` 表示一组生成目录，而不是源码目录。实际工程中可见的 `build`、`build_nomkl`、`build_mkl_probe` 以及多个 `results_*` 目录都属于这一类。

## 顶层文件

- `CMakeLists.txt`：定义 C++17 工程、`bp_fem_core` 静态库、`bp_fem_solver` 可执行程序，以及 `BPFEM_USE_MKL` 求解器开关。
- `README.md`：项目简介、构建命令、运行命令和主要文档入口。
- `wg_bp_filter.aedt`：AEDT 工程输入，提供材料、变量、扫频和波端口定义。
- `current.ngmesh`：NGMesh 网格输入，提供点、表面三角形、facet、body 和四面体单元。

## 源码目录

`include/bpfem/` 存放公开接口，`src/` 存放对应实现。除 `src/main.cpp` 外，源码按相同模块名称成对组织，便于从头文件定位实现。

- `app`：命令行参数解析和端到端运行流程。`Application.cpp` 负责读取输入、构造拓扑、装配矩阵、调用求解器和写出结果。
- `core`：共享数据类型、常量、数学函数、日志和工具函数。`Types.hpp` 是 AEDT、网格、求解结果和 S 参数等核心数据结构的集中定义处。
- `io`：外部输入解析层。`AEDTParser` 将 AEDT 文本转为 `ProjectDefinition`，`NGMeshParser` 将网格文件转为 `Mesh`。
- `fem`：有限元核心。`EdgeTopology` 管理 H(curl) 自由度，`PortModeSolver` 计算端口模态和归一化相关量，`FEMAssembler` 装配频域 curl-curl 系统和边界条件。
- `linalg`：线性代数后端。`SparseMatrix` 存储复数稀疏矩阵，`MklPardisoSolver` 提供 oneMKL PARDISO 路径，`BiCGStabSolver` 提供无 MKL 时的迭代 fallback。
- `post`：后处理输出。`ResultExtractor` 从边自由度投影得到 S 参数，`OutputWriter` 写出 CSV 和 VTU 文件。
- `main.cpp`：轻量入口，只调用 `fem::runApplication(argc, argv)`。

## 文档目录

`docs/` 采用“兼容入口 + 专题目录”的组织方式。根级 Markdown 文件保留为稳定入口，详细内容拆入同名子目录。

- `docs/overview.md` 与 `docs/overview/`：工程阅读入口、运行流、扩展点和安全规则。
- `docs/architecture.md` 与 `docs/architecture/`：模块布局、依赖方向、求解器选择和风险点。
- `docs/primitives.md` 与 `docs/primitives/`：物理、网格、FEM、线性代数、输出和不变量说明。
- `docs/testing.md` 与 `docs/testing/`：构建检查、CLI smoke test、输出文件检查和回归清单。
- `docs/validation.md` 与 `docs/validation/`：残差、端口、S 参数、场可视化和单位验证。
- `docs/task-splitting.md` 与 `docs/task-splitting/`：按任务类型拆分开发上下文。
- `docs/development-workflow.md` 与 `docs/development-workflow/`：编辑前、中、后检查和交接模板。
- `docs/optimization-roadmap.md` 与 `docs/optimization/`：从原型到商业级求解器的优化路线与分期计划。
- `docs/project-tree.md`：当前文件，用于快速识别项目树和目录职责。

## 生成目录

- `build/`：默认 CMake/Visual Studio 构建目录，可能包含 Debug 与 Release 产物。
- `build_nomkl/`：关闭 MKL 或用于无 MKL fallback 路径的构建目录。
- `build_mkl_probe/`：用于探测或验证 MKL/PARDISO 配置的构建目录。
- `results/` 与 `results_*`：求解输出目录，通常包含 `s_parameters.csv`、`field_last.vtu` 和按频点命名的 `field_<frequency>Hz.vtu`。

这些目录可由构建或运行命令重新生成。修改源码或文档时，应避免把生成产物当作维护入口。

## 数据流速览

```text
wg_bp_filter.aedt + current.ngmesh
        |
        v
AEDTParser + NGMeshParser
        |
        v
ProjectDefinition + Mesh
        |
        v
EdgeTopology + PortModeSolver
        |
        v
FEMAssembler -> SparseMatrix + RHS
        |
        v
MklPardisoSolver or BiCGStabSolver
        |
        v
ResultExtractor + OutputWriter
        |
        v
results*/s_parameters.csv + results*/field_*.vtu
```

## 维护约定

- 新增公共数据结构时，优先放入 `include/bpfem/core/Types.hpp`，并同步检查解析、装配、后处理是否需要更新。
- 新增模块时，保持 `include/bpfem/<module>/` 与 `src/<module>/` 的接口/实现对应关系，并在 `CMakeLists.txt` 中登记实现文件。
- 修改端口、边界条件、自由度拓扑或 S 参数投影时，应同时检查 `fem` 与 `post` 模块，避免装配和后处理口径不一致。
- 运行验证时，将输出写入新的 `results_*` 目录，便于和历史结果并排比较。
