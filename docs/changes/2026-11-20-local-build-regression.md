# 本地构建与回归入口

新增 [PowerShell 回归入口](../../scripts/run-local-regression.ps1)，按配置隔离构建目录，运行现有三个 CTest，并可调用独立 CLI 烟测脚本。[操作和依赖配置](../testing/build-checks.md)说明如何选择原生迭代、oneMKL 或 oneMKL＋ZMUMPS。

中文 MSVC 的 `/showIncludes` 输出与 CMake 原先探测后写入 Ninja 的前缀编码不一致。配置现在使用实际编译参数探测前缀并保留原始字节。此前漏记头文件依赖的构建目录需重新建立依赖基线。

`mkl` 和 `mkl-mumps` 配置要求实际找到 oneMKL，手工 CMake 配置仍可允许 fallback。顺序 ZMUMPS 的库、链接目录及 Windows DLL 由本机初始缓存提供。

验证使用中文 MSVC 19.51、Ninja 1.13.2、CMake 4.2.3、oneMKL 2025.3 与顺序 ZMUMPS 5.9.2。原生和 oneMKL 各验证 Debug 与 Release，oneMKL＋ZMUMPS 验证 Release。各组均通过三个 CTest。强制隐藏 oneMKL 时，`mkl` 配置在阶段入口报错。

在隔离源码中改动 `ISparseSolver.hpp`，重复执行 `native / Release`，`GmresBackend.cpp` 和 `SparseSolverFactory.cpp` 重新编译，`NGMeshParser.cpp` 的对象未变，增量构建后再次通过三个 CTest。烟测接口用当前构建的 `bp_fem_solver --help` 验证参数传递和独立输出目录。
