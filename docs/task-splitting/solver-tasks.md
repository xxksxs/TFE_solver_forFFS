# 线性求解任务

## 目标

修改稀疏矩阵存储、直接求解器集成或 fallback 迭代求解器行为。

## 优先阅读文件

- `../primitives/linear-algebra.md`
- `../../src/linalg/SparseMatrix.cpp`
- `../../src/linalg/MklPardisoSolver.cpp`
- `../../src/linalg/BiCGStabSolver.cpp`
- `../../src/app/Application.cpp`

## 可能修改文件

- `../../include/bpfem/linalg/*.hpp`
- `../../src/linalg/*.cpp`
- `../../CMakeLists.txt`

## 验证命令

```powershell
cmake --build build --config Release
.\build\Release\bp_fem_solver.exe --out results_solver_check --max-sweep-points 1
```

涉及 fallback 或 CMake 时，还应运行：

```powershell
cmake -S . -B build_nomkl -DBPFEM_USE_MKL=OFF
cmake --build build_nomkl --config Release
```

## 风险提示

- 保持 `SolveResult` 语义稳定。
- MKL 和非 MKL 两条路径都必须可编译。
- 先看残差，再判断输出差异。
- 不要让求解器依赖应用层或 CLI 逻辑。
