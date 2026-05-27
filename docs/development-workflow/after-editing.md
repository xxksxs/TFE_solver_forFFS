# 编辑后检查

修改完成后应运行与任务范围匹配的最小验证序列。

## 通用检查

```powershell
cmake --build build --config Release
.\build\Release\bp_fem_solver.exe --help
.\build\Release\bp_fem_solver.exe --out results_check --max-sweep-points 1
```

## 求解器或 CMake 修改

涉及求解器选择、MKL、CMake 或 fallback 路径时，还应运行：

```powershell
cmake -S . -B build_nomkl -DBPFEM_USE_MKL=OFF
cmake --build build_nomkl --config Release
```

## 文档修改

纯文档修改不需要编译，但应检查：

- 文件路径是否存在。
- 内部链接是否指向正确位置。
- 是否仍有旧路径或乱码。
- 顶层索引是否同步更新。
