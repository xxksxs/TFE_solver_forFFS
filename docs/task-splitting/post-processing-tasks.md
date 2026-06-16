# 后处理任务

## 目标

修改 S 参数提取或 VTU/CSV 输出。

## 优先阅读文件

- `../primitives/output.md`
- `../primitives/invariants.md`
- `../../src/post/ResultExtractor.cpp`
- `../../src/post/OutputWriter.cpp`
- `../../src/fem/PortModeSolver.cpp`
- `../validation.md`

## 可能修改文件

- `../../include/bpfem/post/*.hpp`
- `../../src/post/*.cpp`

## 验证命令

```powershell
.\build\Release\bp_fem_solver.exe --out result --max-sweep-points 1 --write-all-fields
```

## 风险提示

- S 参数依赖与装配一致的端口模式约定。
- VTU 场来自边自由度重构。
- 不要把边自由度当作节点标量值。
- 输出异常时先检查残差和端口激励。
