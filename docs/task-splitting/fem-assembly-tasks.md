# FEM 装配任务

## 目标

修改体装配、材料处理、PEC 约束或波端口边界项。

## 优先阅读文件

- `../primitives/fem.md`
- `../primitives/invariants.md`
- `../validation.md`
- `../../src/fem/EdgeTopology.cpp`
- `../../src/fem/FEMAssembler.cpp`
- `../../src/fem/PortModeSolver.cpp`

## 可能修改文件

- `../../include/bpfem/fem/*.hpp`
- `../../src/fem/*.cpp`

## 验证命令

```powershell
.\build\Release\bp_fem_solver.exe --out result --max-sweep-points 1
```

## 风险提示

- 保持边方向一致：`localNodes` 应在生成时按全局点号小到大存储，体装配不要单独重排端点。
- 不要对波端口边施加 PEC 零约束。
- 修改端口边界后必须重新验证 S 参数。
- 修改材料项后检查单位和频率相关系数。
