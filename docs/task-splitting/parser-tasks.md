# 解析器任务

## 目标

新增或修复 AEDT、NGMesh 输入文件中的数据解析。

## 优先阅读文件

- `../../include/bpfem/core/Types.hpp`
- `../../src/io/AEDTParser.cpp`
- `../../src/io/NGMeshParser.cpp`
- `../primitives/README.md`
- `../primitives/mesh.md`

## 可能修改文件

- `../../include/bpfem/core/Types.hpp`
- `../../src/io/AEDTParser.cpp`
- `../../src/io/NGMeshParser.cpp`

## 验证命令

```powershell
.\build\Release\bp_fem_solver.exe --max-sweep-points 1
```

## 风险提示

- 不要静默改变单位解释。
- 解析默认值必须明确。
- 共享类型变化后必须更新下游消费者。
- faceId 解析变化会影响端口和 PEC 边界。
