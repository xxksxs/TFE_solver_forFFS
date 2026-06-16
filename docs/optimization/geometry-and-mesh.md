# 几何与网格

本文给出从"依赖外部 NGMesh"演进到"内置几何核 + 网格生成 + 自适应"的方案。

## 1. 现状

- **几何输入**：无 CAD，依赖手动准备 `current.ngmesh`。
- **网格**：仅四面体一阶单元；背景体被 `NGMeshParser` 移除。
- **网格质量**：未计算单元品质（aspect ratio、shape factor、min dihedral angle）。
- **自适应**：无。
- **曲面/曲边**：无（影响曲面共形结构精度）。

## 2. 几何升级

| 项 | 商业级要求 | 实现方案 |
|----|-----------|---------|
| CAD 导入 | STEP / IGES / Parasolid / ACIS / SAT / 3MF | 集成 OpenCascade (OCCT) 提供 STEP/IGES；商业版可考虑 Parasolid SDK |
| 几何修复 | 缝合 / 容差合并 / 自相交检测 | OCCT `ShapeFix_*` |
| 布尔运算 | union / cut / common | OCCT BOPAlgo |
| 拓扑命名 | 持久面 ID | OCCT TNaming 或自管理面 ID 表 |
| 几何参数化 | 设计变量 / 约束 | 内置 expression engine |

模块布局建议：

```
include/bpfem/geometry/
  Geometry.hpp        // 几何顶层接口
  CadImporter.hpp     // STEP / IGES
  Booleans.hpp
  FaceIdRegistry.hpp  // CAD 面 → faceId 持久映射
src/geometry/
  ...
```

## 3. 网格生成升级

| 项 | 商业级要求 | 实现方案 |
|----|-----------|---------|
| 内置 4 面体网格 | 3D Delaunay + advancing-front | 集成 [Netgen](https://github.com/NGSolve/netgen) 库 (LGPL) 或 Gmsh API |
| 表面网格 | 共形 + 曲面采样 | OCCT BRepMesh 或 Netgen surface mesher |
| 曲面 / 曲边单元 | p2 几何 + 等参变换 | 与基函数升阶联动（`numerical-methods.md`） |
| 网格质量度量 | aspect ratio、min/max angle、Jacobian | 新增 `mesh/Quality.{hpp,cpp}` |
| 网格自适应 | h-refine（red-green / longest-edge bisection）+ 平滑 | 新增 `mesh/Refiner.{hpp,cpp}` |
| 非协调网格 | hanging node + H(curl) 约束 | `EdgeTopology` 升级 (`numerical-methods.md`) |
| 网格分区 | METIS / SCOTCH / ParMETIS | DDM 时（Phase 3） |
| 多域共形 | 多体共享面 + 共形细化 | NGMesh 已部分支持，但需要明确多体接口 |

落地策略：

- **Phase 1 末**：包装 Netgen 作为可选依赖；当用户没有提供 `.ngmesh` 时，自动从 STEP 文件做几何 + 网格生成。
- **Phase 2**：曲面单元 + 各单元独立 p（与高阶基函数同时落地）。

## 4. faceId 与 CAD 持久化

- `Mesh` 中现有 `surfaceTriangles[i].id` / `facetId` 是网格级别的 ID。
- 引入 *CAD faceId → 网格 facetId* 双向表；CAD 几何变化时，AEDT 端的 port faceId 仍能匹配网格 surface triangles。
- 这是商业级求解器最容易出错的地方，`primitives/invariants.md` 已经强调"faceId 必须匹配"。

## 5. 网格自适应工作流

```
初始网格 (Netgen / 外部)
   │
   v
装配 + 求解 (FEM)
   │
   v
误差估计 (residual / DWR)  ──┐
   │                          │
   v                          │
标记单元 (Dörfler / Greedy)   │ 直到目标量收敛
   │                          │
   v                          │
h / p / hp 细化               │
   │                          │
   └──────────────────────────┘
```

- 误差估计与细化策略详见 `numerical-methods.md`。
- 自适应必须保证 `EdgeTopology` 与 `PortModeSolver` 在重新装配时一致；建议每轮 adapt 重建 topology。

## 6. 网格诊断

- 新增 `--mesh-report result/mesh_report.json`：包含点数、四面体数、min/max edge length、aspect ratio 直方图、退化单元列表。
- 在 *AEDTParser* 加载完成后立即对比 AEDT port faceId 与网格 facetId，发现不匹配立刻报错（当前可能在求解阶段才被发现）。

## 7. 任务清单

| 任务 | 优先级 | 阶段 | 估算 |
|------|--------|------|------|
| Netgen 网格生成集成 | 高 | P2 | L |
| OCCT STEP/IGES 导入 | 高 | P2 | L |
| CAD faceId ↔ 网格 facetId 持久化 | 高 | P2 | M |
| 网格质量度量 + 报告 | 高 | P1 末 | S |
| h-refine（longest-edge bisection） | 高 | P1 ~ P2 | M |
| 曲面单元 p2 / 曲边 | 中 | P2 | L |
| 非协调网格 + hanging-node 约束 | 中 | P2 | L |
| METIS 网格分区 | 中 | P3 | M |
| 多体共形与接口检测 | 中 | P2 | M |
