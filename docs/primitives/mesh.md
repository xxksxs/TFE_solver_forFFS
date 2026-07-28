# 网格原语

本文说明 NGMesh 输入在工程中的核心网格概念。修改 `NGMeshParser`、`EdgeTopology` 或边界条件前应先阅读本文。

## 点

点是三维坐标，存储在 `Mesh::pointsById` 中，并通过整数 point id 引用。

注意事项：

- 坐标单位必须与频率、波数计算保持一致。
- 不要假设 point id 连续。
- 访问点坐标时应通过 `pointsById` 查找。

## 表面三角形

表面三角形表示边界上的三角面片，通常包含：

- 三个点 ID。
- 所属 faceId。
- 与端口、PEC 或其他边界条件相关的边界语义。

端口模式计算依赖端口 faceId 上的表面三角形集合。

## 四面体

四面体是体网格单元，用于 Nedelec/Whitney 棱元装配。

每个四面体产生：

- 6 条局部边。
- 局部边到全局边的映射。
- 零阶棱元的 6 个局部基函数，或一阶层次棱元的 20 个局部基函数。
- 单元体积、梯度、curl-curl 项和质量项。

NGMesh 中 `body_name background` 对应 HFSS 背景网格。该背景体只用于导出环境包围盒，不参与 FEM 计算；解析结束后会从 `Mesh::tetrahedra`、`Mesh::bodies` 和未引用点中移除。

## faceId

faceId 是连接 AEDT 边界定义与 NGMesh 表面三角形的关键字段。

HFSS 波端口既可以直接引用 Faces(...)，也可以引用单独的 sheet object。
后一种写法由 AEDTParser 保留 objectId，再由 PortFaceResolver 使用 NGMesh
body 包围盒与表面三角网格包围盒做唯一几何匹配，最终仍转换成 FEM 内核统一
使用的 faceId。若没有匹配或出现多个重合面，解析阶段会直接报错，不会带着
无效端口进入矩阵装配。

必须满足：

- AEDT 中定义的端口 faceId 能在 NGMesh 表面三角形中找到。
- 非端口边界可以被 PEC 约束处理。
- 端口边界不能被 PEC 零约束覆盖。

## Facet 信息

Facet 信息是 NGMesh 中关于边界面和体的附加元数据。当前主要用于帮助解析和定位边界，但后续也可扩展为更细粒度的边界分类。

## 网格相关修改清单

- **新增边界类型**：更新共享类型和边界分类逻辑。
- **修改 faceId 解析**：同步验证 AEDT port faceId 是否仍能匹配。
- **修改坐标单位处理**：同步验证截止频率、传播常数和场输出。
- **修改拓扑构建**：同步检查 `EdgeTopology`、`FEMAssembler` 和 `OutputWriter`。
