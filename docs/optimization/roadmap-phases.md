# 分期路线图

本文把所有专题文档中的任务汇总为 4 期分阶段路线图，并给出里程碑、交付与依赖关系。

## 阶段总览

| 阶段 | 周期估算 | 主题 | 关键交付 |
|------|---------|------|---------|
| **Phase 1** | 1 ~ 2 季 | 数值核心可信 + 基本输出 | 二/三阶 hierarchical 基函数、PML/ABC、误差估计 + h-adapt、AMS 预处理、Touchstone/远场/Y-Z 输出、最小基准库与 CI |
| **Phase 2** | 2 ~ 3 季 | 工程通用性 | 多端口/多模 + 集总端口 + 平面波 + 电压源、各向异性 + Debye/Drude/Lorentz、AWE/Krylov MOR、CAD 导入 + 内置 Netgen、OpenMP 装配、Python API |
| **Phase 3** | 3 ~ 4 季 | 大规模并行与多物理 | DDM/FETI-DP、MPI、GPU 装配/求解、特征模、EM-Thermal、瞬态 FFT 反演、参数扫描与全局优化 |
| **Phase 4** | 4 ~ 6 季 | 产品化 | 工程文件 + GUI、检查点续算、伴随梯度优化、DGTD 原型、许可机制、商业级 CI 与基准库 |

## Phase 1：可信核心（M1 ~ M6）

### 目标
任何用户用本工程跑出的 S 参数与 HFSS / CST 在 0.1 dB 以内重合，并自动诊断收敛性。

### 必交付
- **数值方法**：
  - 二阶 H(curl) hierarchical 基函数（Demkowicz / Schöberl 风格）
  - 通用 Gauss 求积
  - 残差型 a posteriori 估计
  - h-adapt 单轮自适应
- **边界与激励**：
  - PMC / 对称面
  - 一阶 ABC
  - UPML 张量层
  - 多模波端口（Top-N 模）
- **求解器**：
  - GMRES + ILU(0)
  - hypre AMS 预处理（首次集成）
  - 自适应频点加密扫频
- **后处理**：
  - Touchstone v1 / v2 输出
  - Y/Z/VSWR/群时延
  - 多端口 N×N S 矩阵
- **可用性**：
  - CLI11 + YAML 配置
  - 结构化 JSON 日志 + summary HTML
- **测试**：
  - Catch2 + 单元测试覆盖 ≥ 60%
  - 5 个基准（矩形腔本征、矩形波导 cutoff、矩形波导 S21、同轴线 TEM、bp_filter 端到端）
  - GitHub Actions CI 矩阵（4 配置）
  - 收敛性回归框架

### 里程碑
- **M1（4 周）**：CLI11 / YAML / Catch2 / GitHub Actions / 单元测试 60%。
- **M2（8 周）**：二阶 H(curl) + 通用 Gauss + 矩形波导 / 矩形腔基准。
- **M3（12 周）**：PML/ABC/PMC + 多模波端口 + 多端口 S 矩阵。
- **M4（16 周）**：GMRES + AMS + 自适应频点。
- **M5（20 周）**：残差估计 + h-adapt 单轮 + Touchstone/Y/Z 输出。
- **M6（24 周）**：收敛回归 CI 全绿 + summary HTML。

### 退出标准
- 5 个基准：本工程 vs 解析或文献误差 < 0.5%。
- bp_filter S21 vs HFSS Direct < 0.1 dB（至少 21 个频点）。
- 收敛回归在 PR 上自动跑前 3 个基准并通过。

## Phase 2：工程通用性（M7 ~ M12）

### 目标
能解决工业实际问题：多端口多激励、典型材料、CAD 直接导入、扫频时间显著缩短。

### 必交付
- **数值方法**：
  - 三阶 H(curl) + 单元独立 p
  - hp-adapt 策略
  - 曲面 / 曲边 (geometry p2)
- **边界与激励**：
  - 集总端口 + Z₀ 重归一
  - 入射平面波 (SF/TF) 激励
  - 阻抗 / SIBC 边界
  - Floquet / 周期主从
- **材料**：
  - 各向异性 ε / μ + 复非对称 PARDISO 路径（mtype 13）
  - Debye / Drude / Lorentz 多极色散
  - AEDT 频率相关 / 张量解析
  - 表面阻抗 / 薄层
- **求解器与扫频**：
  - AWE 单点 Padé
  - Krylov MOR 多点
  - 自适应 + AWE 混合
  - 多 RHS 一次求解
- **几何与网格**：
  - OCCT STEP/IGES 导入
  - Netgen 网格生成集成
  - CAD faceId ↔ 网格 facetId 持久化
  - 网格质量度量 + 报告
- **后处理**：
  - 远场方向图 + 增益 + 轴比
  - RCS 单 / 双站
  - 表面电流 + 体损耗
  - 端口去嵌 + 重归一化
  - HDF5 / XDMF 体场
- **并行**：
  - OpenMP 体装配 + atomic
  - 多线程 PARDISO 调优
  - 跨频点 TBB 并行
- **可用性**：
  - 工程文件格式 + 输出布局
  - Python API (pybind11)
  - 检查点 / 续算

### 退出标准
- 同轴线 / 微带线 / SRR / Yagi 天线四个新基准与文献误差 < 1%。
- 21 频点扫频在 16 核机器上 < 5 min（中等规模 1M DOF 模型）。
- Python 脚本可完成完整 import → solve → S/远场 输出流程。

## Phase 3：大规模并行与多物理（M13 ~ M18）

### 目标
能在多节点集群上跑 10M+ DOF 模型，支持基本多物理与全局优化。

### 必交付
- **求解器**：
  - DDM / FETI-DP 内核
  - MPI 跨节点频点并行
  - cuSPARSE / cuDSS GPU 求解
  - 谐振腔本征 (SLEPc / 自实现)
  - 特征模 (CMA)
- **多物理**：
  - 节点 H¹ 热 FEM 求解器
  - EM-Thermal 松耦合
  - Touchstone-SPICE 联仿
  - Krylov MOR → 状态空间
  - 参数扫描 + 全局优化 (nlopt / pagmo)
- **几何与网格**：
  - 非协调网格 + hanging-node 约束
  - METIS 网格分区
- **后处理**：
  - Q 值 + 谐振检测
  - 传输线 Z_c, γ 提取
  - Tecplot / FFE / FFS 兼容
- **可用性**：
  - C API SDK
  - 作业调度模板
  - Web 仪表盘
- **测试**：
  - 求解器三路一致性回归
  - 性能 baseline + 回归
  - 全量基准 ≥ 15 个

### 退出标准
- 10M DOF 模型在 8 节点 × 64 核完成扫频。
- EM-Thermal 单步收敛在 5 次外迭代以内。
- 全局优化能驱动滤波器调谐到目标 S21 < −30 dB。

## Phase 4：产品化（M19 ~ M30）

### 目标
最终产品形态：完整 GUI + 优化 + 瞬态 + 商业级稳定性。

### 必交付
- **数值方法**：
  - 完整 hp 自适应（光滑度选 h/p）
  - DWR 目标量误差估计
- **多物理与时域**：
  - DGTD 原型（PEC + 平面波 + ABC/PML）
  - TD-FEM 隐式（实验性）
  - 拓扑优化（Density-based）
  - 伴随梯度优化（频域）
- **可用性**：
  - GUI（Qt + VTK / ParaView 集成）
  - 许可机制
  - 商业级 CI + 全量基准 ≥ 30 个

### 退出标准
- 普通工程师从 STEP 文件到 S 参数报告全程在 GUI 内 < 30 min。
- DGTD 在 PEC 球散射上与解析 RCS 误差 < 1 dB。
- 全量基准在 CI 上 8 小时内跑完。

## 依赖关系图

```
                Phase 1                       Phase 2                Phase 3              Phase 4
  Test/CI ──────► ────────────────────► ────────────────────► ─────────► (extend)
  Numerical p2 ──► p3 + hp ──────────────────► hp + DWR ────► DGTD/TD-FEM
  BCs PMC/ABC/PML ──► Lumped/PlaneWave/SIBC ──► (refine)
  Materials const ──► tensor/dispersion ──► gyromagnetic ────► nonlinear
  Solvers PARDISO+GMRES+AMS ──► AWE/MOR ──► DDM/MPI/GPU ───► adjoint/topo
  Geometry external ngmesh ──► OCCT+Netgen ──► non-conforming ──► full GUI
  Post S/Y/Z/Touchstone ──► farfield/RCS/probe ──► circuit/MOR-state ──► full GUI
  Usability CLI/YAML ──► Project/PythonAPI ──► CAPI/Web ──► GUI/License
```

## 风险与里程碑回顾

每月召开里程碑评审：

1. 上月已完成基准 vs 本月新增基准。
2. 性能 baseline 偏离。
3. 跨模块不变量是否被破坏（参考 `primitives/invariants.md`）。
4. 第三方依赖（MKL / hypre / Netgen / OCCT / pybind11）版本与许可。
5. 文档同步（每个新增模块/选项必须有对应专题文档）。
