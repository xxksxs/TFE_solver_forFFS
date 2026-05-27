# 数值方法升级

本文给出把 FEM 数值核心从 0/1 阶 Nedelec + 静态网格演进到商业级"高阶 + hp 自适应 + 数值积分稳健"的完整任务列表。

## 1. 高阶 H(curl) 向量元

**现状**：`EdgeTopology` 提供 `LocalDofKind::EdgeZero / EdgeFirst / FaceFirst0 / FaceFirst1`，每四面体 6（0 阶）或 20（1 阶）局部基函数。

**目标**：

- 至少支持 `p = 0, 1, 2, 3` 的 hierarchical Nedelec 第二族向量元（Demkowicz / Schöberl / Webb）。
- 支持四面体内 *p 可变*（每个单元独立阶数），为 hp 自适应铺路。
- 支持端口三角面元 0 ~ p 阶 H(curl) 投影。
- 引入 *interior bubble* 自由度（cell-based）以支持 p ≥ 2，并在 `EdgeTopology` 中加入 `LocalDofKind::Cell*`。

**关键改动**：

- `EdgeTopology`：扩 `LocalDofKind`，按阶数分配边/面/体 DOF；保留全局点号小到大的边方向不变量。
- `FEMAssembler`：用通用四面体 H(curl) 基函数生成器替换写死的 6/20 路径；要求局部矩阵自动按基函数族大小分配。
- `PortModeSolver`：端口面 2D H(curl) 投影同步升阶。
- `ResultExtractor`、`OutputWriter`：按高阶基函数重构端口投影与点场。

**风险**：

- 高阶基函数对 `localNodes` 排序敏感，必须按 `primitives/invariants.md`「边方向」一次性升级。
- 端口模式必须同步升阶，否则会出现"体内场是 p 阶但端口投影只有 1 阶"的不一致，S 参数会有 0.1 ~ 1 dB 系统误差。

## 2. hp 自适应

**目标**：

- **残差型 a posteriori 误差估计**：实现单元残差 + 面跳跃估计（Beck-Hiptmair-Hoppe-Wohlmuth 风格）。
- **目标量误差估计 (DWR)**：以 S 参数 / 端口耦合等为目标量，求伴随场，计算单元贡献。
- **细化策略**：
  - h 细化：四面体 longest-edge bisection 或 red-green refinement，保持非协调网格的 H(curl) hanging-node constraint。
  - p 细化：单元独立升阶。
  - hp：根据光滑度指示器（Texas 三选一、analyticity estimator）选 h 或 p。

**关键改动**：

- 新增 `fem/ErrorEstimator.{hpp,cpp}` 与 `fem/RefinementStrategy.{hpp,cpp}`。
- `EdgeTopology` 支持 *constraint matrix* 把 hanging-node DOF 映射到无约束基。
- `Application.cpp` 增加 `--adapt-passes N`、`--target-error 1e-3`、`--target-quantity S21|S11|none`。

**验证**：以矩形腔本征频率收敛、矩形波导 cutoff、SRR 等基准检测 h 与 hp 自适应的收敛阶。

## 3. 数值积分稳健化

**现状**：体装配按四面体常梯度 + 解析公式（零阶 Whitney）。

**目标**：

- 通用四面体 *Gauss / Grundmann-Möller / Stroud* 求积，按基函数阶 `2p + 2` 自动选阶。
- 端口三角形求积同步通用化（当前 `PortMode::quadrature` 仅 area/3 节点权重）。
- 对曲边/曲面单元（Phase 2 几何升级）使用等参变换 + Jacobian 求积。

**关键改动**：

- 新增 `fem/Quadrature.{hpp,cpp}` 提供四面体/三角形求积规则表。
- `FEMAssembler` 与 `PortModeSolver` 通过统一接口取求积点。

## 4. 数值条件稳健化

- **basis scaling**：高阶 hierarchical 基函数应做合适缩放（边按边长，面按面积）以避免条件数随 p 增长。
- **质量集中可选**：低频或弱阻抗模式下提供 lumped mass，以加速 PEC 谐振诊断。
- **复对称 vs Hermitian**：当前 PARDISO 用复对称（`mtype = 6`）。当导入有耗 / 各向异性张量时，需要支持 *complex symmetric* 与 *complex non-symmetric* 双路径，并在 `MklPardisoSolver` 中按矩阵类型切换 `mtype`。

## 5. DOF 编号与稀疏图复用

**现状**：`FEMAssembler` 已缓存 `SparsePattern`，扫频时复用 symbolic。

**目标**：

- 加入 *RCM / nested dissection* 重排序（直接交给 PARDISO `iparm[1]=3` 的 OpenMP nested dissection）。
- 高阶时 `SparsePattern` 应在升阶后自动重建，并在 `FEMAssembler::sparsityPattern_` 缓存键中加入阶数与约束 DOF 集合。
- 提供 `--reorder amd|nd|metis` 选项。

## 6. 风险与验收

- 升阶 / 自适应必须保留 `SolveResult` 契约。
- 升阶后必须先在矩形波导 cutoff 误差 < 0.1% 和矩形腔本征频率误差 < 0.05% 上验收，再合并。
- 自适应必须在 *S21 收敛曲线* 上展现单调性；否则视为失败。

## 7. 任务列表（Phase 1 ~ 2）

| 任务 | 优先级 | 估算 | 依赖 |
|------|--------|------|------|
| 高阶 H(curl) 基（p ≤ 2） | 高 | L | `EdgeTopology` 改造 |
| 高阶 H(curl) 基（p = 3 + 可变 p） | 中 | L | 上一项 |
| 通用 Gauss 求积 | 高 | M | 高阶基 |
| Hanging-node H(curl) 约束 | 中 | L | 高阶基 |
| 残差型 a posteriori 估计 | 高 | M | 高阶基 |
| DWR 目标量误差估计 | 中 | L | 残差型 + 伴随求解 |
| h-adapt 流程 | 高 | M | 估计器 |
| hp-adapt 策略 | 中 | L | h-adapt + p 可变 |
| 重排序 / 复对称-Hermitian 双路径 | 中 | M | 求解器层 |
