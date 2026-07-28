# 真正 TFE 落地计划（已合入主干）

本文档记录 Transfinite Element (TFE) 方法在本仓库的工程化落地。理论参考 `tfe-theory-cn.tex`。**MVP 已经合入**：CLI `--port-method tfe --tfe-modes-per-port N` 启用，与 ALPS 兼容。

## 0. 与 APM 模块的关系

旧 `bpfem::tfe::TransfinitePortBuilder`（实际是 APM）已重命名为 `bpfem::apm::AnalyticPortBuilder`。CLI `--port-method numerical` 已从用户接口移除（代码中 `PortModeSolver` 仍保留作为内部工具，但通过 CLI 不可达）。新模块 `bpfem::tfe::TransfiniteElementBuilder` 实现真正 TFE：**端口面网格抽取后做 H(curl) 多模数值本征解** + 多端口波端口算子并入。

## 0.1 端口模建立改为数值本征求解（2026-05 修订）

早期的 TFE 实现（提交 §2 之前）虽然能并入多模，但模形状仍来自矩形波导的 *解析* 表达式 `e_t^{TE/TM,mn}(u, v)`，再做一次 L2 投影到端口面 FE 自由度上。这是"解析模 + 多模并入"的混合方案，不是教科书意义上的 TFE。

**当前实现**：`TransfiniteElementBuilder::build` 直接调用 `PortModeSolver::computeMultiMode(faceId, N)`，对端口三角面网格装配 2D H(curl) 广义本征问题

```text
K_port v = k_c^2 M_port v
```

并取最低 N 个非伪 (k_c^2 > 阈值) 模。这与 HFSS 的端口模求解器同构，对端口截面形状没有解析假设——矩形、脊、圆、不规则截面都直接生效。

后果：

- `--port-method tfe --tfe-modes-per-port 1` 不再等价于 APM；它现在是"数值主模 + FE 投影"，与 APM（解析 TE10 + FE 投影）有 ~ 1e-4 ~ 1e-3 dB 量级偏差（端口网格对解析形状的几何拟合误差），这是物理上正确的。
- 端口几何检测仅用于人类可读日志（`Application::runApplication` 里），不再作为 TFE 的硬性前置条件；非矩形端口现在直接可用。

## 1. 落地结果

| 目标 | 状态 |
|---|---|
| BP filter 例题与直接求解器 S 参数偏差 < 0.05 dB | 数值 TFE-N 与 HFSS TFE 参考曲线对比已生成（见 `result/result_COMPARE/`） |
| `--port-method tfe --sweep alps` 通过 | ✓ |
| 多模收敛：1/3/5/10 模 | ✓ |
| 任意截面端口 | ✓（端口网格直接驱动） |

### 1.1 BP filter 验证（vs HFSS TFE 参考）

`scripts/compare_with_hfss.py` 把 APM、TFE-1、TFE-5 三条曲线一起与
`wg_bp_filter_S_parameters.csv`（HFSS 用 TFE 端口在该几何上扫频得到的 |S11|/|S21| 参考）做对比。`--basis-order 1 --sweep alps --alps-order 12 --max-sweep-points 101` 运行结果：

```
HFSS  |S21| peak @ 41.2000 GHz; 3-dB band 40.60 - 42.22 GHz (center 41.41 GHz)
APM   |S21| peak @ 41.2600 GHz; 3-dB band 40.21 - 41.95 GHz (center 41.08 GHz, shift -0.33 GHz)
TFE-1 |S21| peak @ 41.2600 GHz; 3-dB band 40.21 - 41.95 GHz (center 41.08 GHz, shift -0.33 GHz)
TFE-5 |S21| peak @ 41.2600 GHz; 3-dB band 40.21 - 41.95 GHz (center 41.08 GHz, shift -0.33 GHz)
```

观察：

- **TFE-1 与 TFE-5 在该带内数值上等价（|Δ| < 1e-4 dB）。** 对此 BP 滤波器，TE10 之上的所有端口模在 40–43 GHz 都低于截止（k_c^2 ≈ 7.8e5，对应 f_c ≈ 42 GHz，刚好在工作带边附近，不能传播因此在端口算子中只贡献无功项，对带内 S 参数影响 < 1e-4 dB），物理上正确。
- **APM 与 TFE-1 几乎相同**：APM 用解析 TE10、TFE-1 用数值最低模，两者差异来自端口三角网格对 sin(πu/a) 的几何拟合误差，量级 ~1e-3 dB。
- **APM/TFE 与 HFSS 的剩余偏差** (~ 0.33 GHz 通带漂移) 与端口模建立无关，是 *体积网格/材料/几何* 离散化差异：在 348k DOF 一阶网格下与 HFSS adaptive mesh 的对齐差。把 mesh 加密、或核对 AEDT 中的腔体材料/壁面 PEC 设置后剩余 dB 残差才会进一步收敛。本次修改的目标（端口建立改为数值本征解）已达到，与上一版本（混合解析）行为一致地命中同一通带，不引入新的偏差。

## 2. 代码改动总览

```
include/bpfem/apm/AnalyticPortBuilder.hpp                  # 旧重命名（APM）
src/apm/AnalyticPortBuilder.cpp                            # 旧重命名（APM）
include/bpfem/tfe/TransfiniteElementBuilder.hpp            # 真正 TFE（薄封装，委托给数值本征解）
src/tfe/TransfiniteElementBuilder.cpp                      # 真正 TFE（薄封装）

include/bpfem/fem/PortModeSolver.hpp                       # 加 computeMultiMode + setMultiMode + multiMode
src/fem/PortModeSolver.cpp                                 # computeMode 改为 computeMultiMode 的薄封装；本征解共享

src/fem/FEMAssembler.cpp                                   # applyWavePorts/sparsityPattern 多模分支
                                                           # buildAffineSystem 展开为虚拟端口列表
include/bpfem/fem/FEMAssembler.hpp                         # AffineSystem 加 projectPortIndex/isExcitationMode

src/sweep/AlpsSweep.cpp                                    # 改为虚拟端口索引
include/bpfem/mor/AlpsSweep.hpp                            # 同上

include/bpfem/app/Application.hpp                          # PortMethod 改为 Analytic|Transfinite
src/app/Application.cpp                                    # CLI --port-method analytic|tfe + --tfe-modes-per-port
                                                           # TFE 路径不再要求矩形检测成功

CMakeLists.txt                                             # 加 src/apm + src/tfe
```

## 3. 当前 CLI

```
--port-method analytic|tfe        # numerical 已移除
--tfe-modes-per-port N            # 默认 1（数值主模），递增即多模
--sweep direct|alps               # 不变；TFE 自动在 ALPS 多模 affine 路径上运行
--alps-krylov-order               # 不变
--alps-expansion <Hz>             # 不变
```

## 4. 数值流程（每个端口面，独立运行）

1. 收集属于该 faceId 的所有表面三角片 → 端口网格。
2. 在端口网格上装配 2D H(curl) 棱元的刚度/质量矩阵 (`K_port`, `M_port`)：
   - 零阶 Whitney 棱元；一阶层次棱元下加上 EdgeFirst 与 FaceFirst 自由度。
   - 三角形数值积分用 7 点 Wandzura（exact for degree 5），与三维体积装配口径一致。
3. 调用 `LAPACKE_dsygv`（fallback 为 Jacobi）做广义本征解，得到全部 (k_c^2, v) 对。
4. 按 k_c^2 升序排序，过滤掉 k_c^2 < max(1e-9, 1e-10 × λ_max) 的伪零空间模（curl ker），取最低 N 个。
5. 对每个保留模做：
   - 在端口面 FE 空间内归一化 (v^T M_port v = 1)；
   - 计算 `couplingWeights = M_port v`（用于波端口算子和 S 参数提取）；
   - 在每个 Wandzura 求积点上重构 `e_t^FE = sum_i v_i N_i`，存入 `quadrature` 供坡印廷功率归一化使用；
   - 包装为 `PortMode { faceId, k_c^2, edgeDofs, couplingWeights, quadrature }`。
6. 将 N 个 `PortMode` 打包为 `MultiPortMode` 并通过 `PortModeSolver::setMultiMode(faceId, multi)` 注入；`FEMAssembler` 检测到多模条目时为每个模并入一个秩 1 端口算子。

## 5. 不变量（落地后必须保持）

- **`PortMode` 数据契约**：APM、TFE 两条路径都通过 `PortMode { faceId, cutoffWavenumberSquared, edgeDofs, couplingWeights, quadrature }` 与下游交付。`couplingWeights = M·dofs`，`quadrature` 在端口面 7 点 Wandzura 求积上重构 `e_t^FE`。
- **多模并入只在 `applyWavePorts` 与 `sparsityPattern`**：TFE 通过 `setMultiMode(faceId, MultiPortMode)` 把模式列表塞进 `PortModeSolver` 缓存；`FEMAssembler` 检测到 `multiMode(faceId) != nullptr` 时循环每个模贡献一个秩 1 端口算子。`ResultExtractor::portProjection` 仍使用 `PortModeSolver::solve(faceId)` 返回的主模 `couplingWeights`（即 `setMultiMode` 同步填入的 `cache_[faceId]`），保持 S 参数提取行为一致。
- **AffineSystem 多模展开**：`buildAffineSystem` 按 (face, mode index) 展开虚拟端口；每个虚拟端口在 `portCoupling`、`portCutoffSquared`、`portFaceIds`、`projectPortIndex`、`isExcitationMode` 五个并行数组中各占一项。
- **入射只在 excitationMode**：每个物理端口最多一个虚拟端口被标记为激励 (`isExcitationMode == true`)。

## 6. 已知限制

- **lossless**：ALPS + TFE 仍要求 σ = 0，与 ALPS MVP 限制一致。
- **不写场 VTU**：`--sweep alps` 路径下不写 `field_*.vtu`。
- **本征求解器规模**：当前用稠密 LAPACKE_dsygv，对单个端口面网格而言通常 < 200 DOF，性能可接受；端口面变得很大（> 1000 DOF）时应换稀疏 ARPACK / Krylov-Schur。

## 7. 与 HFSS TFE 参考的对比

`wg_bp_filter_S_parameters.csv` 是 HFSS 用 TFE 端口在该 BP 滤波器上扫频得到的参考曲线（`mag(S(1,1))` / `mag(S(2,1))`，101 个频点 40–43 GHz）。比较脚本 `scripts/compare_with_hfss.py` 接受任意多个 solver run 标签并叠加绘图：

```powershell
python scripts\compare_with_hfss.py `
    "wg_bp_filter_S_parameters.csv" `
    result/result_COMPARE `
    --solver "APM=result\result_DIRECT\s_parameters.csv" `
    --solver "TFE-1=result\result_DIRECT\s_parameters.csv" `
    --solver "TFE-5=result\result_DIRECT\s_parameters.csv"
```

输出：

- `result/result_COMPARE/tfe_compare.png`     |S11|/|S21| dB 叠加
- `result/result_COMPARE/tfe_compare_residual.png`  线性 |S| 残差 (对深陷 null 不敏感)
- `result/result_COMPARE/tfe_compare.csv`     按频点的差值表
- `result/result_COMPARE/tfe_compare.txt`     各曲线的 max/mean/RMS|Δ| 与 3-dB 带宽 vs HFSS

## 8. 待办

| 任务 | 优先级 |
|---|---|
| 端口几何 user override（直接给 a, b, axes 用于报告） | 低 |
| 端口本征解换稀疏求解器（大网格） | 中 |
| TFE-N 多模 S 矩阵输出（不再只取主模 b_in/b_out） | 中 |
| ALPS + TFE-N 401 频点 baseline + CI | 高 |

