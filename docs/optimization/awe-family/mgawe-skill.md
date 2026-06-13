# MGAWE Skill：多点 Galerkin AWE

## 适用目标

MGAWE（Multipoint Galerkin Asymptotic Waveform Evaluation）用于宽带 FEM 快速扫频。它把多个展开点
的局部 AWE/GAWE 向量放入同一个正交化流程，形成一个统一的 Galerkin 降阶空间 `V`。

这点很重要：MGAWE 不是多个局部 Padé 模型的拼接，也不是每个展开点各自求一个 ROM 后再按距离切换。
它的核心是“多个展开点同时贡献同一个 reduced model”。

## 输入

- `AffineModelView model`：来自 `FEMAssembler::buildAffineSystem()` 的仿射 FEM 系统。
- `expansionPoints`：多个展开点 `{omega_l}`，初始可取频带左端、中心、右端。
- `orders`：每个展开点的局部阶数 `{q_l}`。
- `localGenerator`：每点的局部矩/残差生成器，可使用严格 AWE 递推或 GAWE/Krylov 向量生成。
- `orthogonalization`：全局 modified Gram-Schmidt + reorthogonalization，含 deflation 阈值。
- `projectionConfig`：Galerkin 投影约定，当前项目默认与 ALPS 一致使用复对称 bilinear `V^T A V`。
- `reliabilityConfig`：direct 抽检频点、残差阈值、最大 ROM 维度。

## 输出

- `globalBasis V`：所有展开点共同生成并正交化后的统一全局基。
- `reducedPolynomialModel`：投影后的多项式/仿射矩阵系数，例如 `K_r`、`M_r`、`m_{p,r}`。
- `reducedRhsAndOutputs`：投影后的 RHS、端口输出矩阵和 S 参数还原所需信息。
- `sParameters`：每个目标频点的 S 参数。
- `orthogonalityMetrics`：`V^T V` 或选定内积下的正交性指标。
- `coverageDiagnostics`：每个展开点的局部阶数、deflation 数量、残差峰值和频带覆盖情况。

## 算法步骤

1. 选择展开点 `{omega_l}` 和每点阶数 `{q_l}`。
2. 对每个展开点组装并因式分解 `A_l = A(omega_l)`。
3. 在每个展开点生成局部向量：
   - 严格 AWE：由 `PolynomialMomentRecurrence` 生成 `{x_j^(l)}`。
   - GAWE/Krylov：使用 `A_l^{-1}` 驱动局部残差/矩方向。
4. 把所有局部向量按展开点和阶数顺序送入同一个全局正交化器。
5. 对线性相关列执行 deflation，并记录来源展开点、阶数和被剔除原因。
6. 用统一基 `V` 一次性投影所有频率无关矩阵和端口向量：
   `K_r=V^T K V`，`M_r=V^T M V`，`m_{p,r}=V^T m_p`。
7. 在线扫频时只求解小系统：
   `A_r(omega) x_r(omega) = b_r(omega)`。
8. 将 `x_r` 投影到端口输出，生成 S 参数；需要场输出时可用 `x_hat=V x_r`。
9. 在探针频点计算全空间代数残差和 direct 抽检误差；若超限，则添加展开点或提高局部阶数。

## 实现建议

- `MgaweSweep` 的核心产物是 `GalerkinReducedModel`，不应输出多个互相独立的 Padé 分式。
- 正交化器必须知道每个候选列来自哪个展开点，方便诊断宽带失效位置。
- 初版可复用现有 ALPS 的 `V^T A V` 约定；后续若改成 Hermitian 投影，必须同步所有 S 参数基准。
- 保留 `V`，因为 MGAWE 能天然支持近似场重构和 VTU 输出。

## 风险点

- 展开点太少会导致频带中段或边缘残差高。
- 展开点太多会导致 ROM 维度过大，抵消在线加速收益。
- 多展开点向量可能高度相关，必须通过 deflation 防止病态投影。
- 若某个端口模式靠近截止，单个全局基可能不足以覆盖分支变化，需要拆分频带。

## 可靠性验证

- **展开点 direct 对比**：每个展开点处 MGAWE 的 S 参数必须接近 direct。
- **全局残差正交性**：验证 `V^T(A(omega)Vx_r-b(omega))` 在 reduced 空间内接近零。
- **宽带误差对比**：同样展开点预算下，MGAWE 应比单点 AWE/GAWE 覆盖更宽频带。
- **deflation 记录**：输出每个展开点、每个阶数的候选列是否被保留，防止静默降阶失败。
- **BP filter 基准**：40--43 GHz 上 direct 抽检 `max delta |S11|, |S21| <= 0.05 dB`。
