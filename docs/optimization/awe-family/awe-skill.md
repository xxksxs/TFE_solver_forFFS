# AWE Skill：单点矩递推 + Padé 极点留数

## 适用目标

AWE（Asymptotic Waveform Evaluation）用于从一个展开点附近的少量 FEM 矩信息构造低阶有理响应。
论文中的典型形式是把 FEM 系统写成关于波数 `k` 的多项式矩阵方程：

`A(k) x(k) = b(k)`，其中 `A(k) = A0 + k A1 + k^2 A2 + ...`。

在本项目中，等价入口来自 `FEMAssembler::buildAffineSystem()`：

`A(omega) = K - k0^2 M + sum_p Y_p(omega) m_p m_p^T`。

AWE 适合单个展开点附近的窄带扫频、输入阻抗、S 参数、远场或单点场量的快速评估。

## 输入

- `AffineModelView model`：FEM 多项式/仿射矩阵系数，至少包含 `K`、`M`、端口耦合向量和 RHS 规则。
- `expansionPoint`：展开波数 `k0` 或展开频率 `omega0`，需要记录端口 `beta(omega0)` 分支。
- `order`：系统矩阶数与 Padé 阶数，初版默认 `4 <= order <= 12`。
- `rhsSeries`：激励 RHS 的局部级数系数，例如 `{b0, b1, ...}`。
- `observables`：需要近似的输出量，例如端口投影、输入阻抗、`S11/S21` 或某个场自由度。
- `factorizer`：展开点矩阵 `A(expansionPoint)` 的直接因式分解器，优先 PARDISO。
- `padeConfig`：`[L/M]` 阶数选择、Hankel 条件数阈值、极点过滤阈值。

## 输出

- `systemMoments`：解向量矩 `{x0, x1, ...}`，或只保存可重放的递推摘要。
- `responseMoments`：每个观测量的标量/端口传递函数矩 `{m0, m1, ...}`。
- `PadeApproximant`：Padé 分子、分母、极点、留数和拟合条件数。
- `sParameters`：目标频点上的 S 参数。
- `diagnostics`：矩递推误差、Padé 矩匹配误差、Hankel 条件数、spurious pole、direct 抽检误差。

## 算法步骤

1. 选择展开变量。论文中常用波数 `k`；本项目可先使用 `lambda=k0^2` 简化体项，再把端口
   `beta(omega)` 作为可导标量项处理。
2. 在展开点组装并因式分解：
   `A_e = A(expansionPoint)`。
3. 构造局部级数：
   `A(t)=A_e + t A_1 + t^2 A_2 + ...`，
   `b(t)=b_0 + t b_1 + ...`，其中 `t = k-k_e` 或 `omega-omega_e`。
4. 递推系统矩：
   `x_0 = A_e^{-1} b_0`，
   `x_j = A_e^{-1}(b_j - sum_{i=1..j} A_i x_{j-i})`。
   对论文中的二次 FEM 形式，`A1/A2` 项必须同时进入递推。
5. 对观测量做投影，得到响应矩：
   `m_j = observable_j({x_0...x_j}, {C_0...C_j})`。
6. 用响应矩构造 `[L/M]` Padé，使它匹配前 `L+M+1` 阶 Taylor 系数。
7. 对 Padé 做极点-留数分解，用极点位置判断单点有效带宽。
8. 在目标频点评估 Padé 并转成项目 CSV 中的 S 参数格式。
9. 若 Hankel 条件数过高、极点落入目标频带或 direct 抽检误差超限，则缩小频带、降低阶数，
   或切换到 GAWE/MGAWE/WCAWE。

## 实现建议

- `PolynomialMomentRecurrence` 只负责生成系统矩，不直接构造 Padé。
- `PadeApproximant` 只接收响应矩，输出有理函数和极点留数。
- `AweSweep` 不保存巨大矩向量时，应保存展开点、阶数、RHS、观测量和可重放的诊断。
- 初版只支持 lossless，与当前 ALPS 限制一致；有损材料需要先把导电项纳入仿射/多项式系数。

## 风险点

- 传统 AWE 直接使用高阶矩，数值上容易病态。
- Padé Hankel 系统在阶数升高后条件数会快速变差。
- 单点 AWE 的有效带宽由最近极点/主导留数控制，不能默认覆盖整个 40--43 GHz。
- 端口 `beta(omega)=sqrt(k0^2-kc^2)` 靠近截止时不可简单 Taylor 外推。

## 可靠性验证

- **矩递推有限差分**：小矩阵上比较 `x_j` 与 `d^j x/dt^j / j!` 的中心差分结果。
- **Padé 矩匹配**：验证 `[L/M]` 近似与输入响应矩匹配到 `L+M` 阶，误差小于 `1e-10`。
- **极点位置检查**：Padé 极点不得落在目标采样频带内；若靠近采样点，必须标记结果不可靠。
- **单点有效带宽检查**：从展开点向两侧扩大频带，误差应单调或分段可解释地增长。
- **direct FEM 抽检**：在展开点、左右邻近点和频带边缘运行 direct，对 `|S11|/|S21|` 设
  `max delta <= 0.01 dB` 的窄带门槛。
