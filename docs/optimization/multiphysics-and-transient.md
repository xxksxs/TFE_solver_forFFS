# 多物理与瞬态

本文给出从纯频域 FEM 走向多物理 / 时域的扩展计划。

## 1. 现状

当前求解器是**纯频域**，单物理（电磁），不支持瞬态、热、结构、流体或电路联仿。

## 2. 时域求解 (TD-FEM / DGTD)

目的：宽带激励、瞬态响应、非线性材料、瞬时器件。

| 方案 | 特点 | 落地难度 |
|------|------|----------|
| **TD-FEM (Newmark / Bossak)** | 隐式时积分，在每步求解类似频域系统；可复用 PARDISO/AMS | 中 |
| **DGTD (Discontinuous Galerkin Time Domain)** | 显式 Runge-Kutta + 高阶向量元；高度并行；可与 GPU 配合 | 高 |
| **基于 FEM 频域 + IFFT** | 宽带扫频后做 FFT 反演；不能处理非线性 | 低（已具雏形） |

建议路线：

- **Phase 3**：实现 *FFT 反演* 工具（基于已有 frequency sweep）。
- **Phase 3 末 / Phase 4**：实现 *DGTD* 原型（基函数与端口一致）；最先支持 PEC + 平面波 + ABC/PML，再扩展端口。
- **Phase 4**：考虑 TD-FEM 隐式 + 复用直接/AMS 预处理。

## 3. 多物理耦合

商业级电磁求解器（ANSYS Workbench / CST Studio Suite / COMSOL）必备：

| 耦合 | 用途 | 实现要点 |
|------|------|---------|
| **电-热（EM-Thermal）** | 介质损耗 / 导体焦耳热 / 微波加热 / 高功率器件 | EM 体损耗 `P_v = ½ σ \|E\|² + ½ ω ε'' \|E\|²` 作为热源；热场反馈到 ε(T)/σ(T)；松耦合或全耦合 Newton |
| **电-力（EM-Mechanical）** | 加工、谐振位移 | Maxwell 应力张量 → 结构 FEM |
| **电-电路（EM-Circuit）** | 集总元件嵌入 / Touchstone 联仿 | 时域：MNA + 状态空间；频域：S 矩阵嵌入电路网表 |
| **电-流体 / 等离子体** | 高功率微波 / EBG / RF MEMS | 长期目标 |

推荐先做 **EM-Thermal**：

- 体损耗输出：在 `post-processing.md` 已规划。
- 热 FEM 模块：可独立实现节点 H¹ 求解器（标量 PDE，远比 H(curl) 简单）。
- 数据交换：`fem/coupling/EmToThermal.{hpp,cpp}` 把单元损耗映射到热网格（同网格或非协调插值）。
- 控制流：`Application` 增加迭代外层（Picard / Aitken）。

## 4. 电路联仿

- **频域 SPICE 联仿**：把 EM 的 N 端口 S 矩阵输出为 Touchstone（`post-processing.md`），交给 SPICE/Cadence/ADS。
- **时域电路嵌入**：电压/电流源 + 集总 RLC + 受控源；与 TD-FEM 集成。
- **状态空间模型**：从 Krylov MOR 输出降阶 (A,B,C,D)；可被电路求解器读取。

## 5. 优化与参数化

商业级必备的*设计支持*能力：

| 功能 | 用途 | 实现 |
|------|------|------|
| 参数扫描 (DOE) | 拓扑/尺寸研究 | `--param "L1=0.5:0.1:1.5"`，逐次重网格 + 求解 |
| 全局优化 (GA/PSO/CMA-ES) | 滤波器/天线综合 | 调用 `nlopt` / `pagmo` |
| 梯度优化（伴随） | 高维设计变量 | 求解伴随场 → 形状/材料梯度 |
| 拓扑优化 | 任意形状 | Density-based + 滤波 |
| 灵敏度分析 | 公差研究 | 伴随或有限差分 |

伴随求解需要：

- 求解 `A^H λ = ∂Q/∂E*`（Q 为目标量），与正向求解共享分解。
- 梯度积分：`∂Q/∂p = Re[ λ^T (∂A/∂p) E ]`。
- 几何梯度需要 *velocity field* 形状导数。

## 6. 数据契约扩展

```cpp
struct ThermalProperties {
    double thermalConductivity = 0.0;   // W/m/K
    double specificHeat = 0.0;           // J/kg/K
    double density = 0.0;                // kg/m^3
    double convectionH = 0.0;            // W/m^2/K
    double ambientTemperature = 293.15;
};

struct CircuitPort {
    int portId = -1;
    std::string spiceModelPath;
    double referenceImpedance = 50.0;
};

struct OptimizationVariable {
    std::string name;
    double initial = 0.0, lower = 0.0, upper = 1.0;
};
```

## 7. 任务清单

| 任务 | 优先级 | 阶段 | 估算 |
|------|--------|------|------|
| 体损耗 / 热源输出 | 高 | P2 | S |
| 节点 H¹ 热 FEM 求解器 | 中 | P3 | L |
| EM-Thermal 松耦合 | 中 | P3 | M |
| Touchstone-SPICE 联仿 | 高 | P2 | S |
| Krylov MOR → 状态空间 | 中 | P3 | M |
| 参数扫描 + DOE | 高 | P2 | S |
| 全局优化 (nlopt / pagmo) | 中 | P3 | M |
| 伴随梯度 | 中 | P3 | L |
| FFT 反演时域 | 中 | P3 | S |
| DGTD 原型 | 中 | P4 | XL |
| TD-FEM 隐式 | 低 | P4 | XL |
| 拓扑优化 | 低 | P4 | XL |
