# 材料模型扩展

本文给出从"标量、各向同性、频率独立"到商业级材料库的演进路线。

## 1. 现状

`Material` 仅含 `relativePermittivity / relativePermeability / conductivity`，三者皆为标量、实数、频率独立。`FEMAssembler` 通过 `materialForBody(int)` 拿到当前体单元的 ε_r、μ_r、σ。损耗仅靠 σ 体导率项给出。

这种模型只能正确描述**无色散均质介质 + 弱有耗导体内部**，对铁氧体、各向异性陶瓷、薄层金属、Drude/Lorentz 金属、温度敏感材料等都不适用。

## 2. 目标材料族

| 类别 | 典型应用 | 必要参数 |
|------|---------|---------|
| 各向异性张量 | 各向异性介质陶瓷、PML、晶体 | `tensor3<complex>` ε_r、μ_r |
| Debye 单 / 多极 | 生物组织、水 | (ε_∞, ε_s, τ) × N |
| Drude | 金属（光学频段） | (ω_p, γ) |
| Lorentz | 介质共振 | (ω_0, ω_p, γ) × N |
| 旋磁 / Tensor μ | 铁氧体（环行器、隔离器） | (μ_s, ω_0, α, H_dc, M_s) |
| 表面阻抗 / SIBC | 良导体表面 | Z_s(f) 或 σ + 厚度 |
| 薄层 (Two-sided sheet) | 金属化、薄电介质层、石墨烯 | 双侧 Y/Z 矩阵 |
| 温度相关 | 多物理 / 热电 | ε_r(T), σ(T), tan δ(T) |
| 非线性 | 强场 | ε(E), 三阶 χ(3) |
| 复合 / 频率插值 | 实验数据 | (f, ε', ε'', μ', μ'') 表 |

## 3. 数据契约（`include/bpfem/core/Types.hpp`）

```cpp
enum class MaterialModel {
    Constant, Anisotropic,
    DebyePoles, DrudePoles, LorentzPoles,
    Gyromagnetic, ThinSheet, TabulatedFreq, UserCallback
};

struct PoleTerm {
    std::complex<double> a0, a1, a2;     // numerator
    std::complex<double> b0, b1, b2;     // denominator
};

struct Material {
    std::string name;
    MaterialModel model = MaterialModel::Constant;

    // Constant / Anisotropic
    std::array<std::complex<double>, 9> epsilonTensor{};   // 3x3 row-major
    std::array<std::complex<double>, 9> muTensor{};
    double conductivity = 0.0;

    // Pole-based dispersion
    double epsInf = 1.0;
    double muInf = 1.0;
    std::vector<PoleTerm> epsPoles;
    std::vector<PoleTerm> muPoles;

    // Gyromagnetic
    Vec3 dcMagneticField{};
    double saturationMagnetization = 0.0;
    double gyroDamping = 0.0;

    // Thin sheet
    double thickness = 0.0;
    std::array<std::complex<double>, 4> sheetAdmittance{}; // 2x2 Y matrix

    // Tabulated
    std::vector<std::tuple<double, std::complex<double>, std::complex<double>>> table;
};
```

## 4. 装配层影响

- `FEMAssembler` 中所有 `eps * I`、`mu * I` 替换为 `evaluateMaterialTensor(material, frequencyHz)`。
- 各向异性 ε / μ 张量会破坏复对称性 → `MklPardisoSolver` 必须切换到 *complex non-symmetric* 路径（`mtype = 13`）。这是从 Phase 1 末或 Phase 2 起必须支持的双路径。
- 体导率与色散一起出现时，要按 `evaluateComplexEpsilon(material, ω) = ε(ω) − jσ/ω` 合并，避免双重损耗。
- 旋磁材料的 μ 张量会带来 *非对称* 矩阵，需要单独的本征求解器配合圆极化分量。

## 5. 色散扫频策略

- **逐频点装配**：每个频率 `f` 重新装配（已实现），扩展为按色散模型生成 ε(f) / μ(f)。
- **基函数与 ε(f) 解耦**：稀疏图 (`SparsePattern`) 不变，可继续复用 PARDISO symbolic。
- **MOR 友好型展开**：把 ε(ω) 写成关于 ω 的有理函数后，把矩阵展开成 `K + jωC + (jω)² M + Σ_p (1/(jω − s_p)) R_p`，便于后续 Krylov MOR（见 `solvers-and-sweeps.md`）。

## 6. AEDT 解析层影响

`AEDTParser::parseMaterials` 当前只取标量。需要扩展：

- 支持 AEDT 中 *FrequencyDependent dataset*。
- 支持 *Anisotropic* 张量 (`relative_permittivity[1,1] ... [3,3]`)。
- 支持 *Magnetic Saturation / DC bias / Damping*（铁氧体）。
- 支持 *Surface roughness*（Hammerstad / Groisse / Huray 模型）。

## 7. 验收基准

- **Debye 水**：1 GHz ~ 100 GHz ε(f) 与 Cole-Cole 模型重合。
- **Drude 金（光频）**：表面等离激元波数误差 < 1%。
- **铁氧体 Y 结环行器**：3 端口 S 参数与文献基准在 0.3 dB 内。
- **各向异性 PML**：内层反射 < −60 dB。
- **石墨烯薄层**：表面等离激元色散与解析关系 < 1% 误差。

## 8. 任务清单

| 任务 | 优先级 | 估算 |
|------|--------|------|
| 张量 ε / μ + 复非对称 PARDISO 路径 | 高 | M |
| Debye / Drude / Lorentz 多极 | 高 | M |
| AEDT 频率相关 / 张量 / 铁氧体解析 | 高 | M |
| 旋磁 μ + 铁氧体 | 中 | L |
| 薄层 / 两侧 Y/Z 表面 | 中 | M |
| 表面粗糙度 (Huray) | 中 | S |
| 温度相关 / 非线性接口 | 低 | L |
| 频率表插值材料 | 低 | S |
