# 后处理与结果

本文给出从"S11 / S21 + VTU"演进到商业级后处理体系的完整方案。

## 1. 现状

- `ResultExtractor::extract` 输出 `SParameterPoint{ frequencyHz, s11, s21 }`。
- `OutputWriter::writeSParameters` 输出含 `S11_dB / S11_real / S11_imag / S21_*` 的 CSV。
- `OutputWriter::writeVTU` 按 *VTK Lagrange Tetrahedron*（cell type 71）输出 `E_real / E_imag / H_real / H_imag` 四个矢量场。可选 `--field-output-order 1|2|3` 控制每个四面体的 Lagrange 节点数（4 / 10 / 20）。一阶基函数下默认 p_out=2，可看见 EdgeFirst 边中点；显式 p_out=3 可看见 FaceFirst 面贡献。模值 / 切向分量 / dB 等派生量交给 ParaView Calculator 现场计算，不再写入 VTU。

距离商业级：缺 Touchstone、Y/Z 参数、群时延、远场、近场、增益、轴比、RCS、Q 值、表面电流、端口去嵌、端口重归一化、多端口/多模 S 矩阵。

## 2. S 参数与电路兼容输出

| 输出 | 优先级 | 实现要点 |
|------|--------|---------|
| **Touchstone v1 / v2 (`.s2p`, `.s4p`, `.snp`)** | 高 | 每端口/每模一行；规范 `# Hz S MA R Z₀`；多模时使用 v2 `[Reference]` |
| **Y / Z 参数** | 高 | 由 S 矩阵和 Z₀ 直接换算（`Y = (I−S)(I+S)⁻¹ / Z₀`） |
| **VSWR / Return Loss / Insertion Loss** | 高 | 由 S 计算 |
| **群时延 τ_g(f) = −dφ(S21)/dω** | 高 | 在频率插值后做有限差分；与 AWE/MOR 联动得到解析 dφ/dω |
| **多端口 / 多模 S 矩阵** | 高 | `SParameterPoint` 升级为 `SParameterMatrix(numPorts, numModes)` |
| **端口去嵌 (De-embedding)** | 中 | 支持指定参考面平移距离 d，按 `S' = exp(jβd) S exp(jβd)` |
| **端口重归一化** | 中 | 支持任意 Z_ref 列表 |
| **失配损耗 / Mismatch loss** | 低 | 衍生量 |

## 3. 远场 / 增益 / RCS

| 输出 | 实现要点 |
|------|---------|
| **远场方向图 E_θ(θ,φ), E_φ(θ,φ)** | 在外围辐射面（PML 内表面或 ABC 上）取等效电/磁流，做 Stratton-Chu 远场积分 |
| **增益 / 方向性 / 效率** | `D = 4π U(θ,φ) / P_rad`；`G = 4π U / P_in` |
| **极化分量 / 轴比 (AR)** | 圆极化天线分析 |
| **EIRP / SLL / HPBW** | 派生量 |
| **RCS (单/双站)** | 入射平面波 + 远场积分；`σ = lim 4π r² |E_s|² / |E_inc|²` |
| **球面/笛卡尔 / UV 切片** | 输出格式：`farfield.csv`、`farfield.ffe`（FEKO 格式）、`farfield.ffs`（CST 格式） |

模块新增：

```
include/bpfem/post/
  Touchstone.hpp
  NetworkParameters.hpp
  FarField.hpp
  RCS.hpp
src/post/
  ...
```

## 4. 近场 / 体积场 / 其他场量

| 输出 | 实现要点 |
|------|---------|
| **近场 (E, H, S, J_d, J_s)** | 在用户指定的点 / 线 / 面 / 体上采样；`H = (1/jωμ) ∇ × E` |
| **表面电流 J_s = n × H** | PEC / SIBC 表面 |
| **体损耗密度 P_v = ½ σ |E|²** | 用于热-电耦合 RHS |
| **能量密度 / 能流 (Poynting)** | 物理诊断 |
| **场探针 (Field probe)** | 频率点 / 时间点 / 任意位置点列 |

输出格式：

- ParaView：升级 VTU 为二进制 + 压缩；增加 `XDMF + HDF5` 支持长扫频。
- Tecplot：可选 `.plt` ASCII。
- HDF5：统一频率维 + 多场量。

## 5. 谐振 / 滤波器 / 传输线参数

- **Q-factor**：从复本征值 `Q = ω_r / (2 |ω_i|)`；或从 S 参数 3 dB 带宽。
- **谐振频率自动检测**：S21 峰值 + 二次抛物线插值。
- **滤波器参数提取**：耦合矩阵、外部 Q、群时延。
- **传输线参数 (Z, β, γ, α)**：从波端口本征模与多端口 S（双端口 ABCD 反演）得到 `Zc(f)`, `γ(f)`。

## 6. 端口去嵌与参考面控制

商业级求解器允许：

- 指定每个端口参考面相对网格端口面的偏移距离 `d_ref`。
- 在求解后做 `S_meas → S_ref` 的相位补偿。
- 支持参考面跨越色散介质（需要本征模 β(f)）。

## 7. 数据契约

```cpp
struct PortReference {
    int portId = -1;
    double deembedDistance = 0.0;         // m, positive = into device
    double referenceImpedance = 0.0;      // Ω, 0 means use port mode Zc
};

struct SParameterMatrixPoint {
    double frequencyHz = 0.0;
    int numPorts = 0;
    int numModesPerPort = 1;
    std::vector<std::complex<double>> data;   // row-major (numPorts*numModes)^2
};

struct FarFieldPoint {
    double frequencyHz = 0.0;
    double thetaRad = 0.0;
    double phiRad = 0.0;
    std::complex<double> Etheta;
    std::complex<double> Ephi;
    double directivity = 0.0;
    double gain = 0.0;
    double axialRatio = 0.0;
};
```

## 8. CLI / 配置升级

```
--write-touchstone result/device.s2p
--reference-impedance "1=50,2=50"
--deembed "1=2mm,2=2mm"
--farfield "elev:-90:90:181;azim:0:360:361"
--rcs "azim:0:360:361"
--field-probe "0.0,0.5e-3,1.0e-3"
--write-fields-stride 5
```

## 9. 任务清单

| 任务 | 优先级 | 阶段 | 估算 |
|------|--------|------|------|
| Touchstone v1 / v2 输出 | 高 | P1 | S |
| Y / Z / VSWR / 群时延 | 高 | P1 | S |
| 多端口 / 多模 S 矩阵 | 高 | P1 | M |
| 端口去嵌 + 重归一化 | 高 | P2 | M |
| 远场方向图 + 增益 + 轴比 | 高 | P2 | L |
| RCS 单/双站 | 中 | P2 | M |
| 表面电流 + 体损耗 | 中 | P2 | M |
| Q 值 + 谐振检测 | 中 | P2 | S |
| 传输线 Z_c, γ 提取 | 中 | P2 | M |
| HDF5/XDMF 体场 | 中 | P2 | M |
| Tecplot / FFE / FFS 兼容 | 低 | P3 | S |
