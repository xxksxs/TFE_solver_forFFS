# ALPS 快速扫频开发计划

本文档记录 ALPS（自适应 Lanczos-Padé 扫频）模块的工程落地。理论参考 `theory-cn.tex`。

## 1. 现状（MVP 已落地）

第一迭代实现已合入 `src/sweep/AlpsSweep.cpp`，可用 CLI `--sweep alps` 启用。覆盖了 `theory-cn.tex` §5（单点 Padé via Lanczos）但当前工程口径明确为 **ALPS-Galerkin MVP**：使用 Krylov-Galerkin/PRIMA 形式，而不是完整 Lanczos-Padé 三对角化。

- **基底构造**：单展开点 + 块 shift-and-invert Krylov，列式增长（每列做 modified Gram-Schmidt + 一次再正交化），避免显式块三对角矩阵。
- **投影**：复对称 Galerkin（左右用同一个 V，bilinear 形式 `V^T A V`）。
- **降阶系统**：每个频点构造 `q×q` 复对称稠密矩阵并用部分主元 LU 解。
- **端口投影**：`b_p = (V^T m_p)^T x_tilde / s_p`，公式与 `ResultExtractor` 完全一致。
- **lossless 限定**：当前要求所有材料 `σ = 0`；检测到任何有损材料会抛错并提示用户改用 `--sweep direct`。

这与 `theory-cn.tex` 中 ALPS（非对称 Lanczos + Padé via Lanczos）严格意义上不同：MVP 走 Krylov-Galerkin 路径（更稳，不需要 look-ahead Lanczos），并保留 ALPS 的核心特性"少分解 + 多评估"。完整的非对称 Lanczos + 极点-留数显式构造留作下一迭代。

### 实测性能（BP filter 基准，48k tet / 一阶 hierarchical / 348 K DOF）

| 阶段 | 时间 |
|------|------|
| 离线（建 60 维 ROM，`q=30/port`） | ~103 s（含 60 次 PARDISO solve） |
| 在线 101 个频点 | < 1 s（约 8 ms/点） |
| **ALPS 全程** | **~104 s** |
| direct 101 点全程估算 | ~22 min |
| **加速比** | **~12×** |

精度：所有频点上 ALPS 与 direct 的 S 参数差异 < 1e-10（远好于 0.05 dB 的工程容差）。

## 2. 模块布局

```
include/bpfem/mor/AlpsSweep.hpp   — 公开 API（AlpsOptions / AlpsSweep）
src/sweep/AlpsSweep.cpp           — 实现：buildOffline / evaluate / reconstructField
src/fem/FEMAssembler.cpp           — 新增 buildAffineSystem() 提供 K, M, m_p, k_c²
```

`mor` 模块依赖 `fem`、`linalg`、`core`，不向上依赖 `app`，符合 `architecture/dependency-direction.md`。

## 3. 已实现 CLI

```
--sweep direct|alps             # 默认 direct
--alps-krylov-order N           # q per port，默认 30；ROM 维 ≈ N × Np
--alps-expansion <Hz>           # 展开点；省略时取扫频带中心
```

## 4. 数据流

```
ProjectDefinition + Mesh + EdgeTopology + PortModeSolver
        │
        ▼
FEMAssembler::ensureElementCache()                                ← 已有
        │
        ▼
FEMAssembler::buildAffineSystem() → (K, M, m_p, k_c², lossless)   ← 新增
        │
        ▼
AlpsSweep::buildOffline(ω₀)
  ├─ assemble(ω₀) + PARDISO LU                                    ← 一次
  ├─ V = orth([A₀⁻¹ m_1, …, A₀⁻¹ m_Np, A₀⁻¹ M v_1, …])           ← 块 Krylov
  ├─ K̃ = V^T K V (q×q complex symmetric)
  ├─ M̃ = V^T M V (q×q complex symmetric)
  └─ m̃_p = V^T m_p (q-vector each)
        │
        ▼  对每个 ω
AlpsSweep::evaluate(ω)
  ├─ Ã(ω) = K̃ - k₀² M̃ + j Σ β_p(ω) m̃_p m̃_p^T
  ├─ b̃(ω) = Σ 2j β_p s_p √P_W m̃_p
  ├─ x̃ = Ã(ω)⁻¹ b̃(ω)        ← 稠密 LU，O(q³)
  └─ S_pq from m̃_p^T x̃ / s_p
```

## 5. 待办（后续迭代）

| 任务 | 优先级 | 估算 | 收益 |
|------|--------|------|------|
| 多展开点 + 自适应加点（`theory-cn.tex` §7） | 高 | M | 宽带扫频精度提升 |
| 残差驱动停止 + ROM 维度自适应 | 高 | S | 节省离线时间 |
| 有损材料支持（每元 `εc(ω)` 仿射展开） | 中 | M | 工程必需 |
| 非对称 Lanczos + 极点-留数显式输出（`theory-cn.tex` §5.3） | 中 | L | 电路联仿、谐振辨识 |
| ALPS ROM 写出 / 加载（避免重算） | 中 | S | 多次评估同一带 |
| 任意频率全场 VTU（`x = V x̃`） | 中 | M | 用户体验 |
| 集成到验证基准（与 direct 对比图） | 高 | S | 防止退化 |
| LAPACK ZSYTRF/ZSYTRS 替换手卷 LU | 低 | S | 性能微调（q 很小，影响有限） |

## 6. 不变量（落地后必须保持）

- **K, M 实对称**：`buildAffineSystem` 产生的 `K, M` 是上三角 CSR 复矩阵，但虚部为 0；如未来引入张量介质须扩展为多分量缓存。
- **复对称投影**：`bilinear(a, b) = Σ a_i b_i`（不带共轭）才使 `V^T K V` 保留复对称性；`hdot` 仅用于正交化稳定。
- **端口归一化**：`evaluate()` 中 β_p / s_p 与 `FEMAssembler::applyWavePorts` / `ResultExtractor::extract` 共用同一份 `powerNormalizationFactor`，已实现。
- **PEC 约束**：在 `buildAffineSystem` 阶段就把端口耦合向量中受约束的 dof 剔除，确保 `m_p` 与 `K, M` 的稀疏图一致。

## 7. 已知限制

- **lossless 限制**：σ ≠ 0 时直接抛错。理论扩展需要把 `M` 拆为 `M_eps + j ω⁻¹ M_sigma`，并用 3-term 仿射展开重做 Galerkin。
- **单点展开**：宽带（> 1 倍频程）建议增加展开点；当前 1 GHz / 41.5 GHz 单点对 3 GHz 带宽足够。
- **ROM 维度固定**：当前用 `--alps-krylov-order` 显式指定，不做自适应；过高会浪费离线时间，过低会丢精度。
- **不写场 VTU**：`--sweep alps` 路径不写 `field_*.vtu`；如需场，混合用 `--sweep direct` 在感兴趣的频点单独跑。

## 8. 验证

- 与 `--sweep direct` 在 5 / 101 频点上的对比已通过；偏差 < 1e-10。
- 能量守恒在通带内 |S11|² + |S21|² ≈ 0.9999；阻带因端口吸收边界条件而 < 1，与 direct 一致。
- 后续应加入 `validation/sweep-alps.md` 描述容差与失败签名。
