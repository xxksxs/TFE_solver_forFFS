# 快速扫频模块开发计划：Krylov MOR

本文给出"基于 Krylov 子空间矩匹配的快速扫频模块"的实现计划。理论模型见 `theory.tex`。本计划在用户确认前**不会触及代码**。

## 1. 目标

- **核心目标**：把宽带 N_f 频点扫频从"逐频点重装-重分解"改为"建一个降阶模型 (ROM)，在线只算密小矩阵"，目标加速比 ≥ 30×（在 401 点 / 290k DOF / 一阶 hierarchical / PARDISO 路径下）。
- **精度**：与直接扫频相比，所有频点的 |S_pq| 误差 ≤ 0.05 dB；通带内 ≤ 0.01 dB。
- **健壮**：必须自动检测 ROM 不收敛并回退到逐点求解（fail-open，不丢精度）。

## 2. 非目标

- 不替换现有 `MklPardisoSolver` / `BiCGStabSolver`。MOR 是"扫频策略"，正交于"线性求解器"。
- 不引入 hp 自适应、AMS 预处理、DDM 等其他优化路线（见 `optimization/solvers-and-sweeps.md`），它们各走各的里程碑。
- 不处理强非线性介质和带宽内多个尖锐谐振穿插的极端工况（见 `theory.tex` §11）。

## 3. 范围与边界

| 维度 | 在范围内 | 不在范围内（本期） |
|------|---------|------------------|
| 基函数 | 0 阶 + 1 阶 hierarchical | p ≥ 2、hp 自适应 |
| 边界 | PEC + 数值波端口 | PML / 集总端口 / 周期 / 平面波 |
| 材料 | 标量 ε / μ + σ + Debye/Drude/Lorentz（仿射展开） | 张量 / 旋磁 / 强色散经验数据 |
| 扫频 | 单频带 [ω_a, ω_b]，离散 + 自适应加点 | 跨带、多带拼接 |
| 求解器后端 | MKL PARDISO（直接）作为 offline 阶段的 LU 提供者 | 迭代 / GPU / DDM |
| 输出 | S 参数（CSV + Touchstone 后续）；任意频点全场 VTU（按需） | 远场 / RCS（独立路线） |

## 4. 顶层架构

新增模块 `bpfem::mor`，在 `src/mor/` 与 `include/bpfem/mor/`，依赖方向遵循 `architecture/dependency-direction.md`：

```
app
 |-- io
 |-- fem
 |-- linalg
 |-- mor          <-- 新增
 |   |-- core
 |   |-- linalg
 |   |-- fem      （只读：FEMAssembler / PortModeSolver / EdgeTopology）
 |-- post
```

`mor` 只读 `fem` 模块，不反向耦合。这与 `architecture/dependency-direction.md` 中"app 负责组装、底层不调用上层"的原则保持一致。

### 4.1 关键类草案

```
namespace fem::mor {

// 仿射展开：A(omega) = sum_i alpha_i(omega) * A_i
struct AffineOperator {
    SparsePattern pattern;                 // 共享，与体装配相同
    std::vector<SparseMatrix> terms;       // {K, M, m_p m_p^T, ...}
    std::vector<std::function<std::complex<double>(double)>> coefficients;
};

// 输入/输出投影
struct PortIO {
    std::vector<std::vector<std::complex<double>>> rhsVectors;   // 列：每个端口的 m_p
    std::vector<std::vector<std::complex<double>>> outputVectors;// 行：每个端口的 m_p
    std::vector<std::function<double(double)>> inputScales;      // 2*j*beta_p*s_p(omega)
    std::vector<std::function<double(double)>> outputScales;     // 1/s_p(omega)
};

// 抽象的"全空间求解器适配器"，让 MOR 不感知 PARDISO/迭代细节
class FullSpaceSolver {
public:
    virtual void factorize(const SparseMatrix& A) = 0;
    virtual std::vector<std::complex<double>> solve(
        const std::vector<std::complex<double>>& rhs) = 0;
};

// 多点 Krylov ROM
class KrylovMorSweep {
public:
    KrylovMorSweep(const FEMAssembler& assembler,
                   const PortModeSolver& portModeSolver,
                   FullSpaceSolver& linear,
                   MorOptions options);

    void buildOffline(const SweepRequest& request);
    SweepResult evaluateOnline(const std::vector<double>& frequencies);

private:
    void appendKrylovBlock(double omega0, int q);
    void orthonormalize();
    void assembleReducedOperators();
    double residualIndicator(double omega) const;
};

} // namespace fem::mor
```

### 4.2 关键数据结构

| 名称 | 含义 | 维度 |
|------|------|------|
| `V` | ROM 子空间正交基 | n × q |
| `Atilde[i]` | 各项 V*A_i V | q × q（密） |
| `Btilde` | V*B | q × N_p |
| `Ctilde` | C V | N_p × q |
| `expansionPoints` | 展开点 ω_0^(ℓ) | L 个标量 |
| `validationPoints` | 验证频点 | N_v 个标量 |

## 5. 算法落地

### 5.1 离线阶段

1. 调用 `FEMAssembler` 装配仿射 term：`K`、`M`、每个端口的 `m_p`。要求 `FEMAssembler` 暴露这些 term，**而不是只暴露最终复合矩阵**（见 §8 接口契约扩展）。
2. 选择初始展开点 Σ₀（band 端点 + 中点）。
3. 对每个 ω₀^(ℓ)：
   - 装配 A(ω₀^(ℓ))（在 PARDISO 阶段做一次复合，CSR 图复用）。
   - PARDISO 数值分解（symbolic 复用）。
   - Block-Arnoldi q_ℓ 步，得到 `V` 列。
4. 全局 MGS 重正交化 + 去 deflation。
5. 装配 `Atilde[i] = V* A_i V` 与 `Btilde = V* B`。

### 5.2 在线阶段

对每个目标频率 ω：
1. 计算每个 `alpha_i(ω)`。
2. 累加 `Atilde(ω) = sum_i alpha_i(ω) * Atilde[i]`。
3. LDLᵀ 分解 q×q 复对称矩阵，求 `xtilde = Atilde(ω)^{-1} Btilde`。
4. `S(ω)` 从 `Ctilde * xtilde` 加上 `s_p(ω)` 与入射幅度（见 `theory.tex` 式 (3)）。

### 5.3 自适应展开

- 在 [ω_a, ω_b] 上密集探针集 Ω_p（如 200 点）评估 `rho_alg(ω)`（一次全空间 SpMV，不需要分解）。
- 取最大处的 ω* 加入展开点集，重复，直到 `rho_alg ≤ ε_alg` 或达到 q_max。
- 用户可见 tolerance 仍是 `rho_S`（在 N_v 个验证点上跑一次直接求解器）。

### 5.4 收敛失败回退

- 若 `q == q_max` 且 `rho_S > tol_S`：日志 `[warn] MOR did not converge`，自动回退到 `--sweep direct`。
- `SweepResult` 记录每个频点是否走 ROM、对应残差，便于下游分析。

## 6. 模块布局与文件清单

```
include/bpfem/mor/
  AffineOperator.hpp
  PortIO.hpp
  FullSpaceSolver.hpp           // 抽象适配器
  PardisoFullSpaceSolver.hpp    // 包装 MklPardisoSolver
  KrylovMorSweep.hpp
  MorOptions.hpp                // CLI/YAML 入口
  MorTypes.hpp                  // SweepRequest / SweepResult / SweepStrategy

src/mor/
  AffineOperator.cpp
  PortIO.cpp
  PardisoFullSpaceSolver.cpp
  KrylovMorSweep.cpp
  RationalKrylov.cpp            // 多点扩展
  Reorthogonalize.cpp           // MGS + 二次重正交
  ReducedSolve.cpp              // 复对称 LDL^T (LAPACKE_zhetrf 或 zsysv)

tests/unit/
  affine_operator_test.cpp
  krylov_mor_test.cpp
  reduced_ldl_test.cpp

tests/integration/
  bp_filter_mor_sweep_test.cpp
```

## 7. CLI 与配置

新增 CLI 选项（保持 `Options` 结构与现有约定一致）：

```
--sweep direct|adaptive|mor                     # 默认 direct（兼容现状）
--mor-tolerance 1e-3                            # rho_S 上限（dB-equivalent）
--mor-expansion "auto"                          # 或 "0,40e9,41.5e9,43e9"
--mor-q-max 80                                  # 单点 Krylov 最大长度
--mor-validation-points 5                       # N_v
--mor-residual-tolerance 1e-4                   # rho_alg 自适应阈值
--mor-write-rom result/rom.bin                 # 可选：导出 ROM 复用
```

YAML（待 P1 末交付，与 `optimization/usability.md` 计划一致）：

```yaml
sweep:
  strategy: mor
  mor:
    tolerance: 1.0e-3
    expansion: auto
    q_max: 80
    validation_points: 5
    residual_tolerance: 1.0e-4
```

## 8. 对现有代码的接口扩展（只列变更，不改实现）

落地实现期需要的最小契约改动，列在这里供评审：

1. **`FEMAssembler`**：新增 `assembleAffine(const SparsePattern&) -> AffineOperator`，按需返回 `K`、`M`、每个端口 `m_p m_p^T`，复用现有 `sparsityPattern_`。
   - 不改 `assemble(double, ...)`，旧路径完全保留。
2. **`PortModeSolver`**：新增 `couplingVector(int faceId) const -> std::vector<std::complex<double>>` 返回稀疏向量 `m_p`（已经在 `couplingWeights` 中以 `(globalIdx, value)` 形式存在；这里给 MOR 一个紧凑稠密视图）。
   - `s_p(ω) = powerNormalizationFactor` 已经存在，复用即可。
3. **`MklPardisoSolver`**：新增 `factorize(const SparseMatrix&)` 与 `solve(const std::vector<std::complex<double>>&)` 两个接口（或包一层 `PardisoFullSpaceSolver`）。当前 `solve(matrix, rhs)` 把两步合一，对 MOR 不友好。
4. **`Application`**：根据 `--sweep` 选用 `DirectFrequencySweep` 或 `KrylovMorSweep`。两者都返回 `std::vector<SParameterPoint>`，下游 `OutputWriter::writeSParameters` 不变。
5. **`ResultExtractor`**：保持不变。MOR 只是给它一个 `bm x(ω)`，公式（见 `theory.tex` §2.3）完全一致。

任一项契约改动会同步更新 `primitives/invariants.md`、`architecture/risk-points.md`，并添加单元测试覆盖旧路径不退化。

## 9. 验证、测试、CI

| 测试 | 类型 | 验收 |
|------|------|------|
| 空矩形波导 S21 | 单元 | 与解析 e^{-jβL} 误差 ≤ 1e-6 |
| BP filter（本工程） | 集成 | 401 点 vs direct，最大 |S| 偏差 ≤ 0.05 dB |
| 收敛阶 q-curve | 回归 | log10(error) 随 q 单调递减直到饱和 |
| 能量守恒 | 回归 | \|S11\|² + \|S21\|² ≤ 1.001 |
| 回退路径 | 单元 | 故意把 q_max=2，应自动 fallback 到 direct |
| 对称性保持 | 单元 | 装配后 \|Atilde - Atilde^T\|_F / \|Atilde\|_F ≤ 1e-10 |
| Cutoff 频段 | 单元 | β_p 实数 → 虚数过渡，s_p 连续到 0 |

CI：
- 在 `BPFEM_USE_MKL=ON` 与 `OFF` 两个矩阵下都构建。
- MOR 路径在 OFF 时仍要可编译，但 offline 阶段会因为缺 PARDISO 而回退到 BiCGStab + ILU(0)（fallback 的具体策略放到 `theory.tex` §6 的扩展或独立 follow-up）。

## 10. 风险与缓解

| 风险 | 影响 | 缓解 |
|------|------|------|
| ROM 在传输零点附近精度不足 | 通带边界 |S| 偏差大 | 自适应加点 + 验证点强制覆盖零点频段 |
| 端口截止下 β_p → 0 数值奇异 | s_p / β_p 不连续 | `theory.tex` §6 的 sqrt branch；与现有 `powerNormalizationFactor` 行为一致 |
| 一阶 hierarchical mass 病态 | Krylov 很慢 | 偏好 PARDISO LDL^T；监控 \|Atilde\|/\|Atilde^{-1}\| |
| 复对称 vs 复非对称 | 损耗/各向异性时矩阵非对称 | 检测对称性，自动切到非对称 LDU；与未来 P2 材料扩展对接 |
| 内存 V 占用 | 290k × 80 × 16B = 350 MB 量级 | 可接受；记录在 profile |
| 与色散材料 (Phase 2) 的耦合 | AffineOperator 需要更多 term | `terms` 是 vector，扩展不破坏接口 |
| 未来 hp 自适应改基 | V 失效 | 接口预留 `invalidateBasis()`；本期不实现 |

## 11. 分阶段交付（M1 ~ M4）

每阶段必须有单元/集成测试 + 文档更新 + CI 通过。

### M1（约 2 周）：单点 Krylov + 直接评估
- AffineOperator + PardisoFullSpaceSolver + 单点 block Arnoldi + 直接评估 ROM。
- BP filter 中心点单点展开，对中心 ±100 MHz 范围 21 频点：误差 ≤ 0.1 dB。

### M2（约 2 周）：多点 Krylov + 自适应
- 多展开点 + adaptive 加点 + 验证点回退。
- BP filter 401 点：误差 ≤ 0.05 dB；加速 ≥ 20×。

### M3（约 1 周）：CLI / 输出 / 健壮性
- `--sweep mor` 默认体验；JSON 日志含 ROM 维度、展开点列表、各频点残差。
- ROM 序列化（可选）+ 任意频点 VTU 重构。

### M4（约 1 周）：验证 / 文档 / CI
- 全部测试入 CI。
- 更新 `optimization/solvers-and-sweeps.md`、`primitives/linear-algebra.md`、`testing/regression-checklist.md`、`validation/sweep-mor.md`（新增）。

## 12. 文档同步清单（落地后必须更新）

- `optimization/solvers-and-sweeps.md`：把 *AWE* / *Krylov MOR* 段从"计划"改为"已实现"。
- `optimization/roadmap-phases.md`：标记 Phase 2 中相关条目完成。
- `primitives/linear-algebra.md`：新增 *降阶模型* 一节，说明 `Atilde` 与全空间矩阵的语义差。
- `architecture/dependency-direction.md`：补 `mor` 模块依赖图。
- `architecture/solver-selection.md`：明确 `mor` 与 `direct` / `iterative` 的关系。
- `validation/sweep-mor.md`（新建）：MOR 验证流程、容差、典型曲线、失败签名。
- `testing/cli-smoke-tests.md`：补 `--sweep mor` 烟测。
- `testing/regression-checklist.md`：MOR 回归条目。
- `task-splitting/`：可能新增 `mor-tasks.md`，把 `KrylovMorSweep` 升级（多点、p 自适应耦合）作为后续任务派发。

## 13. 评审请求

请重点确认以下项，决定后我再进入实现：

1. **理论范围**：`theory.tex` 中"仿射展开 + 单点/多点 Krylov + 残差指示器"是否覆盖你心中的目标？是否需要加 PRIMA 的 *passivity-by-construction* 证明、或 *AAA / balanced truncation* 备选？
2. **验收基准**：401 点 / 1e-3 容差 / 30× 加速这三条是否合理？是否需要加 *群时延* 或 *谐振 Q* 验证？
3. **接口契约扩展**（§8 的 5 项）：是否同意新增 `assembleAffine`、把 `MklPardisoSolver` 拆 `factorize/solve` 两步？这两项会成为本次开发的最低成本最深耦合点。
4. **CLI 命名**：`--sweep direct|adaptive|mor` 与 `--mor-*` 前缀是否符合现有命名风格（当前都是 `--max-*` / `--no-*`），有没有更偏好的命名？
5. **回退策略**：MOR 失败回退到 `direct` 是否够？还是要严格 fail-closed（直接报错让用户决定）？
6. **顺序**：M1～M4 分期是否可以接受？或者要求"先把多点 + 自适应一次性做完"？

确认后我再启动代码改动。
