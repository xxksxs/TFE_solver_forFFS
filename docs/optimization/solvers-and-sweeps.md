# 线性求解器与扫频

本文给出线性求解器、本征求解器与扫频策略的完整升级方案。

## 1. 现状

- **直接求解**：`MklPardisoSolver` 复对称（`mtype = 6`），symbolic 复用，扫频时只重数值分解。
- **迭代求解**：`BiCGStabSolver` 无预处理；`primitives/linear-algebra.md` 已经指出"高残差时结果不可信"。
- **扫频**：`Application.cpp` 走 `buildFrequencies(...)` 均匀离散；逐频点重装、重分解、重求解。
- **本征**：无谐振 / 端口本征 / 特征模求解。

## 2. 直接求解器升级

| 项 | 现状 | 目标 |
|----|------|------|
| 矩阵类型 | 复对称（mtype 6） | 同时支持复非对称（13）以处理各向异性/旋磁/PML/SIBC |
| 重排序 | 默认 | 可选 METIS / nested dissection / RCM；提供 `--reorder` |
| 多 RHS | 单 RHS | 支持 N 端口 / N 平面波一次求解（PARDISO `nrhs`） |
| Out-of-core | 无 | 大模型时启用 OOC（PARDISO 64-bit + iparm[59]=2） |
| Iterative refinement | 默认 | 控制 `iparm[7]` 与 `iparm[10]` 提升复对称病态系统稳定性 |
| Pivot 控制 | 默认 | 暴露 `iparm[12]` 给 CLI / 配置 |
| 64-bit indexing | 32-bit | 大于 ~2 ⁄ 3 亿非零元时切到 ILP64 |
| 备选直接核 | 无 | 集成 MUMPS、SuperLU_DIST、STRUMPACK 之一作为开源后端 |

## 3. 迭代求解器升级（Phase 1 重点）

`BiCGStabSolver` 在频域 curl-curl 系统上几乎不可用，因为这类矩阵是 *indefinite* 的（`A = K − k₀² M`），普通预处理无法使其变好。商业级求解器普遍使用 *AMS / HX 预处理 + GMRES/FGMRES*：

- **AMS (Auxiliary-space Maxwell Solver, Hiptmair-Xu)**：把 H(curl) 系统通过梯度算子 G、单位插值 P 与节点 H¹ 系统耦合，再分别用 AMG 处理三个 H¹ 子系统。
- **HX 预处理（Hiptmair-Xu 1998 ~ 2007）**：是 AMS 的更通用版本，在 H(curl) 与 H(div) 上都成立。

落地路径：

- 引入 *hypre* 或 *PETSc + ML/AMG/GAMG* 第三方库，复用其 AMS 实现（hypre `HYPRE_AMSCreate`）。
- 自实现轻量版本：依靠 `EdgeTopology` 提供 *离散梯度* `G : H¹ → H(curl)`、*节点位置* 与 AMG 算子；这是中长期工作。
- Krylov 外层使用 GMRES / FGMRES，允许变预处理。

任务列表：

| 任务 | 优先级 | 估算 |
|------|--------|------|
| GMRES + Restart + ILU(0) | 高 | M |
| FGMRES（变预处理） | 高 | M |
| AMS 预处理（hypre 集成） | 高 | L |
| HX 预处理 / 自实现 | 中 | XL |
| 多重网格几何 / 代数 (GMG/AMG) | 中 | L |
| 块 Jacobi + 局部 LU 平滑器 | 中 | M |
| Krylov-Schur / IRA 节点降阶 | 低 | L |

## 4. 区域分解 (DDM)

- **FETI-DP**：把全局问题切成多子域，子域内部 PARDISO，跨界面用 Lagrange 乘子 + 共轭梯度。
- **IETI / 基础 OS-DDM**：商业级求解器（HFSS DDM、CST FD）的主力做法。
- **依赖**：网格分区 (METIS/ParMETIS)、MPI、子域接口构造（界面 H(curl) 模式）。

任务列表：

| 任务 | 优先级 | 估算 |
|------|--------|------|
| 网格分区（METIS） | 中 | M |
| MPI 通信骨架 | 中 | L |
| FETI-DP 内核 | 中 | XL |
| 跨子域端口与 S 参数一致性 | 中 | L |

DDM 实现拨入 Phase 3。

## 5. 本征 / 谐振 / 特征模

商业级求解器至少提供两类本征求解：

- **谐振腔本征 (Eigenmode)**：求 `(K − k² M) e = 0`，用 *Krylov-Schur* 或 *LOBPCG*。
- **端口本征 (Port mode)**：已有 `PortModeSolver`，但仅最低阶；扩展为前 N 个 TE/TEM 模式。
- **特征模分析 (CMA, Characteristic Mode)**：求 `X J = λ R J`，用于天线模式分解；是商业 FEKO/CST 的核心特性。

任务列表：

| 任务 | 优先级 | 估算 |
|------|--------|------|
| 谐振腔本征求解 (Krylov-Schur via ARPACK / SLEPc) | 中 | L |
| 端口高阶本征 + 多模 | 高 | M |
| 特征模 (CMA) | 中 | L |

## 6. 快速扫频 (Frequency Sweep Acceleration)

逐频点重装-重求解是商业级求解器的主要计算成本。优化方法：

- **AWE (Asymptotic Waveform Evaluation)**：在中心频率求一次 + 矩阵泰勒展开，直接 Padé 外推附近频点。
- **多点 AWE / Krylov MOR**：在 N 个频点共建立 Krylov 子空间，构造 *Padé via Lanczos* 或 *PRIMA / SyMPVL*，得到全频段降阶模型。
- **Self-adaptive sweep**：监控 S 参数曲率，按目标精度自动加密频点。
- **离散插值扫频 (Interpolative)**：HFSS 的 *Interpolating Sweep* 思路，对 S(f) 做 Cauchy / Stoer-Bulirsch 插值。
- **Discrete vs Fast (Pole-residue) sweep**：为 MOR 输出 *pole-residue* 表示，可后处理时直接评估任意频点。

**落地状态**：单点 Padé 风格 ALPS（PRIMA/Galerkin 路径）已经合入主干。CLI `--sweep alps` 启用，详见 `alps-sweep/plan.md`；BP filter 基准 101 点扫频加速 ~12×，与 direct 偏差 < 1e-10。

AWE / GAWE / MGAWE / WCAWE 已接入 CLI：

```powershell
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 101 --sweep awe --awe-order 8 --no-write-all-fields --out result
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 101 --sweep gawe --gawe-order 12 --no-write-all-fields --out result
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 101 --sweep mgawe --mgawe-points 3 --mgawe-order 4 --no-write-all-fields --out result
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 1 --max-sweep-points 101 --sweep wcawe --wcawe-order 12 --no-write-all-fields --out result
```

其中 GAWE 是单展开点 Galerkin AWE：先生成局部 AWE 矩向量，再正交化并投影为一个 ROM；MGAWE 使用多个展开点的局部 AWE/GAWE 矩向量，统一正交化为一个 Galerkin ROM。WCAWE 使用单展开点传统 AWE 矩序列生成良条件正交基，并输出 `basis_condition.csv` 记录 AWE/WCAWE 条件曲线。下一步是多点 WCAWE 和自适应展开点。

AWE / GAWE / MGAWE / WCAWE 的模块化算法 skill 与校验方法见 `awe-family/README.md`。这组文档把显式 Padé AWE、单点/多点 Galerkin AWE 和良条件 AWE 拆成可独立实现、可抽离的 fast-sweep 子模块。

当前 fast-sweep 共享内核位于 `include/bpfem/fastsweep/` 与 `src/fastsweep/`：`PolynomialPortMomentBuilder` 统一端口矩递推，`GalerkinReducedModel` 统一 K/M/端口投影与在线 reduced solve，`WellConditionedBasisBuilder` 负责 WCAWE 正交基、上三角系数和 `X≈VR` 诊断。Direct、ALPS、AWE、GAWE、MGAWE、WCAWE 均输出 `diagnostics.json`；WCAWE 额外输出 `basis_condition.csv`。

落地接口：

```cpp
struct SweepStrategy {
    enum Kind { Discrete, Adaptive, AWE, KrylovMOR } kind;
    double tolerance = 1e-3;
    int seedPoints = 3;
    int maxOrder = 30;
};
```

任务列表：

| 任务 | 优先级 | 估算 | 状态 |
|------|--------|------|------|
| 自适应频点加密 | 高 | M | 计划中 |
| AWE 单点 Padé / ALPS Galerkin | 高 | M | **已实现 (`--sweep awe`, `--sweep alps`)** |
| GAWE 单点 Galerkin AWE | 高 | S | **已实现 (`--sweep gawe`)** |
| MGAWE 多点 Galerkin AWE | 高 | M | **已实现 (`--sweep mgawe`)** |
| WCAWE 良条件 AWE | 高 | M | **已实现 (`--sweep wcawe`)** |
| AWE / GAWE / MGAWE / WCAWE 模块化 skill 与校验 | 高 | S | **已规划并落地 (`awe-family/`)** |
| Krylov MOR 多点 | 高 | L | 计划中 |
| Pole-residue 输出 + 任意频评估 | 中 | M |

## 7. 求解器接口契约（Phase 1 已落地）

Phase 1 接口化已合入主干。当前抽象：

```cpp
namespace fem::linalg {
class ISparseSolver {
public:
    virtual SolveResult solve(const SparseMatrix& A,
                              const std::vector<std::complex<double>>& b,
                              const SolverConfig& cfg = {}) = 0;
    virtual void rememberPatternForReuse(bool enable) {}
    virtual const char* name() const = 0;
};

class IPreconditioner {
public:
    virtual void setup(const SparseMatrix& A) = 0;
    virtual void apply(const std::vector<std::complex<double>>& x,
                       std::vector<std::complex<double>>& y) const = 0;
    virtual const char* name() const = 0;
};
}

namespace fem::sweep {
class ISweepStrategy {
public:
    virtual SweepResult run(const std::vector<double>& frequencies,
                            const SweepContext& ctx) = 0;
    virtual const char* name() const = 0;
};
}
```

具体实现 / 路由：

- `linalg::PardisoBackend` / `linalg::BiCGStabBackend` / `linalg::PreconJacobi`
- `sweep::DirectSweep` / `sweep::AlpsSweep`
- 工厂位于 `factory::makeSparseSolver(opts)` 与 `factory::makeSweepStrategy(opts, ...)`
- CLI: `--linear-solver auto|direct|bicgstab`、`--precon none|jacobi`、
  `--sweep direct|alps`

§3 中的 GMRES + AMS 路径在这套接口下要做的事就是新增 `linalg::GmresBackend` +
`linalg::PreconAMS` + 工厂枚举 + CLI 字符串映射；§6 中的多点 ALPS / AAA 插值
扫频是新增 `sweep::AaaSweep` + 工厂枚举。详见
`strategy-interfaces/example.md`。

`MultiRhsSolveResult` 风格的批 RHS 接口尚未引入，留给 §2 直接求解器升级一并
做（PARDISO 的 `nrhs` 字段 + 接口加 `solveBatch(...)`）。

## 8. 验收

- 同一矩阵：PARDISO ↔ GMRES+AMS 解差异 < 1e-6。
- AWE：单点展开 5 阶在频带内误差 < 0.01 dB。
- Krylov MOR：100 频点 sweep 时间 < 10 频点直接 sweep。
- 谐振腔本征：5 个低频本征频率与解析解差异 < 0.05%。
