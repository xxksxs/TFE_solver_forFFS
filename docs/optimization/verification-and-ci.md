# 验证、测试、CI

本文给出商业级求解器需要的验证体系与 CI 框架。

## 1. 现状

- `docs/validation/` 已经规划了最小验证、残差预期、端口验证、S 参数验证、场可视化、单位检查、失败特征、验证流程。
- `docs/testing/` 包含构建检查、CLI smoke test、输出检查、回归清单与"建议补充的单元测试"。
- 工程内只有一个示例 `wg_bp_filter`；没有自动化 CI、没有规范化基准库、没有收敛性回归。

## 2. 单元测试层

引入 *Catch2* 或 *GoogleTest*（开源、跨平台、与 CMake 友好）：

```
tests/
  unit/
    core_types_test.cpp
    sparse_matrix_test.cpp
    edge_topology_test.cpp
    port_mode_test.cpp
    fem_assembler_test.cpp
    solver_test.cpp
    output_writer_test.cpp
  integration/
    rectangular_waveguide_test.cpp
    rectangular_cavity_test.cpp
    bp_filter_test.cpp
  benchmarks/
    convergence_h.cpp
    convergence_p.cpp
    sweep_acceleration.cpp
```

最低覆盖目标：
- `core / linalg / fem`：行覆盖 ≥ 80%。
- `io / post`：行覆盖 ≥ 60%（受 I/O 限制）。

## 3. 数值基准库（Phase 1 强制）

商业级求解器普遍维护数百个内部基准。本工程的最小可信基准库：

| 基准 | 检验目标 |
|------|---------|
| **矩形腔本征频率（PEC，空气）** | 体装配 + 本征求解 |
| **矩形波导 cutoff (TE10/TE20/TM11)** | 端口模式 |
| **矩形波导主模传输 S21** | 体装配 + 端口装配 + S 参数 |
| **同轴线 TEM (Z₀=50Ω)** | TEM 模式归一化 + 集总端口 |
| **介质腔本征 (ε_r=10)** | 介质装配 |
| **WR-90 带通滤波器（本工程示例）** | 端到端 |
| **微带线传输 S21** | 多介质 + 集总端口（待支持） |
| **金属球 RCS（Mie）** | 平面波激励 + 远场（待支持） |
| **Y 结环行器** | 旋磁材料（待支持） |
| **矩形阵列单元 Floquet** | 周期边界（待支持） |
| **WR-28 偶极子辐射方向图** | 远场 + ABC/PML（待支持） |
| **Drude 金薄膜 SPP 色散** | Drude 色散（待支持） |
| **PML 反射 < −60 dB** | PML（待支持） |

每个基准维护：

- 输入文件 (`benchmarks/<name>/`)。
- 期望输出（解析解 / 文献值）。
- 收敛阶预期（h / p）。
- 容差。

## 4. 收敛性回归

- *h-收敛*：固定 p，加密网格，期望 `S21 误差 ∝ h^{p+1}`。
- *p-收敛*：固定 h，升 p，期望指数收敛。
- 每次合并 PR 时跑前 3 ~ 5 个最快基准；每周跑全集。
- 输出 `convergence.json`，与历史 baseline 对比；偏差超过容差自动失败。

## 5. 求解器路径回归

- 同一矩阵：PARDISO ↔ GMRES+AMS ↔ MUMPS 三路必须给出 < 1e-6 一致解。
- 多 RHS vs 单 RHS 一致性。
- 复对称 vs 复非对称路径在小算例上的等价性（输入是对称的）。
- AWE / Krylov MOR 与离散扫频在 N 个频点上的差异 < 0.01 dB。

## 6. 性能回归

- 维护 `bench/baseline.json`：每个基准的 *DOF 数 / 装配时间 / 求解时间 / 内存峰值*。
- CI 在专用 runner 上每周跑一次；偏离 > 20% 触发 review。
- 集成 *Tracy* trace 的关键阶段计时。

## 7. CI 流水线

建议 GitHub Actions / Azure Pipelines / GitLab CI：

- **构建矩阵**：
  - Windows MSVC + MKL ON
  - Windows MSVC + MKL OFF
  - Linux GCC + MKL ON
  - Linux GCC + MKL OFF
  - Linux Clang + MKL OFF
- **阶段**：
  1. 配置 + 编译（warnings as errors）
  2. 单元测试（Catch2/GoogleTest）
  3. 快速基准（前 3 ~ 5 个）
  4. 静态分析（clang-tidy + cppcheck，工作目录优先级警告）
  5. 内存/UB 检查（每周一次的 ASan / UBSan / MSan）
  6. 性能回归（每周一次）

## 8. 静态质量

- *clang-format*：统一代码风格。
- *clang-tidy*：启用 `modernize-*`、`performance-*`、`bugprone-*`、`cppcoreguidelines-*` 子集。
- *include-what-you-use*：减少头文件依赖。
- *cmake-format*：CMake 样式。
- *pre-commit hooks*：本地阻断风格 / 编译错误。

## 9. 验证文档与基准联动

- 每个基准写一份 `docs/validation/<name>.md`，描述：
  - 物理设置、解析解或文献参考、网格说明、扫频范围、容差、可视化预期。
- `optimization/roadmap-phases.md` 中每项落地与基准列表 1:1 对应。
- 引入 `docs/validation/benchmark-index.md` 作为基准库总索引。

## 10. 任务清单

| 任务 | 优先级 | 阶段 | 估算 |
|------|--------|------|------|
| Catch2 / GoogleTest 集成 | 高 | P1 | S |
| 单元测试覆盖 ≥ 60% | 高 | P1 | M |
| 前 5 个数值基准 | 高 | P1 | L |
| GitHub Actions / GitLab CI | 高 | P1 | M |
| 收敛性回归框架 | 高 | P1 末 | M |
| 求解器三路一致性回归 | 高 | P2 | M |
| ASan / UBSan / 静态分析 | 中 | P2 | S |
| 性能 baseline + 回归 | 中 | P2 | M |
| 全量基准库（≥ 30 个） | 中 | P3 ~ P4 | XL |
