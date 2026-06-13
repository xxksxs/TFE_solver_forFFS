# 优化路线专题索引

本目录把"从原型走向商业级 FEM 电磁求解器"的优化任务按主题拆分。顶层入口仍为 `../optimization-roadmap.md`。

## 阅读顺序

1. **总览与差距分析**：`overview.md`
2. **数值方法升级**：`numerical-methods.md`
3. **边界条件与激励**：`boundary-and-excitations.md`
4. **材料模型扩展**：`materials.md`
5. **线性求解器与扫频**：`solvers-and-sweeps.md`
6. **并行化与性能**：`parallelism-and-performance.md`
7. **几何与网格**：`geometry-and-mesh.md`
8. **后处理与结果**：`post-processing.md`
9. **多物理与瞬态**：`multiphysics-and-transient.md`
10. **可用性 / API / UI**：`usability.md`
11. **验证、测试、CI**：`verification-and-ci.md`
12. **分期路线图**：`roadmap-phases.md`

## 文档列表

- **`overview.md`**：当前能力、差距分析与演进哲学。
- **`numerical-methods.md`**：高阶 / hp 自适应 / 误差估计 / 数值积分。
- **`boundary-and-excitations.md`**：PMC / ABC / PML / 周期 / 集总端口 / 平面波 / 电压源。
- **`materials.md`**：各向异性、色散、铁氧体、薄层、非线性、温度相关。
- **`solvers-and-sweeps.md`**：AMS/HX/AMG/DDM/FETI、AWE/Krylov MOR、本征/特征模。
- **`awe-family/README.md`**：AWE、MGAWE、WCAWE 的算法 skill、模块边界与统一校验方法。
- **`parallelism-and-performance.md`**：OpenMP / TBB / MPI / GPU / I/O / cache。
- **`geometry-and-mesh.md`**：CAD 导入、网格生成器、曲面/曲边、非协调网格。
- **`post-processing.md`**：Touchstone、远场、近场、Q、群时延、RCS、表面电流。
- **`multiphysics-and-transient.md`**：TD-FEM / DGTD、热电、电路联仿。
- **`usability.md`**：Python API、GUI、工程文件、调度、许可。
- **`verification-and-ci.md`**：基准库、单元/系统测试、覆盖、回归 CI。
- **`roadmap-phases.md`**：分期目标、里程碑、交付与依赖图。

## 架构演进与已落地的优化

- **`strategy-interfaces/plan.md`**：Phase 1 接口化主计划（IBoundaryCondition /
  ISparseSolver / IPreconditioner / ISweepStrategy + 工厂分发）。**已落地**。
- **`strategy-interfaces/api.md`**：四个抽象接口的契约说明（参数、调用顺序、错误模型）。
- **`strategy-interfaces/example.md`**：加新策略的完整示例（以 ImpedanceBC 为例）。
- **`alps-sweep/plan.md`**：单点 ALPS Krylov MOR。**已落地**。
- **`awe-family/README.md`**：AWE / MGAWE / WCAWE 快速扫频算法族。**设计中**。
- **`analytic-port/plan.md`**：解析端口模 (APM)。**已落地**。
- **`transfinite-element/tfe-plan.md`**：超限元数值多模端口 (TFE)。**已落地**。

## 与现有文档的关系

- 任何修改都遵循 `../primitives/invariants.md` 中的不变量链。
- 验证用 `../validation/` 中的最小验证流程；新增基准追加到 `../validation/` 与 `verification-and-ci.md`。
- 新增模块或扩展点须同步更新 `../overview/extension-points.md`。
- 求解器/CMake 修改须同步 `../architecture/solver-selection.md`。
