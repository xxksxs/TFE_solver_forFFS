# 可用性 / API / UI

本文给出从 CLI 原型走向商业级可用性的完整方案。

## 1. 现状

- 仅 CLI：`bp_fem_solver --aedt ... --mesh ... --out ... --basis-order ... --max-sweep-points ... --write-all-fields`。
- 无脚本 API、无 GUI、无工程文件、无作业管理、无许可机制。
- 日志直接打印到 stdout，结构化困难。

## 2. CLI 升级

- 替换手写参数解析为成熟库（CLI11 / cxxopts）。
- 支持 `--config solver.yaml`，把所有参数集中到 YAML / TOML / JSON。
- 支持 `--dry-run`：仅装配 + 报告稀疏图大小。
- 支持 `--check-mesh`、`--check-ports`、`--validate-aedt`。
- 引入 *exit code* 规范：0 成功，1 输入错误，2 数值错误，3 收敛失败。

## 3. 工程文件格式

- 引入项目根文件 `bpfem_project.json`（或 .bpfemproj）：
  - 引用 AEDT / NGMesh / 自有 STEP。
  - 求解配置（基函数阶、扫频策略、求解器选择）。
  - 后处理任务列表（S 参数、Touchstone、远场、近场）。
- 输出根目录结构：

```
results/
  s_parameters.csv
  s_parameters.s2p
  fields/
    field_<freq>.vtu
  farfield/
    farfield_<freq>.csv
  reports/
    solver_log.json
    convergence.json
    mesh_report.json
    profile.trace.json
```

## 4. Python API

商业级求解器的"可编程性"：HFSS 用 IronPython/PyAEDT、CST 用 VBA/Python、FEKO 用 Lua、COMSOL 用 LiveLink。建议直接选 Python 3：

- 使用 *pybind11* 暴露：
  - `bpfem.Project.fromAedt(path)`
  - `bpfem.Mesh.fromNgmesh(path)`
  - `bpfem.Solver(project, mesh, options)`
  - `solver.run() -> Results`
  - `Results.s_matrix(port_pairs)`、`Results.farfield(...)`、`Results.write_touchstone(...)`
- 支持脚本化参数扫描、批量后处理、基准测试。

模块布局：

```
bindings/python/
  bpfem.cpp
  setup.py
  bpfem/__init__.py
  bpfem/post.py
  bpfem/optimize.py
```

## 5. C API / SDK

- 为非 Python 用户（MATLAB / .NET / Julia）提供 `bpfem_capi.h` ABI 稳定接口。
- 适合分发到第三方框架或集成到现有 EDA 流程。

## 6. GUI

商业级求解器的 GUI 需要：

- **几何视图**：CAD + 网格 + 面/边/点拾取。
- **属性面板**：材料 / 边界 / 端口 / 激励 / 扫频。
- **求解监控**：迭代曲线、残差、内存、剖面。
- **结果浏览**：S 参数表 + 频谱图 + 远场极坐标 + 3D 场可视化。

落地路径：

- **Phase 4**：以 *Qt 6 + VTK*（开源）或集成 *ParaView 客户端* 做最小 GUI。
- 中期可以提供 *Web 仪表盘*（FastAPI + React + plotly）作为轻量结果浏览器。

## 7. 调度与作业管理

- 长扫频或参数扫描需要 HPC / 集群支持：
  - 作业脚本（SLURM / PBS / LSF）模板。
  - 检查点 / 续算（每个频点保存中间状态）。
  - 多用户许可 / 并发控制。
- `Application` 加 `--checkpoint dir` / `--resume dir`。

## 8. 日志与诊断

- 引入 `core/Logger` 升级为 *结构化 JSON 日志*（spdlog + JSON sink）。
- 关键事件 schema：
  - `phase`: project_import / fem_assembly / linear_solve / post_processing
  - `frequency`: f_Hz
  - `timing_ms`: ...
  - `metrics`: residual, iterations, nnz, dof
- 输出 `results/reports/solver_log.json` 与人可读 `solver_log.txt`。
- 扫频结束后自动生成 *summary HTML 报告*（包含 S 参数图 + 远场图 + 剖面）。

## 9. 错误信息与诊断

商业级求解器最被吐槽的不是数学，而是错误信息。建议：

- 把"Unknown option"、"Missing value"等错误改为带建议的形式。
- 端口 / 材料 / 扫频等输入校验前置；提前在 `Application` 阶段失败。
- 网格 / AEDT 一致性检查（faceId 匹配、单位一致、端口面非空）独立成步骤。
- 提供 `--diagnose` 一次性输出问题清单。

## 10. 许可

- 商业版本需要许可锁（FlexLM / RLM / 自有签名 token）。
- 引入 `core/License.hpp` 抽象，开源版本可用 always-on stub。

## 11. 任务清单

| 任务 | 优先级 | 阶段 | 估算 |
|------|--------|------|------|
| CLI11 + YAML 配置 | 高 | P1 | S |
| 结构化 JSON 日志 + summary 报告 | 高 | P1 | M |
| 工程文件格式 + 输出布局 | 高 | P2 | M |
| Python API (pybind11) | 高 | P2 | M |
| 检查点 / 续算 | 中 | P2 | M |
| C API SDK | 中 | P3 | M |
| 作业调度模板 | 中 | P3 | S |
| GUI（Qt + VTK） | 低 | P4 | XL |
| 许可机制 | 低 | P4 | M |
| Web 仪表盘 | 低 | P3 ~ P4 | M |
