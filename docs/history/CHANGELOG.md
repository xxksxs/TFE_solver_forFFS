# 变更日志 (Changelog)

本文件记录本工程的版本演进，遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/) 风格。
版本号采用 `YYYY.MM.DD-描述` 形式以方便追溯（项目尚未发布到外部，故未走 SemVer）。

## [Unreleased]

### Added

### Changed

### Fixed

---

## [2026.05.27-baseline] - 2026-05-27

首次进入 Git 管理；该 commit 即为基线 (`baseline-v1` tag)。

### Added

- **`docs/optimization/build_affine_system.tex`** — 独立 LaTeX 文档，从频域
  Maxwell 方程出发推导 `FEMAssembler::buildAffineSystem()` 的数学模型，
  包含 curl-curl 弱形式、单元参考矩阵 $\bm K^T,\bm M^T$、$\bm A(\omega)$ 仿射拆分
  $\bm A=\bm K-k_0^2\bm M+\jj\sum_p\beta_p\bm m_p\bm m_p^{\!\top}$、端口面 2D
  H(curl) 本征推导，以及代码字段对照表 (8 页 PDF)。
- **`PortMethod::Numerical`** — `--port-method numerical` 端口模 CLI 选项，
  默认值由 `Analytic` 改为 `Numerical`。在端口面三角网上求解 2D H(curl) 广义
  本征 $\bm K_{\text{port}}\bm v=k_c^2\bm M_{\text{port}}\bm v$，取最低非伪本征对作为
  端口主模，适用于任意截面。底层复用既有
  `PortModeSolver::computeMultiMode(faceId, 1)`。
- **`.gitattributes` / `.gitignore`** — Windows / 跨平台行尾约定 + 构建目录、
  扫频结果、IDE 元数据屏蔽规则。
- **`results_npm_alps_101/`**（不进 Git，但在 `compare/` 子目录留下了对比报告）—
  NPM + 1 阶基函数 + ALPS 101 频点扫频与 HFSS TFE 参考曲线对比：全频带
  $\lvert S_{11}\rvert$ max 偏差 0.022 dB，$\lvert S_{21}\rvert$ max 0.004 dB，
  3 dB 中心频率 41.0650 GHz 完全对齐。

### Changed

- `Application.hpp::PortMethod` 增加 `Numerical` 项，默认值改为 `Numerical`。
  `Analytic` 保留为闭式 TE10 参考路径（仅矩形截面），`Transfinite` 多模本征不变。
- `Application.cpp` CLI 解析接受 `numerical|analytic|tfe`，并新增数值本征
  分发分支与对应日志格式。
- `RunReport.cpp` 三选项 `portMethod` 报告。
- `docs/primitives/physical.md`、`docs/testing/cli-smoke-tests.md`、`README.md`
  端口建模章节同步描述：默认走 NPM，APM 降为闭式参考。
- `docs/optimization/build_affine_system.tex` 端口模来源章节改为说明 NPM/APM/TFE
  三条路径，并把端面 2D 本征公式 $\bm K_{\text{port}}\bm v=k_c^2\bm M_{\text{port}}\bm v$
  正式写入。

### Validated

| 配置 | 频点数 | 频域求解时间 | 与 HFSS TFE 偏差 |
| --- | --- | --- | --- |
| NPM + basis-order 1 + ALPS, krylov 30, expansion 41.5 GHz | 101 | 82.8 s | $\lvert S_{11}\rvert$ max 0.022 dB / $\lvert S_{21}\rvert$ max 0.004 dB |

3 dB 通带中心频率 41.0650 GHz，与 HFSS 完全对齐 (shift 0.0000 GHz)。

---

## 维护约定

- 每次有用户可见变化（CLI、API 行为、数学模型修订、性能数字、生成产物布局）
  都至少在 `Unreleased` 下记一条；当一组改动构成里程碑时，发版并打 tag。
- 性能 / 精度数字走 `### Validated` 子节，附上完整命令行与对比指标，便于回归。
- 仅文档措辞修订或注释微调可只写 commit message，不必动 CHANGELOG。
