# 端口与装配烟测

使用已构建的 `bp_fem_solver` 和同一次构建的 `io_parser_tests`，无需新增 CMake 目标。脚本用 Node.js 18 或更新版本运行，从仓库根目录执行。

~~~powershell
node scripts/port_assembly_smoke.mjs --solver <build>/bp_fem_solver.exe --probe <build>/io_parser_tests.exe --benchmark ../../benchmark --out <new-output-directory>
~~~

默认使用 PARDISO 和两个线程。可用 `--linear-solver` 和 `--threads` 更改。多配置构建请传入对应配置目录中的程序。输出路径须为新目录，包含每组 CLI 日志、求解器结果以及 `port-assembly-smoke.json`。报告记录频点、端口和模式编号、截止、约束数、矩阵与 RHS 相对差、可观测功率和真实残差。工程和网格输入记录 SHA-256。

也可以从已配置好 MSVC 和 oneMKL 的 PowerShell 通过[本地回归入口](build-checks.md)运行。环境里有多个 Node.js 时，适配脚本使用 PATH 中排在首位的 Node.js。

~~~powershell
.\scripts\run-local-regression.ps1 -Profile mkl -Configuration Release -SmokeScript .\scripts\run-port-assembly-smoke.ps1
~~~

CTest 通过后，适配脚本选用同一构建目录中的求解器和探针，在 `build/local-regression/mkl-Release/smoke/port-assembly-<id>/` 生成独立报告。设置 `-BuildRoot` 可指定构建和烟测结果根目录。该入口要求可用 PARDISO，可使用 `mkl` 或 `mkl-mumps`。

## 算例与验收

| 算例 | 输入与频点 | 验收 |
|---|---|---|
| 均匀矩形波导 | `tests/fixtures/rectangular-waveguide` 中的配套 AEDT 和 NGMesh。截面 20 × 12 mm，长度 20 mm，真空、外壁 PEC，9–10 GHz 等距 5 点 | 数值单模和旧双模各运行一次。比较 TE10/TE01 的解析截止 $c_0/(2a)$、$c_0/(2b)$，离散相对容差分别为 1% 和 6%。主模反射幅值 ≤ 0.05，传播相位误差模 $\pi$ ≤ 0.10 rad。两个独立求得的实模态向量允许符号相反。 |
| Pillbox 陶瓷窗 | 算例库 `pillbox_window/project/Project4.aedt` + `current.ngmesh`。脚本仅在 AEDT 副本中替换扫频字段，取 5.6–6.0 GHz 共 9 点 | 单模和旧双模下每点 $|S_{21}|\ge0.01$。既有 `io_parser_tests` 额外检验陶瓷赋材、完整电导率字段、内部介质界面切向迹和显式 PEC。 |

每段以第一个端口的主模入射，对同频完整 Direct 的 $A(f)x$、$b(f)$ 与局部 Taylor 的 $A_0x,b_0$ 做对照，同时独立调用实际 `WavePortBC::apply` 核对端口块，避免体刚度掩盖端口差异。使用一般复向量 $x$，误差容差 $10^{-11}$。检查虚拟模数量、耦合向量的 PEC 过滤、被动模无入射，以及 PEC 零阶单位行和一阶齐次行。CLI 要求真实相对线性残差 ≤ $10^{-7}$。

功率和定义为 $\sum_{(p,m)\in\mathcal P(f)}|S_{(p,m),\,\mathrm{in}}|^2$，$\mathcal P$ 为当前频点保留且传播的输出通道。旧 CSV 只输出前两个端口的主模，脚本要求这两个主模恰好覆盖全部保留的传播通道，否则功率检查失败。本例中该和为 $|S_{11}|^2+|S_{21}|^2$。单模允许偏差 $10^{-6}$，旧双模允许 $5\times10^{-4}$，记录当前兼容模型在截止以下使用 $jk_0$ 项时的偏差。

## 适用范围

烟测覆盖固定真空截面的双端口、零阶体装配和旧 `--port-method tfe`。内部介质与 PEC 小测试另覆盖基阶 0 和 1。本入口检查同网格装配与传输，未做 HFSS 曲线对齐、精确截止、跨截止、原生 `tfe-block` 或多传播模输出验收。
