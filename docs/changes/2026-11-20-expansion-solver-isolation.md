# 展开点求解会话与迭代残差

承接 [MGAWE 与 WCAWE 审查](2026-11-20-mgawe-wcawe-audit.md)。在同一后端上，一次只允许一个固定矩阵会话。构造失败会释放会话占用。顺序使用后端时，每个展开点重新准备数值因子或预条件器，PARDISO 仍可复用稀疏结构分析。工厂每次创建独立的预条件器。

并发构造时以原子操作争用同一个后端，未获得占用的构造立即抛出 `std::logic_error`。会话内须顺序调用 `solve` 和 `solveBatch`，会话关闭后可以在另一线程重新打开。GMRES 和 BiCGSTAB 后端保留原有的移动构造能力。

BiCGSTAB 对右端项归一化后迭代，还原解向量后用原始矩阵检查相对残差。GMRES 和 BiCGSTAB 都报告非有限的迭代状态。残差和提前停止状态均通过测试复核。

## 验证

- 并发回归让两个线程争用同一 GMRES 后端，并在启用 oneMKL 的构建中重复验证 PARDISO。每次仅一个会话打开，获胜方用自己的矩阵求解，会话关闭后可重新打开。旧 Ninja 构建中 `GmresBackend.cpp.obj` 记录的头文件依赖数为零，验收使用全新构建目录，避免混用类布局不同的旧对象文件。
- MSVC 19.51 Release，Native 和 oneMKL 2025.3＋MUMPS 5.9.2 两个独立构建的 `ctest --output-on-failure` 均通过 3 项。稀疏求解回归覆盖复数不同 RHS、同稀疏结构的新矩阵、工厂预条件器隔离、准备失败后重开、迭代上限、数值中断和非有限状态。
- 修改前，在同一个 PARDISO 后端构造重叠会话使首个会话同一 RHS 的相对残差从 `4.84e-17` 增至 `6.62e-2`。BiCGSTAB 对缩放 `1e-24` 的有效 RHS 在第 1 次迭代返回 `breakdown`，相对残差为 1。修正后拒绝重叠构造，缩放 RHS 回归达到 `1e-11`。
- `example/bp_filter`，零阶棱元、数值主模，40–43 GHz 共 21 点、2 个 oneMKL 线程。WCAWE 41.5 GHz 单点 16 阶并启用 `--wcawe-check-moments`。PARDISO 一次分析、一次数值分解、32 次 RHS 回代，独立矩投影误差 `1.15e-14`。GMRES、BiCGSTAB 使用 `--precon schwarz --schwarz-subdomains 1 --schwarz-overlap 0 --iterative-spmv mkl --max-iterations 25 --tolerance 1e-11`，均达到 16 阶。独立矩投影误差分别为 `1.19e-14`、`1.31e-14`，最大递推残差为 `9.81e-13`、`1.13e-12`。单子域包含完整 PARDISO 分解，用于验证递推和后端接线。
- 固定点 MGAWE 三点各 6 阶，PARDISO 一次稀疏结构分析、三次数值分解、18 次 RHS 回代。相同配置下，单子域 Schwarz 的 GMRES 与 PARDISO 的复数 S11、S21 最大绝对差分别为 `1.03e-12`、`3.48e-12`。
- 将 WCAWE 线性后端设为 GMRES＋Jacobi 并限制为一步，首个 RHS 以 `iteration_limit` 退出，真残差为 `0.825279`，退出码为 1，无 S 参数 CSV。两个 Schwarz 子域、1 层重叠、350 次迭代的配置在当前 BP 滤波器上仍未收敛。MUMPS BLR 在 `1e-11` 容差下拒绝第 2 个 RHS，其真残差为 `1.24e-11`。采用递推验收容差 `1e-10` 时完成 16 阶，独立矩投影误差为 `8.94e-14`。
