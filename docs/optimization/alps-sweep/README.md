# ALPS 快速扫频专题

本目录维护 Adaptive Lanczos-Padé Sweep（ALPS）的理论、实现约定和验证方法。

## 当前实现

`--sweep alps` 使用中心单展开点的标准非对称双边 Lanczos-Padé：

- `PolynomialPortMomentBuilder` 构造 P1 系统 `A(t)=A0+t*A1`、`b(t)=b0+t*b1`。
- 增广状态保留端口激励的一阶导数，不把 `b1` 当成常量丢弃。
- 固定符号约定 `G=-L0^{-1}L1`、`W^T V=I`、`T=W^T G V`，在线系统为 `(I-tT)y=W^T r0`。
- `LanczosPadeModel` 按 `alpha/beta/gamma` 标准三项递推直接生成三对角 `T`；仅在正交性越界时选择性再双正交化。
- 普通耦合击穿但残差未耗尽时尝试一次最小 2x2 look-ahead；失败会终止并写入诊断。
- oneMKL 构建对 `T` 调用 `LAPACKE_zgeev`，在线默认按极点-留数求和；非 MKL 构建及异常情况回退到三对角 Thomas 求解。
- PARDISO 每个展开点只分解一次，Lanczos 非末步将 `Gv` 与 `G^T w` 两个 RHS 合并到一次 `phase=33`。
- S11 和 S21 顺序构造；只有 S11 保留场重构所需右基，S21 只保留降阶数据。

## 运行

```powershell
.\build_pardiso\Release\bp_fem_solver.exe --basis-order 0 --max-sweep-points 101 --sweep alps --alps-order 12 --linear-solver direct --no-write-all-fields --out result
```

`--alps-order q` 直接表示每个标量 `[q-1/q]` Padé 模型的阶数。`--alps-krylov-order q` 是弃用别名，也直接映射到同一个 `q`；两个参数不能同时指定。省略 `--alps-expansion` 时使用频带中心。

## 诊断

`diagnostics.json` 记录矩匹配、双正交误差、三项递推残差、look-ahead、选择性再正交、极点-留数重构误差和疑似伪极点。`timing.json` / `run.json` 区分 RHS 向量数、`phase=33` 调用数和最大批大小，并拆分端口线性化、Lanczos 算子、正交化、极点分解和在线评估时间。

## 边界

当前只支持无损材料、P1 端口线性化、单展开点和两个独立 SISO 模型。本轮不包含 SyMPVL、block Lanczos、P2、有损材料和自动残差驱动多展开点。极点只做诊断标记，不会被静默删除。