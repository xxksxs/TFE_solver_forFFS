# ALPS 标准双边 Lanczos-Padé 实现说明

本文档记录 `src/sweep/AlpsSweep.cpp` 与 `src/fastsweep/LanczosPadeModel.cpp` 的当前工程约定。理论推导以 `theory-cn.tex` 为准。

## 1. 代数约定

P1 增广系统写为

```text
L0 = [ A0  -b0 ],  L1 = [ A1  -b1 ],  z0 = [ A0^-1 b0 ]
     [  0    1 ]        [  0    0 ]       [      1      ]
G  = -L0^-1 L1
```

因此

```text
h(t) = c^T (I-tG)^-1 z0,
T    = W^T G V,
W^T V = I,
h_q(t) = c^T V (I-tT)^-1 W^T z0.
```

`A1` 同时包含体质量项和端口导纳一阶导数，`b1` 包含激励一阶导数。左右过程使用非共轭转置双线性积；不得把公式改成共轭内积。

## 2. 标准三项递推

对第 `j` 对基向量执行

```text
alpha_j = w_j^T G v_j
r_j = G v_j - alpha_j v_j - gamma_(j-1) v_(j-1)
s_j = G^T w_j - alpha_j w_j - beta_(j-1) w_(j-1)
beta_j gamma_j = s_j^T r_j
```

代码选择 `beta_j=sqrt(s_j^T r_j)`、`gamma_j=(s_j^T r_j)/beta_j`，由此保持 `w_(j+1)^T v_(j+1)=1`。正交性监测超过阈值才做两遍选择性再双正交化。耦合接近零但左右残差不为零时尝试最小 2x2 look-ahead；无法恢复时明确设置 breakdown，不补入数值无效列。

## 3. 阶数与求解次数

`--alps-order q` 是单个 SISO `[q-1/q]` 模型的直接阶数，不再乘虚拟端口数。S11、S21 两个模型顺序构造。固定矩阵只做一次数值分解；每个模型需要 `2q-1` 个算子 RHS，加上展开点场解，总 RHS 向量数为 `1+2(2q-1)=4q-1`。当 `q=6` 时为 23 个 RHS。批量接口将每个非末步的左右 RHS 合成一次 PARDISO `phase=33`。

## 4. 极点-留数

oneMKL 构建对三对角 `T` 求右特征分解 `T=Z Lambda Z^-1`，得到

```text
h_q(t) = sum_i weight_i / (1 - t lambda_i).
```

构造阶段在固定探针上与三对角求解交叉校验；误差超过阈值时禁用极点求和并回退。带内小留数极点仅计入 `spurious_pole_count`，不筛除、不移动。

## 5. 性能与内存

- `solveFactorizedBatch` 用一次 `phase=33` 和 `nrhs>1` 回代多列 RHS。
- 统计分别记录 RHS 向量数与实际回代调用数。
- S21 不保留全空间左右基；S11 只保留右基用于 `field_last.vtu`。
- 在线极点求和复杂度为 `O(q)`，异常时 Thomas 求解同样为 `O(q)`。
- 子阶段计时包括端口线性化、稀疏分解、Lanczos 算子、正交化、极点分解和在线扫频。

## 6. 验证门槛

单元测试检查：

- `W^T V`、`T=W^T G V` 和三项递推残差不高于 `1e-11`。
- 前 `2q` 阶标量矩匹配不高于 `1e-11`。
- P1 的 `A1*x` 与 `b1` 中心有限差分误差不高于 `1e-7`。
- happy breakdown、普通失败和成功 look-ahead 均有确定状态。
- 极点-留数与三对角探针差异不高于 `1e-12`。
- 批量和逐 RHS PARDISO 解的差异不高于 `1e-12`。

大网格基准继续使用零阶 `current.ngmesh`、90--100 GHz、101 点和 `IOStructure_S_parameters.csv`。性能结论必须使用 Release PARDISO、预热后三次中位数；在完成实测前不把目标时间写成已达到事实。

## 7. 当前 30 万网格实测

零阶 `current.ngmesh`（364932 DOF）、90--100 GHz、101 点、PARDISO Release 的阶数筛选结果：

| q | HFSS 相对 L2 S11 | HFSS 相对 L2 S21 |
|---:|---:|---:|
| 6 | 72.79% | 13.53% |
| 8 | 33.23% | 5.54% |
| 10 | 6.46% | 0.483% |
| 12 | 2.91% | 0.622% |

因此默认阶数取 `q=12`，不能用 `q=6` 的时间替代全带精度。最终 `q=12` 与优化前同阶结果的最大复数差异为 S11 `9.93e-13`、S21 `2.80e-15`；1 次数值分解、47 个 RHS 向量、25 次回代调用，峰值约 3101 MB。三次热运行总时间为 26.98/27.14/27.15 s，中位数 27.14 s，相比旧基线 34.01 s 提升约 20.2%。

当前已通过“至少快 20%”、数值一致性、HFSS 精度和内存门槛；`总时间 <=25.5 s` 与 `在线 <=0.005 s` 两项尚未达到（在线中位数约 0.0092 s），后续优化不能以降低阶数破坏精度。

## 8. 后续范围

可靠的全阶残差驱动多展开点、block MIMO、P2 端口校正、有损/色散材料和更完整的 look-ahead block recurrence 保留为后续阶段。本轮不引入 SyMPVL。
