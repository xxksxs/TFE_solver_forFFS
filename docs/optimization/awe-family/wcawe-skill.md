# WCAWE Skill：论文一致的良条件 AWE

## 权威文档

- [中文数学推导 LaTeX](wcawe-theory-cn.tex)
- [中文数学推导 PDF](wcawe-theory-cn.pdf)
- [C++ 代码实现规范](wcawe-implementation-skill.md)
- 论文：R. D. Slone, R. Lee, J.-F. Lee, Well-conditioned asymptotic waveform evaluation for finite elements, IEEE TAP, 2003, DOI 10.1109/TAP.2003.816321。

LaTeX/PDF 是数学定义和代码验收的权威版本；本文件是工程入口。

在本机进入该目录后，使用以下 MiKTeX XeLaTeX 命令运行两次，可重建目录、引用和 PDF：

```powershell
$xelatex = 'C:\Users\30297\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe'
& $xelatex -interaction=nonstopmode -halt-on-error wcawe-theory-cn.tex
& $xelatex -interaction=nonstopmode -halt-on-error wcawe-theory-cn.tex
```

其他机器若已将 XeLaTeX 加入 `PATH`，可直接把 `$xelatex` 替换为 `xelatex`。

## 核心结论

WCAWE 直接处理展开点附近的多项式系统

$$
\left(\sum_{i=0}^{d_A}\sigma^i A_i\right)x(s)
=\sum_{k=0}^{d_b}\sigma^k b_k,\qquad \sigma=s-s_0.
$$

它不对多项式系统做扩维线性化。真正的 WCAWE 也不是先生成完整传统 AWE 矩再做 QR，而是在生成第 n 个候选向量时，使用已经形成的正交基 V 和上三角矩阵 U 构造校正项。

论文的三个核心关系是：

1. 候选矩阵、正交基与上三角矩阵满足 V_n = Vtilde_n U_n^{-1}。
2. 校正项 P_Uw(n,m) 是 U 的连续主子块逆的有序乘积。
3. 第 n 阶候选向量必须按论文式 (7) 读取 U_{n-1} 和 V_{n-1}，随后用 MGS 更新 U 的第 n 列。

当 U=I 时，校正递推退化为传统 AWE；论文实际选择 MGS 系数作为 U。W 表示 well-conditioned，不是频率权重或端口权重。

## 禁止实现

以下路径不能宣称为 Slone 2003 WCAWE：

- 先完整生成传统 AWE 矩，再统一做 QR/MGS。
- 在传统递推中直接使用归一化历史向量，却没有 U 子块逆乘积校正。
- 只用 X≈VR、cond(V)≈1 或基正交性作为论文一致性证明。
- 将 MGS 的 Hermitian 内积与论文复对称 Galerkin 的 transpose 投影混为同一运算。

## 模块输入

- 多项式矩阵系数 A_0...A_dA，按同一 σ 约定构造。
- 激励系数 b_0...b_db，包含频变端口激励导数。
- 展开点 s_0 和目标阶数 q。
- 支持 A_0 一次分解、多 RHS 回代的稀疏求解器。
- MGS、再正交和 breakdown 阈值。
- 明确的复对称 Galerkin 投影约定及端口观测接口。

## 模块输出

- 良条件基 V_q，列长度保持为原 FEM DOF。
- 论文式 (8) 的上三角矩阵 U_q。
- 每阶 P_U1/P_U2 三角回代和校正来源记录。
- 投影后的多项式 reduced model。
- 复数 S11/S21 和按需场重构。
- 式 (7) 递推残差、式 (8) 关系残差、正交性、矩匹配、条件数、时间和峰值内存诊断。

## 推荐模块边界

- PolynomialWcaweRecurrence：按论文式 (7)–(9)逐阶生成候选向量。
- UpperTriangularBlockAction：以三角回代计算校正乘积作用，禁止显式求逆。
- WellConditionedBasisBuilder：执行 MGS、再正交、U 更新和 breakdown。
- GalerkinReducedModel：统一投影、在线求解、S 参数和场重构。
- WcaweDiagnostics：记录论文关系残差和性能数据。

## 当前代码审计结论

当前 `WcaweSweep` 已不再调用 `generateLosslessMoments()` 或保存完整传统 AWE 矩。它通过 `PolynomialWcaweRecurrence` 逐阶组装论文式 (7)，用 `UpperTriangularBlockAction` 对式 (9) 的连续主子块执行逆序三角回代，并在每次 MGS 后立即把系数写入固定步长 `U`，供下一阶候选使用。

当前实现已通过非单位复杂 `U` 的乘积顺序、`U=I` 退化、式 (7)、式 (8)、矩子空间匹配、朴素归一化递推反例和 breakdown 停止测试，可以标记为论文一致的单展开点 P1 WCAWE。当前限制仍是无损材料、P1 端口色散和单展开点；在线阶段使用精确投影后的 K/M/端口频率律，属于论文基上的工程扩展。

30 万四面体网格、90–100 GHz、101 点、q=12、numerical port 的工程结果为：S11/S21 对 HFSS 幅值相对 L2 约 0.774%/0.990%，95 GHz 复数 direct 相对误差约 7.9e-13/9.2e-14，1 次数值分解和 12 次回代。q=20/30/40 均保持递推与正交稳定，但宽带误差没有单调改善，因此该算例尚不足以宣称 WCAWE 相对 GAWE/MGAWE 有精度优势。

## 可靠性验证

- U=I 退化测试：逐列恢复传统 AWE，误差不高于 1e-12。
- Vtilde≈VU：Frobenius 相对误差不高于 1e-12。
- 式 (7) 递推残差不高于 1e-11。
- V^H V≈I：正交误差不高于 1e-11。
- 小型多项式系统至少前 q 个系统矩相对误差不高于 1e-10。
- 必须包含一个能区分论文 WCAWE 与无校正朴素递推的反例。
- 基向量不得扩维；复对称 T 投影与 Hermitian H 投影必须分别验证。
- 工程验收比较 direct FEM 与 HFSS `IOStructure_S_parameters.csv` 的复数 S11/S21，以相对 L2 误差为主，深零点另报绝对幅值误差。
