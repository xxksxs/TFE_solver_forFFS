# Step 3：`IPreconditioner` 抽象 + Jacobi 实现

> Phase 1 / Step 3 of `docs/optimization/strategy-interfaces/plan.md`。可与 Step 2 并行。仅依赖 Step 1。

## 目标

`BiCGStabBackend`（Step 1 引入）接受可选预条件，第一个具体实现 = Jacobi 对角预条件。这一步是**抽象接口的最便宜实战检验**：如果 Jacobi 在 H(curl) 系统上没有正向收益，说明该问题需要更结构化的预条件（AMS / SAI / Schwarz），相关结论会指导后续实现选择。

## 优先阅读文件

- `../optimization/strategy-interfaces/plan.md`（必读）
- `../../src/linalg/BiCGStabSolver.cpp`（重点：迭代主循环里 `r̃ = r` 那段，需要插 `M^-1` apply）
- `../../include/bpfem/linalg/BiCGStabSolver.hpp`
- Step 1 的成果：`include/bpfem/linalg/ISparseSolver.hpp` / `BiCGStabBackend.hpp`

## 可能修改文件

```
include/bpfem/linalg/IPreconditioner.hpp      新增
include/bpfem/linalg/PreconJacobi.hpp         新增
src/linalg/PreconJacobi.cpp                   新增
src/linalg/BiCGStabSolver.cpp                 加 setPreconditioner / 在迭代循环 apply M^-1
include/bpfem/linalg/BiCGStabSolver.hpp       加 setPreconditioner 公共接口
src/linalg/BiCGStabBackend.cpp                构造时根据配置插入预条件
src/factory/SparseSolverFactory.cpp           解析 --precon none|jacobi
include/bpfem/app/Application.hpp             Options 加 PreconditionerKind
src/app/Application.cpp                       parseOptions 加 --precon
```

## 接口签名

```cpp
namespace fem::linalg {

class IPreconditioner {
public:
    virtual ~IPreconditioner() = default;

    // 装配时调用一次，让预条件提取它需要的信息（比如对角线）。
    virtual void setup(const SparseMatrix& A) = 0;

    // 每次迭代调用：y = M^{-1} x。in-place 允许 (&x == &y)。
    virtual void apply(const std::vector<std::complex<double>>& x,
                       std::vector<std::complex<double>>& y) const = 0;
};

}
```

## 验证命令

```powershell
cmake --build build_mkl --config Release

# 单频点 BP filter，BiCGSTAB 无预条件 vs Jacobi 预条件
.\build_mkl\Release\bp_fem_solver.exe --max-sweep-points 1 --linear-solver bicgstab --tolerance 1e-9 --max-iterations 2000 --no-write-all-fields --out results_step3_nopre
.\build_mkl\Release\bp_fem_solver.exe --max-sweep-points 1 --linear-solver bicgstab --precon jacobi --tolerance 1e-9 --max-iterations 2000 --no-write-all-fields --out results_step3_jacobi
```

预期：

- 两次 S 参数 \|Δ\| < 1e-7（同一物理问题、不同求解器精度内的等价）
- Jacobi 收敛迭代次数应 < 无预条件版本（任意正向减少都算成功）
- 若 Jacobi *没有*减少迭代（H(curl) curl-curl 主导时这是常见结论），把这个事实记到 `docs/optimization/numerical-methods.md` 并保留 `--precon jacobi` CLI（占位 + 启发后续 AMS/SAI 实现）

## 不变量 / 风险提示

- **预条件不是必须 SPD**：BiCGSTAB 接受任意非奇异 M^{-1}。Jacobi 在复频域 curl-curl 系统上 diag(A) 可能有近零项（curl-curl 的体单元局部 mass 主导项 `-k0² M_e` 在低 k0 时变小），需要保护：`diag_i^{-1} = 1 / diag_i if |diag_i| > 1e-30 else 0`。
- **数值精度对照不要太严**：BiCGSTAB 是不稳定迭代法，无预条件、有预条件得到的 S 参数末位有可能在 1e-8 量级抖动，验证脚本用 < 1e-7 的 tol。
- **Pardiso 路径不影响**：直接求解器忽略 `--precon`，工厂里要写一行警告"--precon ignored when --linear-solver direct"。

## 估时

0.5 day。

## 输出物

- 1 个 PR
- BP filter 上 BiCGSTAB ± Jacobi 的迭代次数 / 残差对比附在 PR 描述
- 不论收益正负，结论写入 `docs/optimization/numerical-methods.md`
