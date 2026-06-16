# 端口模式任务

## 目标

改进数值端口本征模计算、归一化或投影逻辑。

## 优先阅读文件

- `../primitives/physical.md`
- `../primitives/fem.md`
- `../primitives/invariants.md`
- `../../src/fem/PortModeSolver.cpp`
- `../../src/fem/FEMAssembler.cpp`
- `../../src/post/ResultExtractor.cpp`
- `../validation.md`

## 可能修改文件

- `../../include/bpfem/fem/PortModeSolver.hpp`
- `../../src/fem/PortModeSolver.cpp`
- `../../src/fem/FEMAssembler.cpp`
- `../../src/post/ResultExtractor.cpp`

## 验证命令

```powershell
.\build\Release\bp_fem_solver.exe --out result --max-sweep-points 1 --write-all-fields
```

## 风险提示

- 装配和提取必须使用同一套模态约定。
- 端口三角单元边引用必须与体单元一致：`localNodes` 按全局点号小到大存储。
- faceId 必须匹配。
- L2 模态归一化由 `PortModeSolver` 完成（`∫|e_t|² dS = 1`）。
- 功率归一化由坡印廷面积分 `poyntingPowerIntegral(mode, f)` 给出，磁场由 `H_t = (β/ωμ₀)(n̂×E_t)` 重构。
- `FEMAssembler` 入射幅度和 `ResultExtractor` 投影分母必须用同一个 `powerNormalizationFactor`。
- 两端口截面不同时 `s_p` 不同；S21 模值不会再被等截面巧合掩盖。
- 修改求积点缓存（`PortMode::quadrature`）需要与 stiffness/mass 装配的求积规则保持一致。
- 一阶 hierarchical 路径下，`evaluateBasis` 中的两个面 DOF 必须线性独立（`FaceFirst0 = lambda_c * N_ab`、`FaceFirst1 = lambda_a * N_bc`），且端口三角形求积必须 ≥ degree 4。`PortModeSolver::evaluateBasis` 与 `FEMAssembler::evaluateBasis` 必须用同一份公式。
- 传播常数和截止条件依赖单位体系。
