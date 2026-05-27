# 加新策略：完整可跑示例

本文用一个**虚构但可编译**的"实现 ImpedanceBC"路径，展示 Phase 1 接口化
之后加新边界条件需要触碰的所有代码。同样的步骤套用于加新求解器后端、
新预条件、新扫频策略。

> 范围：仅展示工程接入步骤；surface-impedance 数值实现 (`(jωμ₀/Z_s)
> (n × (n × E)) = ...` 在端口面外的金属壁面三角形上的装配) 是物理工作，
> 与本文无关。

## 起点

仓库已经有 `include/bpfem/bc/ImpedanceBC.hpp` 和 `src/bc/ImpedanceBC.cpp`
作为占位骨架（Step 2 引入）。`apply` / `declareSparsity` 当前抛
"not implemented"。

任务：把它改成能在指定面 ID 上加 `Z_s` 边界，并通过新增 CLI 选项
`--impedance-faces 12,15,17 --impedance-zs 100+50j` 启用。

## Step A：填充 ImpedanceBC 实现

`src/bc/ImpedanceBC.cpp`：

```cpp
#include "bpfem/bc/ImpedanceBC.hpp"
#include "bpfem/core/Constants.hpp"

#include <cmath>

namespace fem::bc {

namespace {

// 收集落在 impedance 面上的三角形的边-边耦合 (row, col) 对。
// 不同 face 上的三角形互斥；triangulations 内 PEC 边由 ctx.constrainedDofs
// 自动过滤掉。
struct ImpedanceTriangle {
    int v0;
    int v1;
    int v2;
    int faceId;
};

std::vector<ImpedanceTriangle> collectFaces(
    const AssemblyContext& ctx,
    const std::unordered_set<int>& impedanceFaceIds) {
    std::vector<ImpedanceTriangle> out;
    for (const auto& tri : ctx.mesh.surfaceTriangles) {
        const auto fit = ctx.mesh.facets.find(tri.facetId);
        if (fit == ctx.mesh.facets.end()) continue;
        const int faceId = fit->second.faceId;
        if (impedanceFaceIds.count(faceId) == 0) continue;
        out.push_back({tri.vertexIds[0], tri.vertexIds[1], tri.vertexIds[2], faceId});
    }
    return out;
}

}  // namespace

void ImpedanceBC::apply(SparseMatrix& matrix,
                        std::vector<std::complex<double>>& rhs,
                        const AssemblyContext& ctx) const {
    (void)rhs;  // Z_s 边界仅贡献 LHS。
    const auto triangles = collectFaces(ctx, faceIds_);
    if (triangles.empty()) return;

    const double omega = 2.0 * pi * ctx.frequencyHz;
    const std::complex<double> coeff = std::complex<double>(0.0, omega * mu0) / surfaceImpedance_;

    for (const auto& tri : triangles) {
        // ... 装配 surface mass 矩阵 (n × N_i) · (n × N_j) ...
        // 对每对边-边自由度 (row, col)，constrainedDofs 已知则 skip，否则
        // matrix.add(row, col, coeff * surfaceMass(i, j))。
    }
}

void ImpedanceBC::declareSparsity(SparsePatternBuilder& builder,
                                  const AssemblyContext& ctx) const {
    const auto triangles = collectFaces(ctx, faceIds_);
    for (const auto& tri : triangles) {
        // 同样遍历每对边-边自由度，对每个 (row, col) 调一次 builder.add(row, col)。
    }
}

}  // namespace fem::bc
```

`include/bpfem/bc/ImpedanceBC.hpp` 改一个字段，让构造函数接受面 ID 集合：

```cpp
class ImpedanceBC : public IBoundaryCondition {
public:
    ImpedanceBC(std::unordered_set<int> faceIds, std::complex<double> Zs)
        : faceIds_(std::move(faceIds)), surfaceImpedance_(Zs) {}
    // ... 其它成员同骨架 ...
private:
    std::unordered_set<int> faceIds_;
    std::complex<double> surfaceImpedance_;
};
```

## Step B：CLI / Options

`include/bpfem/app/Application.hpp::Options`：

```cpp
// Step 5 example: surface-impedance BC.
std::vector<int> impedanceFaceIds;     // 默认空 = 不启用
std::complex<double> impedanceZs{0.0, 0.0};
```

`src/app/Application.cpp::parseOptions`：

```cpp
} else if (arg == "--impedance-faces") {
    const std::string list = requireValue(arg);
    // 解析逗号分隔整数列表，填 options.impedanceFaceIds
} else if (arg == "--impedance-zs") {
    const std::string s = requireValue(arg);
    // 解析 "real+imagj" 格式，填 options.impedanceZs
```

`--help` 文本同步加一行。

## Step C：注册 BC

`src/app/Application.cpp::runApplication`，BC 注册段：

```cpp
std::vector<std::shared_ptr<bc::IBoundaryCondition>> bcs;
bcs.push_back(std::make_shared<bc::WavePortBC>());
if (!options.impedanceFaceIds.empty()) {
    std::unordered_set<int> faces(options.impedanceFaceIds.begin(),
                                  options.impedanceFaceIds.end());
    bcs.push_back(std::make_shared<bc::ImpedanceBC>(
        std::move(faces), options.impedanceZs));
    log.info("Impedance BC active on " + std::to_string(options.impedanceFaceIds.size())
             + " faces, Z_s = " + std::to_string(options.impedanceZs.real())
             + "+" + std::to_string(options.impedanceZs.imag()) + "j Ohm");
}
assembler.setBoundaryConditions(std::move(bcs));
```

## Step D：CMake

`src/bc/ImpedanceBC.cpp` 已经在 `CMakeLists.txt` 的 `bp_fem_core` 源文件
列表中。**不需要修改 CMake**。

## Step E：验证

```powershell
# 默认行为：不传 --impedance-faces 时一切如常
.\build\Release\bp_fem_solver.exe --max-sweep-points 3 --no-write-all-fields

# 启用 ImpedanceBC
.\build\Release\bp_fem_solver.exe --max-sweep-points 3 --no-write-all-fields `
    --impedance-faces 12,15,17 --impedance-zs "377.0+0.0j"
```

第一条 CSV 必须与未引入 `ImpedanceBC` 之前 byte-identical（`Compare-Object`
空输出）。第二条会得到不同的 S 参数（边界吸收损耗）。

## 没有触碰的文件

| 文件 | 原因 |
|---|---|
| `src/fem/FEMAssembler.cpp` | BC 列表分发已在 Step 2 完成 |
| `src/sweep/DirectSweep.cpp` / `AlpsSweep.cpp` | 求解循环不知道 BC 类型 |
| `src/factory/SparseSolverFactory.cpp` | 与 BC 无关 |
| `src/factory/SweepStrategyFactory.cpp` | 与 BC 无关 |
| 其它 BC 类 | 互斥实现 |
| `src/main.cpp` | 总是 fem::runApplication 透传 |

## 同样模式适用于

- **新增 GMRES 求解器**：参考 `BiCGStabBackend.cpp` + `enum LinearSolverKind` 加项 +
  `SparseSolverFactory.cpp` switch case + CLI 字符串。
- **新增 SSOR 预条件**：参考 `PreconJacobi.cpp` + `enum PreconditionerKind` +
  `makePreconditioner` switch + CLI。
- **新增 AAA 有理插值扫频**：参考 `DirectSweep.cpp` + `enum SweepStrategy` +
  `SweepStrategyFactory.cpp` switch + CLI。

每条都遵循"接口实现 + 工厂枚举 + CLI 字符串映射"三步法，**不**修改
`Application::runApplication` 主体或 `FEMAssembler::assemble` 主体。
