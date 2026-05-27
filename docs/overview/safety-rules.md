# 安全规则

本文列出修改工程前应遵守的基础规则。

## 保持数据契约

共享结构变化应先反映到 `../../include/bpfem/core/Types.hpp`，再更新解析器、FEM、求解器和后处理消费者。

## 不要孤立修改端口

端口模式、边界装配和 S 参数投影必须保持一致。涉及端口时至少检查：

- `../../src/fem/PortModeSolver.cpp`
- `../../src/fem/FEMAssembler.cpp`
- `../../src/post/ResultExtractor.cpp`
- `../../include/bpfem/core/Constants.hpp`（μ₀ / ε₀）

修改端口归一化时，`PortModeSolver` 求积点缓存、`poyntingPowerIntegral` 实现、`FEMAssembler` 入射幅度、`ResultExtractor` 投影分母四处必须用同一个 `powerNormalizationFactor`。

## 求解器修改后必须验证

至少运行一个单频点，检查残差和 S 参数。直接求解和 fallback 路径都应保持可构建。

## 保持入口精简

`../../src/main.cpp` 应只作为轻量入口，应用流程放在 `../../src/app/Application.cpp`。

## 文档修改保持索引同步

新增或移动文档后，应同步更新：

- 当前目录的 `README.md`。
- 对应顶层兼容入口。
- `../overview.md` 或 `document-map.md` 中的文档地图。
