# 并行化与性能

本文给出从单线程原型走向多核 / 多节点 / GPU 的并行化与性能优化路线。

## 1. 现状

- **装配**：`FEMAssembler::assemble` 单线程逐四面体串行；端口装配同样串行。每个体单元的 `(K_e, M_e)` 几何项 + 材料 (`invMu`, `epsR`, `σ`) 在第一次 `assemble()` 时按 11 点 Keast 求积一次性预计算成 upper-triangle packed 缓存（见 `FEMAssembler::ensureElementCache`），扫频时只做 `invMu·K_e[a,b] + (-k0² · εc(f))·M_e[a,b]` 的标量组合，每点装配从 ~12 s 降至 ~1.5 s（一阶 hierarchical / 48k tet）。PEC 约束 dof 集合 (`cachedConstrainedDofs`) 与稀疏图 (`sparsityPattern_`) 也按 lazy 模式跨频点复用。
- **求解**：`MklPardisoSolver` 内部 OpenMP 并行（取决于 oneAPI 配置）；`BiCGStabSolver` 全单线程。
- **I/O**：`OutputWriter::writeVTU` 单线程写 ASCII，扫频时大模型 VTU 输出耗时占比可观。
- **内存**：`SparsePattern::entryIndexByRow` 用 `std::unordered_map<int, std::size_t>`，装配阶段哈希查找有不可忽视的开销，但 symbolic 复用后只在第一次频点付出。元素缓存按 `dof*(dof+1)/2` 上三角 packed 存放，每元 0 阶 ~336 B / 1 阶 ~3.4 KB；48k 单元一阶模型约 160 MB。

## 2. 节点内并行（Phase 2 优先）

### 2.1 OpenMP 装配

- 体装配：每个四面体的局部矩阵相互独立；将 `for tet in tetrahedra` 并行化。
- 临界区策略：
  - 选项 A：每线程持有局部行 buffer (`std::vector<std::map>`)，最后归并；适合稀疏图未构建时。
  - 选项 B：稀疏图固定后，按 `pattern_->entryIndexByRow` 直接到全局 `values_[idx]` 上做 *atomic* 累加；推荐。
- 端口装配：端口面三角形数量远小于体单元，可保持串行或单层并行。

### 2.2 多线程 PARDISO

- `mkl_set_num_threads(N)` + `iparm[2] = N`。
- 当矩阵规模 < 100 万自由度时，超线程会得不偿失，建议提供 `--threads N`。

### 2.3 TBB / Task-based

- 大型扫频任务：每个频点是独立任务，`oneapi::tbb::parallel_for` 跨频点并行（前提是每个 task 各自一份 `MklPardisoSolver` 实例，避免共享 `pt`）。
- 适用场景：直接求解 + 大量频点 + 内存允许；与 AWE/MOR 互斥。

### 2.4 SIMD 友好的局部装配

- 把局部矩阵改为 `std::array<std::complex<double>, 6×6>` 或 `12×12` 静态尺寸。
- 局部基函数求积循环用 `#pragma omp simd`。

### 2.5 内存与缓存

- `SparsePattern::entryIndexByRow`：用排序后的 `std::vector<std::pair<int,int>>` 替换 `unordered_map`，缩小内存占用并友好 cache。
- 端口模式的 `quadrature` 缓存按端口 face 局部排布，避免跨多个 face 跳读。

## 3. 多节点并行（Phase 3）

- **MPI 抽象层**：新增 `parallel/Mpi.hpp`，仅在编译期开启。
- **域分解 (DDM/FETI-DP)**：`solvers-and-sweeps.md` 第 4 节。
- **频点并行**：把频点分批分给多个 MPI rank，各自调用 PARDISO；最简单也最有效。
- **集合通信**：S 参数与远场量在 root 汇总。

## 4. GPU 加速（Phase 3 ~ 4）

- **稀疏求解**：使用 NVIDIA cuSOLVER / AMD rocSOLVER 的 sparse direct，或 cuDSS，PARDISO-GPU。
- **稀疏 SpMV**：iterative solver 内的 SpMV 卸载到 cuSPARSE / hipSPARSE。
- **装配**：GPU 装配难度高（原子写突冲 + 不规则数据流），通常仍在 CPU 完成后只把矩阵迁移到 GPU。
- **混合精度**：复对称用 FP64，预处理可用 FP32。
- 验收：≥ 5×（中等规模 1 ~ 10 M DOF）GPU/CPU 加速比。

## 5. 性能基线与剖析

- 引入 `core/Profiler.{hpp,cpp}`：作用域计时器、阶段累积、JSON 报告。
- 验证 baseline：1M DOF 单频点装配 + 求解 < 60 s（16 核），扫频 21 点 < 20 min；超过则视为退化。
- 集成 *Tracy* 或 *Chrome trace* 输出，便于 GUI 调试。

## 6. I/O 优化

- VTU 切换到二进制（`format="appended"` + `header_type="UInt64"`），并使用 zlib / LZ4 压缩。
- 大场可选 *XDMF + HDF5* 输出，便于 ParaView / VisIt 加载。
- 扫频时支持 `--write-fields-stride N`，避免每点都写。

## 7. 任务清单

| 任务 | 优先级 | 阶段 | 估算 |
|------|--------|------|------|
| OpenMP 体装配 + atomic | 高 | P2 | M |
| 多线程 PARDISO 调优 | 高 | P2 | S |
| 跨频点 TBB 并行 | 中 | P2 | M |
| SparsePattern 内存重排 | 中 | P2 | S |
| Profiler + Chrome trace | 高 | P1 末 | S |
| MPI 频点并行 | 中 | P3 | M |
| MPI + DDM/FETI-DP | 中 | P3 | XL |
| GPU 直接求解（cuDSS） | 低 | P4 | L |
| GPU 迭代 (SpMV) | 低 | P4 | L |
| VTU 二进制 + 压缩 + XDMF/HDF5 | 中 | P2 | M |
