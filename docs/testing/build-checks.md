# 构建检查

本文说明不同构建配置的检查命令。

## 默认构建

```powershell
cmake -S . -B build
cmake --build build --config Release
```

## 强制使用迭代 fallback

```powershell
cmake -S . -B build_nomkl -DBPFEM_USE_MKL=OFF
cmake --build build_nomkl --config Release
```

## MKL/PARDISO 构建

```powershell
cmake -S . -B build_mkl -DBPFEM_USE_MKL=ON
cmake --build build_mkl --config Release
```

## 期望行为

- 默认构建能生成 `bp_fem_solver`。
- `-DBPFEM_USE_MKL=OFF` 构建能使用 BiCGSTAB fallback。
- 如果 oneMKL 不可用，CMake 应给出警告并回退到 BiCGSTAB，而不是配置失败。
