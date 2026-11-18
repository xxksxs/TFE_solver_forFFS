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

## 串行 ZMUMPS 构建

选择与当前 C++ 工具链兼容的复双精度串行 ZMUMPS。将 Fortran、BLAS 和 LAPACK 的链接目录加入构建环境，并给出链接库及运行时 DLL。示例中变量为对应安装的绝对路径，多个库使用分号分隔。

```powershell
cmake -S . -B build_combined -DBPFEM_USE_MKL=ON -DBPFEM_USE_MUMPS=ON `
  "-DBPFEM_MUMPS_INCLUDE_DIR=$mumpsInclude" "-DBPFEM_MUMPS_LIBRARIES=$mumpsLibraries" `
  "-DBPFEM_MUMPS_RUNTIME_DLLS=$mumpsRuntimeDlls"
cmake --build build_combined --config Release
ctest --test-dir build_combined -C Release --output-on-failure
```

再用独立目录将 `BPFEM_USE_MKL` 设为 `OFF` 复跑 MUMPS。核对进程实际加载的库版本和架构。禁用应用的 MKL 选项时，ZMUMPS 仍会加载其编译时选用的 BLAS/LAPACK。

## 期望行为

- 默认构建能生成 `bp_fem_solver`。
- `-DBPFEM_USE_MKL=OFF` 构建能使用 BiCGSTAB fallback。
- 如果 oneMKL 不可用，CMake 应给出警告并回退到 BiCGSTAB，而不是配置失败。
