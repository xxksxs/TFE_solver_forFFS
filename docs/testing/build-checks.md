# 本地构建与回归

在仓库根目录运行。Windows 使用已配置为 x64 的 Visual Studio Developer PowerShell，确保 CMake 和 Ninja 可用。如果 Ninja 未加入 PATH，可通过 `-NinjaPath` 指定。

```powershell
.\scripts\run-local-regression.ps1 -Profile native -Configuration Debug,Release
.\scripts\run-local-regression.ps1 -Profile mkl -Configuration Release -MklDir 'C:/oneAPI/mkl/latest/lib/cmake/mkl'
.\scripts\run-local-regression.ps1 -Profile mkl-mumps -Configuration Release -DependencyCache 'D:/bpfem/solver-dependencies.cmake'
```

每个后端组合和配置使用独立目录，例如 `build/local-regression/native-Debug`。通过 `-BuildRoot` 可将构建、CTest 日志和 CLI 烟测输出放到其他位置。重复执行相同命令时复用原构建目录。三个配置均从 CTest 运行 `fastsweep_tests`、`io_parser_tests` 和 `sparse_solver_tests`。

| 配置 | 稀疏求解回归实际执行的后端 |
| --- | --- |
| `native` | GMRES、BiCGSTAB |
| `mkl` | GMRES、BiCGSTAB、oneMKL PARDISO 与稀疏矩阵向量乘 |
| `mkl-mumps` | `mkl` 的覆盖范围，另含 MUMPS 的 in-core、BLR 和 out-of-core |

`mkl` 和 `mkl-mumps` 要求 CMake 实际找到 oneMKL。外部库缺失会导致配置、链接或测试失败。手工调用 CMake 时，仍可通过 `BPFEM_REQUIRE_MKL` 决定 oneMKL 缺失时是否允许回退。

## 外部依赖

oneMKL 可以由已初始化的 oneAPI 环境、`-MklDir` 或传给 `-DependencyCache` 的本机 CMake 初始缓存定位。MUMPS 使用预先编译好的顺序版复数双精度 `ZMUMPS`。例如，本机初始缓存可以按以下形式编写，路径和库名以实际安装为准。

```cmake
set(MKL_DIR "D:/oneAPI/mkl/lib/cmake/mkl" CACHE PATH "")
set(BPFEM_MUMPS_INCLUDE_DIR "D:/mumps/include" CACHE PATH "")
set(BPFEM_MUMPS_LIBRARIES
    "D:/mumps/lib/zmumps.lib;D:/mumps/lib/mumps_common.lib;D:/mumps/lib/pord.lib;D:/mumps/lib/mpiseq.lib"
    CACHE STRING "")
set(BPFEM_MUMPS_LINK_DIRECTORIES "D:/oneAPI/compiler/lib;D:/oneAPI/mkl/lib" CACHE STRING "")
set(BPFEM_MUMPS_RUNTIME_DLLS "D:/oneAPI/compiler/bin/libifcoremd.dll;D:/oneAPI/compiler/bin/libifportmd.dll"
    CACHE STRING "")
```

链接库列表需包含该 ZMUMPS 构建所需的排序、顺序 MPI、BLAS/LAPACK 等依赖，或从链接目录解析相应库。Windows 的运行库列表需按实际 ZMUMPS 构建补齐，CMake 会将这些 DLL 复制到可执行文件目录。更换依赖安装、工具链或链接选项时使用新的 `-BuildRoot`，以便保留每组配置的独立缓存。

## 检查头文件增量依赖

此前受影响的中文 MSVC＋Ninja 构建曾生成乱码的 `msvc_deps_prefix`，导致 Ninja 漏记头文件，混用新旧对象。项目配置现在按编译器输出的原始字节识别前缀。如果旧目录已经丢失依赖记录，首次切换到修复后的配置时用新目录完成一次基线构建，此后可增量复用。

在隔离的源码副本完成一次基线构建，记录 `GmresBackend.cpp.obj` 和 `SparseSolverFactory.cpp.obj` 的更新时间，再给 `include/bpfem/linalg/ISparseSolver.hpp` 增加一行注释。用同一个 `-BuildRoot` 重跑上面的回归命令。确认两个对象都重新编译且三项测试通过。Windows Ninja 可以用以下命令核对对象的头文件依赖。

```powershell
& ninja -C build/local-regression/native-Release -t deps 'CMakeFiles/bp_fem_core.dir/src/linalg/GmresBackend.cpp.obj'
```

## CLI 烟测接入

使用 `-SmokeScript <脚本路径>`，入口会在每组 CTest 通过后调用该 PowerShell 脚本。脚本接收 `-SolverPath`（当前构建的可执行文件）、`-OutputDir`（当前构建下的 `smoke` 目录）和 `-BackendProfile`（`native`、`mkl` 或 `mkl-mumps`）。脚本应将所有产物放入给定目录，失败时抛出异常或返回非零退出码。
