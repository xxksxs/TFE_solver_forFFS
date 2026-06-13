cmake_minimum_required(VERSION 3.20)
project(bp_frequency_fem_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BPFEM_USE_MKL "Use Intel oneMKL PARDISO sparse direct solver" ON)
option(BPFEM_OUTPUT_PARAVIEW
       "Emit ParaView-style VTK Lagrange tetrahedron VTUs (cell type 71). \
When OFF (default), emit VisIt-friendly cells: VTK_TETRA (type 10) at order 1 \
and VTK_QUADRATIC_TETRA (type 24) at order 2. Order 3 is ParaView-only."
       OFF)

add_library(bp_fem_core STATIC
    src/app/Application.cpp
    src/apm/AnalyticPortBuilder.cpp
    src/bc/AbsorbingBC.cpp
    src/bc/FiniteConductorBC.cpp
    src/bc/ImpedanceBC.cpp
    src/bc/WavePortBC.cpp
    src/core/Logger.cpp
    src/core/Math.cpp
    src/core/RunEnvironment.cpp
    src/core/RunReport.cpp
    src/core/Types.cpp
    src/core/Utilities.cpp
    src/factory/SparseSolverFactory.cpp
    src/factory/SweepStrategyFactory.cpp
    src/fastsweep/PadeApproximant.cpp
    src/fastsweep/PolynomialMomentRecurrence.cpp
    src/fem/EdgeTopology.cpp
    src/fem/FEMAssembler.cpp
    src/fem/PortModeSolver.cpp
    src/io/AEDTParser.cpp
    src/io/NGMeshParser.cpp
    src/linalg/BiCGStabBackend.cpp
    src/linalg/BiCGStabSolver.cpp
    src/linalg/GmresBackend.cpp
    src/linalg/GmresSolver.cpp
    src/linalg/MklPardisoSolver.cpp
    src/linalg/PardisoBackend.cpp
    src/linalg/PreconIlu0.cpp
    src/linalg/PreconJacobi.cpp
    src/linalg/SparseMatrix.cpp
    src/post/OutputWriter.cpp
    src/post/ResultExtractor.cpp
    src/sweep/AlpsSweep.cpp
    src/sweep/AweSweep.cpp
    src/sweep/DirectSweep.cpp
    src/sweep/GaweSweep.cpp
    src/sweep/MgaweSweep.cpp
    src/sweep/WcaweSweep.cpp
    src/tfe/TransfiniteElementBuilder.cpp
)

target_include_directories(bp_fem_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(BPFEM_OUTPUT_PARAVIEW)
    target_compile_definitions(bp_fem_core PUBLIC BPFEM_OUTPUT_PARAVIEW)
endif()

if(BPFEM_USE_MKL)
    find_package(MKL CONFIG QUIET)
    if(MKL_FOUND)
        target_compile_definitions(bp_fem_core PUBLIC BPFEM_USE_MKL)
        target_link_libraries(bp_fem_core PUBLIC MKL::MKL)
    else()
        message(WARNING "Intel oneMKL was not found. Falling back to BiCGSTAB. Configure from a oneAPI environment or set MKL_DIR to enable PARDISO.")
    endif()
endif()

add_executable(bp_fem_solver src/main.cpp)
target_link_libraries(bp_fem_solver PRIVATE bp_fem_core)

# Logger uses GetProcessMemoryInfo (Psapi) on Windows for the working-set /
# peak-memory counters in the run summary. Link the import library that
# provides those symbols. Modern Windows merges psapi into kernel32 at
# runtime, but the import library is still required.
if(WIN32)
    target_link_libraries(bp_fem_core PRIVATE psapi)
endif()

if(MSVC)
    target_compile_options(bp_fem_core PRIVATE /W4 /permissive- /EHsc)
    target_compile_options(bp_fem_solver PRIVATE /W4 /permissive- /EHsc)
else()
    target_compile_options(bp_fem_core PRIVATE -Wall -Wextra -pedantic)
    target_compile_options(bp_fem_solver PRIVATE -Wall -Wextra -pedantic)
endif()
