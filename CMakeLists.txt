cmake_minimum_required(VERSION 3.20)
project(bp_frequency_fem_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC" AND CMAKE_GENERATOR MATCHES "^Ninja")
    # Ninja compares the raw /showIncludes bytes. CMake's compiler-ID probe
    # can transcode them using a different code page on localized Windows.
    set(_include_probe_dir "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/bpfem-showincludes")
    file(MAKE_DIRECTORY "${_include_probe_dir}")
    file(WRITE "${_include_probe_dir}/probe.h" "")
    file(WRITE "${_include_probe_dir}/probe.cpp" "#include \"probe.h\"\n")
    execute_process(
        COMMAND "${CMAKE_CXX_COMPILER}" /nologo /utf-8 /showIncludes /Zs probe.cpp
        WORKING_DIRECTORY "${_include_probe_dir}"
        OUTPUT_VARIABLE _include_probe_output
        ERROR_VARIABLE _include_probe_error
        RESULT_VARIABLE _include_probe_status
        ENCODING NONE)
    if(NOT _include_probe_status STREQUAL "0" OR NOT _include_probe_output MATCHES
        "(^|\n)([^\r\n]+[ \t])([A-Za-z]:[/\\\\]|\\./|\\.\\\\)[^\r\n]*probe\\.h")
        message(FATAL_ERROR "Cannot detect the MSVC /showIncludes prefix for Ninja.\n"
            "${_include_probe_output}${_include_probe_error}")
    endif()
    set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_MATCH_2}")
endif()

option(BPFEM_USE_MKL "Use Intel oneMKL PARDISO sparse direct solver" ON)
option(BPFEM_REQUIRE_MKL "Fail configuration if oneMKL is unavailable" OFF)
option(BPFEM_USE_MUMPS "Enable the double-complex MUMPS backend" OFF)
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
    src/fastsweep/FastSweepDiagnostics.cpp
    src/fastsweep/GalerkinReducedModel.cpp
    src/fastsweep/LinearAlgebra.cpp
    src/fastsweep/LanczosPadeModel.cpp
    src/fastsweep/PadeApproximant.cpp
    src/fastsweep/PolynomialMomentRecurrence.cpp
    src/fastsweep/PolynomialPortMomentBuilder.cpp
    src/fastsweep/PolynomialWcaweRecurrence.cpp
    src/fastsweep/PortModeUtilities.cpp
    src/fastsweep/UpperTriangularBlockAction.cpp
    src/fastsweep/WellConditionedBasisBuilder.cpp
    src/fem/EdgeTopology.cpp
    src/fem/FEMAssembler.cpp
    src/fem/PortModeSolver.cpp
    src/io/AEDTParser.cpp
    src/io/NGMeshParser.cpp
    src/io/PortFaceResolver.cpp
    src/linalg/BiCGStabBackend.cpp
    src/linalg/BiCGStabSolver.cpp
    src/linalg/GmresBackend.cpp
    src/linalg/GmresSolver.cpp
    src/linalg/MklPardisoSolver.cpp
    src/linalg/IterativeOperator.cpp
    src/linalg/MumpsBackend.cpp
    src/linalg/PardisoBackend.cpp
    src/linalg/PreconIlu0.cpp
    src/linalg/PreconJacobi.cpp
    src/linalg/PreconSchwarz.cpp
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

find_package(Git QUIET)
set(_provenance_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
add_custom_target(bpfem_build_provenance
    COMMAND "${CMAKE_COMMAND}" "-DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
            "-DOUTPUT_DIR=${_provenance_dir}" "-DGIT=${GIT_EXECUTABLE}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/BuildProvenance.cmake"
    BYPRODUCTS "${_provenance_dir}/BpfemBuildProvenance.hpp"
    VERBATIM)
add_dependencies(bp_fem_core bpfem_build_provenance)
target_include_directories(bp_fem_core PRIVATE "${_provenance_dir}")
if(NOT WIN32)
    find_package(OpenSSL REQUIRED COMPONENTS Crypto)
    target_link_libraries(bp_fem_core PRIVATE OpenSSL::Crypto)
endif()

if(BPFEM_OUTPUT_PARAVIEW)
    target_compile_definitions(bp_fem_core PUBLIC BPFEM_OUTPUT_PARAVIEW)
endif()

set(_solver_backends "GMRES, BiCGSTAB")
if(BPFEM_REQUIRE_MKL AND NOT BPFEM_USE_MKL)
    message(FATAL_ERROR "BPFEM_REQUIRE_MKL requires BPFEM_USE_MKL=ON.")
endif()
if(BPFEM_USE_MKL)
    find_package(MKL CONFIG QUIET)
    if(MKL_FOUND)
        target_compile_definitions(bp_fem_core PUBLIC BPFEM_USE_MKL)
        target_link_libraries(bp_fem_core PUBLIC MKL::MKL)
        string(APPEND _solver_backends ", oneMKL PARDISO")
    elseif(BPFEM_REQUIRE_MKL)
        message(FATAL_ERROR "oneMKL is required for this build. Set MKL_DIR to the directory with MKLConfig.cmake.")
    else()
        message(WARNING "Intel oneMKL was not found. Falling back to BiCGSTAB. Configure from a oneAPI environment or set MKL_DIR to enable PARDISO.")
    endif()
endif()

if(BPFEM_USE_MUMPS)
    find_path(BPFEM_MUMPS_INCLUDE_DIR NAMES zmumps_c.h REQUIRED)
    set(BPFEM_MUMPS_LIBRARIES "" CACHE STRING
        "Link libraries for a sequential ZMUMPS build, including dependencies")
    set(BPFEM_MUMPS_RUNTIME_DLLS "" CACHE STRING
        "Runtime DLL paths for the selected ZMUMPS build on Windows")
    set(BPFEM_MUMPS_LINK_DIRECTORIES "" CACHE STRING
        "Additional directories for dependencies named by the ZMUMPS libraries")
    if(NOT BPFEM_MUMPS_LIBRARIES)
        message(FATAL_ERROR "Set BPFEM_MUMPS_LIBRARIES to the sequential ZMUMPS libraries and dependencies.")
    endif()
    target_compile_definitions(bp_fem_core PUBLIC BPFEM_USE_MUMPS)
    target_include_directories(bp_fem_core PRIVATE "${BPFEM_MUMPS_INCLUDE_DIR}")
    if(BPFEM_MUMPS_LINK_DIRECTORIES)
        target_link_directories(bp_fem_core INTERFACE ${BPFEM_MUMPS_LINK_DIRECTORIES})
    endif()
    target_link_libraries(bp_fem_core PRIVATE ${BPFEM_MUMPS_LIBRARIES})
    string(APPEND _solver_backends ", ZMUMPS")
endif()
message(STATUS "BP-FEM backends: ${_solver_backends}")

add_executable(bp_fem_solver src/main.cpp)
target_link_libraries(bp_fem_solver PRIVATE bp_fem_core)

enable_testing()
add_executable(fastsweep_tests tests/fastsweep_tests.cpp)
target_link_libraries(fastsweep_tests PRIVATE bp_fem_core)
add_test(NAME fastsweep_tests COMMAND fastsweep_tests)
add_executable(io_parser_tests tests/io_parser_tests.cpp)
target_link_libraries(io_parser_tests PRIVATE bp_fem_core)
add_test(NAME io_parser_tests COMMAND io_parser_tests)
add_executable(sparse_solver_tests tests/sparse_solver_tests.cpp)
find_package(Threads REQUIRED)
target_link_libraries(sparse_solver_tests PRIVATE bp_fem_core Threads::Threads)
add_test(NAME sparse_solver_tests COMMAND sparse_solver_tests)

set(_solver_runtime_files)
if(WIN32 AND BPFEM_USE_MKL AND MKL_FOUND AND NOT MKL_LINK STREQUAL "static")
    # MKL loads architecture-specific kernels at runtime. They are not all
    # visible in the executable's import table.
    set(_mkl_runtime_stems
        mkl_core mkl_def mkl_mc3 mkl_avx2 mkl_avx512
        mkl_vml_def mkl_vml_mc3 mkl_vml_avx2 mkl_vml_avx512 mkl_vml_cmpt)
    if(MKL_LINK STREQUAL "sdl")
        list(APPEND _mkl_runtime_stems mkl_rt mkl_sequential mkl_intel_thread mkl_tbb_thread)
    else()
        list(APPEND _mkl_runtime_stems "mkl_${MKL_THREADING}")
    endif()
    foreach(_stem IN LISTS _mkl_runtime_stems)
        file(GLOB _matches "${MKL_DLL_DIR}/${_stem}.*.dll")
        if(_matches)
            list(SORT _matches COMPARE NATURAL ORDER DESCENDING)
            list(GET _matches 0 _runtime)
            list(APPEND _solver_runtime_files "${_runtime}")
        endif()
    endforeach()
    if(OMP_DLL_DIR AND EXISTS "${OMP_DLL_DIR}/libiomp5md.dll")
        list(APPEND _solver_runtime_files "${OMP_DLL_DIR}/libiomp5md.dll")
    endif()
    if(TARGET TBB::tbb)
        list(APPEND _solver_runtime_files "$<TARGET_FILE:TBB::tbb>")
    endif()
endif()
if(WIN32 AND BPFEM_USE_MUMPS)
    list(APPEND _solver_runtime_files ${BPFEM_MUMPS_RUNTIME_DLLS})
endif()
if(_solver_runtime_files)
    add_custom_target(bpfem_solver_runtime
        COMMAND "${CMAKE_COMMAND}" -E make_directory "$<TARGET_FILE_DIR:bp_fem_solver>"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
            ${_solver_runtime_files} "$<TARGET_FILE_DIR:bp_fem_solver>"
        VERBATIM)
    add_dependencies(bp_fem_solver bpfem_solver_runtime)
    add_dependencies(fastsweep_tests bpfem_solver_runtime)
    add_dependencies(io_parser_tests bpfem_solver_runtime)
    add_dependencies(sparse_solver_tests bpfem_solver_runtime)
endif()

# Logger uses GetProcessMemoryInfo (Psapi) on Windows for the working-set /
# peak-memory counters in the run summary. Link the import library that
# provides those symbols. Modern Windows merges psapi into kernel32 at
# runtime, but the import library is still required.
if(WIN32)
    target_link_libraries(bp_fem_core PRIVATE psapi bcrypt)
endif()

if(MSVC)
    target_compile_options(bp_fem_core PRIVATE /W4 /permissive- /EHsc /utf-8)
    target_compile_options(bp_fem_solver PRIVATE /W4 /permissive- /EHsc /utf-8)
    target_compile_options(fastsweep_tests PRIVATE /W4 /permissive- /EHsc /utf-8)
    target_compile_options(io_parser_tests PRIVATE /W4 /permissive- /EHsc /utf-8)
    target_compile_options(sparse_solver_tests PRIVATE /W4 /permissive- /EHsc /utf-8)
else()
    target_compile_options(bp_fem_core PRIVATE -Wall -Wextra -pedantic)
    target_compile_options(bp_fem_solver PRIVATE -Wall -Wextra -pedantic)
    target_compile_options(fastsweep_tests PRIVATE -Wall -Wextra -pedantic)
    target_compile_options(io_parser_tests PRIVATE -Wall -Wextra -pedantic)
    target_compile_options(sparse_solver_tests PRIVATE -Wall -Wextra -pedantic)
endif()
