# Record the sources seen at build time, including dirty and untracked code.
file(GLOB_RECURSE sources RELATIVE "${SOURCE_DIR}" LIST_DIRECTORIES false
    "${SOURCE_DIR}/src/*.cpp" "${SOURCE_DIR}/src/*.hpp" "${SOURCE_DIR}/include/*.h"
    "${SOURCE_DIR}/include/*.hpp" "${SOURCE_DIR}/cmake/*.cmake")
list(APPEND sources CMakeLists.txt)
list(SORT sources)
set(manifest "")
foreach(source IN LISTS sources)
    file(SHA256 "${SOURCE_DIR}/${source}" digest)
    string(APPEND manifest "${source} ${digest}\n")
endforeach()
string(SHA256 source_sha256 "${manifest}")

set(revision unknown)
set(state unknown)
if(GIT AND EXISTS "${GIT}")
    execute_process(COMMAND "${GIT}" --no-optional-locks -c "safe.directory=${SOURCE_DIR}"
        rev-parse HEAD WORKING_DIRECTORY "${SOURCE_DIR}" OUTPUT_VARIABLE head
        OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE head_status ERROR_QUIET)
    if(head_status EQUAL 0 AND head MATCHES "^[a-f0-9]+$")
        set(revision "${head}")
        execute_process(COMMAND "${GIT}" --no-optional-locks -c "safe.directory=${SOURCE_DIR}"
            status --porcelain --untracked-files=normal -- CMakeLists.txt cmake include src
            WORKING_DIRECTORY "${SOURCE_DIR}" OUTPUT_VARIABLE changes
            OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE status ERROR_QUIET)
        if(status EQUAL 0)
            if(changes STREQUAL "")
                set(state clean)
            else()
                set(state dirty)
            endif()
        endif()
    endif()
endif()

file(MAKE_DIRECTORY "${OUTPUT_DIR}")
file(CONFIGURE OUTPUT "${OUTPUT_DIR}/BpfemBuildProvenance.hpp" CONTENT
    "#pragma once\n#define BPFEM_GIT_REVISION \"@revision@\"\n#define BPFEM_GIT_STATE \"@state@\"\n#define BPFEM_SOURCE_SHA256 \"@source_sha256@\"\n" @ONLY)
