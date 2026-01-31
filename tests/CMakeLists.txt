cmake_minimum_required(VERSION 3.16)
project(dinara_tests)

# Repository root (works both when configuring `tests/` directly and when added
# via `add_subdirectory(tests)` from the repo root).
set(DINARA_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")

# C++ dialect
add_definitions(-std=c++23)

# Debug symbols for tests
add_definitions(-ggdb3 -O0)

# Enable test-only hooks in dinara sources.
add_definitions(-DDINARA_TESTING)

# Variant clustering / ClusterGraph / het-site anchors.
if(DINARA_ENABLE_VARIANT_CLUSTERING)
    add_definitions(-DDINARA_ENABLE_VARIANT_CLUSTERING=1)
else()
    add_definitions(-DDINARA_ENABLE_VARIANT_CLUSTERING=0)
endif()

# Include directories
include_directories(${DINARA_ROOT}/src)
include_directories(${DINARA_ROOT}/srcMain)
include_directories(${DINARA_ROOT}/external/catch2)

# Figure out if we are on Intel (as opposed to ARM).
set(X86_64 OFF)
if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
    set(X86_64 ON)
endif()

# 16-byte compare and swap.
if(X86_64)
    add_definitions(-mcx16)
endif()

# Dinara build directory for libraries
set(DINARA_BUILD_DIR $ENV{HOME}/.dinaraBuild)
include_directories(${DINARA_BUILD_DIR}/include)
link_directories(${DINARA_BUILD_DIR}/lib)

# Catch2 main (provides main())
add_library(catch_main OBJECT catch_main.cpp)

# ============================================================================
# INTEGRATION TESTS (link against dinara sources)
# ============================================================================

# Get all dinara source files EXCEPT main.cpp
file(GLOB DINARA_SRC_FILES CONFIGURE_DEPENDS
    ${DINARA_ROOT}/src/*.cpp
    ${DINARA_ROOT}/srcMain/*.cpp
)
# Remove main.cpp (it has its own main())
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*main\\.cpp$")
if(NOT DINARA_ENABLE_VARIANT_CLUSTERING)
    list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/AssemblerVariantClustering\\.cpp$")
    list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/mode3-AnchorFromHetSites\\.cpp$")
    list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/ClusterGraph\\.cpp$")
    list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/AssemblerClusterGraph\\.cpp$")
endif()

# Find required libraries
find_library(ASTARPA_LIB astarpa_c PATHS ${DINARA_BUILD_DIR}/lib NO_DEFAULT_PATH)
find_library(POASTA_LIB poasta_c PATHS ${DINARA_BUILD_DIR}/lib NO_DEFAULT_PATH)
find_library(SIMD_MINIMIZERS_LIB simd_minimizers_c PATHS ${DINARA_BUILD_DIR}/lib NO_DEFAULT_PATH)
find_library(SHASTA2_LIB shasta2 PATHS ${DINARA_BUILD_DIR}/lib NO_DEFAULT_PATH)
find_library(ABPOA_LIB abpoa PATHS ${DINARA_BUILD_DIR}/lib NO_DEFAULT_PATH)
find_package(ZLIB REQUIRED)

# Integration test executable
add_executable(test_integration
    test_integration.cpp
    ${DINARA_SRC_FILES}
    $<TARGET_OBJECTS:catch_main>
)

# Some toolchains in this environment can produce non-executable output files
# (missing +x), which makes CTest fail with "permission denied".
add_custom_command(TARGET test_integration POST_BUILD
    COMMAND chmod +x $<TARGET_FILE:test_integration>
    VERBATIM
)

# Link all required libraries
if(APPLE)
    target_link_libraries(test_integration
        atomic boost_system boost_program_options boost_chrono boost_serialization
        spoa png lapack blas gfortran m pthread
        ${SHASTA2_LIB} ${ABPOA_LIB} ${ASTARPA_LIB} ${POASTA_LIB} ${SIMD_MINIMIZERS_LIB}
        ZLIB::ZLIB
    )
else()
    # Most of our dependencies are static .a libraries, so we need either strict
    # ordering or a link group to resolve transitive symbols.
    target_link_libraries(test_integration
        atomic boost_system boost_program_options boost_chrono boost_serialization
        spoa png lapack blas gfortran m pthread
        -Wl,--start-group
        ${SHASTA2_LIB} ${ABPOA_LIB} ${ASTARPA_LIB} ${POASTA_LIB} ${SIMD_MINIMIZERS_LIB}
        ZLIB::ZLIB
        -Wl,--end-group
    )
endif()

# Ensure the test binary uses the same curated Dinara dependency builds at runtime.
set_target_properties(test_integration PROPERTIES
    BUILD_RPATH "${DINARA_BUILD_DIR}/lib"
)

# Enable CTest
enable_testing()
add_test(NAME IntegrationTests COMMAND test_integration)

# CI/dev environments can have LD_LIBRARY_PATH pointing at /usr/local/lib first.
# Force our curated dependency builds (in ${DINARA_BUILD_DIR}/lib) to take precedence.
set_tests_properties(IntegrationTests PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=${DINARA_BUILD_DIR}/lib"
)
