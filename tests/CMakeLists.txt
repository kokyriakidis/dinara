cmake_minimum_required(VERSION 3.16)
project(dinara_tests)

# C++ dialect
add_definitions(-std=c++23)

# Debug symbols for tests
add_definitions(-ggdb3 -O0)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/../src)
include_directories(${CMAKE_SOURCE_DIR}/../srcMain)
include_directories(${CMAKE_SOURCE_DIR}/../external/catch2)

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
file(GLOB DINARA_SRC_FILES 
    ${CMAKE_SOURCE_DIR}/../src/*.cpp
    ${CMAKE_SOURCE_DIR}/../srcMain/*.cpp
)
# Remove main.cpp (it has its own main())
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*main\\.cpp$")

# Find required libraries
find_library(ASTARPA_LIB astarpa_c HINTS ${DINARA_BUILD_DIR}/lib)
find_library(POASTA_LIB poasta_c HINTS ${DINARA_BUILD_DIR}/lib)
find_library(SIMD_MINIMIZERS_LIB simd_minimizers_c HINTS ${DINARA_BUILD_DIR}/lib)
find_library(SHASTA2_LIB shasta2 HINTS ${DINARA_BUILD_DIR}/lib)
find_library(ABPOA_LIB abpoa HINTS ${DINARA_BUILD_DIR}/lib)

# Integration test executable
add_executable(test_integration
    test_integration.cpp
    ${DINARA_SRC_FILES}
    $<TARGET_OBJECTS:catch_main>
)

# Link all required libraries
target_link_libraries(test_integration
    atomic boost_system boost_program_options boost_chrono boost_serialization 
    spoa png z lapack blas gfortran z m pthread
    ${ASTARPA_LIB} ${POASTA_LIB} ${SIMD_MINIMIZERS_LIB} ${SHASTA2_LIB} ${ABPOA_LIB}
)

# Enable CTest
enable_testing()
add_test(NAME IntegrationTests COMMAND test_integration)
