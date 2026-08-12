cmake_minimum_required(VERSION 3.16)

project(dinaraPythonModule)

# C++ dialect. Scoped to CXX only (add_definitions would also apply it to the
# project's one C source, ksw2_extd2_sse.c, which clang rejects -std=c++23 for).
add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-std=c++23>)

# Compilation warnings.
# Keep actionable warnings enabled by default, but avoid conversion-noise flooding
# from legacy code and external dependency headers.
add_definitions(-Wall -Wextra -Wno-unused-result -Wno-trigraphs -Wno-psabi -Wno-unused-parameter -Wno-unused-function)

# Optimization and debug options.
if(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O0)
else(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O3)
    # NDEBUG is required to turn off SeqAn debug code.
    add_definitions(-DNDEBUG)
endif(BUILD_DEBUG)

# 16-byte compare and swap.
# This is recommended for dset64.hpp/dset64-gccAtomic.hpp".
# It's available only on x86 architectures.
if(X86_64)
    add_definitions(-mcx16)
endif()

# Native build.
if(BUILD_NATIVE)
    add_definitions(-march=native)
endif(BUILD_NATIVE)

# This makes sure the code in PythonModule.cpp is compiled.
add_definitions(-DDINARA_PYTHON_API)

add_definitions(-DDINARA_ENABLE_VARIANT_CLUSTERING=0)

# Sources files.
file(GLOB SOURCES CONFIGURE_DEPENDS ../src/*.cpp ../src/edlib/src/edlib.cpp ../src/ksw2/ksw2_extd2_sse.c)
list(FILTER SOURCES EXCLUDE REGEX ".*/AssemblerVariantClustering\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*/mode3-AnchorFromHetSites\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*/ClusterGraph\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*/AssemblerClusterGraph\\.cpp$")
# Theseus disabled: these files use fork-only theseus APIs (align_from,
# multi-segment constructor) that upstream theseus-lib does not provide.
# Excluded from the build until theseus usage is reworked. Callers are stubbed.
list(FILTER SOURCES EXCLUDE REGEX ".*/AssemblerTheseusMarkerGraphMSA\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*/AssemblerTheseusReadWindowMSA\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*/AssemblerPhasingMSA\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*/AssemblerDirectOverlapMSA\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*/AssemblerMultiSegmentMSA\\.cpp$")

# # Include directory.
# include_directories(../src ~/.dinaraBuild/abpoa/abPOA/include)
include_directories(../src ../src/edlib/include ../src/ksw2 ../external/snarls/vendored ../external/snarls/libhandlegraph/src/include ../external/snarls/structures/src/include)
include_directories(../external/hifiasmCandidates)

# # Additional link directory.
# link_directories(~/.dinaraBuild/abpoa/abPOA/lib)

# Point to local build directory
set(DINARA_BUILD_DIR "$ENV{HOME}/.dinaraBuild")
include_directories(SYSTEM ${DINARA_BUILD_DIR}/include)
link_directories(${DINARA_BUILD_DIR}/lib)
list(APPEND CMAKE_PREFIX_PATH "${DINARA_BUILD_DIR}")

# Define our library.
add_library(dinaraPythonModule SHARED ${SOURCES})
set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/MurmurHash2.cpp
    PROPERTIES COMPILE_OPTIONS "-Wno-implicit-fallthrough")

# Make sure the library is named dinara.so.
set_target_properties(dinaraPythonModule PROPERTIES OUTPUT_NAME "dinara")
set_target_properties(dinaraPythonModule PROPERTIES PREFIX "")
set_target_properties(dinaraPythonModule PROPERTIES DEFINE_SYMBOL "")
set_target_properties(dinaraPythonModule PROPERTIES INSTALL_RPATH "$ORIGIN")



# Include directories and libraries for Python3 and pybind11.
execute_process(COMMAND /usr/bin/python3-config --embed --libs OUTPUT_VARIABLE DINARA_PYTHON_LIBRARIES)
execute_process(COMMAND python3 -m pybind11 --includes OUTPUT_VARIABLE DINARA_PYTHON_INCLUDES)
add_definitions(${DINARA_PYTHON_INCLUDES})
string(STRIP ${DINARA_PYTHON_LIBRARIES} DINARA_PYTHON_LIBRARIES)
SET(CMAKE_LINKER_FLAGS  "${CMAKE_LINKER_FLAGS} ${DINARA_PYTHON_LIBRARIES}")



# Find astarpa-c library
find_library(ASTARPA_LIB astarpa_c)
if(NOT ASTARPA_LIB)
    message(FATAL_ERROR "astarpa_c library not found. Please install it.")
endif()

# Find poasta-c library
find_library(POASTA_LIB poasta_c)
if(NOT POASTA_LIB)
    message(FATAL_ERROR "poasta_c library not found. Please install it.")
endif()

# Find simd-minimizers-c library
find_library(SIMD_MINIMIZERS_LIB simd_minimizers_c)
if(NOT SIMD_MINIMIZERS_LIB)
    message(FATAL_ERROR "simd_minimizers_c library not found. Please install it.")
endif()

# Find shasta2 library
find_library(SHASTA2_LIB shasta2)
if(NOT SHASTA2_LIB)
    message(FATAL_ERROR "shasta2 library not found. Please install it.")
endif()

# Allow multiple definitions to resolve duplicate rust_eh_personality
# from linking both libastarpa_c.a and libpoasta_c.a (each bundles Rust stdlib).
if(NOT APPLE)
    set_target_properties(dinaraPythonModule PROPERTIES LINK_FLAGS "-Wl,--allow-multiple-definition")
endif()

# Libraries to link with.
target_link_libraries(
    dinaraPythonModule
    dinara_snarls hifiasm_overlaps atomic png boost_program_options boost_serialization pthread z m lapack blas ${DINARA_PYTHON_LIBRARIES} ${ASTARPA_LIB} ${POASTA_LIB} ${SIMD_MINIMIZERS_LIB} ${SHASTA2_LIB})

# Install the shared library and the abpoa library it needs into the bin directory.
install(TARGETS dinaraPythonModule DESTINATION dinara-install/bin)
