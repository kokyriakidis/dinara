cmake_minimum_required(VERSION 3.16)
project(dinaraExecutable)

# C++ dialect.
add_definitions(-std=c++23)

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

add_definitions(-DDINARA_ENABLE_VARIANT_CLUSTERING=0)

# Source files
file(GLOB SOURCES CONFIGURE_DEPENDS ../src/*.cpp ../srcMain/*.cpp ../src/edlib/src/edlib.cpp)
list(FILTER SOURCES EXCLUDE REGEX ".*/AssemblerVariantClustering\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*/mode3-AnchorFromHetSites\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*/ClusterGraph\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*/AssemblerClusterGraph\\.cpp$")

# # Include directory.
# include_directories(../src ~/.dinaraBuild/abpoa/abPOA/include)
include_directories(../src ../src/edlib/include)
include_directories(../external/snarls/vendored ../external/snarls/libhandlegraph/src/include ../external/snarls/structures/src/include)

# Additional link directory.
# link_directories(~/.dinaraBuild/abpoa/abPOA/lib)

# Point to local build directory
set(DINARA_BUILD_DIR "$ENV{HOME}/.dinaraBuild")
include_directories(SYSTEM ${DINARA_BUILD_DIR}/include)
link_directories(${DINARA_BUILD_DIR}/lib)
list(APPEND CMAKE_PREFIX_PATH "${DINARA_BUILD_DIR}")

# On Linux, we build a static executable, so we must find static libraries (.a).
if(NOT APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

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

# Find abpoa library
find_library(ABPOA_LIB abpoa)
if(NOT ABPOA_LIB)
    message(FATAL_ERROR "abpoa library not found. Please install it.")
endif()

# Find shasta2 library
find_library(SHASTA2_LIB shasta2)
if(NOT SHASTA2_LIB)
    message(FATAL_ERROR "shasta2 library not found. Please install it.")
endif()

# Find theseus library (optional — enables journey MSA experiment)
find_library(THESEUS_LIB theseus)
if(THESEUS_LIB)
    add_definitions(-DDINARA_HAVE_THESEUS)
    message(STATUS "theseus found: journey MSA experiment enabled.")
else()
    message(WARNING "theseus library not found. Run InstallPrerequisites to enable journey MSA experiment.")
    set(THESEUS_LIB "")
endif()

# Define our executable.
add_executable(dinaraExecutable ${SOURCES})
set_target_properties(dinaraExecutable PROPERTIES OUTPUT_NAME "dinara")
set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/MurmurHash2.cpp
    PROPERTIES COMPILE_OPTIONS "-Wno-implicit-fallthrough")

# Request a static executable (Linux only).
if(NOT APPLE)
    set_target_properties(dinaraExecutable PROPERTIES LINK_FLAGS "-static -Wl,--allow-multiple-definition" )
endif()

# Libraries to link with.
if(APPLE)
    target_link_libraries(
        dinaraExecutable
        dinara_snarls
        atomic boost_system boost_program_options boost_chrono boost_serialization spoa png z
        lapack blas gfortran z m
        ${ASTARPA_LIB} ${POASTA_LIB} ${SIMD_MINIMIZERS_LIB} ${SHASTA2_LIB} ${ABPOA_LIB} ${THESEUS_LIB})
else()
    # For arcane reasons, statically linking with the pthread
    # library on Linux requires "--whole-archive".
        if(X86_64)
        target_link_libraries(
            dinaraExecutable
            dinara_snarls
            atomic boost_system boost_program_options boost_chrono boost_serialization spoa png z
            lapack blas gfortran quadmath z m
            -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
            ${ASTARPA_LIB} ${POASTA_LIB} ${SIMD_MINIMIZERS_LIB} ${SHASTA2_LIB} ${ABPOA_LIB} ${THESEUS_LIB} hts htscodecs deflate lzma bz2 curl z m)
    else(X86_64)
        target_link_libraries(
            dinaraExecutable
            dinara_snarls
            atomic boost_system boost_program_options boost_chrono boost_serialization spoa png z
            lapack blas gfortran z m
            -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
            ${ASTARPA_LIB} ${POASTA_LIB} ${SIMD_MINIMIZERS_LIB} ${SHASTA2_LIB} ${ABPOA_LIB} ${THESEUS_LIB} z m)
    endif(X86_64)
endif()
  
# The static executable goes to the bin directory.
install(TARGETS dinaraExecutable DESTINATION dinara-install/bin)
