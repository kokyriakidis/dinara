cmake_minimum_required(VERSION 3.16)
project(dinaraExecutable)

# C++ dialect.
add_definitions(-std=c++23)

# Compilation warnings.
add_definitions(-Wall -Wconversion -Wno-unused-result -Wno-trigraphs -Wno-psabi -Wunused-parameter -Wno-unused-function)

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

# Source files
file(GLOB SOURCES ../src/*.cpp ../srcMain/*.cpp)

# # Include directory.
# include_directories(../src ~/.dinaraBuild/abpoa/abPOA/include)
include_directories(../src)

# # Additional link directory.
# link_directories(~/.dinaraBuild/abpoa/abPOA/lib)

# Find astarpa-c library
find_library(ASTARPA_LIB astarpa_c)
if(NOT ASTARPA_LIB)
    message(FATAL_ERROR "astarpa_c library not found. Please install it.")
endif()

# Define our executable.
add_executable(dinaraExecutable ${SOURCES})
set_target_properties(dinaraExecutable PROPERTIES OUTPUT_NAME "dinara")

# Request a static executable (Linux only).
if(NOT APPLE)
    set_target_properties(dinaraExecutable PROPERTIES LINK_FLAGS "-static" )
endif()

# Libraries to link with.
if(APPLE)
    target_link_libraries(
        dinaraExecutable
        atomic boost_system boost_program_options boost_chrono boost_serialization spoa png z
        lapack blas gfortran z m
        ${ASTARPA_LIB})
else()
    # For arcane reasons, statically linking with the pthread
    # library on Linux requires "--whole-archive".
    if(X86_64)
        target_link_libraries(
            dinaraExecutable
            atomic boost_system boost_program_options boost_chrono boost_serialization spoa png z
            lapack blas gfortran quadmath z m
            -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
            ${ASTARPA_LIB})
    else(X86_64)
        target_link_libraries(
            dinaraExecutable
            atomic boost_system boost_program_options boost_chrono boost_serialization spoa png z
            lapack blas gfortran z m
            -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
            ${ASTARPA_LIB})
    endif(X86_64) 
endif()
  
# The static executable goes to the bin directory.
install(TARGETS dinaraExecutable DESTINATION dinara-install/bin)