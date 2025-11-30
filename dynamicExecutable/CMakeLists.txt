cmake_minimum_required(VERSION 3.16)
project(dinaraDynamicExecutable)

# C++ dialect.
add_definitions(-std=c++20)

# Compilation warnings.
add_definitions(-Wall -Wconversion -Wno-unused-result -Wno-trigraphs -Wno-psabi -Wunused-parameter)

# Optimization and debug options.
if(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O0)
    # Address sanitizer does not work well with the Python API,
    # so it is commented out here. But it could be useful
    # to debug the dynamic executable. 
    # add_definitions(-fsanitize=address)
else(BUILD_DEBUG)
    add_definitions(-g0)
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

# Build id.
add_definitions(-DBUILD_ID=${BUILD_ID})

# DINARA_LONG_MARKERS.
if(DINARA_LONG_MARKERS)
    add_definitions(-DDINARA_LONG_MARKERS)
endif(DINARA_LONG_MARKERS)

# Definitions needed to eliminate dependency on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# This is needed to avoid some Boost warnings about deprecated messages.
# We canot fix this in Dinara as it is a Boost problem.
# It can be removed if Boost is fixed.
add_definitions(-DBOOST_ALLOW_DEPRECATED_HEADERS)

add_definitions(-DDINARA_PYTHON_API)

# Source files
file(GLOB SOURCES ../srcMain/*.cpp)

# Include directory.
include_directories(../src)

# Define our executable.
add_executable(dinaraDynamicExecutable ${SOURCES})
set_target_properties(dinaraDynamicExecutable PROPERTIES OUTPUT_NAME "dinaraDynamic")

# This is needed to make sure the executable looks for the shared
# library in the same directory.
set_target_properties(dinaraDynamicExecutable PROPERTIES INSTALL_RPATH "\$ORIGIN")

# Libraries to link with.
if(X86_64)
    target_link_libraries(
        dinaraDynamicExecutable
        dinaraDynamicLibrary
        atomic boost_system boost_program_options boost_chrono boost_serialization spoa png z 
        lapack blas quadmath pthread)
else(X86_64)
    target_link_libraries(
        dinaraDynamicExecutable
        dinaraDynamicLibrary
        atomic boost_system boost_program_options boost_chrono boost_serialization spoa png z 
        lapack blas pthread)
endif(X86_64)

if(BUILD_DEBUG) 
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif(BUILD_DEBUG)

# Install the dynamic executable into the bin directory.
install(TARGETS dinaraDynamicExecutable DESTINATION dinara-install/bin)

