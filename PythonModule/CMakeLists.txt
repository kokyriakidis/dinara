cmake_minimum_required(VERSION 3.16)

project(dinaraPythonModule)

# C++ dialect.
add_definitions(-std=c++23)

# Compilation warnings.
# Compilation warnings.
if(APPLE)
    add_definitions(-w)
else()
    add_definitions(-Wall -Wconversion -Wno-unused-result -Wno-trigraphs -Wno-psabi -Wunused-parameter -Wno-unused-function)
endif()

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

# Sources files.
file(GLOB SOURCES ../src/*.cpp)

# Fetch SeqAn 2.4.0
include(FetchContent)
FetchContent_Declare(
    seqan
    URL https://github.com/seqan/seqan/archive/refs/tags/seqan-v2.4.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_GetProperties(seqan)
if(NOT seqan_POPULATED)
    FetchContent_Populate(seqan)
endif()

# # Include directory.
# include_directories(../src ~/.dinaraBuild/abpoa/abPOA/include)
include_directories(../src ${seqan_SOURCE_DIR}/include)

# # Additional link directory.
# link_directories(~/.dinaraBuild/abpoa/abPOA/lib)

if(APPLE)
    # Add Homebrew paths for Apple Silicon and Intel Macs
    include_directories(/opt/homebrew/include /usr/local/include)
    link_directories(/opt/homebrew/lib /usr/local/lib)
endif()

# Define our library.
add_library(dinaraPythonModule SHARED ${SOURCES})

# Make sure the library is named dinara.so.
set_target_properties(dinaraPythonModule PROPERTIES OUTPUT_NAME "dinara")
set_target_properties(dinaraPythonModule PROPERTIES PREFIX "")
set_target_properties(dinaraPythonModule PROPERTIES DEFINE_SYMBOL "")
set_target_properties(dinaraPythonModule PROPERTIES INSTALL_RPATH "$ORIGIN")



# Include directories and libraries for Python3 and pybind11.
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
include_directories(${Python3_INCLUDE_DIRS})

# Find pybind11
execute_process(COMMAND "${Python3_EXECUTABLE}" -m pybind11 --includes OUTPUT_VARIABLE DINARA_PYTHON_INCLUDES)
add_definitions(${DINARA_PYTHON_INCLUDES})



# Add astarpa-c include directory
include_directories(${ASTARPA_DIR}/astarpa-c)

# Libraries to link with.
target_link_libraries(
    dinaraPythonModule 
    atomic png boost_program_options boost_serialization pthread z m lapack blas Python3::Module ${ASTARPA_LIB})

# Install the shared library and the abpoa library it needs into the bin directory.
install(TARGETS dinaraPythonModule DESTINATION dinara-install/bin)
