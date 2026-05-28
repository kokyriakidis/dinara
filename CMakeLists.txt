cmake_minimum_required(VERSION 3.16)
project(dinara)

# This is the top level cmake file for a Dinara build. 

# Before starting a build,
# run dinara/scripts/InstallPrerequisites-Ubuntu.sh

# Use "cmake -D" to modify the options defined below 
# in each "option" line. For example, to turn off
# building the dynamic executable use
# cmake -DBUILD_DYNAMIC_EXECUTABLE=OFF



# Figure out if we are on Intel (as opposed to ARM).
set(X86_64 OFF)
if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
    set(X86_64 ON)
endif()



# Summarize the environment.
message(STATUS "Processor architecture is " ${CMAKE_HOST_SYSTEM_PROCESSOR})
message(STATUS "CMAKE_SYSTEM_NAME is " ${CMAKE_SYSTEM_NAME})



# # Decide what we want to build.
# option(BUILD_STATIC_LIBRARY "Build the static library." ON)
# option(BUILD_STATIC_EXECUTABLE "Build the static executable." ON)
# option(BUILD_DYNAMIC_LIBRARY "Build the shared library." ON)
# option(BUILD_DYNAMIC_EXECUTABLE "Build the dynamic executable." OFF)

# Decide what we want to build.
option(BUILD_EXECUTABLE "Build the executable." ON)
option(BUILD_PYTHON_MODULE "Build the Python module." ON)

# Write out what we will actually build.
message(STATUS "BUILD_EXECUTABLE is " ${BUILD_EXECUTABLE})
message(STATUS "BUILD_PYTHON_MODULE is " ${BUILD_PYTHON_MODULE})



# Option to build with -march=native.
option(BUILD_NATIVE "Build with -march=native." OFF)
message(STATUS "BUILD_NATIVE is " ${BUILD_NATIVE})



# Option to request a debug build.
option(BUILD_DEBUG "Make a debuggable build." OFF)
message(STATUS "BUILD_DEBUG is " ${BUILD_DEBUG})


# Discover htslib and libcurl transitive dependencies for static linking.
# pkg-config --static reports all transitive deps, but not all may have
# static archives (.a) available. We filter to only keep libs that exist.
find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(HTS QUIET htslib)
    pkg_check_modules(CURL QUIET libcurl)
endif()

# Helper: filter a list of library names, keeping only those whose
# static archive (lib<name>.a) can be found by the compiler.
function(filter_available_static_libs INPUT_LIBS OUTPUT_VAR)
    set(_result "")
    foreach(_lib ${INPUT_LIBS})
        find_library(_found_${_lib} NAMES "lib${_lib}.a" PATHS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib)
        if(_found_${_lib})
            list(APPEND _result ${_lib})
        endif()
        unset(_found_${_lib} CACHE)
    endforeach()
    set(${OUTPUT_VAR} ${_result} PARENT_SCOPE)
endfunction()

if(HTS_STATIC_LIBRARIES)
    filter_available_static_libs("${HTS_STATIC_LIBRARIES}" HTS_LINK_LIBS)
else()
    set(HTS_LINK_LIBS hts htscodecs deflate lzma bz2)
endif()
if(CURL_STATIC_LIBRARIES)
    filter_available_static_libs("${CURL_STATIC_LIBRARIES}" CURL_LINK_LIBS)
else()
    set(CURL_LINK_LIBS curl)
endif()

# Snarl finder (vendored from vg, stripped of protobuf/SnarlManager).
add_subdirectory(external/snarls)

# Add the subdirectories we need.
if(BUILD_EXECUTABLE)
    add_subdirectory(Executable)
endif(BUILD_EXECUTABLE)

if(BUILD_PYTHON_MODULE)
    add_subdirectory(PythonModule)
endif(BUILD_PYTHON_MODULE)

# Tests
include(CTest) # defines BUILD_TESTING option and enables CTest integration
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()


# Install to the dinara-install directory.
set(CMAKE_INSTALL_PREFIX .)

# Install the scripts.
file(GLOB SCRIPTS scripts/*.py scripts/*.sh)
install(PROGRAMS ${SCRIPTS} DESTINATION dinara-install/bin)

# Install the configuration files.
install(DIRECTORY conf DESTINATION dinara-install)

# Install the docs directory.
install(DIRECTORY docs DESTINATION dinara-install)
