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


# Discover htslib and libcurl link flags for static linking.
# We run pkg-config --libs --static to get the full transitive dependency
# list and pass the raw flags to the linker.
find_package(PkgConfig)
set(HTS_LINK_FLAGS "")
set(CURL_LINK_FLAGS "")
if(PKG_CONFIG_FOUND)
    execute_process(
        COMMAND ${PKG_CONFIG_EXECUTABLE} --libs --static htslib
        OUTPUT_VARIABLE HTS_LINK_FLAGS
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _hts_result)
    if(NOT _hts_result EQUAL 0)
        set(HTS_LINK_FLAGS "")
    endif()

    execute_process(
        COMMAND ${PKG_CONFIG_EXECUTABLE} --libs --static libcurl
        OUTPUT_VARIABLE CURL_LINK_FLAGS
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _curl_result)
    if(NOT _curl_result EQUAL 0)
        set(CURL_LINK_FLAGS "")
    endif()
endif()

# Fallback if pkg-config is unavailable or failed.
if(NOT HTS_LINK_FLAGS)
    set(HTS_LINK_FLAGS "-lhts -lhtscodecs -ldeflate -llzma -lbz2")
endif()

# Extract -l and -L flags, then filter out libs without a static archive.
# This is needed because pkg-config --static may list libraries that only
# ship as shared (.so) on the system, which breaks -static linking.
function(extract_static_link_flags INPUT_STRING OUTPUT_VAR)
    separate_arguments(_tokens UNIX_COMMAND "${INPUT_STRING}")
    set(_lib_dirs "")
    set(_result "")
    # First pass: collect -L paths.
    foreach(_tok ${_tokens})
        if(_tok MATCHES "^-L(.+)")
            list(APPEND _lib_dirs "${CMAKE_MATCH_1}")
            list(APPEND _result ${_tok})
        endif()
    endforeach()
    # Standard search paths.
    list(APPEND _lib_dirs "/usr/lib" "/usr/lib/x86_64-linux-gnu" "/usr/local/lib"
                          "/usr/lib/aarch64-linux-gnu")
    # Second pass: keep -l flags only if lib<name>.a exists somewhere.
    foreach(_tok ${_tokens})
        if(_tok MATCHES "^-l(.+)")
            set(_name "${CMAKE_MATCH_1}")
            set(_found FALSE)
            foreach(_dir ${_lib_dirs})
                if(EXISTS "${_dir}/lib${_name}.a")
                    set(_found TRUE)
                    break()
                endif()
            endforeach()
            if(_found)
                list(APPEND _result ${_tok})
            else()
                message(STATUS "  Skipping -l${_name} (no static archive found)")
            endif()
        endif()
    endforeach()
    set(${OUTPUT_VAR} ${_result} PARENT_SCOPE)
endfunction()

extract_static_link_flags("${HTS_LINK_FLAGS}" HTS_LINK_FLAGS)

# Only link curl if htslib actually depends on it.
if(HTS_LINK_FLAGS MATCHES "-lcurl")
    extract_static_link_flags("${CURL_LINK_FLAGS}" CURL_LINK_FLAGS)
else()
    set(CURL_LINK_FLAGS "")
endif()

message(STATUS "HTS link flags: ${HTS_LINK_FLAGS}")
message(STATUS "CURL link flags: ${CURL_LINK_FLAGS}")

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
