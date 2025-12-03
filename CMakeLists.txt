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


# Add the subdirectories we need.
if(BUILD_EXECUTABLE)
    add_subdirectory(Executable)
endif(BUILD_EXECUTABLE)

if(BUILD_PYTHON_MODULE)
    add_subdirectory(PythonModule)
endif(BUILD_PYTHON_MODULE)


# Install to the dinara-install directory.
set(CMAKE_INSTALL_PREFIX .)

# Install the scripts.
file(GLOB SCRIPTS scripts/*.py scripts/*.sh)
install(PROGRAMS ${SCRIPTS} DESTINATION dinara-install/bin)

# Install the configuration files.
install(DIRECTORY conf DESTINATION dinara-install)

# Install the docs directory.
install(DIRECTORY docs DESTINATION dinara-install)


