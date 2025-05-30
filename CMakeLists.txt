cmake_minimum_required(VERSION 3.5)
project(numav LANGUAGES CXX)

# CMake policies
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

# Global settings
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(Julia REQUIRED)

# Output directories
if(MSVC)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/bin")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/bin")
else()
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
endif()

# RPath settings
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${Julia_LIBRARY_DIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)

# Common compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-parameter -Wextra -Wreorder -fPIC")
endif()

set(CMAKE_CXX_STANDARD 20)
add_definitions(-DJULIA_ENABLE_THREADING)

# Platform-specific flags
if(WIN32)
    if(DEFINED ENV{MSYSTEM} AND $ENV{MSYSTEM} STREQUAL "MINGW32")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=pentium4 -static-libstdc++")
    endif()
    if(MINGW)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj -O2")
    endif()
    if(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
        set(CMAKE_SHARED_LIBRARY_PREFIX "lib")
        set(CMAKE_STATIC_LIBRARY_PREFIX "lib")
    endif()
endif()

# Include required modules
include(CheckCXXCompilerFlag)
include(CMakePackageConfigHelpers)
include(GenerateExportHeader)
include(GNUInstallDirs)
include(InstallRequiredSystemLibraries)

# Build components
add_subdirectory(src)              # Main library
add_subdirectory(third-party/cxxwrap)  # CxxWrap dependency
option(BUILD_JULIA_BINDINGS "Build Julia bindings" ON)
if(BUILD_JULIA_BINDINGS)
    add_subdirectory(bindings/julia)   # Julia bindings
endif()