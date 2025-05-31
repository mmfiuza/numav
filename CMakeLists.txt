cmake_minimum_required(VERSION 3.5)
project(numav LANGUAGES CXX)

# CMake policies
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

# Global settings
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# -------------------------------------------------------------------
# Cross-compilation setup MUST COME FIRST
# -------------------------------------------------------------------
if(CMAKE_CROSSCOMPILING)
    message(STATUS "Configuring for Windows cross-compile")
    if(NOT Julia_ROOT)
        message(FATAL_ERROR "Julia_ROOT must be set to Windows Julia path")
    endif()
    
    # Check for DLL name variations
    if(EXISTS "${Julia_ROOT}/bin/julia.dll")
        set(Julia_DLL_NAME "julia.dll")
    elseif(EXISTS "${Julia_ROOT}/bin/libjulia.dll")
        set(Julia_DLL_NAME "libjulia.dll")
    else()
        message(FATAL_ERROR "Could not find julia.dll or libjulia.dll in ${Julia_ROOT}/bin")
    endif()

    # Set paths explicitly
    set(Julia_INCLUDE_DIRS "${Julia_ROOT}/include/julia")
    set(Julia_LIBRARY "${Julia_ROOT}/lib/${Julia_DLL_NAME}.a")
    set(Julia_LIBRARY_DIR "${Julia_ROOT}/lib")
    set(Julia_MAIN_HEADER "${Julia_INCLUDE_DIRS}/julia.h")
    set(JULIA_HOME "${Julia_ROOT}/bin")
    
    # Set as cache variables to prevent re-detection
    set(Julia_FOUND TRUE CACHE BOOL "Julia found" FORCE)
    set(Julia_INCLUDE_DIRS ${Julia_INCLUDE_DIRS} CACHE PATH "Julia include directories" FORCE)
    set(Julia_LIBRARY ${Julia_LIBRARY} CACHE FILEPATH "Julia library" FORCE)
    set(Julia_LIBRARY_DIR ${Julia_LIBRARY_DIR} CACHE PATH "Julia library directory" FORCE)
    set(Julia_MAIN_HEADER ${Julia_MAIN_HEADER} CACHE FILEPATH "Julia main header" FORCE)
    
    # Verify paths exist
    if(NOT EXISTS "${Julia_LIBRARY}")
        message(FATAL_ERROR "Julia import library missing: ${Julia_LIBRARY}")
    endif()
    if(NOT EXISTS "${Julia_MAIN_HEADER}")
        message(FATAL_ERROR "Julia header missing: ${Julia_MAIN_HEADER}")
    endif()
else()
    # Native build
    find_package(Julia REQUIRED)
endif()

include_directories(${Julia_INCLUDE_DIRS})
# -------------------------------------------------------------------

# Output directories
if(MSVC OR CMAKE_CROSSCOMPILING)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/bin")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/bin")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
else()
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
endif()

# RPath settings - skip for Windows/cross-compile
if(NOT CMAKE_CROSSCOMPILING AND NOT WIN32)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${Julia_LIBRARY_DIR}")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    set(CMAKE_MACOSX_RPATH 1)
endif()

# Common compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-parameter -Wextra -Wreorder -fPIC")
endif()

set(CMAKE_CXX_STANDARD 20)
add_definitions(-DJULIA_ENABLE_THREADING)

# Platform-specific flags - FIXED SECTION
if(CMAKE_CROSSCOMPILING OR MINGW)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc -Wa,-mbig-obj")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
    set(CMAKE_SHARED_LIBRARY_PREFIX "lib")
    set(CMAKE_STATIC_LIBRARY_PREFIX "lib")
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