# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

cmake_minimum_required(VERSION 3.5)
project(numav LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_TESTS "Build tests" OFF)
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_SHARED_LIBS "Build shared library" OFF)
option(GDB "Build with debug symbols for gdb" OFF)

# find dependencies
find_package(Eigen3 5.0.1 REQUIRED)
find_package(spdlog REQUIRED)

# solver-specific configuration
set(SOLVER "ONEMKL" CACHE STRING "The linear system solver")
set_property(CACHE SOLVER PROPERTY STRINGS "INTERNAL" "EIGEN" "ONEMKL" "MUMPS")
if(SOLVER STREQUAL "INTERNAL")
    set(SELECTED_SOLVER_DEFINITION NUMAV_USE_INTERNAL_SOLVER)
    set(SOLVER_SPECIFIC_LINK_LIBRARIES ldlt_solver)
    add_subdirectory(third-party/ldlt-solver)
elseif(SOLVER STREQUAL "EIGEN")
    set(SELECTED_SOLVER_DEFINITION NUMAV_USE_EIGEN_SOLVER)
    set(SOLVER_SPECIFIC_LINK_LIBRARIES -lblas)
    set(SOLVER_SPECIFIC_DEFINITIONS EIGEN_USE_BLAS)
elseif(SOLVER STREQUAL "ONEMKL")
    set(SELECTED_SOLVER_DEFINITION NUMAV_USE_ONEMKL_SOLVER)
    set(SOLVER_SPECIFIC_INCLUDE_DIRECTORIES
        /opt/intel/oneapi/mkl/2025.2/include
    )
    set(SOLVER_SPECIFIC_LINK_LIBRARIES
        -Wl,--start-group
        /opt/intel/oneapi/mkl/2025.2/lib/libmkl_gf_ilp64.a
        /opt/intel/oneapi/mkl/2025.2/lib/libmkl_gnu_thread.a
        /opt/intel/oneapi/mkl/2025.2/lib/libmkl_core.a
        -Wl,--end-group
        -lgomp
    )
    set(SOLVER_SPECIFIC_DEFINITIONS MKL_ILP64)
elseif(SOLVER STREQUAL "MUMPS")
    set(SELECTED_SOLVER_DEFINITION NUMAV_USE_MUMPS_SOLVER)
    set(SOLVER_SPECIFIC_INCLUDE_DIRECTORIES /usr/include)
    set(SOLVER_SPECIFIC_LINK_LIBRARIES -lzmumps_seq)
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/lib")

# add numav
file(GLOB_RECURSE NUMAV_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE NUMAV_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
add_library(numav ${NUMAV_SOURCES} ${NUMAV_HEADERS})

target_include_directories(numav 
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/third-party/indicators/include
        ${CMAKE_CURRENT_SOURCE_DIR}/third-party/safe-ptr/include
        ${SOLVER_SPECIFIC_INCLUDE_DIRECTORIES}
)

target_compile_definitions(numav PUBLIC 
    FMT_HEADER_ONLY
    SAFE_PTR_DEBUG
    EIGEN_DEFAULT_DENSE_INDEX_TYPE=int64_t
    EIGEN_NO_DEBUG
    EIGEN_NO_AUTOMATIC_RESIZING
    EIGEN_FAST_MATH
    ${SELECTED_SOLVER_DEFINITION}
    ${SOLVER_SPECIFIC_DEFINITIONS}
)

target_link_libraries(numav PRIVATE
    Eigen3::Eigen
    spdlog::spdlog_header_only
    ${SOLVER_SPECIFIC_LINK_LIBRARIES}
)

if(NOT GDB)
    set(COMPILER_FLAGS -Ofast -ffast-math)
    set(LINKER_FLAGS -flto)
else()
    set(COMPILER_FLAGS -g -O0 -fno-inline)
endif()

target_compile_options(numav PRIVATE -Wall -Wextra -m64 ${COMPILER_FLAGS})
target_link_options(numav PRIVATE -m64 ${LINKER_FLAGS})

# add subdirectories
if(BUILD_TESTS)
    add_subdirectory(tests)
endif()
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
