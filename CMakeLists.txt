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
find_package(Eigen3 3.4 REQUIRED)
find_package(spdlog REQUIRED)

# Set MKL paths
set(MKLROOT "/opt/intel/oneapi/mkl/2025.2")
set(MKL_COMPILER_ROOT "/opt/intel/oneapi/compiler/2025.2")

# add third-party subdirectories
add_subdirectory(third-party/ldlt-solver)

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
        ${MKLROOT}/include
)

target_compile_definitions(numav PUBLIC 
    MKL_ILP64
    EIGEN_DEFAULT_DENSE_INDEX_TYPE=int64_t
    FMT_HEADER_ONLY
)

if(NOT GDB)
    target_compile_options(numav PRIVATE -Wall -Wextra -O3 -m64)
    target_link_options(numav PRIVATE -flto -m64)
else()
    target_compile_options(numav PRIVATE -Wall -Wextra -g -O0 -fno-inline -m64)
    target_link_options(numav PRIVATE -m64)
endif()

target_link_libraries(numav PRIVATE
    Eigen3::Eigen
    spdlog::spdlog_header_only
    ldlt_solver
    -Wl,--start-group
    ${MKLROOT}/lib/libmkl_gf_ilp64.a 
    ${MKLROOT}/lib/libmkl_gnu_thread.a
    ${MKLROOT}/lib/libmkl_core.a
    -Wl,--end-group
    -lgomp
    -lpthread -lm -ldl
)

# add subdirectories
if(BUILD_TESTS)
    add_subdirectory(tests)
endif()
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
