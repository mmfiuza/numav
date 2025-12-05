# Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

cmake_minimum_required(VERSION 3.5)
project(numav LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

# add subdirectories
if(BUILD_TESTS)
    add_subdirectory(tests)
endif()
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# find dependencies
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(spdlog REQUIRED)

# Set MKL paths (adjust versions if needed)
set(MKLROOT "/opt/intel/oneapi/mkl/2025.2")
set(MKL_COMPILER_ROOT "/opt/intel/oneapi/compiler/2025.2")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/lib")

# add numav
file(GLOB_RECURSE NUMAV_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE NUMAV_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
if(NOT DYNAMIC)
    add_library(numav STATIC ${NUMAV_SOURCES} ${NUMAV_HEADERS})
else()
    add_library(numav SHARED ${NUMAV_SOURCES} ${NUMAV_HEADERS})
endif()

target_include_directories(numav PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/third-party/safe-ptr/include
    /usr/include/eigen3
    ${MKLROOT}/include
)

target_compile_definitions(numav PUBLIC MKL_ILP64)
if(NOT GDB)
    target_compile_options(numav PRIVATE -Wall -Wextra -O3 -m64)
    target_link_options(numav PRIVATE -flto -m64)
else()
    target_compile_options(numav PRIVATE -Wall -Wextra -g -O0 -fno-inline -m64)
    target_link_options(numav PRIVATE -m64)
endif()

target_link_libraries(numav PRIVATE
    -Wl,--start-group
    ${MKLROOT}/lib/libmkl_intel_ilp64.a
    ${MKLROOT}/lib/libmkl_intel_thread.a
    ${MKLROOT}/lib/libmkl_core.a
    -Wl,--end-group
    ${MKL_COMPILER_ROOT}/lib/libiomp5.a
    -lpthread -lm -ldl
)
target_link_libraries(numav PRIVATE spdlog::spdlog_header_only)
