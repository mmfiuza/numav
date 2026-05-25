# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

cmake_minimum_required(VERSION 3.5)
project(numav LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_TESTS "Build tests" OFF)
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_SHARED_LIBS "Build shared library" OFF)
option(GDB "Build with debug symbols for gdb" OFF)

# solver-specific configuration
set(SOLVER "ONEMKL" CACHE STRING "The linear system solver")
set_property(CACHE SOLVER PROPERTY STRINGS "INTERNAL" "EIGEN" "ONEMKL" "MUMPS")
if(SOLVER STREQUAL "INTERNAL")
    set(SELECTED_SOLVER_DEFINITION NUMAV_USE_INTERNAL_SOLVER)
    set(SOLVER_SPECIFIC_LINK_LIBRARIES ldlt_solver)
    add_subdirectory(third-party/ldlt-solver)
elseif(SOLVER STREQUAL "EIGEN")
    set(SELECTED_SOLVER_DEFINITION NUMAV_USE_EIGEN_SOLVER)
    # set(SOLVER_SPECIFIC_LINK_LIBRARIES -lblas)
    # set(SOLVER_SPECIFIC_DEFINITIONS EIGEN_USE_BLAS)
elseif(SOLVER STREQUAL "ONEMKL")
    set(SELECTED_SOLVER_DEFINITION NUMAV_USE_ONEMKL_SOLVER)
    set(MKL_LINK dynamic)
    find_package(MKL REQUIRED)
    set(SOLVER_SPECIFIC_LINK_LIBRARIES
        MKL::MKL
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
        ${CMAKE_CURRENT_SOURCE_DIR}/third-party/eigen3/include
        ${CMAKE_CURRENT_SOURCE_DIR}/third-party/spdlog/include
        ${CMAKE_CURRENT_SOURCE_DIR}/third-party/fmt/include
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
    ${SOLVER_SPECIFIC_LINK_LIBRARIES}
)

# set numav's compile and link flags
set(NUMAV_COMPILER_FLAGS -Wall -Wextra)
set(NUMAV_LINKER_FLAGS)
if(NOT GDB)
    set(NUMAV_COMPILER_FLAGS ${NUMAV_COMPILER_FLAGS} -O3)
    set(NUMAV_LINKER_FLAGS ${NUMAV_LINKER_FLAGS} -flto)
else()
    set(NUMAV_COMPILER_FLAGS ${NUMAV_COMPILER_FLAGS} -g -O0 -fno-inline)
endif()
if(BIND_JULIA)
    set(NUMAV_COMPILER_FLAGS ${NUMAV_COMPILER_FLAGS} -fPIC)
endif()
target_compile_options(numav PRIVATE ${NUMAV_COMPILER_FLAGS})
target_link_options(numav PRIVATE ${NUMAV_LINKER_FLAGS})

if(BIND_JULIA)
    find_package(JlCxx REQUIRED)
    get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
    get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
    
    add_library(numav_jl SHARED
        ${CMAKE_CURRENT_SOURCE_DIR}/julia-bindings/numav_jl.cpp
    )
    target_link_libraries(numav_jl
        numav
        JlCxx::cxxwrap_julia
        JlCxx::cxxwrap_julia_stl
    )
    install(TARGETS numav_jl LIBRARY DESTINATION lib)
endif()

# add subdirectories
if(BUILD_TESTS)
    add_subdirectory(tests)
endif()
