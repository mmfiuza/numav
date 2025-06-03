project(numav LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
cmake_minimum_required(VERSION 3.5)

# Global settings
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/lib")
set(CMAKE_MODULE_PATH 
    ${CMAKE_CURRENT_SOURCE_DIR}/third-party/libcxxwrap-julia/install/lib/cmake/JlCxx)
set(CMAKE_PREFIX_PATH
    ${CMAKE_PREFIX_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/third-party/libcxxwrap-julia/install)

find_package(Julia REQUIRED)
find_package(JlCxx REQUIRED)

# Common compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-parameter -Wextra -Wreorder -fPIC")

# Build components
add_subdirectory(src) # Main library
option(BUILD_JULIA_BINDINGS "Build Julia bindings" ON)
if(BUILD_JULIA_BINDINGS)
    add_subdirectory(bindings/julia) # Julia bindings
endif()