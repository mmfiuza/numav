cmake_minimum_required(VERSION 3.5)
project(numav LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

# Global settings
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}/lib")

find_package(Julia REQUIRED)
find_package(JlCxx REQUIRED)

# Common compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-parameter -Wextra -Wreorder -fPIC")

add_definitions(-DJULIA_ENABLE_THREADING)

# Build components
add_subdirectory(src) # Main library
option(BUILD_JULIA_BINDINGS "Build Julia bindings" ON)
if(BUILD_JULIA_BINDINGS)
    add_subdirectory(third-party/cxxwrap) # CxxWrap dependency
    add_subdirectory(bindings/julia) # Julia bindings
endif()