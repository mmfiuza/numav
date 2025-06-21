# Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

project(numav LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
cmake_minimum_required(VERSION 3.5)

# find Eigen
find_package (Eigen3 3.4 REQUIRED NO_MODULE)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/lib")

# add numav
set(NUMAV_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/numav.cpp)
set(NUMAV_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/numav.hpp)
add_library(numav STATIC ${NUMAV_SOURCES} ${NUMAV_HEADERS})
target_include_directories(numav PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include /usr/include/eigen3)

# set install path
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

if(BIND_JULIA)
    
    # find cxxwrap_julia
    find_package(JlCxx REQUIRED)
    get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
    get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
    
    # add numav_jl
    add_library(numav_jl SHARED ${CMAKE_CURRENT_SOURCE_DIR}/bindings/julia/numav_jl.cpp)

    # define link targets of numav_jl
    target_link_libraries(numav_jl
        numav
        JlCxx::cxxwrap_julia
        JlCxx::cxxwrap_julia_stl
    )

    # install numav_jl
    install(TARGETS numav_jl
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
    )

else() # default c++ build

    # install numav
    install(TARGETS numav
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
    )

endif()
