# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

project(ldlt_solver LANGUAGES CXX)

add_library(ldlt_solver STATIC
    src/ldlt-solver.cpp
)

target_include_directories(ldlt_solver
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_options(ldlt_solver PRIVATE -Wall -Wextra -Ofast -ffast-math)
