# Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

set(EXAMPLES_BIN_DIR "${CMAKE_BINARY_DIR}/examples_bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXAMPLES_BIN_DIR})

add_executable(example1 example1.cpp)
target_link_directories(example1 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../build)
target_link_libraries(example1 PRIVATE numav)
target_compile_options(example1
    PRIVATE -g -O0 -fno-inline
)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh1.bdf
    ${EXAMPLES_BIN_DIR}/mesh1.bdf
    COPYONLY
)
