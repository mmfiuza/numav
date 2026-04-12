# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

set(TESTS_BIN_DIR "${CMAKE_BINARY_DIR}/tests_bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TESTS_BIN_DIR})

add_executable(test1 test1.cpp)
target_link_directories(test1 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../build)
target_link_libraries(test1 PRIVATE numav)
target_compile_options(test1
    PRIVATE -g -O0 -fno-inline
)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/test1.bdf
    ${TESTS_BIN_DIR}/test1.bdf
    COPYONLY
)
