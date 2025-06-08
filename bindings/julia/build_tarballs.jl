# Copyright (c) 2025 Matheus Machado Fiuza <matheus.fiuza@eac.ufsm.br>

using BinaryBuilder

run(`rm -rf build`)
run(`rm -rf /tmp/numav`)
run(`cp -r . /tmp/numav`)

name = "libnumav_jl"
version = v"0.0.1"

sources = [
    DirectorySource("/tmp/numav"),
]

script = raw"""

    mkdir build && cd build

    cmake \
        -D CMAKE_BUILD_TYPE=Release \
        -D BIND_JULIA=TRUE \
        -D Julia_PREFIX=${prefix} \
        -D CMAKE_FIND_ROOT_PATH=${prefix} \
        -D JlCxx_DIR=${prefix}/lib/cmake/JlCxx \
        -D CMAKE_INSTALL_PREFIX=${prefix} \
        -D CMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
        ..

    cmake --build . --config Release --target install -- -j${nproc}

"""

platforms = [
    Platform("x86_64", "linux")
]
platforms = expand_cxxstring_abis(platforms)

products = [
    LibraryProduct("libnumav_jl", :libnumav_jl),
]

dependencies = [
    BuildDependency("libjulia_jll"),
    Dependency("libcxxwrap_julia_jll")
]

build_tarballs(
    ARGS, name, version, sources, script, platforms, products,
    dependencies; preferred_gcc_version=v"12"
)
