# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

using BinaryBuilder, Pkg

run(
    `sh -c "
        rm -rf build &&
        rm -rf install &&
        rm -rf products &&
        rm -rf /tmp/binary_builder &&
        cp -r . /tmp/binary_builder
    "`
)

name = "libnumav_jl"
version = v"0.1.0"

sources = [
    DirectorySource("/tmp/binary_builder")
]

script = raw"""
    mkdir build && cd build
    cmake \
        -D CMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
        -D CMAKE_BUILD_TYPE=Release \
        -D SOLVER=EIGEN \
        -D BIND_JULIA=TRUE \
        -D Julia_PREFIX=${prefix} \
        -D JlCxx_DIR=${prefix}/lib/cmake/JlCxx \
        -D CMAKE_FIND_ROOT_PATH=${prefix} \
        -D CMAKE_INSTALL_PREFIX=${prefix} \
        ..
    cmake --build . --config Release --target install -- -j${nproc}
"""

platforms = [Platform("x86_64", "linux"; libc=:glibc, julia_version=v"1.11")]
platforms = expand_cxxstring_abis(platforms)

products = [LibraryProduct("libnumav_jl", :libnumav_jl)]

dependencies = [
    Dependency("libcxxwrap_julia_jll", compat="0.14.9"),
    BuildDependency(PackageSpec(;name="libjulia_jll", version="1.11.1"))
]

build_tarballs(
    ARGS, name, version, sources, script, platforms, products, dependencies;
    preferred_gcc_version=v"12", julia_compat="~1.11"
)
