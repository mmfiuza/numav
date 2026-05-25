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

sources = [ DirectorySource("/tmp/binary_builder") ]

script = raw"""
    if [[
        ${target} == x86_64-linux-gnu* ||
        ${target} == x86_64-w64* ||
        ${target} == x86_64-apple*
    ]]; then
        SOLVER=ONEMKL
    else
        SOLVER=EIGEN
    fi

    mkdir build && cd build
    cmake \
        -D CMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
        -D CMAKE_BUILD_TYPE=Release \
        -D SOLVER=${SOLVER} \
        -D BIND_JULIA=TRUE \
        -D Julia_PREFIX=${prefix} \
        -D JlCxx_DIR=${prefix}/lib/cmake/JlCxx \
        -D CMAKE_FIND_ROOT_PATH=${prefix} \
        -D CMAKE_INSTALL_PREFIX=${prefix} \
        ..
    cmake --build . --config Release --target install -- -j${nproc}
"""

julia_version = "1.11"
jl_ver_num = VersionNumber(julia_version)

platforms = [
    Platform("x86_64" , "linux"  ; libc=:glibc, julia_version=jl_ver_num),
    Platform("aarch64", "linux"  ; libc=:glibc, julia_version=jl_ver_num),
    Platform("x86_64" , "linux"  ; libc=:musl , julia_version=jl_ver_num),
    Platform("aarch64", "linux"  ; libc=:musl , julia_version=jl_ver_num),
    Platform("x86_64" , "windows";              julia_version=jl_ver_num),
    Platform("x86_64" , "macos"  ;              julia_version=jl_ver_num),
    Platform("aarch64", "macos"  ;              julia_version=jl_ver_num),
]
platforms = expand_cxxstring_abis(platforms)

products = [ LibraryProduct("libnumav_jl", :libnumav_jl) ]

mkl_linux_windows = filter(p ->
    (arch(p) == "x86_64" && Sys.islinux(p) && libc(p) == "glibc") ||
    (arch(p) == "x86_64" && Sys.iswindows(p)),
    platforms
)
mkl_macos = filter(p -> arch(p) == "x86_64" && Sys.isapple(p), platforms)

dependencies = [
    Dependency("libcxxwrap_julia_jll", compat="0.14.9"),
    BuildDependency(PackageSpec(;name="libjulia_jll", version="1.11.1")),

    # oneMKL for x86_64 Linux(glibc) and Windows
    Dependency("MKL_jll", compat="=2025.2.0"; platforms=mkl_linux_windows),
    BuildDependency(
        PackageSpec(;name="MKL_Headers_jll", version="2025.2.0");
        platforms=mkl_linux_windows
    ),

    # oneMKL for x86_64 MacOS is 2023.2.0 (last supported version)
    Dependency("MKL_jll", compat="=2023.2.0"; platforms=mkl_macos),
    BuildDependency(
        PackageSpec(;name="MKL_Headers_jll", version="2023.2.0");
        platforms=mkl_macos
    ),
]

build_tarballs(
    ARGS, name, version, sources, script, platforms, products, dependencies;
    preferred_gcc_version=v"12", julia_compat="~"*julia_version
)
