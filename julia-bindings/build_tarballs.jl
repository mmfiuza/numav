# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

using BinaryBuilder, Pkg

name = "numav_julia"
version = v"0.1.0"

sources = [ 
    GitSource(
        "https://github.com/mmfiuza/numav.git", 
        "09d7c407a396ebf7018b2b7f757ff6647abe57b8"
    )
]

script = raw"""
    # oneMKL is picked for the available platforms, otherwise Eigen is picked
    if [[
        ${target} == x86_64-linux-gnu* ||
        ${target} == x86_64-w64* ||
        ${target} == x86_64-apple*
    ]]; then
        SOLVER=ONEMKL
    else
        SOLVER=EIGEN
    fi
    
    cd ${WORKSPACE}/srcdir/numav && mkdir build && cd build
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

julia_ver = "1.11" # Julia version that the user can use numav_julia_jll
jl_ver_num = VersionNumber(julia_ver)

platforms = [
    Platform("x86_64" , "linux"  ; libc=:glibc, julia_version=jl_ver_num),
    # Platform("aarch64", "linux"  ; libc=:glibc, julia_version=jl_ver_num),
    # Platform("x86_64" , "linux"  ; libc=:musl , julia_version=jl_ver_num),
    # Platform("aarch64", "linux"  ; libc=:musl , julia_version=jl_ver_num),
    # Platform("x86_64" , "windows";              julia_version=jl_ver_num),
    # Platform("x86_64" , "macos"  ;              julia_version=jl_ver_num),
    # Platform("aarch64", "macos"  ;              julia_version=jl_ver_num),
]
platforms = expand_cxxstring_abis(platforms)

products = [ LibraryProduct("numav_julia", :numav_julia) ]

mkl_linux_windows = filter(p ->
    (arch(p) == "x86_64" && Sys.islinux(p) && libc(p) == "glibc") ||
    (arch(p) == "x86_64" && Sys.iswindows(p)),
    platforms
)
mkl_macos = filter(p -> arch(p) == "x86_64" && Sys.isapple(p), platforms)

dependencies = [
    Dependency("libcxxwrap_julia_jll", compat="0.14.9"),
    BuildDependency(PackageSpec(;name="libjulia_jll", version=julia_ver*".1")),
    BuildDependency(PackageSpec(;name="Eigen_jll", version="5.0.1")),
    BuildDependency(PackageSpec(;name="spdlog_jll", version="1.15.0")),

    # oneMKL for x86_64 Linux(glibc) and Windows
    Dependency("MKL_jll", compat="=2025.2.0"; platforms=mkl_linux_windows),
    BuildDependency(
        PackageSpec(;name="MKL_Headers_jll", version="2025.2.0");
        platforms=mkl_linux_windows
    ),

    # oneMKL for x86_64 MacOS is 2023.2 (last supported version for MacOS)
    Dependency("MKL_jll", compat="=2023.2.0"; platforms=mkl_macos),
    BuildDependency(
        PackageSpec(;name="MKL_Headers_jll", version="2023.2.0");
        platforms=mkl_macos
    ),
]

build_tarballs(
    ARGS, name, version, sources, script, platforms, products, dependencies;
    preferred_gcc_version=v"10", julia_compat="~"*julia_ver
)
