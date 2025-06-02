using Pkg
Pkg.activate(".")

using BinaryBuilder

name = "Numav"
version = v"0.0.1"

run(`rm -rf build`)
run(`rm -rf /tmp/repo_copy`)
run(`cp -r . /tmp/repo_copy`)

sources = [
    DirectorySource("/tmp/repo_copy"),
]

script = raw"""

cd $WORKSPACE/srcdir

cmake -B build \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_BUILD_TYPE=Release

cmake --build build --parallel ${nproc}

cmake --install build

"""

platforms = Platform("x86_64", "linux", libc=:musl)
platforms = expand_cxxstring_abis(platforms)

products = [
    LibraryProduct("libnumav", :libnumav),
    FileProduct("share/numav/Numav.jl", :Numav_jl)
]

dependencies = [
    Dependency("libcxxwrap_julia_jll")
]

build_tarballs(
    ARGS, name, version, sources, script, platforms, products, dependencies;
    julia_compat="1.7",
    preferred_gcc_version=v"9"
)