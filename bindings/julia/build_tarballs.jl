using BinaryBuilder, Pkg

name = "libnumav"
version = v"0.0.1"

sources = [
    DirectorySource("src")
]

script = raw"""
cd ${WORKSPACE}/srcdir/src

# Build with proper dependency paths
cmake -B build \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_PREFIX_PATH=${prefix} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON

cmake --build build --config Release -- -j${nproc}

# Install library to standard location
cmake --install build --component Runtime

# Install Julia bindings to proper location
mkdir -p ${prefix}/share/numav
install -Dvm 644 bindings/julia/Numav.jl ${prefix}/share/numav/Numav.jl
"""

platforms = Platform("x86_64", "linux", libc=:musl)
platforms = expand_cxxstring_abis(platforms)

# Corrected products section
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