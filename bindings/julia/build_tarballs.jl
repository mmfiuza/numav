using BinaryBuilder, Pkg

name = "Numav"
version = v"0.0.1"

run(`rm -rf build`)
run(`rm -rf /tmp/numav`)
run(`cp -r . /tmp/numav`)

sources = [
    DirectorySource("/tmp/numav"),
]

script = raw"""
ls
    Julia_PREFIX=$prefix

    mkdir build && cd build
    cmake \
        -DJulia_PREFIX=$Julia_PREFIX \
        -DCMAKE_FIND_ROOT_PATH=$prefix \
        -DJlCxx_DIR=$prefix/lib/cmake/JlCxx \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
        -DCMAKE_CXX_STANDARD=17 \
        $macos_extra_flags \
        -DCMAKE_BUILD_TYPE=Release ..
    VERBOSE=ON cmake --build . --config Release --target install -- -j${nproc}

"""

platforms = [Platform("x86_64", "linux")]
# platforms = expand_cxxstring_abis(platforms)

products = [
    LibraryProduct("libnumav", :libnumav),
]

# Dependencies that must be installed before this package can be built
dependencies = [
    BuildDependency("libjulia_jll"),
    Dependency("libcxxwrap_julia_jll")
]

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(
    ARGS, name, version, sources, script, platforms, products,
    dependencies; preferred_gcc_version = v"12"
)
