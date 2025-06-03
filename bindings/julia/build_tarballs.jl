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
    # Override compiler ID to silence the horrible "No features found" cmake error
    if [[ $target == *"apple-darwin"* ]]; then
        macos_extra_flags="-DCMAKE_CXX_COMPILER_ID=AppleClang -DCMAKE_CXX_COMPILER_VERSION=10.0.0 -DCMAKE_CXX_STANDARD_COMPUTED_DEFAULT=11"
    fi

    Julia_PREFIX=$prefix

    mkdir build && cd build
    cmake \
        -DJulia_PREFIX=$Julia_PREFIX \
        -DCMAKE_FIND_ROOT_PATH=$prefix \
        -DJlCxx_DIR=$prefix/lib/cmake/JlCxx \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
        $macos_extra_flags \
        -DCMAKE_BUILD_TYPE=Release ..
    VERBOSE=ON cmake --build . --config Release --target install -- -j${nproc}

"""

platforms = [
    Linux(:armv7l; libc=:glibc, compiler_abi=CompilerABI(cxxstring_abi=:cxx11)),
    Linux(:x86_64; libc=:glibc, compiler_abi=CompilerABI(cxxstring_abi=:cxx11)),
    MacOS(:x86_64; compiler_abi=CompilerABI(cxxstring_abi=:cxx11)),
    Windows(:x86_64; compiler_abi=CompilerABI(cxxstring_abi=:cxx11)),
]
platforms = expand_cxxstring_abis(platforms)

products = [
    LibraryProduct("libnumav", :libnumav),
]

# Dependencies that must be installed before this package can be built
dependencies = [
    Dependency("libcxxwrap_julia_jll"),
    BuildDependency(PackageSpec(name="Julia_jll", version=v"1.4.2"))
]

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies; preferred_gcc_version = v"7")
