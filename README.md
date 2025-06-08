# Numav

Secret project!

# How to build libnumav
```
rm -rf build
rm -rf install
cmake -B build \
    -D CMAKE_BUILD_TYPE=Release \
    -D CMAKE_INSTALL_PREFIX=./install \
    -D CMAKE_PREFIX_PATH=./libcxxwrap-julia/install
cmake --build build --parallel ${nproc}
cmake --install build
```

# How to build libnumav_jl
```
rm -rf build
rm -rf install
cmake -B build \
    -D CMAKE_BUILD_TYPE=Release \
    -D BIND_JULIA=TRUE \
    -D CMAKE_INSTALL_PREFIX=./install \
    -D CMAKE_PREFIX_PATH=./libcxxwrap-julia/install
cmake --build build --parallel ${nproc}
cmake --install build
```

# How to build libcxxwrap-julia for local testing
```
git clone https://github.com/JuliaInterop/libcxxwrap-julia.git
cd libcxxwrap-julia
rm -rf build
rm -rf install
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../install ..
cmake --build . --config Release -j $(nproc)
cmake --install .
cd ../..
```

# How to build the JLL
```
juliaup default 1.7
julia --project --color=yes ./bindings/julia/build_tarballs.jl --deploy-jll=local --verbose
juliaup default release
```

# TODO
get libcxxwrap-julia from libcxxwrap_julia_jll intead of building it