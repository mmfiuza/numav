# Numav

Secret project!

# How to build
```
rm -rf build
rm -rf install
cmake -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=./install \
    -DCMAKE_PREFIX_PATH=./libcxxwrap-julia/install
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
Create script to generate overrides
Make the cmakelists either install libnumav_jl.so or libnumav.a
