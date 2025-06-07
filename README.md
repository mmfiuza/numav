# Numav

Secret project!

# How to build
```
rm -rf build
cmake -B build -DCMAKE_INSTALL_PREFIX=./install -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel ${nproc}
cmake --install build
```

# How to build the JLL
```
julia --project --color=yes ./bindings/julia/build_tarballs.jl --deploy-jll=local --verbose
```
