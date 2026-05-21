# Numav Solver

The Numav Solver is a C++ library to perform acoustics and vibrations simulations.

## How to start the development container
The development container can be started in two ways:
- directly in a linux terminal environment;
- through the Dev Containers VS Code extension.

### Linux terminal
To start the container directly in a Linux terminal environment, go to the root directory of this repository and run:
```
./docker/start-container.sh <path-to-your-github-ssh-keys>
```
`<path-to-your-github-ssh-keys>` is an optional argument. Its default value is `~/.ssh`.

### VS Code
To start the container through the Dev Containers VS Code extension, install it and open the root folder of this repository. Then, click the remote indicator in the bottom-left corner and select **Reopen in Container**.

## How to build libnumav (static)
Run:
```
rm -rf build &&
cmake -B build -D CMAKE_BUILD_TYPE=Release &&
cmake --build build --parallel ${nproc}
```

## How to build libnumav (dynamic)
Run:
```
rm -rf build &&
cmake -B build -D CMAKE_BUILD_TYPE=Release -D BUILD_SHARED_LIBS=TRUE &&
cmake --build build --parallel ${nproc}
```

## How to build and run tests
Run:
```
rm -rf build &&
cmake -B build -D CMAKE_BUILD_TYPE=Release -D BUILD_TESTS=TRUE &&
cmake --build build --parallel ${nproc} &&
./build/tests_bin/test1
```

## How to link a code against libnumav.a (Static Numav)
Run:
```
g++ -o my_simulation my_simulation.cpp -Iinclude -L./build/lib -L/opt/intel/oneapi/mkl/2025.2/lib -Wl,--start-group -l:libmkl_core.a -l:libmkl_gf_ilp64.a -l:libmkl_gnu_thread.a -Wl,--end-group -lgomp -lpthread -lm -ldl -lnumav -m64 -flto
```

## How to link a code against libnumav.so (Dynamic Numav)
Run:
```
g++ -o my_simulation my_simulation.cpp -I include -L ./build/lib -l numav -Wl,-rpath,./build/lib
```

## How to build libnumav_jl
```
rm -rf build &&
rm -rf julia_bindings/override/lib &&
cmake -B build \
    -D CMAKE_BUILD_TYPE=Release \
    -D SOLVER=EIGEN \
    -D BIND_JULIA=TRUE \
    -D CMAKE_INSTALL_PREFIX=julia_bindings/override \
    -D CMAKE_PREFIX_PATH=/usr/local/share/julia/artifacts/662f181f562225e139306f1e6e383c70bc9255f9 &&
cmake --build build --parallel ${nproc} &&
cmake --install ./build
```

## How to build the JLL
```
julia +1.7.3 /workspace/julia_bindings/build_tarballs.jl --deploy-jll=local

```

## How to generate the libnumav_jl_jll override
```
rm -rf /workspace/julia_bindings/override
mkdir /workspace/julia_bindings/override
tar -xzvf /workspace/products/libnumav_jl.v0.1.0.x86_64-linux-gnu-cxx11-julia_version+1.11.0.tar.gz -C /workspace/julia_bindings/override
julia /workspace/julia_bindings/generate_override.jl
```

## Dev the Julia packages
```
] dev /usr/local/share/julia/dev/libnumav_jl_jll
] dev /workspace/julia_bindings/Numav
```
