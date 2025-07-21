# Numav (work in progress)

Numav (Numerical Acoustics and Vibrations) is a C++ library with Julia bindings for acoustics and vibrations simulations.

## How to build libnumav (static)
```
rm -rf build &&
rm -rf install &&
cmake -B build \
    -D CMAKE_BUILD_TYPE=Release \
    -D CMAKE_INSTALL_PREFIX=./install &&
cmake --build build --parallel ${nproc} &&
cmake --install build
```

## How to build libnumav (dynamic)
```
rm -rf build &&
rm -rf install &&
cmake -B build \
    -D CMAKE_BUILD_TYPE=Release \
    -D DYNAMIC=TRUE \
    -D CMAKE_INSTALL_PREFIX=./install &&
cmake --build build --parallel ${nproc} &&
cmake --install build
```

## How to build and run tests
```
rm -rf build &&
rm -rf install &&
cmake -B build \
    -D CMAKE_BUILD_TYPE=Release \
    -D BUILD_TESTS=TRUE \
    -D CMAKE_INSTALL_PREFIX=./install &&
cmake --build build --parallel ${nproc} &&
cmake --install build &&
./build/tests_bin/test1
```

## How to build examples
```
rm -rf build &&
rm -rf install &&
cmake -B build \
    -D CMAKE_BUILD_TYPE=Release \
    -D BUILD_EXAMPLES=TRUE \
    -D CMAKE_INSTALL_PREFIX=./install &&
cmake --build build --parallel ${nproc} &&
cmake --install build
```

## How to build libnumav_jl
```
rm -rf build
rm -rf install
cmake -B build \
    -D CMAKE_BUILD_TYPE=Release \
    -D BIND_JULIA=TRUE \
    -D CMAKE_INSTALL_PREFIX=./install \
    -D CMAKE_PREFIX_PATH=/root/.julia/artifacts/68b3990305fd468fa06214ccfafd50120554feca/
cmake --build build --parallel ${nproc}
cmake --install build
```

## How to build the JLL
```
julia +1.7 --color=yes ./bindings/julia/build_tarballs.jl --deploy-jll=local --verbose
```

## How to generate the libnumav_jl_jll override
```
julia ./bindings/julia/generate_override.jl
```
Before using Numav.jl locally, build it first.

## TODO
