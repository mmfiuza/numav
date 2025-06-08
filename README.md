# Numav

Secret project!

## How to build libnumav
```
rm -rf build
rm -rf install
cmake -B build \
    -D CMAKE_BUILD_TYPE=Release \
    -D CMAKE_INSTALL_PREFIX=./install \
    -D CMAKE_PREFIX_PATH=/root/.julia/artifacts/68b3990305fd468fa06214ccfafd50120554feca/
cmake --build build --parallel ${nproc}
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

## TODO
