# Numav (work in progress)

Numav (Numerical Acoustics and Vibrations) is a C++ library for acoustics and vibrations simulations.

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
