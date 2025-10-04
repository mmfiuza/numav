# Numav Solver

The Numav Solver is a C++ library to perform acoustics and vibrations simulations.

## How to start the development container
Run:
```
./start_container.sh <path-to-github-ssh-keys>
```
`<path-to-github-ssh-keys>` is an optional argument. Its default value is `~/.ssh`.

## How to build libnumav (static)
Run:
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
Run:
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
Run:
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
Run:
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
