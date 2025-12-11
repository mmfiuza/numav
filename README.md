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
cmake --build build --parallel ${nproc}
./build/tests_bin/test1
```

## How to build examples
Run:
```
rm -rf build &&
cmake -B build -D CMAKE_BUILD_TYPE=Release -D BUILD_EXAMPLES=TRUE &&
cmake --build build --parallel ${nproc}
```

## How to link a code against libnumav.so (Static Numav)
Run:
```
g++ -o my_simulation my_simulation.cpp -Iinclude -L./build/lib -L/opt/intel/oneapi/mkl/2025.2/lib -Wl,--start-group -l:libmkl_core.a -l:libmkl_gf_ilp64.a -l:libmkl_gnu_thread.a -Wl,--end-group -lgomp -lpthread -lm -ldl -lnumav -m64 -flto
```

## How to link a code against libnumav.so (Dynamic Numav)
Run:
```
g++ -o my_simulation my_simulation.cpp -I include -L ./build/lib -l numav -Wl,-rpath,./build/lib
```
