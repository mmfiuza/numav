# Numav

Secret project!

# How to build
```
rm -rf build
cmake -B build -DCMAKE_INSTALL_PREFIX=./install -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel ${nproc}
cmake --install build
```

# How to build libcxxwrap-julia
```
cd third-party/libcxxwrap-julia
rm -rf build
rm -rf install
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../install ..
cmake --build . --config Release -j $(nproc)
cmake --install .
cd ../../../
```