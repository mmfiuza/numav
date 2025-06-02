# Numav

Secret project!

# How to build
```
cmake -B build -DCMAKE_INSTALL_PREFIX=./install -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel ${nproc}
cmake --install build
```

git clone https://github.com/JuliaInterop/libcxxwrap-julia.git
cd libcxxwrap-julia
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../../third-party/jlcxx ..
cmake --build . --config Release -j $(nproc)
cmake --install .
cd ../../
rm -rf libcxxwrap-julia