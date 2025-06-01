# Numav

Secret project!

# How to build
```
cmake -S . -B build
cmake --build build --config Release -j $(nproc)
```

git clone https://github.com/JuliaInterop/libcxxwrap-julia.git
cd libcxxwrap-julia
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../../third-party/jlcxx ..
cmake --build . --config Release -j $(nproc)
cmake --install .
cd ../../
rm -rf libcxxwrap-julia