# Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# start from debian:bookworm-slim
FROM debian@sha256:df52e55e3361a81ac1bead266f3373ee55d29aa50cf0975d440c2be3483d8ed3

# install dependencies
RUN apt-get update && apt-get install -y \
    libc6-dev=2.36-9+deb12u13 \
    git=1:2.39.5-0+deb12u2 \
    gcc=4:12.2.0-3 \
    g++=4:12.2.0-3 \
    gdb=13.1-3 \
    make=4.3-4.1 \
    cmake=3.25.1-1 \
    curl=7.88.1-10+deb12u14 \
    squashfs-tools=1:4.5.1-1 \
    libeigen3-dev=3.4.0-4 \
    libspdlog-dev=1:1.10.0+ds-0.4 \
    wget=1.21.3-1+deb12u1

# install intel oneMKL
RUN wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/47c7d946-fca1-441a-b0df-b094e3f045ea/intel-onemkl-2025.2.0.629_offline.sh
RUN sh ./intel-onemkl-2025.2.0.629_offline.sh -a --silent --eula accept
RUN rm intel-onemkl-2025.2.0.629_offline.sh

# install Julia
RUN curl -fsSL https://install.julialang.org | sh -s -- --yes
ENV PATH="/root/.juliaup/bin:$PATH"
RUN juliaup add 1.11
RUN juliaup add 1.7
RUN juliaup default 1.11

# install Julia packages
RUN julia +1.7 -e 'using Pkg; Pkg.add([ \
    PackageSpec(name="CxxWrap", version="0.17"), \
    PackageSpec(name="BinaryBuilder", version="0.6"), \
    PackageSpec(name="libjulia_jll", version="1.10"), \
    PackageSpec(name="libcxxwrap_julia_jll", version="0.14") \
])'

# set privileged mode environment variable
ENV BINARYBUILDER_RUNNER=privileged

# agree to the Xcode license terms
ENV BINARYBUILDER_AUTOMATIC_APPLE=true
