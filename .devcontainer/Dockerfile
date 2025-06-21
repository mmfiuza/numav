# Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

FROM debian:bookworm-slim

# install dependencies
RUN apt-get update && apt-get install -y \
    libc-dev git gcc g++ make cmake curl squashfs-tools libeigen3-dev

# install Julia
RUN curl -fsSL https://install.julialang.org | sh -s -- --yes
ENV PATH="/root/.juliaup/bin:$PATH"
RUN juliaup add 1.7

# install Julia packages
RUN julia +1.7 -e \
    'using Pkg; Pkg.add(["CxxWrap", "BinaryBuilder", "libjulia_jll", "libcxxwrap_julia_jll"])'

# set privileged mode environment variable
ENV BINARYBUILDER_RUNNER=privileged

# agree to the Xcode license terms
ENV BINARYBUILDER_AUTOMATIC_APPLE=true