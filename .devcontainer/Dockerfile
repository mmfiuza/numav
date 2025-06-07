FROM debian:bookworm-slim

# Install dependencies
RUN apt-get update && apt-get install -y \
    libc-dev git gcc g++ make cmake curl squashfs-tools

# Install Julia
RUN curl -fsSL https://install.julialang.org | sh -s -- --yes
ENV PATH="/root/.juliaup/bin:$PATH"
RUN juliaup add 1.7 && juliaup default 1.7

# Install Julia packages
RUN julia -e \
'using Pkg; Pkg.add(["CxxWrap", "BinaryBuilder", "libjulia_jll", "libcxxwrap_julia_jll"])'

# Set privileged mode environment variable
ENV BINARYBUILDER_RUNNER=privileged