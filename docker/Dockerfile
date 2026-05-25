# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# start from debian:bookworm-slim
ARG DEB_DIGEST=df52e55e3361a81ac1bead266f3373ee55d29aa50cf0975d440c2be3483d8ed3
FROM debian@sha256:${DEB_DIGEST}

# install apt dependencies
RUN apt-get update && apt-get install -y \
    libc6-dev=2.36-9+deb12u13 \
    git=1:2.39.5-0+deb12u2 \
    gcc=4:12.2.0-3 \
    g++=4:12.2.0-3 \
    gdb=13.1-3 \
    clang=1:14.0-55.7~deb12u1 \
    make=4.3-4.1 \
    cmake=3.25.1-1 \
    curl=7.88.1-10+deb12u14 \
    squashfs-tools=1:4.5.1-1 \
    libspdlog-dev=1:1.10.0+ds-0.4 \
    wget=1.21.3-1+deb12u1 \
    gosu=1.14-1+b10 \
    tree=2.1.0-1 \
    sudo=1.9.13p3-1+deb12u3 \
    libblas-dev=3.11.0-2 \
    liblapack-dev=3.11.0-2 \
    libmumps-seq-dev=5.5.1-1 \
    patchelf=0.14.3-1+b1

# install Eigen
ARG EIGEN_VERSION="5.0.1"
RUN wget https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/\
eigen-${EIGEN_VERSION}.tar.gz
RUN tar -xzf eigen-${EIGEN_VERSION}.tar.gz
RUN rm eigen-${EIGEN_VERSION}.tar.gz
RUN cd eigen-${EIGEN_VERSION} && \
    mkdir build && \
    cd build && \
    cmake \
        -D CMAKE_INSTALL_PREFIX=/opt \
        -D EIGEN_BUILD_BLAS=OFF \
        -D EIGEN_BUILD_LAPACK=OFF .. && \
    cmake --install . && \
    cd ../.. && \
    rm -rf eigen-${EIGEN_VERSION}

# install Intel oneMKL
RUN wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/\
47c7d946-fca1-441a-b0df-b094e3f045ea/intel-onemkl-2025.2.0.629_offline.sh
RUN sh ./intel-onemkl-2025.2.0.629_offline.sh -a --silent --eula accept
RUN rm intel-onemkl-2025.2.0.629_offline.sh

# install Julia
ARG JULIA_VERSION="1.11.9"
ARG BINARY_BUILDER_JULIA_VERSION="1.7.3"
ARG JULIAUP_VERSION="1.20.1"
ENV JULIA_DEPOT_PATH=/usr/local/share/julia
ENV JULIAUP_DEPOT_PATH=/usr/local/share
RUN case "$(uname -m)" in \
        x86_64)    JULIAUP_ARCH="x86_64" ;; \
        aarch64)   JULIAUP_ARCH="aarch64" ;; \
        *)         echo "Unsupported arch: $(uname -m)" && exit 1 ;; \
    esac && \
curl -fsSL "https://github.com/JuliaLang/juliaup/releases/download/\
v${JULIAUP_VERSION}/juliaup-${JULIAUP_VERSION}-${JULIAUP_ARCH}\
-unknown-linux-musl-portable.tar.gz" | tar -xz -C /usr/local/bin
RUN juliaup add ${JULIA_VERSION}
RUN juliaup add ${BINARY_BUILDER_JULIA_VERSION}
RUN juliaup default ${JULIA_VERSION}
RUN chmod -R 777 ${JULIAUP_DEPOT_PATH}/juliaup

# install Julia packages
RUN julia +${JULIA_VERSION} -e 'using Pkg; Pkg.add([ \
    PackageSpec(name="CxxWrap", version="0.17.5"), \
    PackageSpec(name="libcxxwrap_julia_jll", version="0.14.9"), \
    PackageSpec(name="Plots"), \
    PackageSpec(name="Symbolics"), \
    PackageSpec(name="Revise") \
])'
RUN julia +${BINARY_BUILDER_JULIA_VERSION} -e 'using Pkg; Pkg.add([ \
    PackageSpec(name="BinaryBuilder", version="0.6.6") \
])'
RUN chmod -R 777 ${JULIA_DEPOT_PATH}

# set privileged mode environment variable
ENV BINARYBUILDER_RUNNER=privileged

# agree to the Xcode license terms
ENV BINARYBUILDER_AUTOMATIC_APPLE=true

# set the entrypoint
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]
