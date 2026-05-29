# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# start from debian:bookworm-slim
ARG DEB_DIGEST=2477d9ee0ead4370c778ce3aa42258a0b07684d1a84ded8f4af518383fbc3f2d
FROM debian@sha256:${DEB_DIGEST}

# install apt dependencies
RUN apt-get update && apt-get install -y \
    libc6-dev \
    git \
    gcc \
    g++ \
    gdb \
    clang \
    make \
    cmake \
    curl \
    squashfs-tools \
    libspdlog-dev \
    wget \
    gosu \
    tree \
    sudo \
    libblas-dev \
    liblapack-dev \
    libmumps-seq-dev \
    patchelf

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
    PackageSpec(name="Revise"), \
    PackageSpec(name="HDF5") \
])'
RUN julia +${BINARY_BUILDER_JULIA_VERSION} -e 'using Pkg; Pkg.add([ \
    PackageSpec(name="BinaryBuilder", version="0.6.6") \
])'
RUN chmod -R 777 ${JULIA_DEPOT_PATH}

# set privileged mode environment variable
ENV BINARYBUILDER_RUNNER=privileged

# agree to the Xcode license terms
ENV BINARYBUILDER_AUTOMATIC_APPLE=true

# install HDF5
RUN wget https://github.com/HDFGroup/hdf5/releases/\
download/2.1.1/hdf5-2.1.1-ubuntu-2404_gcc.deb
RUN apt install ./hdf5-2.1.1-ubuntu-2404_gcc.deb
RUN rm hdf5-2.1.1-ubuntu-2404_gcc.deb

# create devuser
ARG DEFAULT_UID=1000
ARG DEFAULT_GID=1000
RUN groupadd -g ${DEFAULT_GID} devgroup && \
    useradd -m -u ${DEFAULT_UID} -g devgroup -G sudo -s /bin/bash devuser && \
    echo "devuser ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/devuser && \
    chmod 440 /etc/sudoers.d/devuser

# pre-create SSH agent socket target so bind-mount works
RUN touch /ssh-agent

# set the entrypoint
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]
