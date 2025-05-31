FROM julia:1.7.3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libatomic1 \
    python3 \
    gfortran \
    perl \
    wget \
    m4 \
    cmake \
    pkg-config \
    git

# Initialize Julia package manager and install packages
RUN julia -e 'using Pkg; Pkg.add("CxxWrap")' && \
    julia -e 'using Pkg; Pkg.add("BinaryBuilder")'