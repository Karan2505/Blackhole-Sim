# =============================================================================
# BlackHole-Sim Production Docker Image
# 
# Multi-stage build for minimal production image with CUDA support
# =============================================================================

# -----------------------------------------------------------------------------
# Stage 1: Build Environment
# -----------------------------------------------------------------------------
FROM nvidia/cuda:12.2.0-devel-ubuntu22.04 AS builder

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libhdf5-dev \
    libopenmpi-dev \
    openmpi-bin \
    python3 \
    python3-dev \
    python3-pip \
    python3-numpy \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install pybind11
RUN pip3 install pybind11[global]

# Set working directory
WORKDIR /build

# Copy source files
COPY CMakeLists.txt .
COPY core/ core/
COPY cuda/ cuda/
COPY python/ python/
COPY examples/ examples/

# Build
RUN mkdir -p build && cd build && \
    cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DENABLE_CUDA=ON \
        -DENABLE_MPI=ON \
        -DENABLE_OPENMP=ON \
        -DENABLE_PYTHON=ON \
        -DENABLE_TESTS=OFF \
        -DCMAKE_INSTALL_PREFIX=/install && \
    make -j$(nproc) && \
    make install

# -----------------------------------------------------------------------------
# Stage 2: Runtime Environment
# -----------------------------------------------------------------------------
FROM nvidia/cuda:12.2.0-runtime-ubuntu22.04 AS runtime

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libhdf5-103 \
    libopenmpi3 \
    python3 \
    python3-pip \
    python3-numpy \
    python3-matplotlib \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install --no-cache-dir \
    numpy \
    matplotlib \
    h5py \
    scipy \
    jupyter

# Copy built artifacts
COPY --from=builder /install /usr/local
COPY --from=builder /build/build/python/_bhsim_core*.so /usr/local/lib/python3/dist-packages/bhsim/

# Copy Python package
COPY python/bhsim /usr/local/lib/python3/dist-packages/bhsim

# Copy example notebooks and scripts
COPY examples/ /app/examples/
COPY python/notebooks/ /app/notebooks/

# Set library path
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
ENV PYTHONPATH=/usr/local/lib/python3/dist-packages:$PYTHONPATH

WORKDIR /app

# Default command: start Jupyter
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]

# Expose Jupyter port
EXPOSE 8888

# -----------------------------------------------------------------------------
# Stage 3: Web Demo Image
# -----------------------------------------------------------------------------
FROM nginx:alpine AS webdemo

# Copy static web files
COPY index.html /usr/share/nginx/html/
COPY docs/ /usr/share/nginx/html/docs/

# Custom nginx config for SPA
RUN echo 'server { \
    listen 80; \
    root /usr/share/nginx/html; \
    index index.html; \
    location / { \
        try_files $uri $uri/ /index.html; \
    } \
}' > /etc/nginx/conf.d/default.conf

EXPOSE 80

CMD ["nginx", "-g", "daemon off;"]