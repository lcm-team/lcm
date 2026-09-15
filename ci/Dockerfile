# Reproducible toolchain for the LCM CMake/CTest gate.
FROM debian:bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive \
    OMPI_MCA_rmaps_base_oversubscribe=1

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        ca-certificates \
        cmake \
        gfortran \
        libopenmpi-dev \
        ninja-build \
        openmpi-bin \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
