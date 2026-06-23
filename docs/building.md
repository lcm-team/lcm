# Building LCM

This repository is being moved toward a standalone LCM library. The current
CMake build is an early compilation target, not a complete scientific driver.

## Requirements

- CMake 3.20 or newer
- A Fortran compiler
- MPI Fortran compiler support

On macOS with Homebrew, this usually means `gcc`, `open-mpi`, and `cmake`.

## Default Build

Configure and build:

```sh
cmake -S . -B build
cmake --build build
```

The default target builds a static `lcm` library from the core files that can
currently compile outside SAM with temporary compatibility stubs.

## Temporary Compatibility Stubs

The `stubs/sam_compat.f90` file provides placeholder SAM modules such as
`grid`, `params`, `vars`, timers, and output hooks. These stubs are only a
compilation aid. They are not a scientific coupling layer and should be replaced
as the standalone LCM interface is designed.

## Full SAM-Derived Source Probe

To ask CMake to try every copied SAM-derived source file:

```sh
cmake -S . -B build -DLCM_BUILD_ALL_SAM_SOURCES=ON
cmake --build build
```

This is expected to expose remaining SAM coupling assumptions, including
precision choices, thermodynamic fields, SGS tendencies, precipitation arrays,
MPI/domain decomposition, and output interfaces.
