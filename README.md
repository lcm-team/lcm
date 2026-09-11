# LCM

LCM is a Lagrangian Cloud Model being prepared as a standalone Fortran library.

The current development goal is to separate the inherited LCM source from its
original SAM-LES integration and define a clean interface that can later be used
from DALES.

## Current Status

This repository currently contains the LCM source inherited from the SAM-LCM
implementation. A first standalone CMake build is available for a core subset of
the code using temporary compatibility stubs.

The compatibility stubs are only a build aid while SAM-specific dependencies are
identified and moved behind a cleaner host-model interface.

## Repository Layout

- `src/`: LCM Fortran source files
- `stubs/`: temporary compatibility stubs for standalone builds
- `docs/`: build notes and legacy/reference documentation

## Building

See [`docs/building.md`](docs/building.md) for the current standalone build
instructions.

The accepted ownership and MPI boundary for standalone particle exchange is
documented in
[`docs/particle_exchange_architecture.md`](docs/particle_exchange_architecture.md).

## Development Checks

This repository uses the Python `pre-commit` framework for lightweight,
shared checks before each commit.

Install it once in your Python environment:

```sh
pip install pre-commit
```

Enable the repository hooks once per clone:

```sh
pre-commit install
```

Run the checks manually on all files:

```sh
pre-commit run --all-files
```

## Legacy Notes

The original imported SAM-LCM readme is kept for reference at
[`docs/legacy_readme.txt`](docs/legacy_readme.txt).
