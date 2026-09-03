# Build and Run DALES–LCM on LISE

This guide builds standalone LCM, links DALES to it, and runs an LCM
experiment. The build and run directories are separate from source checkouts.

## 1. Choose source and build directories

The paths below are examples. Choose the locations that suit your work; source
and build directories do not need to share a parent directory. CMake creates
the build directories when it configures the project. Do not point `LCM_ROOT`
at the LCM source directory: it must point at the completed LCM *build*
directory.

```bash
export LCM_SRC="/home/${USER}/src/lcm"
export DALES_SRC="/home/${USER}/src/dales"
export LCM_BUILD="/scratch/usr/${USER}/build/lcm-smoke"
export DALES_BUILD="/scratch/usr/${USER}/build/dales-smoke"
```

Use fresh, uniquely named build directories for a different compiler, source
revision, or build type. Never reuse an Intel CMake build directory for GNU,
or vice versa.

## 2. Prepare the DALES checkout

DALES uses Git submodules, including `external/ecbuild`, which CMake requires
during configuration. Run this once after cloning DALES, or whenever the
submodules have not been initialized:

```bash
git -C "${DALES_SRC}" submodule update --init --recursive
```

## 3. Intel + Intel MPI build

Load one consistent Intel MPI dependency stack:

```bash
module purge
module load NHRZIBenv
module load sw.clx.el9
module load slurm
module load intel/2025.3
module load impi/2021.15
module load cmake/3.31.9
module load hdf5-parallel/impi/intel/1.14.4
module load netcdf-parallel/impi/intel/4.9.2
module load fftw3/intel/3.3.10

export FC=mpiifx
export NetCDF_ROOT="$(nf-config --prefix)"
```

Configure and compile LCM:

```bash
cmake -S "${LCM_SRC}" -B "${LCM_BUILD}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=mpiifx \
  -DCMAKE_Fortran_FLAGS:STRING="-fpp"

cmake --build "${LCM_BUILD}" --parallel 8
```

The `-fpp` flag tells Intel Fortran to process LCM's conditional `#ifdef`
source blocks before compilation. Without it, Intel may compile both branches
and report duplicate declarations.

Configure and compile DALES linked to this LCM build:

```bash
cmake -S "${DALES_SRC}" -B "${DALES_BUILD}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=mpiifx \
  -DENABLE_LCM=ON \
  -DLCM_ROOT="${LCM_BUILD}" \
  -DNetCDF_ROOT="${NetCDF_ROOT}"

cmake --build "${DALES_BUILD}" --parallel 8
```

The executable is `${DALES_BUILD}/bin/dales`.

## 4. GNU + Open MPI build

Alternatively, use GNU Fortran with Open MPI instead of the Intel compiler
stack. Load the following modules; do not mix them with the Intel MPI modules
above.

```bash
module purge
module load NHRZIBenv
module load sw.clx.el9
module load slurm
module load gcc/13.3.0
module load openmpi/gcc/5.0.3
module load cmake/3.31.9
module load hdf5-parallel/ompi/gcc/1.14.4
module load netcdf-parallel/ompi/gcc/4.9.2
module load fftw3/gcc/3.3.10

export FC=mpif90
export NetCDF_ROOT="$(nf-config --prefix)"
```

Use GNU build directories separate from the Intel build directories:

```bash
export LCM_BUILD="/scratch/usr/${USER}/build/lcm-smoke-gfortran"
export DALES_BUILD="/scratch/usr/${USER}/build/dales-smoke-gfortran"

cmake -S "${LCM_SRC}" -B "${LCM_BUILD}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=mpif90
cmake --build "${LCM_BUILD}" --parallel 8

cmake -S "${DALES_SRC}" -B "${DALES_BUILD}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DENABLE_LCM=ON \
  -DLCM_ROOT="${LCM_BUILD}" \
  -DNetCDF_ROOT="${NetCDF_ROOT}"
cmake --build "${DALES_BUILD}" --parallel 8
```

LCM's CMake configuration enables GNU preprocessing (`-cpp`) automatically.

## 4. Run LCM tests with CTest

CMake registers LCM's unit and integration tests during configuration. `ctest`
reads that test list from the LCM build directory and runs it. Some tests use
MPI, so run the suite in a compute-node allocation with at least four tasks:

```bash
salloc --partition=cpu-clx:test --nodes=1 --ntasks=4 --time=00:10:00
ctest --test-dir "${LCM_BUILD}" --output-on-failure
exit
```

To list the registered tests without running them:

```bash
ctest --test-dir "${LCM_BUILD}" -N
```

To run tests whose names match a pattern, for example particle exchange tests:

```bash
ctest --test-dir "${LCM_BUILD}" -R particle_exchange --output-on-failure
```

## 5. Run a DALES–LCM experiment

Choose a run directory containing the DALES input files and edit
`namoptions.001` for the scientific case. To activate LCM, it must
contain at least:

```fortran
&NAMMICROPHYSICS
imicro = 12
/

&NAMLCM
np_xy = 50
nzl = 10
init_aero = 'bin'
unit_aero = '1/m3'
switch_ice = .false.
/
```

The MPI layout in the namelist must agree with the Slurm task count. The smoke
example below uses `nprocx = 2`, `nprocy = 2`, and four tasks. Change those
values together when using another decomposition.

Create `submit_lcm_smoke.sbatch` in the run directory. This is the Intel
smoke-test template:

```bash
#!/bin/bash
#SBATCH --job-name=lcm-smoke4
#SBATCH --partition=cpu-clx:test
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH --time=00:10:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

module purge
module load NHRZIBenv
module load sw.clx.el9
module load slurm
module load intel/2025.3
module load impi/2021.15
module load hdf5-parallel/impi/intel/1.14.4
module load netcdf-parallel/impi/intel/4.9.2
module load fftw3/intel/3.3.10

export LD_LIBRARY_PATH="$(nf-config --prefix)/lib:${LD_LIBRARY_PATH:-}"
export DALES=/scratch/usr/bekmah93/lcm_particle_scalability/build/dales-smoke/bin/dales

srun --nodes=1 --ntasks=4 --ntasks-per-node=4 --cpus-per-task=1 \
  --cpu-bind=cores "${DALES}" namoptions.001 > dales.out 2> dales.err
```

For a GNU build, use the same Slurm resources and `srun` line, but load the
same GNU/Open-MPI modules used at build time and set `DALES` to the GNU binary.

Submit and monitor:

```bash
sbatch submit_lcm_smoke.sbatch
squeue -u "$USER"
```

After completion, inspect `dales.out` and `dales.err`.
