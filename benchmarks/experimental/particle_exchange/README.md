# Experimental particle-exchange benchmark

This opt-in benchmark measures the existing standalone LCM count and typed
particle-payload transport. It is intentionally separate from the correctness
test suite: timing results depend on the MPI implementation, compiler, machine,
rank placement, and system load, so the executable has no performance pass/fail
threshold and is not registered with `ctest`.

Configure and build it with:

```sh
cmake -S . -B build-bench -DLCM_BUILD_EXPERIMENTAL_BENCHMARKS=ON
cmake --build build-bench --target benchmark_particle_exchange
```

Run with two or four MPI ranks:

```sh
mpiexec -n 2 build-bench/benchmarks/experimental/particle_exchange/benchmark_particle_exchange
mpiexec -n 4 build-bench/benchmarks/experimental/particle_exchange/benchmark_particle_exchange
```

Save machine-specific output in the ignored `results` directory:

```sh
mpiexec -n 4 \
  build-bench/benchmarks/experimental/particle_exchange/benchmark_particle_exchange \
  20 10000 \
  > benchmarks/experimental/particle_exchange/results/particle_exchange_4ranks.csv
```

The `results` directory retains only its `.gitignore`; generated CSV files and
other local benchmark output are not tracked by Git.

The first optional argument is the number of measured iterations. The second
is the maximum number of messages sent in each direction. Defaults are 20
iterations and 10,000 messages per direction:

```sh
mpiexec -n 4 build-bench/benchmarks/experimental/particle_exchange/benchmark_particle_exchange 50 100000
```

The driver exercises both cardinal directions in each x or y phase and calls
the production `lcm_exchange_particle_counts()` and
`lcm_exchange_particle_payloads()` routines. Directional-buffer preparation,
MPI barriers, and correctness checks are outside the timed regions. Output is
CSV with the maximum elapsed time across ranks averaged over all measured
iterations. `message_payload_bytes` is the number of bytes MPI transfers for
one derived-datatype message, while `message_extent_bytes` is its stride in the
contiguous directional buffer. Four ranks form a 2-by-2 Cartesian topology, so
the x and y phases target different neighboring ranks. With two ranks, MPI
creates a 2-by-1 topology and the y phase communicates with the same rank; use
the four-rank result when assessing both dimensions.
