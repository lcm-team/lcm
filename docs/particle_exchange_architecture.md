# Standalone Particle-Exchange Architecture

Status: accepted design for Issue #40
Date: 2026-08-11

## Objective

Standalone LCM must provide one reusable particle-exchange implementation for
MPI-based LES hosts. DALES, SAM, and future host models must not each implement
LCM particle transport themselves.

Standalone in this context means independent of host-model internals. It does
not mean that the MPI-enabled LCM build is forbidden from using MPI.

## Ownership decision

LCM owns the complete particle-exchange algorithm and its MPI implementation.
The host adapter supplies the host's existing Cartesian communicator and
domain data, but it does not implement LCM communication loops.

### LCM owns

- outgoing-particle classification;
- the two cardinal exchange phases (x followed by y);
- temporary directional send and receive buffers;
- temporary destination-cell metadata;
- periodic global-coordinate wrapping;
- directional count exchange;
- particle-payload exchange;
- MPI tags, requests, statuses, and error handling;
- creation and cleanup of the particle-message MPI datatype;
- removal of outgoing particles;
- insertion of received particles into cell-owned storage; and
- safe serial behavior when parallel exchange is not configured.

### Every LES adapter supplies

- the MPI communicator already used for its horizontal decomposition;
- local and global grid geometry;
- attached atmospheric fields;
- LCM configuration; and
- the timestep used by the host advance sequence.

### The DALES adapter supplies

- `comm3d` as the existing two-dimensional Cartesian communicator;
- DALES grid dimensions, bounds, spacing, and vertical coordinates;
- DALES atmospheric field arrays; and
- DALES configuration and timestep values.

The DALES adapter must not own LCM particle buffers, MPI tags, particle-message
datatypes, or particle send/receive loops.

## Communicator contract

LCM never creates a replacement domain decomposition and never assumes
`MPI_COMM_WORLD`. The host passes its existing two-dimensional Cartesian
communicator through the public parallel interface:

```fortran
call lcm_init_parallel(lcm_parallel, host_cartesian_communicator)
call lcm_set_parallel(lcm_parallel)
```

For DALES, `host_cartesian_communicator` is `comm3d`.

LCM queries the communicator for rank, size, processor coordinates, processor
grid dimensions, periodicity, and the four cardinal neighbor ranks. Deriving
these ranks with `MPI_CART_SHIFT` preserves correctness when MPI rank
reordering is enabled. Diagonal ranks are not part of the LCM topology state.

MPI execution state belongs in `lcm_parallel_t` or its private state module.
Physical grid geometry belongs in `lcm_grid_t`. MPI topology data must not be
duplicated in both types without a documented validation requirement.

## Direction and buffer contract

Each rank exchanges particles with four cardinal neighbors in two ordered
phases:

| Phase | Send directions |
| --- | --- |
| x | east, west |
| y | south, north |

LCM owns four send buffers and four receive buffers. A buffer contains all
particles moving in one cardinal direction during the current exchange step.
LCM does not normally issue one MPI operation per particle.

The x phase is completed first: classify east/west particles, exchange counts
and payloads, and insert received particles. The y phase then repeats those
steps for south/north movement. A corner-crossing particle therefore travels
through two cardinal transfers. For example, a northeast particle first moves
east and, after insertion on the receiving rank, moves north. No diagonal rank,
buffer, tag, or direct diagonal MPI message is required.

Some directions can map to `MPI_PROC_NULL` at non-periodic physical edges.
Multiple cardinal directions can map to the same rank in small or periodic
processor grids. Direction-specific tags keep those transfers distinct.

## Particle-message contract

Permanent particle storage remains cell-owned:

```fortran
cell_particles(i,j,k)%particles(n)
```

Communication uses a temporary typed message:

```fortran
type :: lcm_particle_mpi_message_t
   integer :: global_i
   integer :: global_j
   integer :: k
   type(lcm_particle_t) :: particle
end type lcm_particle_mpi_message_t
```

`global_i`, `global_j`, and `k` identify the receiving cell. They are temporary
communication metadata and must not become permanent particle properties. The
buffer direction identifies the destination rank.

## Exchange sequence

Particle exchange runs after physical particle-boundary handling and before
local cell-ownership reassignment:

```fortran
call advance_lcm_particles(time_step)
call apply_lcm_particle_boundaries()
call exchange_lcm_particles()
call update_lcm_particle_cell_ownership()
```

`exchange_lcm_particles()` performs this sequence:

1. Run the x phase for east/west classification, dimension-specific periodic
   wrapping, count exchange, payload exchange, and received-particle insertion.
2. Run the y phase for south/north classification, dimension-specific periodic
   wrapping, count exchange, payload exchange, and received-particle insertion.

Within each phase LCM clears only that phase's two buffers, classifies and
stably compacts cell-owned particles, exchanges one count per direction,
reserves receive storage, posts payload receives and sends, waits for
completion, and inserts received particles. Completing insertion between the
phases is what makes corner movement correct.

Vertical physical-boundary behavior and local microphysics are explicitly
outside this routine.

## MPI implementation boundary

Direct MPI use is restricted to small, low-level standalone LCM modules, such
as:

- `lcm_parallel_state.f90` for communicator and topology initialization; and
- `lcm_particle_count_exchange.f90` for count transport;
- `lcm_particle_mpi_datatype.f90` for the typed message representation; and
- `lcm_particle_mpi_transport.f90` for payload transport.

These modules use `mpi_f08`. Higher-level grid, particle, buffer, boundary,
ownership, and advance modules must not call MPI directly.

The payload implementation uses an MPI derived datatype for
`lcm_particle_mpi_message_t`. It is constructed from actual component
addresses and resized to the correct array-element extent. A multi-rank test
transfers more than one message and verifies every particle field. This
avoids a manually maintained homogeneous wire array while keeping all
representation-specific logic inside one LCM transport module.

The transport posts nonblocking receives before nonblocking sends and completes
each phase with `MPI_Waitall`. It accesses the active contiguous directional
buffers directly, so it does not allocate or copy a second transport-only
payload array. This follows the useful receive-before-send behavior present in
SAM-LCM and the nonblocking communication pattern used by DALES, without
depending on DALES `modmpi`.

After transport, the x-phase receiver inserts a corner particle temporarily in
the nearest local y-edge cell because its y destination is not local yet. The
particle's physical coordinate makes the subsequent y phase extract it. A
y-phase receive requires both horizontal destination indices to be local.

## Serial behavior

If LCM parallel state has not been configured, particle exchange is a no-op.
LCM must return before classifying or removing outgoing particles. This prevents
particle loss in serial hosts and in tests that do not configure MPI exchange.

An optional build that removes the MPI dependency completely can be considered
separately. It is not required to change the ownership decision for the
MPI-enabled standalone library.

## Why LCM does not use DALES `modmpi`

Calling `modmpi` from LCM would make the standalone library depend on DALES and
would require another implementation for every future LES host. DALES's
current point-to-point wrappers also support selected homogeneous intrinsic
arrays, not the mixed typed LCM particle message.

LCM therefore uses the same underlying MPI library and the communicator
created by the host, while keeping its reusable particle-transport code inside
LCM.

## Alternatives rejected for Issue #40

### DALES-owned LCM particle transport

Rejected because it makes DALES responsible for LCM buffers, datatypes, tags,
and communication loops. Every additional LES host would have to duplicate
that implementation.

### Mandatory host transport callbacks

Rejected as the primary MPI design because they also require every host to
implement the same LCM particle-exchange operations. A pluggable non-MPI
backend could be added later without replacing the standard LCM MPI backend.

### DALES `modmpi` dependency inside LCM

Rejected because it violates host independence.

### Manual homogeneous wire arrays

Not selected for the first implementation because they duplicate every
particle field in pack and unpack indexes and are easy to desynchronize from
`lcm_particle_t`.

## Acceptance criteria

The design is satisfied when:

- DALES only passes `comm3d` and host data through its adapter;
- another MPI LES can use the same public LCM interface;
- LCM owns count and payload transport;
- LCM supports all four cardinal directions and routes periodic corners through
  the ordered x and y phases;
- unavailable neighbors are safe;
- multiple messages per direction preserve every particle field;
- outgoing particles are removed exactly once;
- received particles are inserted exactly once;
- serial execution does not remove particles; and
- the complete advance ordering remains unchanged.

## Implementation order

1. Modernize the isolated LCM MPI boundary to `mpi_f08`. **Done.**
2. Separate MPI topology state from physical grid geometry without duplicates.
   **Done.**
3. Make the DALES adapter pass only `comm3d` for parallel initialization.
   **Done in the intentionally uncommitted DALES adapter work.**
4. Implement the LCM-owned particle-message MPI datatype. **Done.**
5. Implement LCM-owned nonblocking payload exchange. **Done.**
6. Insert received messages into cell-owned storage between phases. **Done.**
7. Add standalone multi-rank LCM tests. **Done for typed arrays and all four
   corner routes on a 2x2 Cartesian communicator.**
8. Add a thin DALES integration test that verifies communicator handoff.
