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
- the eight exchange directions;
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
grid dimensions, periodicity, and neighbor ranks. Deriving neighbor ranks from
the communicator preserves correctness when MPI rank reordering is enabled.

MPI execution state belongs in `lcm_parallel_t` or its private state module.
Physical grid geometry belongs in `lcm_grid_t`. MPI topology data must not be
duplicated in both types without a documented validation requirement.

## Direction and buffer contract

Each rank can exchange particles in eight horizontal directions:

| Direction | Opposite direction |
| --- | --- |
| east | west |
| west | east |
| north | south |
| south | north |
| northeast | southwest |
| northwest | southeast |
| southeast | northwest |
| southwest | northeast |

LCM owns eight send buffers and eight receive buffers. A buffer contains all
particles moving in one direction during the current exchange step. LCM does
not normally issue one MPI operation per particle.

Some directions can map to `MPI_PROC_NULL` at non-periodic physical edges.
Multiple directions can map to the same rank in small or periodic processor
grids. Direction-specific tags keep those transfers distinct.

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

`exchange_lcm_particles()` performs these phases:

1. Clear temporary directional buffers.
2. Classify and remove particles leaving the local horizontal subdomain.
3. Apply global periodic wrapping only when crossing a global periodic edge.
4. Pack each outgoing particle into its directional send buffer.
5. Exchange one particle count per direction.
6. Reserve directional receive-buffer storage.
7. Post particle-payload receives.
8. Post particle-payload sends.
9. Wait for communication to complete.
10. Store received typed messages.
11. Insert received particles into local cell-owned storage.

Vertical physical-boundary behavior and local microphysics are explicitly
outside this routine.

## MPI implementation boundary

Direct MPI use is restricted to small, low-level standalone LCM modules, such
as:

- `lcm_parallel_state.f90` for communicator and topology initialization; and
- `lcm_particle_mpi_transport.f90` for count and payload transport.

These modules use `mpi_f08`. Higher-level grid, particle, buffer, boundary,
ownership, and advance modules must not call MPI directly.

The first payload implementation uses an MPI derived datatype for
`lcm_particle_mpi_message_t`. It must be constructed from actual component
addresses and resized to the correct array-element extent. A multi-rank test
must transfer more than one message and verify every particle field. This
avoids a manually maintained homogeneous wire array while keeping all
representation-specific logic inside one LCM transport module.

The transport should post nonblocking receives and sends and complete them with
`MPI_Waitall`. This follows the useful receive-before-send behavior present in
SAM-LCM and the nonblocking communication pattern used by DALES, without
depending on DALES `modmpi`.

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
- LCM supports all eight directions and periodic corners;
- unavailable neighbors are safe;
- multiple messages per direction preserve every particle field;
- outgoing particles are removed exactly once;
- received particles are inserted exactly once;
- serial execution does not remove particles; and
- the complete advance ordering remains unchanged.

## Implementation order

1. Modernize the isolated LCM MPI boundary to `mpi_f08`.
2. Separate MPI topology state from physical grid geometry without duplicates.
3. Make the DALES adapter pass only `comm3d` for parallel initialization.
4. Implement the LCM-owned particle-message MPI datatype.
5. Implement LCM-owned nonblocking payload exchange.
6. Insert received messages into cell-owned storage.
7. Add standalone multi-rank LCM tests.
8. Add a thin DALES integration test that verifies communicator handoff.
