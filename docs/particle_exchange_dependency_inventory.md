# Particle Exchange Dependency Inventory

The accepted standalone ownership and MPI design is documented in
[`particle_exchange_architecture.md`](particle_exchange_architecture.md).
This inventory remains the detailed comparison with the inherited SAM-LCM
implementation.

This note inventories SAM-LCM's particle-exchange implementation before
implementing standalone `exchange_lcm_particles()` (Issue #40). The reference
is `MICRO_LAGRANGE/micro_exchange.f90`, especially
`micro_exchange_particles` and `task_exchange_particle`.

## Scope and advance contract

Particle exchange runs after physical boundary handling and before the final
local cell-ownership update:

```fortran
call advance_lcm_particles(time_step)
call apply_lcm_particle_boundaries()
call exchange_lcm_particles()
call update_lcm_particle_cell_ownership()
```

Exchange owns horizontal inter-rank transfer only. It must not apply vertical
top/bottom rules, microphysics, or DALES feedback. A serial host must be able
to call it as a no-op.

## SAM-LCM algorithm inventory

SAM's `micro_exchange_particles` has the following phases:

1. Classify every flat-array particle as local or outgoing in one of eight
   horizontal directions: east, west, southeast, northwest, south, north,
   southwest, and northeast.
2. Apply global horizontal periodic wrapping to coordinates before transfer.
3. Build one linked list per outgoing direction and count its particles.
4. Exchange directional counts with neighbor ranks.
5. Allocate flat send/receive buffers using the maximum incoming or outgoing
   directional count.
6. Pack each linked list into a `particle_array` buffer, exchange payloads,
   and collect received particles in a flat `particle_in(:)` array.
7. Rebuild the legacy flat `particle(:)` array from surviving local particles
   plus received particles.

`task_exchange_particle` posts nonblocking receives before it performs
blocking sends. It uses tags to map an outgoing direction to the corresponding
incoming direction. If periodicity makes a neighbor equal to the local rank,
SAM copies the outgoing buffer locally instead of using MPI.

## Direction mapping used by SAM

| SAM buffer index | Outgoing direction | Opposite incoming direction |
| --- | --- | --- |
| 1 | east | west |
| 2 | west | east |
| 3 | southeast | northwest |
| 4 | northwest | southeast |
| 5 | south | north |
| 6 | north | south |
| 7 | southwest | northeast |
| 8 | northeast | southwest |

Standalone LCM intentionally does not reproduce these eight direct-neighbor
identifiers. It owns four cardinal directions and two ordered phases:

| Standalone phase | Directions |
| --- | --- |
| x | east, west |
| y | south, north |

A corner crossing is routed through both phases. The x payload is received and
inserted before y classification, so a northeast crossing becomes east then
north. This follows the cell-owned, cardinal-pass approach used by PALM and
removes diagonal topology, buffers, and transport operations from standalone
LCM.

### Local subdomain-bound contract

LCM particle cells are always locally indexed as `1:nx_local` and
`1:ny_local`.  The inclusive global bounds used to locate and classify a
particle are derived from the grid state:

```fortran
global_i_first = global_i_start
global_i_last  = global_i_start + nx_local - 1
global_j_first = global_j_start
global_j_last  = global_j_start + ny_local - 1
```

`i_start:i_end` and `j_start:j_end` remain host-array interior indices (for
DALES, usually `2:i1` and `2:j1`); they must not be used as LCM particle-cell
indices.  This avoids duplicating the derived global end indices and keeps the
DALES array layout separate from LCM cell-owned storage.

### Directional message buffers

LCM owns temporary, flat send and receive buffers for each of the four cardinal
directions. Each entry is an `lcm_particle_mpi_message_t`; the buffer holds a
count and grows as needed.  Buffers are cleared at the start of a later
exchange phase but retain their allocation for reuse. Clearing one phase does
not disturb the other phase's transport state. They are only
temporary MPI-transfer state: they do not change cell-owned particle storage.

After classification removes an outgoing particle from its local cell, LCM
calculates its exact global horizontal cell indices from position, packs those
indices and `k` with the particle payload, and appends the message to the send
buffer identified by its exchange direction. The MPI transport sends these
buffers directly without building another outgoing payload copy.
At a global periodic edge, the outgoing copy's destination index and matching
physical coordinate are normalized into the opposite side of the global
domain. Internal MPI-subdomain crossings are not wrapped.

### MPI count exchange

Before payload exchange, LCM exchanges one integer message count per outgoing
direction. For each cardinal phase it posts both nonblocking count receives,
then both nonblocking count sends, and completes all four requests with one
`MPI_WAITALL`. A count sent toward direction `D` uses `D` as its MPI tag and is
received from the opposite neighbor with that same tag. It is stored as the
expected count for the receive buffer whose direction denotes the source
neighbor. Receive storage is reserved from this count, but the buffer's actual
message count remains zero until payload transfer completes.
At non-periodic edges, each unavailable source or destination is mapped
independently to `MPI_PROC_NULL`, so a rank can still communicate through its
valid side.  A missing source produces an expected receive count of zero.
The MPI tests include a four-rank `2x2` Cartesian communicator that checks all
four LCM cardinal neighbor lookups and verifies directional count delivery in
both phases.

### MPI payload exchange and insertion

After count exchange, LCM posts typed nonblocking receives before typed
nonblocking sends for the active phase and completes them together with
`MPI_WAITALL`. The MPI derived datatype represents
`lcm_particle_mpi_message_t`, including destination metadata and every particle
field. Transport obtains contiguous views of the directional storage, so the
send path does not duplicate its payload in a second flat array.

Completed receive storage is inserted into cell-owned particle containers. An
x-phase corner receive may still have an out-of-range y destination; it is
placed temporarily in the nearest y-edge cell and is classified again in the
y phase. A four-rank non-periodic `2x2` test routes NE, NW, SE, and SW particles
through these two cardinal transfers and verifies final ownership and payload
fields. A one-rank periodic test exercises the same production transport path
for wrapping and self-neighbor communication.

When no LCM parallel state has been configured, the complete exchange step is
an explicit no-op: it does not classify or remove particles, clear directional
buffers, or call MPI.  This protects serial hosts and prevents particle loss if
no communicator is available.

### MPI interface boundary

Direct MPI interfaces and calls in standalone LCM are restricted to low-level
topology and transport modules:

- `lcm_parallel_state.f90` queries the host-provided Cartesian topology.
- `lcm_particle_count_exchange.f90` performs directional count transport.

Higher-level grid, particle, boundary, buffer, and advance modules operate on
LCM types and APIs only.  A CMake test scans standalone `lcm_*.f90` sources and
fails if `mpif.h`, `use mpi`, `use mpi_f08`, or a direct `MPI_*` call appears
outside the declared boundary.  This keeps a later migration from `mpif.h` to
a modern MPI module local to the transport/topology layer.

SAM performs all eight directions in a 3-D horizontal run and only east/west
in its two-dimensional mode. Standalone LCM instead selects ordered x and y
cardinal phases. Particles that cross both local horizontal boundaries remain
defined because the first phase's received payload is inserted before the
second phase classifies it.

## SAM dependencies and standalone decisions

| SAM dependency | Role in SAM | Standalone LCM decision |
| --- | --- | --- |
| `MPI_COMM_WORLD`, `rank`, `rankee`, `rankww`, `rankss`, `ranknn`, and diagonal ranks | Communicator and neighbor topology. | The host supplies its Cartesian communicator through `lcm_parallel_t`. LCM derives only cardinal ranks with `MPI_CART_SHIFT`, so it handles rank reordering without storing diagonal topology. LCM must never assume `MPI_COMM_WORLD`. |
| `gpxs`, `gpys`, `xwidth`, `ywidth` | Classify local vs outgoing by coordinate and wrap periodic global coordinates. | Receive local x/y bounds and periodicity through grid/parallel state. Keep wrapping explicitly separate from physical vertical boundaries. |
| `particle(:)`, linked-list heads/tails, `particle_in(:)` | Separate outgoing/local particles and rebuild legacy storage. | Traverse `cell_particles(:,:,:)`; use flat temporary communication buffers only. Classify and stably compact each cell in one pass, then append incoming particles without rebuilding all permanent storage. |
| `n_out(:)`, `n_in(:)`, `np_out`, `np_in` | Directional count exchange and allocation. | Retain directional counts as temporary exchange metadata. |
| `particle_array`, `micro_particle_array2mpi_datatype` | Sends arrays of SAM integer/real particle fields through a custom MPI datatype. | Define an explicit standalone packed message. It must include temporary destination `i`, `j`, `k` metadata plus named `lcm_particle_t` payload fields. |
| `micro_task_recv_particle`, `micro_task_bsend_particle` | Posts receives and sends a payload with SAM's datatype. | Isolate typed MPI operations in a standalone exchange module; use the host communicator. Start with a clear count-then-payload exchange. |
| `RUN3D` | Selects two versus eight directions. | Standalone defines four cardinal directions and ordered x/y phases explicitly; it does not inherit SAM's mode flag. |
| `particle(n)%next` | Linked-list storage topology. | No standalone counterpart. Cell ownership and flat temporary buffers replace it. |

## Important compatibility observations

- SAM classifies from continuous x/y coordinates, not permanent cell indices.
- SAM wraps only at global periodic limits (`xwidth`, `ywidth`), not at every
  local subdomain boundary. Crossing a local bound selects a neighbor.
- The existing standalone particle type already keeps `x`, `y`, and `z`, so
  coordinate-based classification remains possible without reintroducing SAM's
  flat field indices.
- Directional `lcm_particle_mpi_message_t` buffers are the sole temporary
  storage for particles removed from local cell ownership. Keeping a second
  flat outgoing-particle array would duplicate every payload and add avoidable
  memory traffic.
- SAM's custom MPI datatype assumes a particular in-memory arrangement of its
  legacy `particle_array`. Standalone packing must not make equivalent layout
  assumptions about Fortran derived types.

## Required standalone separation

Permanent particle ownership remains cell-owned:

```fortran
cell_particles(i,j,k)%particles(n)
```

Temporary exchange messages require destination metadata because a receiving
rank must know where to insert the payload:

```fortran
type :: lcm_particle_mpi_message_t
   integer :: global_i
   integer :: global_j
   integer :: k
   type(lcm_particle_t) :: particle
end type lcm_particle_mpi_message_t
```

The horizontal message indices are global destination indices, so sender and
receiver never confuse their local cell numbering. All message indices are
communication metadata only. They must not become permanent physical
properties of `lcm_particle_t`.

## Follow-up design decisions

1. Define the neutral `lcm_parallel_t` state and setter. **Done:** it contains
   only the host communicator; grid decomposition remains in `lcm_grid_t`.
2. Decide the first implementation's direction set. **Done:** use four cardinal
   directions in ordered x/y phases. LCM derives cardinal ranks from the host
   Cartesian communicator and does not store diagonal ranks.
3. Define local x/y bounds used for classification and their relationship to
   the host grid's interior/halo cells.
4. Define portable typed packing/unpacking and count/payload exchange semantics.
5. Update `lcm_advance` to call exchange before ownership. **Done:** exchange
   preparation extracts outgoing particles before local ownership reassignment.
