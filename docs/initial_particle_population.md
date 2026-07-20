# Standalone LCM Initial Particle Population

This note tracks the first implementation task for issue #43: auditing the
SAM-LCM `micro_place_particle` routine before replacing it with an LCM-owned
particle population path.

The standalone direction is to keep the scientific initialization behavior as
close as possible to SAM-LCM while changing the storage model from the legacy
flat `particle(:)` array to cell-owned `cell_particles(i,j,k)` containers.

## SAM-LCM Entry Point

SAM-LCM creates the initial particle population through:

```fortran
call micro_place_particle
```

This call happens in the LCM initialization block of `micro_proc`, after:

- particle IDs are initialized with a rank-dependent offset
- temperature fields are copied into LCM working fields
- `micro_coupling` exchanges LES fields into SAM-LCM internal arrays

In standalone LCM, the equivalent entry point should be called from the fresh
initialization path of `lcm_initialize`, after config, grid, atmospheric field
attachments, and cell-particle storage allocation are complete.

## `micro_place_particle` Responsibilities

`micro_place_particle` currently owns these responsibilities:

| SAM-LCM responsibility | Standalone LCM direction |
| --- | --- |
| Copy `n_init` into `n_particles`. | Use the expected particle count already derived from `lcm_config_t` and `lcm_grid_state`. |
| Compute `n_particles_alloc` with extra growth capacity. | Avoid a process-wide flat allocation. Reserve per-cell capacity through `lcm_allocate_particle_cells(initial_capacity_per_cell=config_state%np_xy)`. |
| Allocate the flat `particle(:)` array. | Use allocated `cell_particles(:,:,:)`; do not reintroduce flat storage as the main container. |
| Allocate temporary `inv_cdf(4,0:n_inv_cdf)` for aerosol initialization. | Keep this tied to aerosol initialization, not to particle placement itself. The later aerosol initializer can own temporary lookup storage. |
| Reset `n_start(:,:,:)` and `n_end(:,:,:)`. | Replace with `cell_particles(i,j,k)%count`; each cell owns its own active particle count. |
| Loop over local cells `k=1:nzl`, `j=1:ny_LCM`, `i=1:nx_LCM`. | Preserve this as the first compatibility baseline, using `grid_state%grid%nx_local`, `ny_local`, and the configured initialization depth. |
| Create `np_xy` particles per initialized cell. | Preserve this baseline using `config_state%np_xy`. |
| Store within-cell particle index `iLEM=m`. | Decide whether this is still needed. If needed for compatibility, it should become explicit standalone state rather than an `ivar` index. |
| Store global cell indices `iilst`, `ijlst`, `iklst`. | The permanent particle should not need to store cell ownership because ownership comes from `cell_particles(i,j,k)`. If global indices are needed for diagnostics or intermediate advection, keep them separate from physical particle properties. |
| Initialize sprayed-particle marker `iSPR=-999`. | Leave to a future particle-source/spray design unless required by the first fresh-run baseline. |
| Optionally initialize tracking marker `iTRACK`. | Leave to a future diagnostics/tracking design. It should not block the first fresh-run particle population. |
| Assign unique particle ID from `n_ID`, then increment `n_ID`. | Implement an LCM-owned ID generator. It should preserve uniqueness across ranks without relying directly on SAM globals. |
| Initialize particle age/time `itime=0`. | Add only if the standalone particle type keeps particle age. The current `lcm_particle_t` does not yet include this field. |
| Initialize continuous position `iX`, `iY`, `iZ` using random offsets inside the cell. | Decide as part of the motion/location task. For the cell-owned design, cell membership is the primary location. Continuous position may still be needed temporarily for the velocity-based advection baseline. |
| Initialize resolved velocities `iU`, `iV`, `iW` to zero. | Keep out of physical `lcm_particle_t` unless the intermediate advection path requires particle-carried velocity state. |
| Initialize SGS velocities `iUsgs`, `iVsgs`, `iWsgs` to zero. | Same as resolved velocities: separate motion state from permanent physical particle properties. |
| Initialize sedimentation velocity `iWsedi=0`. | Covered by `lcm_particle_t%w_sedi`, which defaults to zero and should be explicitly set by the particle initializer. |
| Call `micro_init_particles(i,j,k)` after the particles for a cell are placed. | Replace with a standalone microphysical-property initializer that fills `lcm_particle_t` fields for each cell. |
| Deallocate `inv_cdf`. | Keep temporary aerosol lookup cleanup inside the future aerosol initializer. |

## Related Setup Outside `micro_place_particle`

SAM-LCM computes the initial particle count in `micro_set`:

```fortran
n_init = nx_gl / nsubdomains_x * ny_gl / nsubdomains_y * np_xy * nzl
n_total = nsubdomains * n_init
```

Issue #42 already moved this responsibility into the standalone initialization
lifecycle as expected local and global particle counts. Issue #43 should use
that same count model when creating actual particles.

SAM-LCM initializes unique particle IDs in `micro_proc`:

```fortran
n_ID = n_step_ID * rank
```

The standalone implementation should not use SAM globals directly. The next
design step should define an LCM-owned ID strategy using whatever parallel or
domain metadata the standalone LCM exposes.

## LCM-Owned Particle ID Generation

Standalone LCM uses deterministic IDs derived from the particle's initial
global grid cell and its within-cell particle index:

```text
global_cell_id =
  ((k - 1) * ny_global + (global_j - 1)) * nx_global + global_i

particle_id =
  (global_cell_id - 1) * particles_per_cell + cell_particle_index
```

where:

- `global_i = grid_state%grid%global_i_start + i - 1`
- `global_j = grid_state%grid%global_j_start + j - 1`
- `cell_particle_index` is the `m` loop index from `1:np_xy`
- `particles_per_cell` is the configured initial particles per cell

This replaces SAM-LCM's `n_ID = n_step_ID * rank` strategy. It does not require
a SAM rank or a host-specific process numbering scheme, but it still preserves
unique IDs across decomposed domains as long as the host supplies correct global
cell offsets.

The code entry point is:

```fortran
id = lcm_particle_id_for_cell_member(i, j, k, m, config_state%np_xy)
```

## Current Standalone Particle API

The current standalone API already provides the storage operations needed for a
first implementation slice:

- `cell_particles(:,:,:)`
- `lcm_allocate_particle_cells(initial_capacity_per_cell)`
- `lcm_add_particle(i,j,k, particle)`
- `lcm_reset_particle_cells()`
- `lcm_particle_id_for_cell_member(i,j,k,cell_particle_index,particles_per_cell)`

The current `lcm_particle_t` contains named physical fields such as:

- `id`
- `weight`
- `particle_type`
- `mass_water`
- `mass_aero`
- `w_sedi`
- `temperature`
- `qv`
- ice-shape and ice-state fields

It does not currently contain the full set of SAM-LCM motion/location fields.
That is consistent with the design goal that cell ownership should come from
the cell container, not from permanent particle coordinates.

## `micro_init_particles` Property Audit

`micro_init_particles(i,j,k)` is called after `micro_place_particle` has created
all `np_xy` particles for one grid cell. It fills the particle microphysical
properties and adjusts the cell's vapor and temperature state toward an initial
equilibrium.

The routine is more than a field setter. It also owns aerosol sampling,
ice-nucleating particle assignment, and a local thermodynamic adjustment.

| SAM-LCM responsibility | Standalone LCM direction |
| --- | --- |
| Return immediately when `np_xy == 0`. | Preserve this guard in the standalone fresh-population path. |
| Split `np_xy` particles across four aerosol modes using `shares_aero`. | Preserve as the compatibility baseline, but source values from `lcm_config_t` rather than `micro_vars`. |
| Compute aerosol radius bounds for each mode from the lognormal CDF. | Move into a standalone aerosol-initialization helper. This is scientific behavior and should not be hidden inside placement. |
| Support `init_aero='bin'`. | Preserve, because it defines deterministic aerosol bins and multiplicities. |
| Support `init_aero='random'`. | Preserve, using standalone random-number access rather than direct SAM-LCM globals. |
| Support `init_aero='USM'`. | Preserve as a named aerosol initialization mode. |
| Interpret `unit_aero` as either `1/m3` or `1/kg`. | Preserve and validate through configuration. The conversion needs grid-cell volume and, for `1/kg`, air density. |
| Initialize particle multiplicity `ivar(iA)`. | Map to `lcm_particle_t%weight`. |
| Initialize aerosol mass `var(iN)`. | Map to `lcm_particle_t%mass_aero`. |
| Store aerosol composition portions in `var(iAero_portions)`. | The current `lcm_particle_t` has one `mass_aero` field only. If composition remains scientifically required, the particle type needs additional named composition state or a separate aerosol payload. |
| Reassign zero-multiplicity aerosol bins by splitting the most populated bin. | Preserve in the aerosol helper to avoid zero-weight computational particles. |
| Assign default liquid particle type `ivar(iTYPE)=1`. | Map to `lcm_particle_t%particle_type`. |
| When ice is enabled, randomly convert a fraction of particles to IN-containing type `2`. | Preserve under an ice-enabled config path. This depends on `frac_ice`, IN concentrations, and random sampling. |
| For IN-containing particles, sample IN species and IN size, store IN surface area, and compute freezing temperature `iTFZ`. | Current `lcm_particle_t` has `freezing_temperature`, but not per-species IN surface-area fields. The standalone particle type or an ice payload must be extended before full ice parity. |
| For non-IN particles, zero IN surface area and set `iTFZ=233.15`. | Preserve when ice fields exist; otherwise document as omitted in a liquid-only first slice. |
| Iterate five times to compute equilibrium liquid mass from supersaturation. | Preserve eventually. It couples particle liquid mass, vapor, and temperature and should be handled by a dedicated microphysical initialization helper. |
| Update `micro_field(i,j,k,imv)` and `tabs_LCM(i,j,k)` during equilibrium adjustment. | Standalone LCM should not copy host fields into SAM-style working arrays. It needs a clear decision: either initialize particle mass from attached host fields without mutating host fields in the first one-way baseline, or later write feedback through explicit tendency/state design. |
| Initialize particle liquid mass `var(iM)`. | Map to `lcm_particle_t%mass_water`. |
| Initialize particle temperature `var(iT)` from adjusted cell temperature plus `gamaz(k)`. | Map to `lcm_particle_t%temperature`, using attached temperature/static-energy state rather than SAM globals. |
| Initialize particle vapor `var(iq)`. | Map to `lcm_particle_t%qv`. |
| Initialize absolute supersaturation `var(ieta)`. | Current `lcm_particle_t` has no supersaturation field. Decide whether this should be derived on demand from attached fields or stored explicitly. |
| Initialize monomer count `ivar(iNmono)=1`. | Current `lcm_particle_t` has no monomer-count field. Add only if needed for the preserved microphysics path. |
| Call `micro_user_init_particles(i,j,k)`. | Do not include in the first standalone baseline. Decide separately whether standalone LCM supports user hooks. The current SAM hook is empty by default. |

## Required Field Decisions From The Audit

The current `lcm_particle_t` is enough for a liquid-only skeleton with IDs,
weights, type, aerosol mass, water mass, sedimentation velocity, temperature,
and vapor.

Full parity with `micro_init_particles` still needs decisions for:

- aerosol composition portions across modes
- IN surface-area state per ice-nucleating species
- stored supersaturation versus derived supersaturation
- monomer count
- optional particle age/time
- temporary continuous position and velocity state for the intermediate
  advection baseline
- user particle initialization hooks

## `lcm_particle_t` Field Decision

The first standalone population baseline should not expand `lcm_particle_t`
yet. It should use the fields that already exist and only initialize values
that have a clear SAM-LCM counterpart.

| Field or concept | Current standalone location | Decision |
| --- | --- | --- |
| Particle ID | `lcm_particle_t%id` | Required for the first baseline. The ID generator is a separate task. |
| Multiplicity / weighting factor | `lcm_particle_t%weight` | Required for the first baseline. It maps from SAM-LCM `ivar(iA)`. |
| Particle type | `lcm_particle_t%particle_type` | Required for the first baseline. Liquid-only initialization can set the SAM-LCM-compatible default type. |
| Liquid water mass | `lcm_particle_t%mass_water` | Required, but a simple first baseline may initialize it to zero until the aerosol/equilibrium initializer is ported. |
| Aerosol mass | `lcm_particle_t%mass_aero` | Required, but one scalar only supports a reduced aerosol representation. Full aerosol-mode parity needs a follow-up decision. |
| Sedimentation velocity | `lcm_particle_t%w_sedi` | Required for the first baseline and initialized to zero, matching `micro_place_particle`. |
| Particle temperature | `lcm_particle_t%temperature` | Required for the first baseline. It should come from attached host fields through `lcm_field_state`. |
| Particle vapor mixing ratio | `lcm_particle_t%qv` | Required for the first baseline. It should come from attached host fields through `lcm_field_state`. |
| Ice axis fields | `axis_a`, `axis_c` | Keep existing fields. They are not required for a liquid-only first baseline. |
| Particle density | `density` | Keep existing field. It is not required for the first placement/count baseline. |
| Rimed mass | `rimed_mass` | Keep existing field. It is not required until ice/riming behavior is ported. |
| Freezing temperature | `freezing_temperature` | Keep existing field. Full ice parity needs IN surface-area state before this can match SAM-LCM. |
| Aerosol composition portions | Missing | Do not add in the first baseline. Add only when porting SAM-LCM aerosol-mode initialization. |
| IN surface area per species | Missing | Do not add in the first baseline. This belongs with ice initialization design. |
| Supersaturation / absolute excess vapor | Missing | Do not add in the first baseline. Prefer deriving it from attached atmospheric state unless later microphysics needs stored per-particle state. |
| Monomer count | Missing | Do not add in the first baseline. Add only if the preserved microphysics path requires it. |
| Particle age / time | Missing | Do not add in the first baseline. It is initialized by SAM-LCM but not needed for count/cell placement. |
| Continuous position `x/y/z` | Missing from `lcm_particle_t` | Do not add as permanent physical particle state here. Decide separately with the motion/location task. |
| Resolved and SGS particle velocities | Missing from `lcm_particle_t` | Do not add as permanent physical particle state here. Keep with motion/advection design if needed. |
| Cell ownership indices | `cell_particles(i,j,k)` container | Do not store in `lcm_particle_t` for the first baseline. Cell ownership is represented by the container. |

This decision means the first implementation can create valid
`lcm_particle_t` objects without changing the type definition. The first
baseline should focus on count, IDs, cell placement, default type,
sedimentation velocity, and simple thermodynamic values from attached fields.
Scientific parity for aerosol composition, ice, supersaturation, and
motion/location should remain explicit follow-up tasks.

## Recommended Next Slice

After the placement and property audits, the smallest implementation slice
should be:

1. add a fresh-run population routine that loops over local cells and creates
   `config_state%np_xy` default particles per initialized cell
2. assign IDs, default type, zero sedimentation velocity, and cell-local counts
3. add tests for count and ID uniqueness

The aerosol, ice, tracking, particle-source, and continuous motion/location
details should remain explicit follow-up slices so they can be compared against
`micro_init_particles` carefully.
