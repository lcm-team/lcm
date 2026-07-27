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

Standalone LCM follows the SAM-LCM block strategy, but derives the subdomain
identity from `lcm_grid_t` instead of reading SAM globals directly:

```text
n_subdomains =
  n_subdomains_x * n_subdomains_y

subdomain_id =
  subdomain_y_index * n_subdomains_x + subdomain_x_index

n_step_id =
  huge(1_int64) / n_subdomains
```

This is the standalone counterpart of:

```fortran
n_step_ID = huge(1_8) / nsubdomains
n_ID = n_step_ID * rank
```

The standalone code keeps the same counter behavior: initialization seeds the
LCM-owned `n_id` counter from the subdomain block, and each newly created
particle takes the current value before the counter is incremented:

```text
n_id =
  n_step_id * subdomain_id

particle_id =
  n_id

n_id =
  n_id + 1
```

where:

- `subdomain_x_index`, `subdomain_y_index`, `n_subdomains_x`, and
  `n_subdomains_y` come from `lcm_grid_t`

This keeps the same broad ID layout as SAM-LCM while keeping the implementation
host-facing through the LCM grid metadata. If a future host needs exact MPI-rank
semantics rather than a regular Cartesian subdomain index, `lcm_grid_t` can be
extended to carry explicit rank metadata.

The code entry point is:

```fortran
call lcm_initialize_particle_id_counter()
id = lcm_next_particle_id()
```

## Initial Placement Baseline

The first standalone implementation preserves the SAM-LCM `np_xy` and `nzl`
placement structure while storing particles in cell-owned containers:

```fortran
call lcm_allocate_particle_cells(initial_capacity_per_cell=config_state%np_xy)
call lcm_initialize_computational_particles()
```

`lcm_initialize_computational_particles` is the standalone counterpart of the
placement part of SAM-LCM `micro_place_particle`. It loops over the initialized
vertical levels and local horizontal cells:

```fortran
do k = 1, min(config_state%nzl, grid_state%grid%nz)
   do j = 1, grid_state%grid%ny_local
      do i = 1, grid_state%grid%nx_local
         do n = 1, config_state%np_xy
            particle = lcm_particle_t()
            particle%id = lcm_next_particle_id()
            call init_particle_properties(i, j, k, n, particle)
            call lcm_add_particle(i, j, k, particle)
         enddo
         call lcm_init_particles(i, j, k)
      enddo
   enddo
enddo
```

This creates the same number of initial computational particles as the
SAM-LCM compatibility baseline, initializes the placement-time particle
metadata, and then calls `lcm_init_particles(i,j,k)` for the first
microphysical initialization slice.

## Initial Physical Properties

`init_particle_properties` handles the physical fields that have a clear
first-slice mapping from SAM-LCM `micro_place_particle` and
`micro_init_particles`:

| SAM-LCM field | Standalone field | First baseline |
| --- | --- | --- |
| `ivar(iTYPE)` | `particle_type` | Set to `1`, the SAM-LCM default liquid particle type. |
| `var(iM)` | `mass_water` | Set by the five-iteration equilibrium liquid initialization after aerosol mass is known. |
| `var(iN)` | `mass_aero` | Set by the aerosol distribution initializer. |
| `var(iWsedi)` | `w_sedi` | Set explicitly to zero, matching `micro_place_particle`. |
| `var(iT)` | `temperature` | Read from the attached host temperature field at `(i,j,k)`. |
| `var(iq)` | `qv` | Read from the attached host water-vapor field at `(i,j,k)`. |
| `ivar(iLEM)` | `cell_particle_index` | Set to the cell-local particle number `n`, matching the old within-cell tag. |
| `ivar(iilst:iklst)` | `last_i`, `last_j`, `last_k` | Set from the grid cell indices. |
| `ivar(iSPR)` | `spray_marker` | Set to the SAM-LCM default unsprayed marker `-999`. |
| `ivar(iTRACK)` | `track_marker` | Set to the SAM-LCM default untracked marker `-999`. |

The SAM-LCM multiplicity field `ivar(iA)` is not initialized in this first
physical-property slice because it is produced by the aerosol distribution
logic. In the current implementation, `lcm_init_particles(i,j,k)` initializes
that field for the SAM-LCM-compatible `init_aero='bin'` path.

SAM-LCM preprocesses aerosol mode configuration before this routine, including
zeroing shares for modes with zero concentration. Standalone LCM keeps the
default in that preprocessed form and stops if a zero-concentration mode is
assigned particles.

## Current Standalone Particle API

The current standalone API already provides the storage operations needed for a
first implementation slice:

- `cell_particles(:,:,:)`
- `lcm_allocate_particle_cells(initial_capacity_per_cell)`
- `lcm_add_particle(i,j,k, particle)`
- `lcm_reset_particle_cells()`
- `lcm_initialize_particle_id_counter()`
- `lcm_next_particle_id()`
- `lcm_initialize_computational_particles()`

The current `lcm_particle_t` intentionally contains named counterparts for the
SAM-LCM computational-particle state. This is a compatibility choice: during
the port, the trusted SAM-LCM routines should be translated by replacing
`particle(n)%var(...)` and `particle(n)%ivar(...)` access with named fields,
not by changing the scientific sequence.

The named fields include:

- `id`
- `weight`
- `particle_type`
- motion fields `x`, `y`, `z`, `u`, `v`, `w`, `u_sgs`, `v_sgs`, `w_sgs`
- `mass_water`
- `mass_aero`
- `w_sedi`
- `temperature`
- `qv`
- `eta`
- ice-shape and ice-state fields
- aerosol composition portions
- IN surface area per species
- cloudy age
- cell-local and tracking metadata

Cell ownership still comes from the `cell_particles(i,j,k)` container. The
motion/location fields are kept because the medium-term compatibility target
still advects particles with continuous position and velocity state. If a later
finite-particle-flux advection path removes that need, those fields can be
deleted after the replacement is verified.

## Aerosol Initialization Port Audit

The current standalone `lcm_init_particles(i,j,k)` keeps the liquid-only aerosol
initialization path close to SAM-LCM `micro_init_particles(i,j,k)`.

### Direct Field Mapping

| SAM-LCM expression | Standalone expression |
| --- | --- |
| `particle(n)%ivar(iA)` | `cell_particles(i,j,k)%particles(n)%weight` |
| `particle(n)%var(iN)` | `cell_particles(i,j,k)%particles(n)%mass_aero` |
| `particle(n)%var(iAero_portions)` | `cell_particles(i,j,k)%particles(n)%aerosol_portions` |
| `particle(n)%ivar(iTYPE)` | `cell_particles(i,j,k)%particles(n)%particle_type` |
| `particle(n)%var(in_startind+1:in_startind+4)` | `cell_particles(i,j,k)%particles(n)%in_surface_area` |
| `particle(n)%var(iTFZ)` | `cell_particles(i,j,k)%particles(n)%freezing_temperature` |
| `dV(k)` | `grid_state%dV(k)` |
| `rho_air(k)` | `field_state%density(k)` |
| aerosol namelist state | fields of `config_state` |

### Ported Liquid-Only Branches

`init_aero='bin'` preserves the SAM-LCM sequence:

1. split `np_xy` over aerosol modes using `shares_aero`
2. compute `r_aero_min` and `r_aero_max` from lognormal CDF percentiles
3. divide each mode into log-radius bins
4. compute multiplicity from `log_normal_pdf`, `dV`, and optionally `rho_air`
5. compute aerosol mass from multiplicity, radius, and aerosol density
6. set aerosol portions and liquid particle type
7. split the most populated bin when another bin receives zero multiplicity

`init_aero='USM'` preserves the SAM-LCM sequence:

1. split `np_xy` over aerosol modes using `shares_aero`
2. compute mode radius bounds with the same regula-falsi percentile search
3. use `p_USM = 1.0 / np_mode`
4. initialize particles `1:np_mode-1` from successive percentile intervals
5. use equal per-mode multiplicity for the interior samples
6. apply SAM-LCM's special largest-particle treatment for `np_done+np_mode`

`init_aero='random'` preserves the SAM-LCM sequence:

1. split `np_xy` over aerosol modes using `shares_aero`
2. use the same inverse-CDF percentile index rule based on `n_inv_cdf`
3. draw random radii using `randomu`
4. keep the SAM-LCM accept/retry loop that checks the generated radius against
   the lognormal PDF
5. use equal per-mode multiplicity
6. set aerosol mass, aerosol portions, and liquid particle type

### Intentional Differences From SAM-LCM

These differences are intentional and should remain visible during review:

- SAM-LCM stores particles in a flat `particle(:)` array and uses
  `n_start/n_end`. Standalone LCM writes directly to
  `cell_particles(i,j,k)%particles(n)`.
- SAM-LCM stores `inv_cdf(4,0:n_inv_cdf)` for `init_aero='random'`, where
  `n_inv_cdf` is 50,000,000. Standalone LCM computes the selected inverse-CDF
  entry on demand. This preserves the same percentile rule without allocating a
  very large table.
- SAM-LCM uses `randomu(rank+1)`. Standalone LCM derives the random seed ID from
  `lcm_grid_t` subdomain metadata. This keeps the host-facing source of
  rank-like information inside the LCM grid structure.
- Ice-enabled IN assignment is ported for the initial particle population.
  Later ice microphysical evolution remains outside this slice.
- The five-iteration equilibrium liquid mass, vapor, and temperature adjustment
  is ported for particle initialization. SAM-LCM mutates its internal
  `micro_field(i,j,k,imv)` and `tabs_LCM(i,j,k)` working fields; standalone LCM
  keeps the same local relaxation but stores the adjusted qv/T only on the
  initialized particles because attached host fields are read-only in this
  lifecycle.
- SAM-LCM's default `micro_user_init_particles(i,j,k)` hook is empty but is
  called after standard particle initialization. Standalone LCM preserves this
  compile-time extension concept through `lcm_user_init_particles`, passing the
  cell-owned particle container instead of SAM's flat `particle(:)` plus
  `n_start/n_end` slice.

### Cell-Local Particle Counts

SAM-LCM needs `n_start(i,j,k)` and `n_end(i,j,k)` because all computational
particles live in one flat `particle(:)` array. Those arrays identify the
contiguous slice of the flat array that belongs to each cell.

Standalone LCM does not preserve that indexing model. Cell ownership is encoded
directly by the container:

```fortran
cell_particles(i,j,k)%particles(1:cell_particles(i,j,k)%count)
```

The active particle range for a cell is therefore `1:count`, and no standalone
initialization loop should reconstruct or depend on `n_start/n_end`.

The only remaining `n_start/n_end` usage in this PR is inside the SAM-LCM source
files and the SAM comparison probe. That probe intentionally allocates
`n_start=1` and `n_end=np_xy` so it can call the original SAM-LCM
`micro_init_particles(i,j,k)` and compare the standalone result against SAM's
actual routine.

### Test Coverage

`tests/test_initial_particle_population.f90` now covers:

- default aerosol configuration from `lcm_init_config`
- one active bin aerosol mode
- multiple bin aerosol modes
- bin zero-weight splitting
- deterministic `USM` aerosol initialization
- structural `random` aerosol initialization
- equilibrium liquid mass, final particle qv/temperature, eta, and monomer
  count

The deterministic `bin` and `USM` cases compare particle weights, aerosol mass,
and aerosol portions against test-local reference implementations following the
SAM-LCM formulas. The `random` case checks the SAM-LCM-compatible invariants that
do not depend on a fixed random radius sequence: particle count, cell placement,
unique IDs, equal per-mode multiplicity, positive aerosol mass, particle type,
and aerosol portions.

`tests/test_sam_particle_population_comparison.f90` also runs the original
SAM-LCM `micro_init_particles(i,j,k)` in a tiny probe harness and compares
standalone particle weights, aerosol mass/composition, IN state, equilibrium
liquid mass, final qv/temperature, eta, and monomer count against that direct
SAM output for representative `bin`, `USM`, and ice-enabled cases.

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
| Compute aerosol radius bounds for each mode from the lognormal CDF. | Ported inside `lcm_init_particles(i,j,k)`, preserving the SAM-LCM regula-falsi sequence. |
| Support `init_aero='bin'`. | Ported for the liquid-only path, preserving deterministic aerosol bins, multiplicities, aerosol mass, aerosol portions, and zero-weight bin splitting. |
| Support `init_aero='random'`. | Ported for the liquid-only path, preserving the SAM-LCM random sampling structure and equal per-mode multiplicity. |
| Support `init_aero='USM'`. | Ported for the liquid-only path, preserving the Shima et al. uniform sampling structure, equal per-mode multiplicity for interior samples, and the SAM-LCM special handling of the largest sampled particle. |
| Interpret `unit_aero` as either `1/m3` or `1/kg`. | Preserve and validate through configuration. The conversion needs grid-cell volume and, for `1/kg`, air density. |
| Initialize particle multiplicity `ivar(iA)`. | Mapped to `lcm_particle_t%weight` for liquid-only `bin`, `USM`, and `random`. |
| Initialize aerosol mass `var(iN)`. | Mapped to `lcm_particle_t%mass_aero` for liquid-only `bin`, `USM`, and `random`. |
| Store aerosol composition portions in `var(iAero_portions)`. | Mapped to `lcm_particle_t%aerosol_portions` for liquid-only `bin`, `USM`, and `random`. |
| Reassign zero-multiplicity aerosol bins by splitting the most populated bin. | Ported for `init_aero='bin'`. |
| Assign default liquid particle type `ivar(iTYPE)=1`. | Map to `lcm_particle_t%particle_type`. |
| When ice is enabled, randomly convert a fraction of particles to IN-containing type `2`. | Ported for the initial population path. This depends on `frac_ice`, IN concentrations, and random sampling. |
| For IN-containing particles, sample IN species and IN size, store IN surface area, and compute freezing temperature `iTFZ`. | Ported with IN surface area mapped to `lcm_particle_t%in_surface_area` and freezing temperature mapped to `freezing_temperature`. |
| For non-IN particles, zero IN surface area and set `iTFZ=233.15`. | Ported for the initial population path. |
| Iterate five times to compute equilibrium liquid mass from supersaturation. | Ported for particle initialization, preserving the SAM-LCM relaxation loop. |
| Update `micro_field(i,j,k,imv)` and `tabs_LCM(i,j,k)` during equilibrium adjustment. | Ported as local `cell_qv` and `cell_temperature` relaxation. Standalone LCM does not mutate attached host fields during initialization. |
| Initialize particle liquid mass `var(iM)`. | Mapped to `lcm_particle_t%mass_water` after the five-iteration equilibrium relaxation. |
| Initialize particle temperature `var(iT)` from adjusted cell temperature plus `gamaz(k)`. | Mapped to `lcm_particle_t%temperature` from the locally adjusted cell temperature. |
| Initialize particle vapor `var(iq)`. | Mapped to `lcm_particle_t%qv` from the locally adjusted cell qv. |
| Initialize absolute supersaturation `var(ieta)`. | Mapped to `lcm_particle_t%eta`. |
| Initialize monomer count `ivar(iNmono)=1`. | Mapped to `lcm_particle_t%n_monomer`. |
| Call `micro_user_init_particles(i,j,k)`. | Ported as `lcm_user_init_particles(i,j,k, particle_cell)`, called after standard initialization with the cell-owned particle container. The default hook is empty, matching SAM-LCM. |

## `lcm_particle_t` Field Decision

`lcm_particle_t` is expanded to cover the SAM-LCM computational-particle state
with named fields. This keeps the physics port close to SAM-LCM while removing
the generic `var(:)` and `ivar(:)` indexing.

| Field or concept | Current standalone location | Decision |
| --- | --- | --- |
| Particle ID | `lcm_particle_t%id` | Required for the first baseline. The ID generator is a separate task. |
| Multiplicity / weighting factor | `lcm_particle_t%weight` | Required for the first baseline. It maps from SAM-LCM `ivar(iA)`. |
| Particle type | `lcm_particle_t%particle_type` | Required for the first baseline. Liquid-only initialization can set the SAM-LCM-compatible default type. |
| Continuous position `x/y/z` | `x`, `y`, `z` | Keep for the medium-term velocity/location advection baseline. |
| Resolved particle velocities | `u`, `v`, `w` | Keep for SAM-LCM-compatible advection. |
| SGS particle velocities | `u_sgs`, `v_sgs`, `w_sgs` | Keep for SAM-LCM-compatible SGS velocity behavior. |
| Liquid water mass | `lcm_particle_t%mass_water` | Required for `var(iM)`. |
| Aerosol mass | `lcm_particle_t%mass_aero` | Required for `var(iN)`. |
| Sedimentation velocity | `lcm_particle_t%w_sedi` | Required for the first baseline and initialized to zero, matching `micro_place_particle`. |
| Particle temperature | `lcm_particle_t%temperature` | Required for the first baseline. It should come from attached host fields through `lcm_field_state`. |
| Particle vapor mixing ratio | `lcm_particle_t%qv` | Required for the first baseline. It should come from attached host fields through `lcm_field_state`. |
| Absolute supersaturation / excess vapor | `eta` | Required for `var(ieta)`. |
| Ice axis fields | `axis_a`, `axis_c` | Keep existing fields. They are not required for a liquid-only first baseline. |
| Particle density | `density` | Keep existing field. It is not required for the first placement/count baseline. |
| Rimed mass | `rimed_mass` | Keep existing field. It is not required until ice/riming behavior is ported. |
| Freezing temperature | `freezing_temperature` | Required for `var(iTFZ)`. |
| Aerosol composition portions | `aerosol_portions(4)` | Required for `var(iAero_portions)`. |
| IN surface area per species | `in_surface_area(4)` | Required for `var(iSC:iSP)`. |
| Monomer count | `n_monomer` | Required for `ivar(iNmono)`. |
| Particle age / time | `cloudy_age` | Required for `var(itime)`. |
| Cell-local particle index | `cell_particle_index` | Counterpart of `ivar(iLEM)`. |
| Last grid indices | `last_i`, `last_j`, `last_k` | Counterpart of `ivar(iilst:iklst)`. Cell ownership still comes from the container. |
| Spray marker | `spray_marker` | Counterpart of `ivar(iSPR)`. |
| Tracking marker | `track_marker` | Counterpart of `ivar(iTRACK)`. |

This decision does not mean all fields are initialized in the first slice. It
means the data structure is ready for a faithful, near line-by-line port of
SAM-LCM routines. Each later port should preserve SAM-LCM control flow and
replace only the storage expression.

## Recommended Next Slices

The initial particle population now covers the fresh-run `micro_place_particle`
and `micro_init_particles` responsibilities in scope for issue #43. Later work
should continue with focused PRs for timestep physics and source/restart paths,
so each slice can still be compared against SAM-LCM carefully.
