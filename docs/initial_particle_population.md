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

## Current Standalone Particle API

The current standalone API already provides the storage operations needed for a
first implementation slice:

- `cell_particles(:,:,:)`
- `lcm_allocate_particle_cells(initial_capacity_per_cell)`
- `lcm_add_particle(i,j,k, particle)`
- `lcm_reset_particle_cells()`

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

## Recommended Next Slice

After this audit, the smallest implementation slice should be:

1. define an LCM-owned initial particle ID strategy
2. add a fresh-run population routine that loops over local cells and creates
   `config_state%np_xy` default particles per initialized cell
3. assign IDs, default type, zero sedimentation velocity, and cell-local counts
4. add tests for count and ID uniqueness

The aerosol, ice, tracking, particle-source, and continuous motion/location
details should remain explicit follow-up slices so they can be compared against
`micro_init_particles` carefully.
