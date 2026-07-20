# Standalone LCM Initialization Lifecycle

This note tracks how the old SAM-LCM `micro_set` responsibilities map onto
the standalone LCM host-facing lifecycle.

The intended host call order is:

```fortran
call lcm_set_config(...)
call lcm_set_grid(...)
call lcm_attach_fields(...)
call lcm_initialize()
call lcm_advance(dt)
```

Fresh initialization is the default. A host can request restart initialization
through `lcm_config_t%start_from_restart`, which should be filled from the host
namelist. That path currently stops with a clear error until the standalone
restart lifecycle is implemented.

## `micro_set` Responsibility Map

| SAM-LCM `micro_set` responsibility | Standalone LCM location | Status |
| --- | --- | --- |
| Validate configuration consistency | `lcm_set_config` / future config validation work | Partly covered. Basic config validation exists in `lcm_set_config`. Broader physics-mode compatibility checks from SAM-LCM `micro_set` should not live in the initialization lifecycle unless they are required to allocate or initialize LCM-owned state. |
| Compute initial particle counts `n_init`, `n_total` | `lcm_initialize` | Covered as expected counts through `lcm_initial_particle_count_local` and `lcm_initial_particle_count_global`. Actual particle creation remains follow-up work because the new design stores particles by cell instead of in a global flat array. |
| Compute domain width `xwidth`, `ywidth` | `lcm_set_grid` | Covered. |
| Compute inverse grid spacing `idx`, `idy` | `lcm_set_grid` | Covered. |
| Compute grid-cell volume `dV` | `lcm_set_grid` | Covered. |
| Compute vertical bounds `ztop`, `zbottom` | `lcm_set_grid` | Covered. |
| Store LCM grid spacings `dx_LCM`, `dy_LCM` | `lcm_set_grid` through `grid_state%grid%dx`, `grid_state%grid%dy` | Covered by the host grid state, not by legacy scalar globals. |
| Compute local/global subdomain offsets | `lcm_set_grid` | Covered through `global_i_start`, `global_j_start`, `nx_sb`, and `ny_sb`. |
| Print SAM-LCM initialization summary | Not yet assigned | Not implemented. If needed, this should become standalone diagnostics/logging rather than SAM-style `PRINT*` blocks. |
| Create MPI datatype for legacy flat particle array | Future MPI particle exchange interface | Not implemented in #42. The standalone design uses cell-owned particles and will need explicit exchange buffers instead of the old flat particle datatype. |
| Fill grid side/center arrays `gpxs`, `gpys`, `gpzs`, `gpxc`, `gpyc`, `gpzc` | `lcm_set_grid` | Covered. |
| Allocate velocity and atmospheric working arrays | Future atmospheric/internal state refactor | Not implemented in #42. The standalone direction is to use attached host fields and only allocate LCM-owned derived state where needed. |
| Initialize output/history state | Future diagnostics/output work | Not implemented in #42. |
| Choose fresh initialization vs restart read | `lcm_config_t` / `lcm_initialize` | Covered for lifecycle selection. Fresh initialization is the default; `config_state%start_from_restart` selects restart mode and stops with a clear error until restart lifecycle work exists. |
| Set `flag_LCM_init` | `lcm_is_initialized` | Partly covered. `lcm_is_initialized` records that host setup and LCM-owned storage initialization are complete before `lcm_advance`. |

## Derived Grid/Helper Audit

SAM-LCM computes several grid helper values in `micro_set`. Standalone LCM
keeps these responsibilities in `lcm_set_grid`, because they are derived from
the host grid and should exist before `lcm_initialize`.

| SAM-LCM helper from `micro_set` | Standalone LCM state | Status |
| --- | --- | --- |
| `xwidth = nx_gl * dx` | `grid_state%xwidth` | Covered in `lcm_set_grid`. |
| `ywidth = ny_gl * dy` | `grid_state%ywidth` | Covered in `lcm_set_grid`. |
| `idx = 1.0 / dx` | `grid_state%idx` | Covered in `lcm_set_grid`. |
| `idy = 1.0 / dy` | `grid_state%idy` | Covered in `lcm_set_grid`. |
| `dV(k) = dx * dy * (zi(k+1)-zi(k))` | `grid_state%dV(:)` | Covered using host `dx`, `dy`, and `dz(:)`. |
| `dz_LCM(k) = zi(k+1)-zi(k)` | `grid_state%grid%dz(:)` | Covered by the host grid type; no duplicate `dz_LCM` scalar array is needed. |
| `ztop = zi(nz) - boundary_distance` | `grid_state%ztop` | Covered with the current standalone boundary offset. Boundary policy remains part of particle-boundary work. |
| `zbottom = boundary_distance` | `grid_state%zbottom` | Covered with the current standalone boundary offset. Boundary policy remains part of particle-boundary work. |
| `dx_LCM = dx` | `grid_state%grid%dx` | Covered by storing the host grid. |
| `dy_LCM = dy` | `grid_state%grid%dy` | Covered by storing the host grid. |
| `task_rank_to_index(rank, i_global_LCM, j_global_LCM)` | `grid_state%i_global_LCM`, `grid_state%j_global_LCM` | Covered from host-supplied local/global index metadata. |
| `nx_sb`, `ny_sb` | `grid_state%nx_sb`, `grid_state%ny_sb` | Covered from `global_i_start` and `global_j_start`. |
| `gpxs(0:nx+1)` | `grid_state%gpxs(:)` | Covered in `lcm_set_grid`. |
| `gpys(0:ny+1)` | `grid_state%gpys(:)` | Covered in `lcm_set_grid`. |
| `gpzs(1:nz)` | `grid_state%gpzs(:)` | Covered from host `z_face(:)`. |
| `gpxc(-1:nx+1)` | `grid_state%gpxc(:)` | Covered in `lcm_set_grid`. |
| `gpyc(-1:ny+1)` | `grid_state%gpyc(:)` | Covered in `lcm_set_grid`. |
| `gpzc(1:nzm)` | `grid_state%gpzc(:)` | Covered from host `z_center(:)`. |
| `x1u/x2u`, `y1u/y2u`, etc. | `grid_state%x1u`, `grid_state%x2u`, etc. | Covered as derived halo/field bounds. |

The audit does not identify additional grid-helper state that must be
initialized in `lcm_initialize`. The one area that still needs a later design
decision is boundary policy: SAM-LCM used `boundary_distance` when setting
`ztop` and `zbottom`; standalone LCM currently stores bounds in `grid_state`,
but the exact particle deletion/reflection behavior belongs to the particle
boundary issue.

## What #42 Implements

Issue #42 introduces the explicit lifecycle boundary:

- `lcm_initialize` validates that config, grid, and required atmospheric fields
  have been supplied.
- `lcm_initialize` allocates LCM-owned cell particle containers with
  `config_state%np_xy` as the initial per-cell capacity.
- `lcm_initialize` records expected local and global initial particle counts.
- `lcm_initialize` fails explicitly when
  `config_state%start_from_restart=.true.` because restart initialization is not
  implemented yet.
- `lcm_is_initialized` records that this setup has completed.
- `lcm_advance` rejects calls before `lcm_initialize`.
- `lcm_host_interface` exports `lcm_initialize`.

## What Remains Outside #42

The following pieces are intentionally left to focused follow-up work:

- actual initial particle population from config, replacing
  `micro_place_particle` and `micro_init_particles`
- broader SAM-LCM physics-mode compatibility checks from `micro_set`
- restart initialization and output/restart file design
- diagnostics and feedback accumulator initialization
- MPI particle exchange buffers and particle communication topology
- replacement of legacy atmospheric working arrays with `lcm_field_state`
