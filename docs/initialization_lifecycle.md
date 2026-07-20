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
| Choose fresh initialization vs restart read | `lcm_initialize` | Partly covered. `lcm_initialize` supports fresh initialization and rejects `from_restart=.true.` with a clear error until restart lifecycle work exists. |
| Set `flag_LCM_init` | `lcm_is_initialized` | Partly covered. `lcm_is_initialized` records that host setup and LCM-owned storage initialization are complete before `lcm_advance`. |

## What #42 Implements

Issue #42 introduces the explicit lifecycle boundary:

- `lcm_initialize` validates that config, grid, and required atmospheric fields
  have been supplied.
- `lcm_initialize` allocates LCM-owned cell particle containers with
  `config_state%np_xy` as the initial per-cell capacity.
- `lcm_initialize` records expected local and global initial particle counts.
- `lcm_initialize(from_restart=.true.)` fails explicitly because restart
  initialization is not implemented yet.
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
