# Standalone LCM Restart Initialization Cycle

This note records how restart initialization works in SAM-LCM and how that
responsibility should connect to the standalone LCM initialization lifecycle.

## SAM-LCM Behavior

SAM-LCM handles the restart decision inside `micro_set`.

Near the end of `micro_set`, after configuration checks, grid helper setup, MPI
datatype setup, and working-array allocation, SAM-LCM decides whether LCM has
already started:

```fortran
if (nstep < nstep_LCM_init) then
   flag_LCM_init = .false.
else if (nstep >= nstep_LCM_init) then
   call micro_restart_read
   flag_LCM_init = .true.
endif
```

The restart read itself is implemented in `micro_restart_read`. That routine:

- opens the LCM restart file, either rank-separated or shared
- reads `n_particles`, `history_length`, `filename_LCM`, `n_particles_alloc`,
  and `n_ID`
- allocates the legacy flat `particle(:)` storage
- reads a temporary particle array from file
- copies temporary particle data into the live `particle(:)` array
- calls `micro_sort_by_grid_box`
- calls `micro_coupling`

The important design point is that SAM-LCM does not create a fresh initial
particle population during restart. It restores particle state from file, sorts
the restored particles into grid-box access structures, couples atmospheric
fields, and then marks LCM initialized.

Fresh particle placement happens later in the timestep path when
`flag_LCM_init` is false and `nstep == nstep_LCM_init`; that path calls
`micro_place_particle`.

## Standalone LCM Direction

Standalone LCM should keep the same scientific lifecycle but make the branch
explicit:

```fortran
call lcm_configure_runtime()
call validate_initialization_fields()

if (config_state%start_from_restart) then
   call initialize_lcm_from_restart()
else
   call initialize_lcm_from_config()
endif

lcm_is_initialized = .true.
```

Fresh initialization should:

- run after config/grid preparation and atmospheric field attachment
- allocate LCM-owned particle cell storage
- record expected local and global particle counts
- later call the initial particle population path from issue #43
- set `lcm_initialization_mode = "fresh"`

Restart initialization should eventually:

- read standalone LCM restart metadata and particle state
- reconstruct `cell_particles(i,j,k)` directly, or read a flat restart buffer and
  insert particles into the correct cell containers
- restore particle IDs, particle counts, and any diagnostic/restart metadata
- restore local cell ownership, equivalent in role to `micro_sort_by_grid_box`
- ensure atmospheric field pointers are attached and ready, equivalent in role
  to the `micro_coupling` dependency after SAM-LCM restart read
- set `lcm_initialization_mode = "restart"`

Until the restart file format and read path are implemented, the standalone
restart branch should fail explicitly. This prevents a restart run from
silently falling back to fresh particle initialization.

## Configuration Source

The restart selection should come from host configuration, not from an extra
argument to `lcm_initialize`.

For DALES, this means the namelist should eventually set a restart option, the
DALES adapter should pack that option into `lcm_config_t%start_from_restart`,
and `lcm_initialize()` should choose the fresh or restart path from
`config_state%start_from_restart`.
