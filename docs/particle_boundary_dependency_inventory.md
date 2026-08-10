# Particle Boundary Dependency Inventory

This note inventories the dependencies of SAM-LCM's
`MICRO_LAGRANGE/micro_boundaries.f90` before implementing standalone
`apply_lcm_particle_boundaries()` (Issue #39).

## Scope and advance contract

The standalone routine applies **physical particle-boundary rules only** after
particle positions have been advanced.  It must not perform MPI communication
or assign the final local cell owner.  The intended sequence is:

```fortran
call advance_lcm_particles(time_step)
call apply_lcm_particle_boundaries()
call exchange_lcm_particles()
call update_lcm_particle_cell_ownership()
```

For DALES, MPI decomposition is horizontal.  Consequently, vertical
out-of-domain particles must be handled here before any cell-location or MPI
work can process them.

## SAM behavior inventory

`micro_boundaries` has four distinct behaviors:

| Behavior | SAM rule | Standalone status |
| --- | --- | --- |
| Vertical-domain test | A particle is inside when `zbottom <= z <= ztop`. | Required. |
| Default loss/reflection | A particle below `zbottom` with radius above 1 micrometre is removed; all other vertical escapes are reflected. | Required policy to preserve, using named fields. |
| Pi-chamber walls | Particles outside configured x/y chamber walls are absorbed or reflected. | Out of scope for the first DALES implementation. |
| Steady aerosol/chamber reinjection | A bottom-escaping large particle is repositioned randomly in an injection layer and its aerosol/water state is reconstructed. | Out of scope; requires atmospheric thermodynamics and a random-stream contract. |

SAM rebuilds a temporary flat `particle(:)` array while applying those rules.
Standalone LCM must instead traverse `cell_particles(:,:,:)` and use
`lcm_remove_particle` only for particles that are actually deleted.  A
reflected particle remains in its existing container until the later ownership
step moves it, if necessary.

## Dependency map

| SAM dependency | Purpose in `micro_boundaries` | Standalone equivalent / decision |
| --- | --- | --- |
| `particle(:)`, `n_particles`, `particle_temp`, `n_particles_alloc` | Flat storage and compaction after deletion or reinjection. | Replace with `cell_particles(:,:,:)` and `lcm_remove_particle`; no flat-array rewrite or allocation counter. |
| `particle%var(iZ)`, `ztop`, `zbottom` | Test and mirror vertical position. | `lcm_particle_t%z`, `grid_state%ztop`, and `grid_state%zbottom`. Already available. |
| `boundary_distance` | SAM uses a 0.001 m inset when deriving limits. | Already represented by `grid_state%zbottom = z_face(1)+0.001` and `grid_state%ztop = z_face(nz+1)-0.001`. Keep the policy centralized in grid setup. |
| `get_r_particle`, `particle%var(iM)`, `particle%ivar(iA)`, `rho_liq` | Determine whether a particle is larger than 1 micrometre at the bottom. | `lcm_particle_radius` uses named fields and preserves SAM's volume approximation. |
| `particle%var(iWsgs)` | On reflection, force subgrid vertical velocity away from the wall. | `lcm_particle_t%w_sgs`; required if reflection is included. The resolved `w` is not changed by SAM. |
| `np_deleted` | Number of particles removed in the call. | Add an LCM-owned per-step boundary diagnostic if needed; do not retain a SAM global. |
| `total_water_prec`, `dx_LCM`, `dy_LCM` | SAM records the change in total liquid-particle mass per horizontal area across its compaction pass. | Needs an explicit LCM diagnostics/feedback design. Do not update DALES or SAM globals from the boundary routine. |
| Pi-chamber wall limits and switches | `wall_x_*_pichamber`, `wall_y_*_pichamber`, `switch_walls_pichamber`, `switch_reflective_walls_pichamber`. | Do not add to the first DALES boundary configuration. A future generic wall-policy interface can model this without Pi-chamber names. |
| Reinjection switches and bounds | `switch_steady_aerosol`, `switch_steady_pichamber`, `inject_bottom`, `inject_top`. | Postpone. These are experiment-specific source/recycling policies, not a default atmospheric boundary. |
| `randomu(rank+1)` | Select reinjection height. | Postpone with reinjection; future use must go through standalone LCM RNG, never an MPI-rank-indexed SAM call. |
| `micro_field`, `tabs_LCM`, `pres_LCM`, `qsatw`, `r_equi`, `get_r_aero`, `rho_aero`, `pi` | Reconstruct mass and aerosol state for reinjected particles. | Postpone with reinjection; this is microphysical initialization, not a basic boundary rule. |
| `idx`, `idy`, `nx_sb`, `ny_sb`, `nx_LCM`, `ny_LCM`, `gpzs`, `search_sorted_array`, `RUN3D` | Locate a reinjected particle's cell. | Postpone with reinjection. Normal ownership belongs to `update_lcm_particle_cell_ownership()`. |
| `gamaz`, `dz_LCM`, and several legacy particle indices | Imported by SAM but unused in this routine's executable statements. | No standalone dependency. |
| `t_startf`, `t_stopf` | SAM performance timers. | No functional dependency; future profiling should use standalone instrumentation. |

## Required named particle fields

The initial vertical-boundary routine needs `z` for every particle.  It also
needs `w_sgs` to reproduce SAM reflection.  Exact preservation of SAM's
large-droplet bottom-loss criterion additionally needs a standalone physical
radius helper derived from `mass_water` and `weight` (and, for relevant
particle types, the existing particle material state).

SAM-LCM's `get_r_particle(n)` in `micro_vars.f90` is the reference for this
helper. It calculates aqueous volume from `iM / iA / rho_liq` and includes
aerosol volume only when `iN / iM > 0.1`; the source comment describes this as
an approximation whose radius error is small when aerosol mass is lower. The
standalone `lcm_particle_radius` helper preserves this same condition using
`mass_aero / mass_water > 0.1`.

## DALES integration decision

Issue #39 adopts the following initial DALES boundary-mode policy.  This is a
deliberate baseline, not an attempt to expose every legacy SAM-LCM mode.

| Mode | Decision | Reason |
| --- | --- | --- |
| Physical vertical bounds | Enabled. Use `grid_state%zbottom` and `grid_state%ztop`. | DALES has a vertically bounded physical domain and these limits are already derived from the host grid. |
| Bottom loss | Enabled for particles with physical radius greater than 1 micrometre. | Preserves the default SAM-LCM precipitation/removal distinction. |
| Bottom reflection | Enabled for particles at or below the 1-micrometre threshold. | Preserves the SAM-LCM aerosol-like particle behavior. |
| Top reflection | Enabled for every particle crossing `ztop`. | Preserves the SAM-LCM default behavior. |
| Pi-chamber x/y walls | Disabled and unsupported. | They are experiment-specific and do not describe a DALES atmospheric-domain boundary. |
| Steady aerosol/chamber reinjection | Disabled and unsupported. | It is a source/recycling policy requiring thermodynamic reconstruction and reproducible random sampling. |
| Horizontal periodic wrapping | Not handled here. | It belongs to the horizontal MPI-exchange design (Issue #40). |
| MPI send/receive | Not handled here. | This routine owns physical rules, not inter-rank transfer. |
| Final local cell ownership | Not handled here. | This follows particle exchange in Issue #41. |
| DALES precipitation/global diagnostics | Not updated directly. | LCM must first expose an LCM-owned diagnostic contract. |

The sole baseline boundary configuration field is
`lcm_config_t%bottom_removal_radius`, which defaults to `1.0e-6 m` and must be
positive. `ztop` and `zbottom` remain derived grid state. This keeps the
preserved SAM threshold named and testable without exposing experiment-specific
wall or reinjection settings.

## Vertical-grid boundary contract

`lcm_grid_t%z_face(:)` is the authoritative host-provided vertical geometry.
For `nz` local vertical cells, the physical domain is bounded by
`z_face(1)` and `z_face(nz + 1)`; vertical MPI decomposition is not assumed.
When `lcm_set_grid` initializes `grid_state`, it derives the particle-boundary
limits directly from those `real64` endpoints:

```fortran
zbottom = z_face(1)      + lcm_vertical_boundary_inset
ztop    = z_face(nz + 1) - lcm_vertical_boundary_inset
```

`lcm_vertical_boundary_inset` is the named standalone equivalent of SAM-LCM's
`boundary_distance` and is currently `0.001 m`. A particle exactly at either
derived limit is in-domain; only `z < zbottom` or `z > ztop` invokes a
vertical-boundary rule. The inset keeps boundary-processed particles inside
the raw face interval that the later cell-location routine uses.

## Cell-owned boundary mutation contract

SAM-LCM copies surviving particles into a temporary flat array and replaces
`particle(:)` after the loop. Standalone LCM does not rebuild a flat array.
Each particle stays in its current `cell_particles(i,j,k)` container while
boundary rules are evaluated:

- an in-domain or reflected particle is retained in the current container;
- a deleted particle is removed with `lcm_remove_particle(i,j,k,n)`;
- removal preserves the order of surviving particles by shifting later
  entries down; and
- the caller must **not** increment `n` after removal, because the next
  unprocessed particle has shifted into that index.

The intended traversal is therefore:

```fortran
n = 1
do while (n <= cell_particles(i,j,k)%count)
   if (particle_is_deleted) then
      call lcm_remove_particle(i, j, k, n)
   else
      n = n + 1
   endif
enddo
```

This is the standalone replacement for SAM's flat-array compaction. It does
not move reflected particles to their final cell and does not transfer them to
another rank; those responsibilities stay with the later ownership and MPI
steps.

## Implemented top-boundary rule

The standalone top rule preserves the SAM-LCM `micro_boundaries` behavior for
`z > ztop` using named fields:

```fortran
particle%z = ztop + (ztop - particle%z)
particle%w_sgs = -abs(particle%w_sgs)
```

It retains the reflected particle in its current cell container. The later
ownership step is responsible for any local cell reassignment, and MPI is not
involved. The bottom rule is intentionally not applied by this implementation
until its dedicated task is complete.

### Required observable behavior

For the initial DALES scope:

1. Use the already-derived `grid_state%zbottom` and `grid_state%ztop`.
2. Keep in-domain particles unchanged.
3. Preserve the SAM default distinction: delete bottom-escaping particles
   larger than 1 micrometre and reflect the remaining vertical escapes.
4. Mirror `z` about the violated limit and set `w_sgs` away from it on
   reflection.
5. Remove deleted particles in place from their current cell container.
6. Leave horizontal wrapping, MPI transfer, and final cell reassignment to
   their respective later routines.

The implementation tests for this policy must demonstrate that a particle
below the bottom is either removed or mirrored according to its radius, a
particle above the top is mirrored, the reflected `w_sgs` points away from the
violated boundary, and neither path performs MPI or ownership updates.

The following must not be silently included in this first implementation:

- Pi-chamber x/y walls;
- steady aerosol or chamber reinjection;
- SAM global precipitation accounting;
- local grid-cell reassignment; or
- MPI sends, receives, and buffers.

## Follow-up decisions exposed by this audit

- Define an LCM-owned boundary diagnostics state for deleted count and removed
  water mass, then let the DALES adapter consume it deliberately.
- Decide whether the 1-micrometre bottom criterion is the desired DALES
  default, and specify its treatment for ice particles.
- Add a generic boundary-policy configuration only when more than the default
  atmospheric behavior is required; do not expose Pi-chamber flags in the
  baseline API.
- Define reinjection separately as a particle-source/recycling feature, with
  its atmospheric-field and reproducibility requirements.
