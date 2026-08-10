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
| `np_deleted` | Number of particles removed in the call. | `boundary_diagnostics%deleted_particle_count`, reset for each boundary call. |
| `total_water_prec`, `dx_LCM`, `dy_LCM` | SAM records the change in total liquid-particle mass per horizontal area across its compaction pass. | `boundary_diagnostics%deleted_water_mass` retains raw local mass. A later host adapter owns reduction and flux conversion; the boundary routine does not update DALES or SAM globals. |
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

## Deferred SAM-LCM Pi-chamber wall behavior

SAM-LCM implements a separate, optional horizontal wall model in
`MICRO_LAGRANGE/micro_boundaries.f90`. It is enabled by
`switch_walls_pichamber`; its wall positions are configured with
`wall_x_w_pichamber`, `wall_x_e_pichamber`, `wall_y_s_pichamber`, and
`wall_y_n_pichamber`. The mode is inactive by default.

The wall test is applied only after SAM has established that a particle is
inside the normal vertical interval `zbottom <= z <= ztop`. A particle inside
the configured x/y rectangle survives unchanged. For a particle outside that
rectangle, SAM selects one of two behaviors using
`switch_reflective_walls_pichamber`:

| Wall policy | SAM-LCM behavior |
| --- | --- |
| Absorbing walls (the default) | Delete the particle while rebuilding the legacy flat `particle(:)` array. |
| Reflective walls | Keep the particle, mirror each coordinate that crossed a wall, and reverse the corresponding SGS velocity component. |

For example, a west-wall crossing becomes
`x = 2 * wall_x_w_pichamber - x` and `u_sgs = -u_sgs`; an east-wall crossing
uses the east-wall coordinate. South and north crossings use the equivalent
`y`/`v_sgs` formulas. A particle that crosses a corner is reflected in both
coordinates. SAM changes the SGS components only; it does not modify the
resolved horizontal velocity in this branch.

This is deliberately **out of scope for Issue #39**. Pi-chamber walls describe
a particular laboratory geometry, rather than the normal DALES atmospheric
horizontal domain. In the DALES integration, horizontal-domain crossings are
handled by the later periodic/MPI exchange design, not as particle loss at a
physical x/y wall. Introducing SAM's Pi-chamber switches directly would also
make the standalone LCM configuration experiment-specific.

If a future host requires physical horizontal walls, introduce a generic,
host-configured wall-policy interface with explicit planar limits and
absorb/reflect behavior. It can reproduce the SAM Pi-chamber case without
embedding Pi-chamber names or assumptions in the standalone LCM core.

## Deferred SAM-LCM steady aerosol/chamber reinjection

SAM-LCM provides two optional recycling modes:
`switch_steady_aerosol` and `switch_steady_pichamber`. Both are disabled by
default. Either switch changes the treatment of a particle that has crossed
the bottom and has radius greater than `1 micrometre`; without either switch,
SAM deletes that particle.

With either switch enabled, SAM does not simply move the same droplet upward.
It copies the particle, then performs the following reinitialization in
`micro_boundaries.f90`:

1. Draw a uniform value from `randomu(rank + 1)` and place the copy uniformly
   between `inject_bottom` and `inject_top`.
2. Set `w_sgs` to zero.
3. Determine the particle's host-grid location and obtain the local
   temperature, pressure, and water-vapor state.
4. Compute a bounded subsaturation (`-0.20 <= SS_gb <= -0.01`), obtain the
   equilibrium/aerosol radius, and rebuild the particle water mass from that
   radius and liquid-water density.
5. Convert the particle type to its non-ice sign, reset monomer count to one,
   and reset rime mass to zero.

Thus, SAM describes this mode as reinjecting absorbed droplets **as aerosols**:
it is a particle source/recycling and microphysical-initialization policy, not
a simple physical boundary condition. The implementation depends on SAM's
random stream, rank, host thermodynamic arrays, saturation calculation,
equilibrium-radius routine, aerosol radius routine, and legacy cell-location
logic. In the current SAM code, the vertical index lookup uses the source
particle's pre-reinjection `z` value even though the copied particle has
already received its new injection height; any future design should specify
the intended sampling location explicitly rather than inheriting this detail.

This behavior is deliberately **out of scope for Issue #39**. Adding it now
would require standalone contracts for reproducible random-number generation,
an injection-layer configuration, host atmospheric sampling, particle
reinitialization, diagnostics, and ownership after insertion. Those decisions
belong in a separate particle-source/recycling feature, after the basic
physical-boundary routine is stable.

A future implementation should expose a generic source policy rather than
SAM-specific `steady_aerosol` or `steady_pichamber` flags. Its interface must
state the injection distribution, RNG/restart guarantees, host fields sampled,
microphysical reset rules, mass/number diagnostics, and whether insertion
occurs before or after MPI exchange and local ownership assignment.

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

## Implemented bottom-boundary rule

For `z < zbottom`, standalone LCM preserves SAM-LCM's default split using the
named radius helper and `lcm_config_t%bottom_removal_radius`:

- if `lcm_particle_radius(particle, rho_aero) > bottom_removal_radius`, remove
  the particle from its current cell container; otherwise
- mirror the particle with `zbottom + (zbottom - z)` and set
  `w_sgs = abs(w_sgs)`.

The removal loop deliberately does not advance its particle index after a
deletion, so compacted consecutive out-of-domain particles are all processed.
This does not include reinjection, precipitation diagnostics, MPI exchange, or
final cell ownership updates.

## Implemented boundary-loss diagnostics

`lcm_boundary_diagnostics` owns the losses from the most recent
`lcm_apply_particle_boundaries()` call. The routine resets this state before it
applies physical rules, so the values are per-call rather than implicitly
cumulative or restart-persistent:

| Field | Meaning | Units |
| --- | --- | --- |
| `deleted_particle_count` | Number of removed computational particles. | count |
| `deleted_water_mass` | Sum of `mass_water` for removed particles. | kg |
| `deleted_aerosol_mass` | Sum of `mass_aero` for removed particles. | kg |

The current vertical policy only removes large particles crossing the bottom,
but the diagnostic names are intentionally general enough for a future
physical-wall policy. Values are local to the rank; they are deliberately not
divided by horizontal area or timestep and are not MPI-reduced. This preserves
the raw physical loss without assuming how a host interprets it. A later DALES
feedback adapter may reduce the values, normalize them, and expose an
appropriate precipitation or mass-budget diagnostic.

This differs from SAM's `total_water_prec` update, which calculates a
before/after flat-array mass difference divided by `dx_LCM * dy_LCM` and writes
directly to a SAM global. The standalone bottom-boundary integration test
verifies that a deleted particle contributes its water and aerosol mass while
reflected particles contribute neither.

## SAM reflection/deletion compatibility

The baseline rules intentionally preserve the executable comparisons in
SAM-LCM `micro_boundaries`:

| Condition | Standalone behavior | SAM behavior |
| --- | --- | --- |
| `z > ztop` | Mirror about `ztop`; set `w_sgs = -abs(w_sgs)`. | Same top-reflection formula. |
| `z < zbottom` and `radius > bottom_removal_radius` | Remove the particle. | SAM's `radius > 1.0E-6` deletion branch. |
| `z < zbottom` and `radius <= bottom_removal_radius` | Mirror about `zbottom`; set `w_sgs = abs(w_sgs)`. | SAM's remaining-particle reflection branch. |
| `z == zbottom` or `z == ztop` | Leave unchanged. | SAM treats both limits as in-domain. |

The strict greater-than comparison is significant: a particle exactly at the
configured removal radius reflects rather than being deleted. The minimal-host
bottom-boundary integration test covers this edge case alongside the smaller
reflected and larger deleted cases.

### Required observable behavior

For the initial DALES scope:

1. Use the already-derived `grid_state%zbottom` and `grid_state%ztop`.
2. Keep in-domain particles unchanged.
3. Preserve the SAM default distinction: delete bottom-escaping particles
   larger than 1 micrometre and reflect the remaining vertical escapes.
4. Mirror `z` about the violated limit and set `w_sgs` away from it on
   reflection.
5. Remove deleted particles in place from their current cell container.
6. Record each removed particle's count, water mass, and aerosol mass in the
   LCM-owned per-call boundary diagnostics.
7. Leave horizontal wrapping, MPI transfer, and final cell reassignment to
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

- Define the DALES feedback adapter that consumes local boundary diagnostics,
  including MPI reduction, area/timestep normalization, and restart/output
  semantics.
- Decide whether the 1-micrometre bottom criterion is the desired DALES
  default, and specify its treatment for ice particles.
- Add a generic boundary-policy configuration only when more than the default
  atmospheric behavior is required; do not expose Pi-chamber flags in the
  baseline API.
- Define reinjection separately as a particle-source/recycling feature, with
  its atmospheric-field and reproducibility requirements.
