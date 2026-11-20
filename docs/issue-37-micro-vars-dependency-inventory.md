# Issue 37: `micro_vars` dependency inventory for `micro_condensation`

This note records the first Issue 37 task: inventory the `micro_vars` dependencies used by the SAM-LCM `micro_condensation` routine. It is intended as a working reference for the standalone refactor; it does not change or prescribe condensation physics.

## Source inspected

- SAM-LCM routine: `/Users/mahnoosh/SAM6.10.10.LCM/SRC/MICRO_LAGRANGE/micro_cond.f90`, `micro_condensation`, lines 2–729.
- Its `USE micro_vars, ONLY:` list is at lines 8–27.
- Particle field index definitions are in `/Users/mahnoosh/SAM6.10.10.LCM/SRC/MICRO_LAGRANGE/micro_vars.f90`, around lines 394–454.

The inventory follows the module import list and direct references in this routine. The routine also imports other modules (`microphysics`, `micro_ds`, `params`, `vars`, `grid`, `rad`, and conditionally `chemistry`); those dependencies are outside this first task.

## Directly used `micro_vars` dependencies

### Particle storage and particle fields

| Dependency | Role in this routine | Standalone mapping to investigate |
|---|---|---|
| `n_start`, `n_end` | Define the flat particle-index range belonging to cell `(i,j,k)`; loops use the range throughout. It is also used to divide cell volume among particles in LEM mode. | Replace with `cell_particles(i,j,k)` iteration and define any per-particle/cell volume convention explicitly. |
| `particle` | Flat particle array; reads and mutates real (`var`) and integer (`ivar`) slots. | Typed particle state/accessors. |
| `iA` | Particle multiplicity/weighting factor; scales mass and diagnostics, and changes during nucleation, fragmentation, and sublimation. | Preserve as a typed weighting field with the existing count semantics. |
| `iM` | Particle total water mass; updated during liquid/ice growth and used in phase changes. | Typed total-water-mass field. |
| `iN` | Aerosol mass used by the film-forming-compound path (`switch_FFC_micro`). | Typed aerosol/composition field if that option remains supported. |
| `iT` | Particle temperature; updated for latent heating/cooling in LEM mode. | Typed particle temperature, subject to the standalone thermodynamic design. |
| `iq` | Particle water-vapor/humidity state used by LEM supersaturation and updated during growth. | Typed particle vapor state only if the standalone formulation retains particle-local vapor. |
| `ieta` | Stored absolute supersaturation from the previous step; used in LEM supersaturation evolution and rewritten after condensation. | Typed particle history/state field or an explicitly designed alternative. |
| `iZ` | Particle vertical coordinate used to derive particle temperature from stored static-energy-like state. | Location/trajectory state; likely belongs outside purely physical particle fields. |
| `iWsedi` | Particle terminal velocity, used for ventilation factors. | Motion property/accessor, with ownership coordinated with sedimentation. |
| `iAXS`, `iCXS` | Ice-crystal axis dimensions; used for capacitance/ventilation and updated during freezing and fragmentation. | Typed ice-shape fields or a phase-specific state representation. |
| `iDNS` | Apparent ice density; reset to ice density during freezing/fragmentation. | Typed ice property if needed by the preserved habit/shape model. |
| `iTYPE` | Encodes phase/sign and IN membership; controls liquid/ice paths and transitions. | Typed phase/type enum plus explicit IN membership semantics. |
| `iTFZ` | Particle freezing threshold used by immersion freezing. | Typed particle nucleation property. |
| `iMrime` | Rimed-mass state, reset during melting/freezing/sublimation transitions. | Typed particle state if retained by the standalone ice processes. |
| `iNmono` | Ice monomer count/state, reset to one during phase changes. | Typed ice state if retained. |
| `iAero_portions` | Per-mode aerosol composition fractions; selects aerosol mode when restoring multiplicity after complete sublimation. | Typed composition array or a replacement composition abstraction. |
| `in_startind` | Offset locating per-mode ice-nucleus composition slots in the particle array; used to select the dominant IN mode after sublimation. | Replace slot arithmetic with named composition fields/accessors. |

### Atmospheric state, geometry, and timestep

| Dependency | Role in this routine | Standalone mapping to investigate |
|---|---|---|
| `tabs_LCM(i,j,k)` | Cell air temperature for saturation, diffusional growth, radiation, nucleation/freezing/melting; in SAM it also receives latent heat feedback after the process. | Read-only host atmospheric input for this sprint. Do not update host temperature/static energy in the one-way implementation; defer that feedback path to a later sprint. |
| `pres_LCM(k)` | Pressure for saturation and supersaturation calculations. | Host atmospheric accessor. |
| `rho_air(k)` | Air density used in cell mixing-ratio conversions and latent heating. | Host atmospheric accessor. |
| `dV(k)` | Cell volume used to convert particle mass/multiplicity to local water amounts and supersaturation terms. | LCM grid/cell geometry input. |
| `dtn_LCM` | Timestep used in analytical supersaturation evolution, liquid/ice growth, and homogeneous-freezing probability. | Public `time_step` argument. |

### Configuration switches and thresholds

| Dependency | Role in this routine |
|---|---|
| `switch_LEM_micro` | Selects particle-local thermodynamic/supersaturation treatment and latent particle temperature/vapor updates versus a cell-mean path. |
| `switch_analytical_micro` | Selects analytical versus Euler liquid condensational growth. |
| `switch_koehler_micro` | Enables Köhler curvature/solute effects. |
| `switch_FFC_micro` | Enables film-forming-compound correction to diffusional growth. |
| `switch_radiation_micro` | Enables radiative enhancement of liquid droplet growth. |
| `switch_ice` | Enables ice supersaturation/nucleation processes. |
| `start_ice` | First timestep at which the ice nucleation block may run. |
| `switch_habit_ice` | Enables/controls ice habit evolution in the growth path. |
| `switch_ds`, `ds_param` | Select droplet-shattering model and its parameterization during freezing. |
| `nstep_cut_supersat`, `supersat_threshold` | Bound liquid supersaturation for initial configured timesteps. |

### Physical constants and helper procedures

| Dependency | Role in this routine |
|---|---|
| `rho_liq`, `rho_ice` | Liquid/ice density in radius/mass conversions and phase transitions. |
| `sigma` | Surface tension in Köhler curvature. |
| `molecular_weight_water`, `rho_aero`, `vanthoff_aero`, `molecular_weight_aero` | Imported but not directly referenced in this routine body. |
| `N_Avogadro` | Converts FFC surface density to molecules per droplet. |
| `k_ice`, `A_ice`, `B_ice` | Imported but not directly referenced in this routine body. |
| `nv`, `boltzk`, `planckh`, `av`, `bv`, `dw`, `srfthick` | Homogeneous-freezing rate calculation. |
| `habitparam` | Selects ice habit handling, including DISKICE setup. |
| `get_bfactor`, `get_r_aero` | Helper functions used for aerosol/Köhler terms and aerosol radius. |
| `np_xy`, `n_ice`, `frac_ice`, `unit_aero`, `n_aero` | Reconstruct particle weighting after complete ice sublimation. |
| `ds_param` | Selects the droplet-shattering helper/parameterization. |

### LCM profile accumulators

These are per-level process accumulators updated by `micro_condensation`; they are not particle fields:

`cnd_prof`, `evp_prof`, `act_prof`, `dea_prof`, `imm_prof`, `hom_prof`, `nuc_prof`, `nml_prof`, `sub_prof`, `dep_prof`, `frz_prof`, `mel_prof`, `nsub_prof`, and `ds_prof`.

They capture condensation/evaporation, activation/deactivation, nucleation/freezing, melting, deposition/sublimation, and droplet-shattering counts or masses. Their reset/ownership lifecycle should be mapped when doing the separate `microphysics`/`micro_field` dependency inventory.

## Imported but not directly referenced

The following names occur in the routine's `USE micro_vars, ONLY:` list but have no direct reference in the routine body as inspected:

`rho_aero`, `molecular_weight_water`, `vanthoff_aero`, `qnmin_LCM`, `k_ice`, `A_ice`, `B_ice`, `hzn_prof`, `molecular_weight_aero`, `dz_LCM`, and `iID`.

Keep this list as a cleanup candidate only. Before removing imports, check compiler/preprocessor configurations and whether a dependency is intentionally retained for conditional builds or nearby refactors.

## Refactor-relevant observations

- This routine combines particle state mutation, cell atmospheric state, configuration, grid geometry, and diagnostics in one 728-line procedure.
- Particle phase and composition are encoded through integer slot indices and magic values (`iTYPE` sign and magnitude); these need an explicit typed representation before changing behavior.
- The LEM path uses per-particle vapor/supersaturation state and particle temperature; the non-LEM path uses cell-mean atmospheric state. The two paths should be preserved distinctly during the port.
- The SAM reference writes `tabs_LCM` directly after accumulating latent heat. That host mutation is intentionally outside the current one-way scope; the later two-way design must define its energy feedback contract before adding it.
- The full routine has conditional `SAMCHEM` arguments and branches, even though those names do not come from `micro_vars`. Handling those hooks is a separate Issue 37 task.

## Scope boundary after task 1

At the time of the first inventory, this document covered only `micro_vars`; the `microphysics`/`micro_field` inventory below was added for Issue 37 task 2. Other module dependencies, dependency classification, implementation, and tests remain separate work.

## `microphysics` and `micro_field` dependencies (Issue 37 task 2)

Source inspected:

- `micro_condensation` imports these names in `micro_cond.f90`, lines 29–35.
- The `micro_field` allocation and diagnostic index definitions are in `microphysics.f90`, lines 50–156.
- `micro_field` is an allocatable 4-D array; its fourth dimension selects prognostic fields and diagnostics. `nmicro_fields` and `nmicro_fields_diag` determine the allocated fourth-dimension extent.

### Prognostic and forcing fields read or written

| Index | Declared meaning | Use in `micro_condensation` |
|---|---|---|
| `imv` | Water vapor mixing ratio | Read to calculate cell-mean supersaturation/forcing; reduced by net liquid plus ice growth at the end of the routine. The optional SAMCHEM call also receives this field. |
| `ilastT` | Advected/previous temperature used for forcing | Read with `ilastp` to calculate the previous-step saturation state; passed to optional chemistry; overwritten with current `tabs_LCM` at the end. |
| `ilastp` | Advected/previous pressure used for forcing | Read with `ilastT` to calculate the previous-step saturation state; passed to optional chemistry; overwritten with current `pres_LCM` at the end. |
| `iss` | Current absolute supersaturation (`qv - qsat`) | Written after vapor and temperature feedback are updated. |
| `isLES` | Relative/LES supersaturation diagnostic | Written as `iss / qsatw(...)` after condensation/phase change. |

These fields are part of the cell atmospheric/forcing state and SAM diagnostic array contract. The standalone one-way path should read host atmospheric inputs and maintain LCM-owned diagnostics; it should not copy SAM's direct host vapor/temperature updates. Two-way feedback is deferred.

### Process diagnostic slots updated

The routine adds weighted particle counts or mass increments to these `micro_field(i,j,k,index)` entries:

| Index | Declared meaning in `microphysics.f90` | Update in `micro_condensation` |
|---|---|---|
| `icon` | Condensation rate | Adds liquid mass growth (`dm * multiplicity`) when radius grows. |
| `ievp` | Evaporation rate | Adds liquid mass change when radius shrinks (the signed `dm` is retained). |
| `idep` | Deposition rate | Adds ice mass growth (`dm * multiplicity`) for positive ice mass change. |
| `isub` | Sublimation rate | Adds ice mass change for negative `dm`. |
| `iACT` | Activation rate | Adds multiplicity when radius crosses above `r_act`. |
| `iDEA` | Deactivation rate | Adds multiplicity when radius crosses back below `r_act`. |
| `ids` | Droplet shattering rate | Adds `prior multiplicity * fragment count` when shattering creates fragments. |
| `inuc` | Total ice nucleation rate | Adds multiplicity for immersion and homogeneous nucleation. |
| `iimm` | Immersion freezing rate | Adds multiplicity for heterogeneous immersion freezing. |
| `ihom` | Homogeneous nucleation rate | Adds multiplicity for homogeneous freezing. |
| `ifrz` | Total freezing rate (including riming in the model-wide definition) | Adds particle mass for immersion and homogeneous freezing events. |
| `imel` | Total melting rate | Adds particle mass when an ice particle melts. |

The corresponding profile accumulators in `micro_vars` are updated in parallel (listed above). The array indices and metadata (`mkname`, `mklongname`, `mkunits`, output scaling) are defined in `microphysics.f90`; preserve the intended diagnostics and units, but do not assume the raw additions alone are the final reported rate because the surrounding model controls reset, scaling, averaging, and output.

### Optional chemistry fields (`SAMCHEM` only)

When compiled with `SAMCHEM`, the routine initializes `ipH` and `idelR` to zero for the current cell, then passes these fields by reference to `chem_proc` along with `iSbyO3t`, `iSbyO3m`, `iSbyH2O2t`, and `iSbyH2O2m`. These are optional chemistry coupling/diagnostic state and need an explicit standalone support decision; they are not needed by the non-chemistry build.

### Imported name not directly used

`ihzn` (haze/deposition nucleation rate) is in the `microphysics` import list but is not directly referenced in the `micro_condensation` body inspected. It is a cleanup candidate, subject to the same conditional-build check as unused `micro_vars` imports.

### Lifecycle/design implications

- `micro_field` mixes prognostic atmospheric values, previous-step forcing state, and process diagnostics in one SAM array. A standalone design should separate these responsibilities rather than recreate a monolithic array by default.
- The process diagnostics are incremented inside the particle process. The parent `micro_proc` resets `iACT` and `iDEA` each active call (around lines 1606–1610); the profile accumulators are reset in that orchestration path (around lines 1613–1642). Determine reset and accumulation scope before transferring these diagnostics.
- In SAM, `micro_condensation` writes `imv`, then updates `tabs_LCM` using the summed liquid, ice, and fusion latent-heat terms, then recomputes `iss` and `isLES`. This is reference behavior only; the current one-way standalone implementation will not update host `qv` or temperature/static energy.
- Under `SAMCHEM`, `chem_proc` is invoked before the radius-growth loops and receives atmospheric/forcing values and mutable chemistry diagnostics. This branch should remain a separately chosen optional integration, rather than silently dropping state.

The inventory tasks are complete. Dependency classification is recorded below; implementation remains separate work.

## Candidate ownership map for `lcm_particle_t` (Issue 37 design task)

This map compares the SAM fields touched by condensation with the existing standalone types in `src/lcm_types.f90`. It is a design reference for a future `lcm_condensation` module; it does not implement the port or change the legacy modules.

### Per-particle state

| SAM field/dependency | Meaning in condensation | Existing standalone representation | Candidate ownership / decision |
|---|---|---|---|
| `ivar(iA)` | Super-particle multiplicity; scales represented mass/number and is sometimes changed by freezing, fragmentation, or sublimation. | `lcm_particle_t%weight` (`integer(int64)`). | Keep on `lcm_particle_t`. This is computational weight rather than a material property, but it is intrinsic to the super-particle state. |
| `var(iM)` | Total water mass represented by the super-particle. | `%mass_water`. | Keep on `lcm_particle_t`; this is the primary evolving condensate/ice mass. |
| `var(iN)` | Total aerosol mass, including FFC-related use. | `%mass_aero`. | Keep on `lcm_particle_t`; needed for aerosol volume and optional FFC behavior. |
| `var(iT)` | Particle temperature in the LEM path; latent heat changes it. | `%temperature`. | Keep as particle-carried thermodynamic state for LEM. It is not a material property; its exact relationship to host static energy needs to be maintained deliberately. |
| `var(iq)` | Particle-local vapor/humidity in the LEM path; reduced by condensation/deposition. | `%qv`. | Keep as particle-carried state only for the LEM formulation. In non-LEM mode the routine uses cell vapor instead. The standalone interface must make this mode distinction explicit. |
| `var(ieta)` | Previous absolute supersaturation in the LEM path. | `%eta`. | Keep as per-particle integration/history state if the LEM path is preserved; do not treat it as a material property. |
| `ivar(iTYPE)` | SAM encodes current phase in the sign (`>0` liquid, `<0` ice) and IN membership in the magnitude (`1` without IN, `2` with IN); freezing and melting change the sign while retaining the magnitude. | `%particle_type` currently has values `1` (liquid/no IN) and `2` (liquid/with IN) assigned during initialization and source creation. The particle also has `%in_surface_area(4)` for four IN modes and `%aerosol_portions(4)` for aerosol composition. | Follow SAM's four-state encoding in `%particle_type`, using named constants for liquid/ice with/without IN. Keep detailed IN-mode composition in its existing fields. |
| `var(iTFZ)` | Per-particle freezing threshold used by immersion freezing. | `%freezing_temperature`. | Keep on `lcm_particle_t`; it is particle-associated nucleation state. |
| `var(iAXS)`, `var(iCXS)` | Ice crystal semi-axes; updated by ice growth and phase changes. | `%axis_a`, `%axis_c`. | Keep on `lcm_particle_t` for the habit/shape formulation. |
| `var(iDNS)` | Apparent ice density; updated during nucleation and growth. | `%density`. | Keep on `lcm_particle_t` while the habit/ice formulation uses it. |
| `var(iMrime)` | Rimed water mass; changed during ice sublimation and reset on phase transitions. | `%rimed_mass`. | Keep on `lcm_particle_t` if the standalone ice representation retains riming. |
| `ivar(iNmono)` | Number of ice monomers; reset during phase changes. | `%n_monomer` (`integer(int64)`). | Keep on `lcm_particle_t` if used by standalone ice aggregation/shape physics. |
| `var(iSC:iSP)` | Surface area for four ice-nucleating aerosol modes; used by nucleation/composition logic and related ice physics. | `%in_surface_area(4)`. | Keep as per-particle composition/nucleation state. |
| `var(iAero_portions)` | Aerosol mass fractions by mode; used to restore weighting after complete sublimation and in aerosol-radius calculations. | `%aerosol_portions(4)`. | Keep on `lcm_particle_t` as composition state. |

The existing standalone type already has named fields for these dependencies. For Issue 37, `%particle_type` will retain SAM's combined sign/magnitude phase and IN-membership encoding, expressed through named constants. No additional particle field is required for phase.

### Dependencies that should not be classified as physical particle properties

| SAM dependency | Proper owner or access path | Reason |
|---|---|---|
| `var(iX:iZ)` and velocities, including `var(iWsedi)` | Particle motion/location state; currently carried by `lcm_particle_t` (`x/y/z`, velocity components, `w_sedi`) for the active advection path. | `iZ` is used to convert particle temperature to a height-adjusted temperature; `iWsedi` is read for ventilation. These are needed by the algorithm but describe position/motion, not the particle's physical composition. Keep their ownership distinct in the design. |
| `n_start`, `n_end` | `cell_particles(i,j,k)` and its cell-owned particle list. | These are legacy flat-array bounds, not per-particle properties. |
| `tabs_LCM`, `pres_LCM`, `rho_air`, `micro_field(...,imv)` | Read host atmospheric fields through `lcm_field_state` (`temperature` or `static_energy`, `qv`, `pressure`, `density`). | For the current one-way coupling scope, these are read-only inputs to LCM. LCM evolves particle state and its own diagnostics without modifying host vapor or energy. Design of two-way host feedback is deferred to a later sprint. |
| `dV(k)`, grid coordinates and dimensions | `lcm_grid_state` geometry (`dV`, cell centers/faces, spacing, grid extents). | These belong to grid/cell state and are shared across particles in the cell. |
| `dtn_LCM` | `time_step` argument passed through `lcm_advance` to `apply_lcm_microphysics`. | This is call/timestep context, not particle state. |
| `switch_*`, thresholds, aerosol/ice-mode parameters, habit and shattering selection | LCM configuration (`lcm_config_state` / `lcm_config_t`) or immutable physical parameters (`lcm_parameters`). | These configure the process or define constants; they are not properties of a single particle. |
| `micro_field` slots, profile accumulators | LCM-owned process diagnostics, with reset and reduction lifecycle defined outside the particle. | They aggregate over a cell or vertical profile and should not be stored on every particle. Host feedback accumulators are deferred with two-way coupling. |
| `radius`, ventilation factors, supersaturation arrays, temporary growth increments | Local scratch/derived values in `lcm_condensation` or helpers. | They are recomputable from particle, atmospheric, and configuration state for the current call; persistent storage is not justified by this inventory. |
| `rho_liq`, `rho_ice`, surface tension, gas constants, latent heats, Avogadro/Boltzmann/Planck constants | Shared physical parameter module (`lcm_parameters`) or process-specific parameter module. | They are model constants shared by many particles. |
| `iID` | Particle identity metadata (`lcm_particle_t%id`), though condensation does not directly use it. | It is identity, not a material property or a condensation input. |

### Existing standalone interfaces relevant to the map

- `lcm_particle_t` already contains particle identity, particle-carried thermodynamic state, water/aerosol mass, ice shape and nucleation properties, composition, motion, and age. The current type intentionally has a broader scope than physical properties alone.
- `lcm_field_state_t` provides attached host pointers for temperature, static energy, `qv`, pressure, and density (plus velocity and cloud liquid). The current `lcm_condensation` design should read required inputs through this state; host feedback is deferred.
- `lcm_grid_state` owns cell geometry including `dV`; `cell_particles(i,j,k)` owns cell membership and particle iteration.
- `lcm_config_state` holds the relevant process switches and mode parameters. `lcm_parameters` holds shared thermodynamic constants.
- No standalone microphysics diagnostic state type was found in the inspected `src` modules. The process-diagnostic ownership/reset/reduction interface still needs design before implementation.

### Particle phase and ice-nucleating identity map (Issue 37 design task)

The code inspection resolves the particle-state mapping needed before implementing condensation:

| Concept | SAM-LCM representation | Existing standalone representation | Design decision for `lcm_condensation` |
|---|---|---|---|
| Current phase | Sign of `ivar(iTYPE)`: positive is liquid, negative is ice. Freezing/melting changes the sign. | `%particle_type` encodes liquid with positive values and ice with negative values. | Follow the SAM phase encoding in `%particle_type`; use named constants rather than bare signed integers. |
| IN membership | Magnitude of `iTYPE` is `1` without IN and `2` with IN; detailed mode state is also stored in mode-specific aerosol/IN fields. | `%particle_type` magnitude records with/without IN; `%in_surface_area(4)` stores four IN-mode surface areas and `%aerosol_portions(4)` stores aerosol mode fractions. | Preserve the SAM combined state encoding and existing detailed mode state. Do not add a separate boolean or phase field. |
| Former initialization/source category | Combined with phase and IN status in `iTYPE`; this is current state, not historical particle provenance. | Existing `%particle_type` values are assigned at initialization/source creation as liquid without IN (`1`) or liquid with IN (`2`). | Use `%particle_type` as current SAM-compatible phase/IN state. `micro_condensation` does not require a separate historical source-category field; detailed IN mode remains in `%in_surface_area(4)`. |

This matches SAM's compact state representation while making each legal value explicit by name. It preserves IN membership when phase changes by switching between the corresponding liquid and ice constants. The detailed nucleating mode remains in `%in_surface_area(4)`. SAM also has a special complete-sublimation path that restores aerosol multiplicity/composition; that behavior must be preserved separately when ported.

Evidence inspected: `src/lcm_types.f90` (`lcm_particle_t` fields), `src/lcm_particle_state.f90` (initialization and IN-mode population), `src/lcm_particle_source.f90` (source categories), `src/micro_vars.f90` (`iTYPE` definition), and `src/micro_cond.f90` (freezing, melting, and complete-sublimation transitions).

### SAM-LCM phase/IN transition map for `micro_condensation` (Issue 37 follow-on task)

The combined `%particle_type` uses named values for SAM's four legal states: `liquid_no_in=1`, `liquid_with_in=2`, `ice_no_in=-1`, and `ice_with_in=-2`. The table below traces only phase-changing paths inside SAM `micro_condensation`; collision and boundary routines also touch `iTYPE`, but are outside this Issue 37 condensation slice.

| Process and condition in SAM `micro_condensation` | SAM state transition | Standalone mapping and adjacent state updates to preserve |
|---|---|---|
| Immersion freezing with IN: liquid particle, `iTYPE == 2`, positive mass, temperature at/below its freezing threshold, positive liquid supersaturation; ice enabled and start step reached (`micro_cond.f90`, around line 549). | Sets `iTYPE = -2`: liquid-with-IN → ice-with-IN. | Map `liquid_with_in` → `ice_with_in`. Before the assignment, the selected shattering parameterization may alter multiplicity and water mass; freezing also initializes ice axes/density, monomer count, and rimed mass, applies the LEM latent-heat update, and increments nucleation/freezing diagnostics. |
| Homogeneous freezing: liquid particle, sufficiently cold cell, positive water mass, supersaturation threshold met, and random probability test succeeds (`micro_cond.f90`, around line 641). | Sets `iTYPE = -ABS(iTYPE)`: `+1 → -1` (liquid without IN → ice without IN), or `+2 → -2` (liquid with IN → ice with IN). | Map each liquid constant to its matching ice constant. Preserve the magnitude/IN membership. Initialize ice axes/density, reset monomer count and rimed mass, apply the LEM latent-heat update, and increment homogeneous nucleation/freezing diagnostics. |
| Instantaneous melting: ice particle and cell temperature above 273.15 K (`micro_cond.f90`, around line 664). | Sets `iTYPE = ABS(iTYPE)`: `-1 → +1` or `-2 → +2`. | Map each ice constant to its matching liquid constant. Preserve IN membership; reset monomer count and rimed mass, apply the LEM fusion-energy update, and increment melting diagnostics. Water mass and particle weight are not reset by this block. |
| Complete sublimation: ice mass per represented particle at/below the tiny-mass threshold (`micro_cond.f90`, around line 683). | Sets `iTYPE = ABS(iTYPE)`, resets monomer/rimed-mass state, then recalculates multiplicity. If the resulting type is `2`, weight is restored from the dominant IN mode; otherwise it is restored from the dominant aerosol portion. The configured `1/kg` unit applies an air-density scale. | Map `-1 → +1` or `-2 → +2` first, then preserve the matching legacy weight-restoration branch, unit conversion, LEM temperature update, and sublimation diagnostics. Do not treat this as melting: it also changes multiplicity. |
| Particle initialization and source creation (context for valid starting values) | Aerosol initialization uses `iTYPE=1` without IN and `iTYPE=2` with IN; both are liquid. | Initialize as `liquid_no_in` or `liquid_with_in`. The ice constants are reached through freezing, not by merely adding IN material. |

Agreed SAM-compatible contract: `%particle_type` uses the four named states declared in `lcm_types`; its default and current creation paths represent liquid with or without IN. Freezing and melting switch between matching named states, while complete sublimation additionally restores multiplicity. The constants and initial state values are implemented; the process transitions remain for the `lcm_condensation` implementation. The existing MPI integer field already transports `%particle_type`, so no MPI schema change is needed.

Scope note: this table is for the `micro_condensation` counterpart only. SAM collision and boundary routines also read or change `iTYPE`; those process ports should preserve the same state encoding when their respective issues are addressed, but are not part of this transition map.

### Decisions to carry into `lcm_condensation` design

1. Follow SAM's combined phase/IN state in `%particle_type` using named constants for the four legal values; preserve detailed IN identity and mode in `%in_surface_area(4)` and related composition fields.
2. Keep LEM particle-local thermodynamic/history fields (`temperature`, `qv`, `eta`) distinct from non-LEM cell-mean host fields.
3. For the current one-way coupling scope, read host thermodynamic fields and do not apply latent heating or vapor depletion to the host. LCM particle-local state and LCM-owned diagnostics may still evolve. Defer the two-way feedback contract, conservation accounting, and host field updates to a later sprint.
4. Keep location, velocity, timestep, grid geometry, configuration, constants, and aggregate diagnostics out of the particle's physical-property set, even where the existing `lcm_particle_t` currently carries motion alongside microphysical state.
5. Preserve the existing fields already present in `lcm_particle_t` unless later source tracing demonstrates a required state is missing; this task is classification, not a schema edit.

This completes the particle-property ownership map. It is a proposed design map for review and does not implement `lcm_condensation` or modify any Fortran module.

## Motion and location ownership map (Issue 37 follow-on task)

This section records the motion/location dependencies used by `micro_condensation` and distinguishes them from physical particle state. The current `lcm_particle_t` stores both categories in one Fortran derived type, so this is a logical ownership distinction; it does not require splitting the type during this task.

| SAM field/dependency | How condensation uses it | Existing standalone representation | Ownership classification |
|---|---|---|---|
| `var(iZ)` | Used with particle temperature to form a height-adjusted particle temperature in the LEM path (`particle_tabs = particle%temperature - particle%z * ggr / cp`). | `lcm_particle_t%z`. | Particle location. Keep distinct from physical/thermodynamic state in the design. |
| `var(iWsedi)` | Read for liquid ventilation and passed to the ice-ventilation calculation. `micro_condensation` does not compute this velocity. | `lcm_particle_t%w_sedi`. | Particle motion/diagnostic velocity supplied by the sedimentation-velocity stage. Treat as read-only input to condensation. |
| `var(iX)`, `var(iY)` | Not directly read by `micro_condensation`; used by other motion, boundaries, and ownership code. | `lcm_particle_t%x`, `%y`. | Particle position, not physical condensate state. |
| `var(iU:iWsgs)` | Not directly read by `micro_condensation`; used by the particle-advance path. | `%u`, `%v`, `%w`, `%u_sgs`, `%v_sgs`, `%w_sgs`. | Particle motion state. Keep available to the advance path without treating it as condensation physics state. |
| `n_start`, `n_end` | Bound the particle loops and are used to derive the LEM per-particle cell-volume allocation. | `cell_particles(i,j,k)%particles` and `%count`. | Cell membership/storage metadata, not particle properties. Any replacement for the LEM volume allocation must be resolved separately; do not infer it from `cell%count` without checking the intended physical volume semantics. |
| `dV(k)` and cell-center height | Used in cell-to-mixing-ratio conversions and temperature/saturation calculations. | `lcm_grid_state%dV`, `%gpzc` and grid geometry. | Grid/cell geometry, shared by all particles in a cell. |

### One-way coupling scope for the current sprint

- `lcm_condensation` may read host `temperature` or `static_energy`, `qv`, `pressure`, and `density` through the attached field state; the specific thermodynamic reconstruction must match the existing standalone field conventions.
- Host fields remain inputs. Condensation changes each LCM particle's mass/phase and any retained LEM particle-local thermodynamic state, and accumulates LCM-owned diagnostics.
- The current one-way implementation must not subtract condensate from host `qv` or add latent heating to host static energy/temperature. Those two-way feedback rules, including conservation and diagnostic refresh, are deferred to a later sprint.
- `w_sedi` is consumed by ventilation as an existing motion input; condensation should not take ownership of its calculation.
- Motion/location fields already share `lcm_particle_t` with physical fields for transport and exchange. This task classifies ownership only; any physical split into nested types or a separate motion component is optional future design and should preserve the current exchange/advance interfaces.

This completes the motion/location classification for the current one-way coupling scope. That classification required no changes to motion/location fields; the particle phase/IN encoding is covered separately above.

## LCM-owned diagnostic classification (Issue 37 follow-on task)

For the current one-way scope, these are LCM-owned process outputs. They describe what happened to LCM particles; they are not tendencies to apply to host `qv`, temperature, or static energy. The SAM reference accumulates raw process amounts/counts in `micro_field`, then normalizes many rate fields later in `micro_diagnose` using cell volume, air density, and timestep.

### Candidate cell-level process diagnostics

| SAM diagnostic | Raw quantity accumulated by `micro_condensation` | Suggested LCM diagnostic meaning | Note |
|---|---|---|---|
| `icon` | Weighted positive liquid mass change (`dm * weight`). | Liquid condensation mass increment. | Count only growth events; SAM later reports a normalized rate. |
| `ievp` | Weighted signed liquid mass change when radius decreases. | Liquid evaporation mass increment, retaining the SAM sign convention or documenting a positive-loss convention. | Choose and document the sign convention explicitly. |
| `idep` | Weighted positive ice mass change. | Ice deposition mass increment. | SAM later reports a normalized rate. |
| `isub` | Weighted signed ice mass change for sublimation. | Ice sublimation mass increment, retaining the SAM sign convention or documenting a positive-loss convention. | Choose and document the sign convention explicitly. |
| `iACT` | Super-particle weight crossing above the activation radius. | Activated represented-particle count/weight. | SAM normalizes by cell air mass and timestep for rate output. |
| `iDEA` | Super-particle weight crossing back below the activation radius. | Deactivated represented-particle count/weight. | Same normalization caveat as activation. |
| `inuc` | Weight of particles nucleated by immersion or homogeneous freezing. | Total ice-nucleation represented count/weight. | The sum of the two mechanism-specific diagnostics below. |
| `iimm` | Weight of particles undergoing immersion freezing. | Immersion-freezing represented count/weight. | Mechanism-specific count diagnostic. |
| `ihom` | Weight of particles undergoing homogeneous freezing. | Homogeneous-freezing represented count/weight. | Mechanism-specific count diagnostic. |
| `ifrz` | Mass of particles undergoing the two freezing transitions. | Freezing event mass. | This is converted phase mass, not a host water-vapor sink in one-way coupling. |
| `imel` | Mass of particles undergoing melting. | Melting event mass. | This is converted phase mass, not host feedback. |
| `ids` | Prior weight multiplied by number of fragments produced by droplet shattering. | Droplet-shattering fragment count/weight. | SAM labels and normalizes this as a rate. |

These raw accumulations should remain distinguishable from reported rates. In SAM, the output path divides the listed process fields by cell volume, air density, and `dtn_LCM` before reporting (see `microphysics.f90`, around lines 2500–2522). The diagnostic implementation should decide whether its stored values are per-step increments, normalized cell rates, or both, and apply conversion exactly once.

### Vertical profile accumulators

`micro_condensation` also accumulates these per-level `micro_vars` profiles: `cnd_prof`, `evp_prof`, `act_prof`, `dea_prof`, `dep_prof`, `sub_prof`, `ds_prof`, `nuc_prof`, `imm_prof`, `hom_prof`, `frz_prof`, `mel_prof`, `nml_prof`, and `nsub_prof`. `nml_prof` is the represented count/weight of particles that melt; `nsub_prof` is the represented count/weight of ice particles that fully sublimate.

They are horizontal/profile summaries of overlapping process information, not independent particle properties. For standalone LCM, prefer deriving profile summaries by reducing cell-level diagnostics unless a consumer requires a separate vertically accumulated output. The SAM parent `micro_proc` zeros these profiles before its per-step work (around lines 1613–1642); the standalone reset scope should be explicit, likely once per `apply_lcm_microphysics` call, with longer output-window accumulation handled separately.

### State diagnostics and deferred feedback

| SAM value | Classification for one-way LCM |
|---|---|
| `iss` (absolute supersaturation) and `isLES` (relative supersaturation) | Derived atmospheric diagnostics, not process accumulators. They may be calculated from read-only host fields if requested, but are not required as condensation feedback. |
| `imv` (water-vapor mixing ratio) | Host atmospheric input. Do not modify it or create a host `qv` tendency in the current one-way implementation. |
| `tabs_LCM`, static energy, and `ilastT`/`ilastp` forcing history | Host input or SAM history state, not process diagnostics. Any LCM-owned history needed by the chosen supersaturation formulation must be designed separately. |
| `ipH`, `idelR`, `iSbyO3t`, `iSbyO3m`, `iSbyH2O2t`, `iSbyH2O2m` | Optional `SAMCHEM` chemistry outputs; exclude from the first standalone diagnostic contract unless chemistry is explicitly retained. |
| Latent-heat and vapor-depletion feedback increments | Not part of the current one-way diagnostic contract. Defer host tendency/feedback outputs and conservation accounting to a later sprint. |

### Ownership recommendation

- Keep condensation diagnostics in an LCM-owned diagnostic state, indexed by local cell or returned as per-cell process increments; do not add them to `lcm_particle_t`.
- Keep diagnostic state separate from attached host fields in `lcm_field_state_t`.
- Reset per-step process increments at the start of `apply_lcm_microphysics`; do not erase longer-window output accumulators there unless the output design explicitly requires it.
- Preserve process identity and sign information (condensation versus evaporation, deposition versus sublimation, and nucleation mechanism) so later output can report useful rates without losing the raw mass/count changes.
- Two-way host feedback remains deferred and must not be inferred from these diagnostics.

No standalone diagnostics module exists in the inspected source today. This classification informs a later `lcm_condensation`/diagnostics interface design; it does not implement one.

## Proposed one-way `lcm_condensation` interface

This is the proposed contract for the first standalone implementation. It uses a new module (`src/lcm_condensation.f90`) and treats `src/micro_cond.f90` only as the SAM-LCM scientific reference. It does not require edits to the legacy routine.

### Public entry point and ownership

```fortran
module lcm_condensation
   implicit none
   private
   public :: lcm_apply_condensation
contains
   subroutine lcm_apply_condensation(time_step)
      ! Apply local condensation and phase changes to LCM-owned particles.
   end subroutine lcm_apply_condensation
end module lcm_condensation
```

The public routine should take `time_step` explicitly. It should iterate local cells and dispatch to a private per-cell routine. To match the current repository pattern, it can access initialized `cell_particles`, `grid_state`, `field_state`, and `config_state` through their owning modules; it should not expose raw SAM arrays or `n_start`/`n_end` through the interface. The private per-cell routine should receive cell indices and the timestep explicitly, and mutate only its cell's particle state and LCM diagnostic state. `apply_lcm_microphysics(time_step)` in `lcm_advance_mod` remains the orchestration point and calls `lcm_apply_condensation(time_step)` once after particle movement, boundary handling, migration/exchange, and source handling.

### Inputs, state changes, and outputs

| Contract element | One-way behavior |
|---|---|
| Timestep | Read `time_step`; never read `dtn_LCM`. |
| Particle storage | Iterate `cell_particles(i,j,k)%particles(1:count)`; read/write named fields on each `lcm_particle_t`. |
| Grid/cell state | Read cell geometry and `dV` from `grid_state`; do not store cell geometry in particles. |
| Atmospheric inputs | Read required attached host fields through `field_state`: `qv`, pressure, density, and temperature/static energy according to the selected SAM-equivalent thermodynamic calculation. The attached host fields are read-only. |
| Configuration/constants | Read process switches and mode parameters from `config_state`; use shared physical constants from `lcm_parameters` or focused helper modules. |
| Particle outputs | Update particle water mass, phase, and relevant particle-carried LEM/ice state using named fields. Motion/location values, including `z` and `w_sedi`, are read-only inputs to condensation. |
| Diagnostics | Accumulate per-step process increments/counts in LCM-owned diagnostic state, separate from particles and host fields. |
| Host outputs | None in this sprint. Do not change host `qv`, temperature, static energy, pressure, or density; do not prepare host tendencies. |

### Per-call sequence

1. Validate the required runtime state and the atmospheric inputs needed by the selected configuration.
2. Reset only the per-step condensation diagnostic increments.
3. Loop over local `(i,j,k)` cells and skip empty particle cells.
4. Read the cell's atmosphere and geometry once for that cell; keep these values read-only.
5. Apply the chosen condensational growth and phase-change processes to particles in the cell, preserving SAM process ordering where applicable.
6. Record per-cell process increments and counts in LCM-owned diagnostics; do not reduce or apply them to host fields here.

### Decisions required before implementing the physics body

- Define a distinct liquid/ice phase value that does not overload the existing `particle_type` initialization category. Preserve the SAM distinction for IN membership separately.
- Choose the supported first thermodynamic path. `switch_lem_micro` is enabled by default and the particle type already carries `temperature`, `qv`, and `eta`; the non-LEM path instead uses cell-mean state and SAM's previous-step forcing history. The implementation must state which path(s) are in the initial port and how the other path is handled.
- Define how `field_state%temperature` and `%static_energy` are used for saturation and freezing calculations. Both are attached and required today, but they are not interchangeable quantities.
- Define the cell volume allocation for particle-local LEM state. SAM divides `dV` by the number of particles in the legacy cell range; do not carry this expression over until its scientific meaning in cell-owned storage is confirmed.
- Confirm diagnostic names, signs, and whether the first implementation stores raw per-step increments only or also computes normalized rates. Never normalize twice.
- Optional `SAMCHEM` hooks are excluded from the first one-way interface unless separately selected for support.

This interface keeps the first port local to LCM: host atmosphere in, particle state and LCM diagnostics out. Two-way moisture/energy feedback and its conservation contract remain explicitly deferred to a later sprint.
