# Issue 37: SAM-LCM condensation port map

This is the source-to-source checklist for the new standalone `lcm_condensation`
module. `src/micro_cond.f90` is the scientific source of truth. The port must
retain the SAM equations, conditions, and process order for every supported
branch. This map is a guard against silently dropping a branch or adding a
replacement process that SAM does not contain.

## Scope decisions already made

- Implement in a new standalone module; leave the SAM-derived `micro_cond.f90`
  unchanged.
- Keep one-way coupling in this sprint: host fields are read-only, particle
  state and LCM-owned diagnostics are updated, and no host vapor or thermal
  feedback is applied. Define two-way feedback in a later sprint.
- Preserve both `switch_LEM_micro` modes. Keep their supersaturation histories
  and calculations distinct as in SAM.
- Preserve the SAM equal allocation `dV / number_of_particles_in_cell` wherever
  the SAM LEM path uses it. Here the number is the cell's active particle count,
  including liquid and ice particles.
- Use the dry static energy provided through `field_state` to derive the SAM
  temperature used by condensation, following the existing LCM field-state
  convention. Do not substitute dry static energy directly into saturation
  equations.
- Preserve SAM's signed `particle_type` encoding: sign is liquid/ice phase and
  magnitude is IN membership. Use the named standalone constants.
- Preserve SAM's process diagnostics as raw per-step increments/counts first;
  do not normalize them in the process routine. SAM's output normalization is
  performed later by its diagnostic path.
- Do not port `SAMCHEM` in the first one-way interface. Keep its omission
  explicit in the code and design notes so it can be added as a separately
  scoped interface later.

## `micro_condensation` process order and parity checklist

Source ranges below refer to `src/micro_cond.f90` in this repository. A checked
item means the standalone code must have a corresponding implementation or a
documented scope exclusion. It does not mean that a simpler substitute is
acceptable.

| Order | SAM source | SAM behavior that must be retained | Standalone mapping / boundary |
|---:|---|---|---|
| 1 | 118–133 | Count liquid and ice particles from signed `iTYPE`, build their index lists, start the routine timer. | Iterate typed particles in the cell; classify by the signed SAM-compatible type. Keep timer behavior if the standalone timer API is used. |
| 2 | 135–179 | Initialize net liquid/ice/fusion changes; optional DISKICE critical supersaturation; choose variable or constant thermal conductivity and diffusivity; compute water and ice saturation/pressure and growth coefficients; set Köhler curvature term or disable both Köhler terms. | Preserve each configured coefficient choice and formula. Use standalone field/parameter inputs and exact SAM saturation functions. |
| 3 | 181–200 | For every particle, compute aerosol radius and LEM temperature; for liquid compute radius, ventilation, and cell liquid mixing ratio; for ice compute shape capacitance and ice ventilation. | Use named particle fields and cell geometry. Keep motion values such as `w_sedi` read-only. |
| 4 | 202–230 | If radiation is enabled, compute droplet absorption `Q_a`, area-weighted radiation allocation, liquid share of radiative heating, and zero it below SAM's liquid-water threshold; otherwise set radiation and `Q_a` to zero. | Preserve the entire branch and its configuration restrictions. Add an explicit host radiation input if the current field interface lacks `qrad`; do not replace it with zero or an invented estimate. |
| 5 | 232–252 | If FFC is enabled, compute droplet area, FFC molecules from aerosol mass and multiplicity, choose the accommodation coefficient using SAM's molecular surface threshold, and derive `r0`; otherwise use the fixed accommodation coefficient and the same `r0` formula. | Preserve constants, formula, threshold, and option. Requires an aerosol-radius/composition calculation equivalent to SAM `get_r_aero`/`get_bfactor`. |
| 6 | 254–338 | Compute average supersaturation. Non-LEM: form liquid/ice depletion sums and timescale, use current and previous cell vapor/temperature/pressure history, include the ice/water saturation difference in forcing, integrate with SAM `eta_new`/`eta_int`, bound the integral, and apply the initial-step supersaturation cap. LEM: repeat the calculation per particle using SAM's equal cell-volume allocation, particle vapor, temperature, and stored `eta`; include the ice forcing when ice is enabled. Derive liquid supersaturation and, when enabled, ice mixing-ratio and vapor-pressure supersaturation. | Preserve both modes and their exact histories, bounds, and formulas. Cell-mean history must be LCM-owned. Host history/fields remain read-only. Do not use an ice-free shortcut when `switch_ice` is true. |
| 7 | 339–343 | For non-LEM, calculate DISKICE diffusion parameters; calculate the fitting coefficients used by the DISKICE deposition formulation. | Keep the exact calls when the corresponding ice/habit branch is enabled; avoid computing unused ice-only values in a liquid-only configuration. |
| 8 | 345–369 (`SAMCHEM`) | Initialize chemistry diagnostics, conditionally call `chem_proc`, and allow its radii to replace the ordinary liquid-growth solution. | Explicit first-port exclusion by prior decision. No chemistry placeholder may alter radius or state. Retain a clear TODO/interface note only. |
| 9 | 371–430 | Integrate liquid radius for every liquid particle using either SAM's Euler/Newton solve or analytic solution; in Euler mode compute particle Köhler solute term when enabled. Update particle water mass; accumulate signed weighted mass change; in LEM mode update particle temperature and vapor; separately accumulate condensation or evaporation based on radius direction; count activation/deactivation threshold crossings. | Preserve both integrators, the full Euler helper, all options, mass arithmetic, latent/vapor changes to particle-local LEM state, and raw signed diagnostics. Host fields are not changed. |
| 10 | 432–487 | If habit is `Chen`, grow/sublimate ice with `mass_ice_euler`; update rimed mass consistently during sublimation and cap it at total ice mass; update habit dimensions/density through `acd_update` or reset to spherical ice when habit evolution is off; accumulate deposition/sublimation; update LEM temperature/vapor. | Required when ice is enabled. Preserve helper equations and guards, typed axes/density/rime fields, and diagnostic semantics. |
| 11 | 488–542 | If habit is `DISKICE`, use `diffparams` in LEM mode, call `ice_growth_diskice`, update mass/rime/habit state, accumulate deposition/sublimation, and update LEM temperature/vapor. | Required for this configured SAM branch. Preserve DISKICE size threshold, fitting, deposition coefficient, and shape update helpers exactly. |
| 12 | 544–624 | Once ice start time is reached, perform heterogeneous immersion freezing only for liquid particles with IN, positive water mass, temperature at/below their freezing threshold, and positive liquid supersaturation. Preserve each enabled droplet-shattering parameterization (`Phillips`, `Lachapelle`, `Lauber`), fragment matching/reassignment and its multiplicity changes, then apply the SAM immersion-freezing transition/state initialization and diagnostics. | Preserve the predicates, each selectable shattering branch, transition updates, and counts/mass. Any random draw must use a reproducible standalone random source with the same sampling semantics. |
| 13 | 627–660 | Perform homogeneous freezing below SAM's temperature threshold and supersaturation condition, using the Kuhn rate and timestep probability; use SAM's random comparison; initialize ice axes/density, preserve IN membership while changing phase sign, reset monomer/rime state, and accumulate fusion and nucleation diagnostics. | Preserve equations, thresholds, random decision, state transition, and raw diagnostics. |
| 14 | 663–680 | Instantly melt ice above 273.15 K; change only the phase sign while preserving IN magnitude, reset monomer/rime state, account for fusion, and accumulate melt diagnostics. | Preserve condition, state, and diagnostic semantics. |
| 15 | 682–707 | On complete sublimation below SAM's per-particle ice mass threshold, return ice to liquid while preserving IN membership, reset monomer/rime state, restore multiplicity from the dominant IN or aerosol mode, apply the `1/kg` density conversion, and account for fusion and fully-sublimated count. | Preserve dominant-mode selection, multiplicity arithmetic, unit branch, and transition. Do not silently clamp/alter the restored count. |
| 16 | 709–725 | SAM applies net vapor depletion and latent heating to host fields, refreshes supersaturation and previous cell thermodynamic history, then updates each LEM particle's absolute supersaturation from its post-process vapor and temperature. | Do not mutate host vapor, temperature, static energy, or host history in one-way coupling. Preserve the LEM particle `eta` update. Store any non-LEM cell history in LCM-owned state for the next call. |

## Helper procedures that belong to the port

The following helpers are called by, or define equations used by, the routine.
They are part of source parity and must not be omitted when their branch is
supported:

| Helper | Source range | Use |
|---|---:|---|
| `venteef` | 731–758 | Liquid ventilation for droplets larger than 40 micrometers. |
| `eta_int` | 760–769 | Integrated absolute supersaturation for the non-LEM and LEM TFL calculation. |
| `eta_new` | 771–780 | End-of-step absolute supersaturation for the TFL calculation. |
| `radius_liquid_euler` | 782–842 | SAM Newton-Raphson Euler growth with ventilation, Köhler terms, radiation, molecular length scale, and aerosol-radius lower bound. |
| `mass_ice_euler` | 844–866 | Ice mass growth/deposition/sublimation for the Chen habit path. |
| `ice_growth_diskice` | 868–952 | DISKICE deposition/sublimation and axis/density evolution. |
| `capacitance` | 953–982 | Spheroid capacitance for ice. |
| `gamma_ice` | 984–1024 | Temperature-indexed Chen-Lamb inherent growth ratio. |
| `getscrit` | 1026–1086 | DISKICE critical supersaturation fit. |
| `diffparams` | 1088–1121 | DISKICE diffusion/kinetic parameters. |
| `fittingparams` | 1137–1171 | DISKICE surface supersaturation fit coefficients. |
| `depcoeffparam` and its continuation | 1174 onward | DISKICE deposition coefficients and local surface supersaturation solution. |
| `acd_update` | 1313 onward | Chen habit dimensions and apparent density update. |
| `ice_ventilation` | 1377 onward | Ice ventilation factors used by both ice habits. |
| `micro_ds` shattering routines | imported at line 48 | Phillips, Lachapelle, and Lauber droplet-shattering implementations selected by the freezing branch. |
| `get_r_aero`, `get_bfactor` | imported from `micro_vars` | Aerosol radius and composition-dependent Köhler solute term. Port their actual algorithms and composition choices, not approximations. |
| `get_V_aero` | `micro_vars.f90:569–581` | Sums aerosol volume by mode using aerosol portions, aerosol mass per multiplicity, and each mode's aerosol density. |
| `get_bfactor` | `micro_vars.f90:513–532` | Uses the one-mode cached coefficient or mode-weighted coefficients, then scales by aerosol mass per multiplicity. |
| `esatw`, `qsatw`, `dtqsatw` | SAM's `vars` module; compatibility counterpart in `stubs/sam_compat.f90:135 onward` | The compatibility implementation currently matches the standalone `qsatw_liquid` expression, but the new module should call a shared exact implementation rather than maintain a divergent duplicate. Confirm precision/host build parity before moving it. |
| `esati`, `qsati`, `dtqsati` | SAM's `vars` module; compatibility counterpart in `stubs/sam_compat.f90` | Required by enabled ice supersaturation and deposition branches; retain the same functions and finite-difference derivative. |
| `sigma_air_liq`, `sigma_ice_liq` | `micro_cond.f90:1406–1444` | Legacy definitions are local functions and are not called by `micro_condensation` in this source version. Do not add them to the port unless another retained counterpart call is found. |

Before coding a helper, trace its complete definition and downstream helpers in
the source. This table is not permission to replace an unavailable dependency
with a new approximation.

## Deliberate one-way exclusions

These are the only known differences from the SAM routine at this stage:

1. `SAMCHEM`'s `chem_proc` branch is deferred by the agreed scope.
2. The final SAM writes to host vapor, temperature, supersaturation, and
   previous-step history are not performed. LCM retains particle-local LEM
   updates and stores any required non-LEM history internally.

Every other physical branch and helper in the checklist remains in the Issue
#37 port scope. If a required standalone input or output is missing, add the
matching interface/state or stop and report the blocker; do not silently omit
the SAM branch or synthesize replacement physics.

## Current implementation status

- Source map: all process branches and directly referenced helpers are listed.
  The full transitive dependency and interface trace remains to be completed
  before claiming implementation parity.
- `src/lcm_condensation.f90`: first code slice implements the SAM liquid
  condensational growth/evaporation path, including both LEM and non-LEM
  supersaturation calculations, analytical and Euler radius updates, liquid
  ventilation, Köhler terms, variable transport coefficients, radiation input,
  particle-local LEM temperature/vapor/eta updates, and raw liquid diagnostics.
  Host fields remain read-only. Ice-enabled configurations and ice particles
  stop with an explicit error until the corresponding SAM branches are ported.
- Still to implement from the map: ice supersaturation and deposition for both
  habits, freezing and shattering, melting, complete-sublimation transitions,
  their helper routines and diagnostics, and corresponding branch-by-branch
  source verification. SAMCHEM and host feedback remain the agreed exclusions.
- Tests: not added or run.
