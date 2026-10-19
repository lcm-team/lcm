# Issue 37 progress

Last updated: 2026-10-19.

This is the short working log for Issue #37, **“refactor: port
micro_condensation to standalone LCM.”** It records decisions, completed work,
current implementation status, and the next work item. The detailed SAM source
inventory and line-by-line parity checklist remain in
[`issue-37-micro-vars-dependency-inventory.md`](issue-37-micro-vars-dependency-inventory.md)
and [`issue-37-sam-condensation-port-map.md`](issue-37-sam-condensation-port-map.md).

## Agreed direction

- Implement condensation in standalone `src/lcm_condensation.f90`; leave the
  SAM legacy modules unchanged.
- Preserve SAM-LCM equations, process order, conditions, and state transitions.
  Do not add physics or helper behavior without a SAM counterpart.
- Use one-way coupling this sprint: attached DALES fields are read-only;
  condensation updates LCM particles and LCM-owned process diagnostics. Design
  host moisture and energy feedback in a later sprint.
- Keep both SAM thermodynamic paths: LEM particle-local supersaturation and
  non-LEM cell-mean supersaturation/history.
- Preserve SAM's `dV / number_of_particles_in_cell` allocation where its LEM
  calculation uses it.
- Preserve SAM's signed particle type convention: sign represents liquid/ice
  phase and magnitude represents IN membership.
- Store raw per-step diagnostic increments/counts; apply any reporting
  normalization once in the later diagnostic/output path.
- Defer `SAMCHEM` for the first one-way port.

## Work completed

- Inventoried the direct `micro_vars` dependencies and the
  `microphysics`/`micro_field` dependencies used by SAM `micro_condensation`.
- Created a SAM process-order and helper parity map so each branch has an
  explicit implementation target or an agreed scope exclusion.
- Added standalone particle phase constants and aligned phase/IN interpretation
  with SAM-LCM (the work is on this feature branch; PR #76 was opened for it).
- Added the initial `lcm_condensation` implementation for liquid growth and
  evaporation. It includes the LEM and non-LEM supersaturation paths, analytical
  and Euler radius solvers, ventilation, Köhler/FFC handling, transport
  coefficient options, radiation input, particle-local LEM updates, and raw
  liquid diagnostics.
- Added the optional attached `qrad` field required by SAM's radiation branch.
- Wired condensation into the LCM advance sequence. DALES calls `lcm_advance`
  through its post-RK `lcm_microphysics` hook, after final atmospheric field
  updates, and passes the full `rdt` timestep.
- Moved fixed SAM physical values (including `r_act` and FFC constants) to
  `lcm_parameters`; kept runtime switches in the configuration type. Exposed
  diagnostic increments as individual allocatable arrays rather than an extra
  diagnostics wrapper type.

## Current state

- Branch: `refactor/issue-37-micro-condensation`.
- The full CMake build succeeds. The local CTest suite reports 47 tests: 37
  passed and 10 failed. Failures are in ice-enabled advance/source cases that
  hit the deliberate unported-ice stop, transport cases with default
  zero-multiplicity particles, and two particle-source temperature assertions
  affected by the newly connected liquid-growth call. These regressions need
  resolution as the SAM branches are completed and the integration is refined.
- The current physics implementation is a partial liquid-only slice. It stops
  explicitly when ice is enabled or an ice particle is encountered, so ice
  physics is not silently omitted.
- Issue #37 is not complete. The exact SAM ice deposition/growth paths for both
  Chen and DISKICE habits, freezing and shattering, melting, complete
  sublimation transitions, their helper procedures and diagnostics, and a
  full source-parity audit remain to be implemented.
- The current one-way interface intentionally does not update DALES `qv`,
  temperature, or static energy. `SAMCHEM` remains deferred.

## Next work

Continue from the parity map: trace the SAM ice-growth dependencies and port the
ice deposition/growth paths and their helpers without replacing their
algorithms. Then port freezing/shattering, melting, and complete sublimation
transitions. Keep each completed source range and any implementation blockers
updated in the parity map; add a brief dated entry here when a meaningful
milestone is reached.
