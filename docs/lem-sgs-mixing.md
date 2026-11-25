# SAM-derived LEM SGS mixing

The full standalone advance sequence with `switch_lem_micro=true` requires an attached LEM environment. It transports/exchanges particles, applies sources, mixes particle-associated air, then computes condensation. Direct calls to `lcm_apply_condensation` remain process-kernel calls and do not perform SGS mixing.

## Host interface

Call `lcm_attach_lem_environment` after grid setup and keep all six TARGET arrays alive, shape `(nx_local,ny_local,nz)`. Refresh their values before each full advance:

| Argument | Meaning and units |
|---|---|
| `eddy_diffusivity` | Scalar SGS diffusivity, m²/s |
| `eddy_length` | SGS mixing length, m; largest LEM eddy is limited to the local cell height |
| `total_water_dry` | Resolved total water, kg water/kg dry air |
| `liquid_static_temperature` | Resolved `T+gz/cp-(L/cp)*r_l`, K |
| `static_temperature_tendency` | External dry-static-temperature forcing including radiative/surface heating, K/s |
| `vapor_mixing_ratio_tendency` | External vapor mass source per current dry-air mass, kg/kg/s |

The supplied mean targets must use the same thermodynamic conventions as LCM. They are independent of the instantaneous LCM particle liquid. In SAM the targets are reconstructed from its coupled mean vapor, temperature and LCM liquid. The host supplies consistent resolved targets in the same convention. This mixing module reads the attached host arrays; atmospheric feedback remains the host coupling’s responsibility.

## Algorithm and scope

This ports `micro_sgs_mixing` and its triplet rearrangement from SAM-LCM. With more than ten computational particles in a cell: apply external forcing, order parcels by their previous cell and LEM position, apply periodic explicit diffusion and stochastic triplet rearrangements, then relax total water, liquid-static energy and supersaturation history to the resolved environment. The SAM relaxation combines 900 s with its sedimentation timescale. Ten or fewer particles use the SAM-style immediate environmental reconstruction, including the vapor correction for each particle's liquid.

Liquid mass and multiplicity are not changed by mixing. Parcel metadata `cell_particle_index,last_i,last_j,last_k` are refreshed and already belong to the existing MPI payload. Particle motion must not skip more than one grid cell per full step. Sparse fallback and relaxation exchange thermodynamics with the host target; they are not a closed global particle-air transport budget.

Warm liquid is supported. Invalid/missing environment, negative vapor, nonfinite states, nonpositive absolute temperature, and an explicit relaxation timestep exceeding its timescale fail with an error. There is no silent vapor clipping. The outer timestep and original Euler growth solver are unchanged. Diffusion uses SAM-style internal stability substeps; the outer timestep remains fixed during this operator. Per-cell vertical extent is used.

## Deliberate differences from the legacy source

The first periodic temperature-diffusion point uses its own old value. The legacy line uses the last old value as the base, violating the periodic diffusion sum for nonuniform temperature. Insertion-sort bounds are checked before indexing. Zero turbulent diffusivity skips stochastic events without dividing by zero. Calculations use real64, exact rational exponents, and the standalone RNG rather than SAM module state. A fresh shuffled origin permutation is drawn per cell; exact random trajectories are not claimed identical to SAM.

The existing SAM event probability/substep rule and the 900 s relaxation constant are retained, not fitted to a DALES liquid curve. No ice processes, FFC or particle-radiative growth branch are added. External radiative heating of particle-associated air is separate from the deferred particle-radiative growth branch.

## Tests

Warm-cloud process tests verify periodic diffusion and triplet-map conservation, sparse fallback, source rates, conserved-state relaxation and unchanged liquid mass. The full advance tests attach an explicit environment and exercise the transport/mixing/condensation sequence.
