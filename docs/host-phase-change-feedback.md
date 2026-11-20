# Host-side phase-change feedback

`lcm_get_liquid_water_mass(mass)` copies represented liquid water in kg per
local grid cell. Multiplicity is already included in the stored water mass;
do not multiply by weight again.

`lcm_get_condensation_mass_change(dm)` copies the net liquid-water increment
from the last condensation call in kg/cell. Condensation is positive,
evaporation negative, and particle transport is excluded. Values are zero
following initialization and are reset on each condensation call.

A host whose humidity normalization depends on transported particle loading
can call `lcm_advance(dt, refresh_atmosphere=refresh)` with an optional
no-argument procedure. It runs after particle boundaries, migration, MPI
exchange and sources, immediately before microphysics. It should refresh
already-attached host atmospheric fields in place. Hosts omitting the
callback retain the previous call sequence.

The host owns its atmospheric prognostic fields and applies the returned
phase change once to its vapor/energy budget. The getters do not modify
host fields or perform feedback. The DALES piggyback adapter uses them on
its private thermodynamic state, leaving the bulk driver unchanged.
