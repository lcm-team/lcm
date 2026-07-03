# LCM Timestep Interface

Host models advance LCM through the public `lcm_advance` routine:

```fortran
use iso_fortran_env, only: real64
use lcm_host_interface, only: lcm_advance

call lcm_advance(real(host_time_step, kind=real64))
```

The host remains responsible for deciding where this call belongs in its own
time loop. For DALES, the natural first integration point is the selected
microphysics branch, where the current model timestep is available as `rdt`.

The initial facade validates that the timestep is positive and that LCM has
received its configuration and grid state. The actual particle advancement will
be wired behind this routine after the internal `micro_vars` and coupling
routines are refactored for the standalone host interface.
