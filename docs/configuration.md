# LCM Configuration Interface

LCM configuration is passed through the host interface as an `lcm_config_t`
object. The host model owns its own runtime input format, reads the options it
wants to expose, fills the config object, and passes it to LCM with
`lcm_set_config`.

This follows the same separation used by COMIN and ICON: the host model decides
how configuration is read, while the coupled component receives a typed runtime
configuration state. LCM should not require DALES, SAM, or another host to
create a separate LCM namelist file.

Typical host-side use is:

```fortran
use lcm_host_interface, only: lcm_config_t, lcm_init_config, lcm_set_config

type(lcm_config_t) :: lcm_config

call lcm_init_config(lcm_config)

! Host namelist values can override selected defaults here.
! For example:
! lcm_config%np_xy = host_lcm_np_xy
! lcm_config%switch_collection = host_lcm_collection

call lcm_set_config(lcm_config)
```

For DALES, this means the DALES namelist should eventually expose only the LCM
settings that are meaningful for a DALES run. The adapter then maps those values
to `lcm_config_t`. SAM can make a different choice without changing the LCM
library interface.
