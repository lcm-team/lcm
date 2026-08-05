!> Public interface for host models coupled to LCM.
module lcm_host_interface
   use lcm_types, only: lcm_grid_t, lcm_config_t, lcm_init_grid, lcm_init_config
   use lcm_config_state, only: lcm_set_config
   use lcm_grid_state, only: lcm_set_grid
   use lcm_initialization_state, only: lcm_configure_runtime, lcm_initialize
   use lcm_advance_mod, only: lcm_advance
   use lcm_field_state, only: lcm_attach_fields

   implicit none

   private

   public :: lcm_grid_t
   public :: lcm_config_t
   public :: lcm_init_grid
   public :: lcm_init_config
   public :: lcm_set_config
   public :: lcm_set_grid
   public :: lcm_configure_runtime
   public :: lcm_initialize
   public :: lcm_advance
   public :: lcm_attach_fields

end module lcm_host_interface
