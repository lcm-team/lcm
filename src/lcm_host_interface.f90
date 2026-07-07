!> Public interface for host models coupled to LCM.
module lcm_host_interface
   use lcm_types, only: lcm_grid_t, lcm_config_t, lcm_init_grid, lcm_init_config
   use lcm_config_state, only: lcm_set_config
   use lcm_grid_state, only: lcm_set_grid

   implicit none

   private

   public :: lcm_grid_t
   public :: lcm_config_t
   public :: lcm_init_grid
   public :: lcm_init_config
   public :: lcm_set_config
   public :: lcm_set_grid

end module lcm_host_interface
