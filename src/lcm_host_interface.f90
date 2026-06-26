!> Public interface for host models coupled to LCM.
module lcm_host_interface
   use lcm_types, only: lcm_grid_t, lcm_init_grid

   implicit none

   private

   public :: lcm_grid_t
   public :: lcm_init_grid

end module lcm_host_interface
