!> Internal state owned by LCM.
module lcm_state
   use lcm_types, only: lcm_grid_t

   implicit none

   private

   public :: lcm_state_t
   public :: lcm_state_data
   public :: lcm_set_grid

   type :: lcm_state_t
      type(lcm_grid_t) :: grid
      logical :: grid_is_set = .false.
   end type lcm_state_t

   type(lcm_state_t), save :: lcm_state_data

contains

   subroutine lcm_set_grid(grid)
      type(lcm_grid_t), intent(in) :: grid

      lcm_state_data%grid = grid
      lcm_state_data%grid_is_set = .true.
   end subroutine lcm_set_grid

end module lcm_state
