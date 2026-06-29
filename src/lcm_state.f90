!> Internal state owned by LCM.
module lcm_state
   use lcm_types, only: lcm_grid_t

   implicit none

   private

   public :: lcm_set_grid

   type :: t_lcm_state
      type(lcm_grid_t) :: grid
      logical :: grid_is_set = .false.
   end type t_lcm_state

   type(t_lcm_state), save :: state

contains

   subroutine lcm_set_grid(grid)
      type(lcm_grid_t), intent(in) :: grid

      state%grid = grid
      state%grid_is_set = .true.
   end subroutine lcm_set_grid

end module lcm_state
