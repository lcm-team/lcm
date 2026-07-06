!> Public timestep facade for host-driven LCM advancement.
module lcm_timestep_state
   use iso_fortran_env, only: real64
   use lcm_config_state, only: config_is_set
   use lcm_grid_state, only: grid_state

   implicit none

   private

   public :: lcm_advance

contains

   subroutine lcm_advance(time_step)
      real(real64), intent(in) :: time_step

      call validate_step_state(time_step)

      ! TODO: call the internal LCM advance routine here once micro_vars and
      ! micro_coupling have been refactored behind the host interface.
   end subroutine lcm_advance

   subroutine validate_step_state(time_step)
      real(real64), intent(in) :: time_step

      if (.not. config_is_set) then
         error stop 'LCM config must be set before lcm_advance'
      endif

      if (.not. grid_state%is_set) then
         error stop 'LCM grid must be set before lcm_advance'
      endif
   end subroutine validate_step_state

end module lcm_timestep_state
