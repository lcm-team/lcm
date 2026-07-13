!> Public advance facade for host-driven LCM advancement.
module lcm_advance_mod
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
      call run_lcm_advance_sequence(time_step)
   end subroutine lcm_advance

   subroutine run_lcm_advance_sequence(time_step)
      real(real64), intent(in) :: time_step

      call initialize_lcm_state()
      call update_lcm_atmospheric_state()
      call advance_lcm_particles(time_step)
      call apply_lcm_microphysics(time_step)
      call exchange_lcm_particles()
      call prepare_lcm_feedback()
   end subroutine run_lcm_advance_sequence

   subroutine validate_step_state(time_step)
      real(real64), intent(in) :: time_step

      if (time_step <= 0.0_real64) then
         error stop 'LCM time_step must be positive'
      endif

      if (.not. config_is_set) then
         error stop 'LCM config must be set before lcm_advance'
      endif

      if (.not. grid_state%is_set) then
         error stop 'LCM grid must be set before lcm_advance'
      endif
   end subroutine validate_step_state

   subroutine initialize_lcm_state()
      ! Future internal initialization point for particle state, restart state,
      ! and other LCM-owned state that depends on config/grid setup.
   end subroutine initialize_lcm_state

   subroutine update_lcm_atmospheric_state()
      ! Future bridge from attached host fields into LCM-owned atmospheric state.
   end subroutine update_lcm_atmospheric_state

   subroutine advance_lcm_particles(time_step)
      real(real64), intent(in) :: time_step

      ! Future particle advection step.
   end subroutine advance_lcm_particles

   subroutine apply_lcm_microphysics(time_step)
      real(real64), intent(in) :: time_step

      ! Future particle microphysics step.
   end subroutine apply_lcm_microphysics

   subroutine exchange_lcm_particles()
      ! Future particle exchange step for particles crossing rank boundaries.
   end subroutine exchange_lcm_particles

   subroutine prepare_lcm_feedback()
      ! Future feedback/tendency preparation step for the host model.
   end subroutine prepare_lcm_feedback

end module lcm_advance_mod
