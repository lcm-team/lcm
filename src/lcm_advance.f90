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

      call prepare_lcm_atmosphere()
      call advance_lcm_particles(time_step)
      call apply_lcm_particle_boundaries()
      call exchange_lcm_particles()
      call update_lcm_particle_cell_ownership()
      call apply_lcm_microphysics(time_step)
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

   subroutine prepare_lcm_atmosphere()
      ! Future preparation/validation step for attached host atmospheric fields.
   end subroutine prepare_lcm_atmosphere

   subroutine advance_lcm_particles(time_step)
      real(real64), intent(in) :: time_step

      ! Future particle advection step.
   end subroutine advance_lcm_particles

   subroutine apply_lcm_particle_boundaries()
      ! Future boundary handling for particles outside the local physical domain.
   end subroutine apply_lcm_particle_boundaries

   subroutine apply_lcm_microphysics(time_step)
      real(real64), intent(in) :: time_step

      ! Future particle microphysics step.
   end subroutine apply_lcm_microphysics

   subroutine exchange_lcm_particles()
      ! Future particle exchange step for particles crossing rank boundaries.
   end subroutine exchange_lcm_particles

   subroutine update_lcm_particle_cell_ownership()
      ! Future local cell ownership update after particle motion and exchange.
   end subroutine update_lcm_particle_cell_ownership

   subroutine prepare_lcm_feedback()
      ! Future feedback/tendency preparation step for the host model.
   end subroutine prepare_lcm_feedback

end module lcm_advance_mod
