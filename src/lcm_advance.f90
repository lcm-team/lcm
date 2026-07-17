!> Public advance facade for host-driven LCM advancement.
module lcm_advance_mod
   use iso_fortran_env, only: real64
   use lcm_config_state, only: config_is_set
   use lcm_grid_state, only: grid_state
   use lcm_particle_state, only: cell_particles, lcm_allocate_particle_cells

   implicit none

   private

   public :: lcm_advance

contains

   subroutine lcm_advance(time_step)
      real(real64), intent(in) :: time_step

      call validate_step_state(time_step)
      call initialize_lcm_particles()

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

   subroutine initialize_lcm_particles()
      if (allocated(cell_particles)) return

      ! This is where LCM-owned particle initialization starts. The current
      ! facade allocates the per-cell containers; future work will populate
      ! them with computational particles.
      call lcm_allocate_particle_cells()
   end subroutine initialize_lcm_particles

end module lcm_advance_mod
