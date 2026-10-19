!> Public advance facade for host-driven LCM advancement.
module lcm_advance_mod
   use iso_fortran_env, only: real64
   use lcm_config_state, only: config_is_set, config_state
   use lcm_condensation, only: lcm_apply_condensation
   use lcm_grid_state, only: grid_state
   use lcm_exchange_directions, only: lcm_exchange_x_phase, &
                                      lcm_exchange_y_phase
   use lcm_initialization_state, only: lcm_is_initialized
   use lcm_parallel_state, only: parallel_state
   use lcm_particle_boundaries, only: lcm_apply_particle_boundaries
   use lcm_particle_exchange, only: lcm_exchange_particle_phase
   use lcm_particle_state, only: lcm_prepare_particle_migration, &
                                 lcm_prepare_particle_velocities, &
                                 lcm_refine_particle_velocities, &
                                 lcm_update_particle_positions
   use lcm_particle_source, only: lcm_apply_particle_sources

   implicit none

   private

   public :: lcm_advance

   integer, save :: current_lcm_step = 0
   logical, save :: current_lcm_step_is_initialized = .false.

contains

   subroutine lcm_advance(time_step)
      real(real64), intent(in) :: time_step

      call validate_step_state(time_step)
      call initialize_lcm_step_counter()
      call run_lcm_advance_sequence(time_step)
      current_lcm_step = current_lcm_step + 1
   end subroutine lcm_advance

   subroutine run_lcm_advance_sequence(time_step)
      real(real64), intent(in) :: time_step

      call prepare_lcm_atmosphere()
      call advance_lcm_particles(time_step)
      call apply_lcm_particle_boundaries()
      call prepare_lcm_particle_migration()
      call exchange_lcm_particles()
      call lcm_apply_particle_sources(time_step, current_lcm_step)
      call apply_lcm_microphysics(time_step)
      call prepare_lcm_feedback()
   end subroutine run_lcm_advance_sequence

   subroutine initialize_lcm_step_counter()
      if (.not. current_lcm_step_is_initialized) then
         current_lcm_step = config_state%nstep_lcm_init
         current_lcm_step_is_initialized = .true.
      endif
   end subroutine initialize_lcm_step_counter

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

      if (.not. lcm_is_initialized) then
         error stop 'LCM must be initialized before lcm_advance'
      endif
   end subroutine validate_step_state

   subroutine prepare_lcm_atmosphere()
      ! Future preparation/validation step for attached host atmospheric fields.
   end subroutine prepare_lcm_atmosphere

   subroutine advance_lcm_particles(time_step)
      real(real64), intent(in) :: time_step

      call lcm_prepare_particle_velocities()
      call lcm_refine_particle_velocities(time_step)
      call lcm_update_particle_positions(time_step, current_lcm_step)
   end subroutine advance_lcm_particles

   subroutine apply_lcm_particle_boundaries()
      call lcm_apply_particle_boundaries()
   end subroutine apply_lcm_particle_boundaries

   subroutine prepare_lcm_particle_migration()
      call lcm_prepare_particle_migration()
   end subroutine prepare_lcm_particle_migration

   subroutine apply_lcm_microphysics(time_step)
      real(real64), intent(in) :: time_step

      call lcm_apply_condensation(time_step, current_lcm_step)
   end subroutine apply_lcm_microphysics

   subroutine exchange_lcm_particles()
      ! Migration preparation has already completed local moves and populated
      ! directional buffers. These phases now perform transport and insertion.
      if (.not. parallel_state%is_set) return

      call lcm_exchange_particle_phase(lcm_exchange_x_phase)
      call lcm_exchange_particle_phase(lcm_exchange_y_phase)
   end subroutine exchange_lcm_particles

   subroutine prepare_lcm_feedback()
      ! Future feedback/tendency preparation step for the host model.
   end subroutine prepare_lcm_feedback

end module lcm_advance_mod
