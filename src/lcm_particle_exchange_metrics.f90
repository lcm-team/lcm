!> Optional, rank-local performance metrics for particle exchange.
!!
!! Metrics are disabled by default so unmeasured runs do not pay clock or
!! counter costs unless a host explicitly enables measurement. Values
!! accumulate until reset and are kept separate for the x and y exchange
!! phases. Cross-rank reductions and output formatting belong to the host or
!! a later reporting layer. This module is isolated so the Stage 1
!! instrumentation can be removed after the scalability work is complete.
module lcm_particle_exchange_metrics
   use iso_fortran_env, only: int64, real64
   use lcm_exchange_directions, only: lcm_exchange_x_phase, &
                                      lcm_exchange_y_phase, &
                                      lcm_require_valid_exchange_phase

   implicit none

   private

   type, public :: lcm_particle_exchange_metrics_t
      real(real64) :: x_outgoing_preparation_seconds = 0.0_real64
      real(real64) :: y_outgoing_preparation_seconds = 0.0_real64
      integer(int64) :: x_particles_examined = 0_int64
      integer(int64) :: y_particles_examined = 0_int64
      integer(int64) :: x_particles_outgoing = 0_int64
      integer(int64) :: y_particles_outgoing = 0_int64
   end type lcm_particle_exchange_metrics_t

   logical, save :: metrics_enabled = .false.
   type(lcm_particle_exchange_metrics_t), save :: particle_exchange_metrics

   public :: lcm_get_particle_exchange_metrics
   public :: lcm_particle_exchange_metrics_are_enabled
   public :: lcm_record_particle_outgoing_preparation
   public :: lcm_reset_particle_exchange_metrics
   public :: lcm_set_particle_exchange_metrics_enabled

contains

   subroutine lcm_set_particle_exchange_metrics_enabled(enabled)
      logical, intent(in) :: enabled

      metrics_enabled = enabled
   end subroutine lcm_set_particle_exchange_metrics_enabled

   logical function lcm_particle_exchange_metrics_are_enabled() result(enabled)
      enabled = metrics_enabled
   end function lcm_particle_exchange_metrics_are_enabled

   subroutine lcm_reset_particle_exchange_metrics()
      particle_exchange_metrics = lcm_particle_exchange_metrics_t()
   end subroutine lcm_reset_particle_exchange_metrics

   subroutine lcm_get_particle_exchange_metrics(metrics)
      type(lcm_particle_exchange_metrics_t), intent(out) :: metrics

      metrics = particle_exchange_metrics
   end subroutine lcm_get_particle_exchange_metrics

   subroutine lcm_record_particle_outgoing_preparation(phase, elapsed_seconds, &
                                                       particles_examined, &
                                                       particles_outgoing)
      integer, intent(in) :: phase
      real(real64), intent(in) :: elapsed_seconds
      integer(int64), intent(in) :: particles_examined
      integer(int64), intent(in) :: particles_outgoing

      call lcm_require_valid_exchange_phase(phase)
      if (.not. metrics_enabled) return

      select case (phase)
      case (lcm_exchange_x_phase)
         particle_exchange_metrics%x_outgoing_preparation_seconds = &
            particle_exchange_metrics%x_outgoing_preparation_seconds + elapsed_seconds
         particle_exchange_metrics%x_particles_examined = &
            particle_exchange_metrics%x_particles_examined + particles_examined
         particle_exchange_metrics%x_particles_outgoing = &
            particle_exchange_metrics%x_particles_outgoing + particles_outgoing
      case (lcm_exchange_y_phase)
         particle_exchange_metrics%y_outgoing_preparation_seconds = &
            particle_exchange_metrics%y_outgoing_preparation_seconds + elapsed_seconds
         particle_exchange_metrics%y_particles_examined = &
            particle_exchange_metrics%y_particles_examined + particles_examined
         particle_exchange_metrics%y_particles_outgoing = &
            particle_exchange_metrics%y_particles_outgoing + particles_outgoing
      end select
   end subroutine lcm_record_particle_outgoing_preparation

end module lcm_particle_exchange_metrics
