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
      real(real64) :: x_count_exchange_seconds = 0.0_real64
      real(real64) :: y_count_exchange_seconds = 0.0_real64
      real(real64) :: x_payload_exchange_seconds = 0.0_real64
      real(real64) :: y_payload_exchange_seconds = 0.0_real64
      real(real64) :: x_received_insertion_seconds = 0.0_real64
      real(real64) :: y_received_insertion_seconds = 0.0_real64
      real(real64) :: ownership_update_seconds = 0.0_real64
      integer(int64) :: x_particles_examined = 0_int64
      integer(int64) :: y_particles_examined = 0_int64
      integer(int64) :: x_particles_outgoing = 0_int64
      integer(int64) :: y_particles_outgoing = 0_int64
      integer(int64) :: x_count_exchange_calls = 0_int64
      integer(int64) :: y_count_exchange_calls = 0_int64
      integer(int64) :: x_payload_exchange_calls = 0_int64
      integer(int64) :: y_payload_exchange_calls = 0_int64
      integer(int64) :: x_received_insertion_calls = 0_int64
      integer(int64) :: y_received_insertion_calls = 0_int64
      integer(int64) :: x_particles_inserted = 0_int64
      integer(int64) :: y_particles_inserted = 0_int64
      integer(int64) :: ownership_update_calls = 0_int64
      integer(int64) :: ownership_particle_examinations = 0_int64
      integer(int64) :: ownership_particles_moved = 0_int64
   end type lcm_particle_exchange_metrics_t

   type, public :: lcm_particle_exchange_metrics_timer_t
      private
      integer(int64) :: count_start = 0_int64
      integer(int64) :: count_rate = 0_int64
      integer(int64) :: count_max = 0_int64
   end type lcm_particle_exchange_metrics_timer_t

   logical, save :: metrics_enabled = .false.
   type(lcm_particle_exchange_metrics_t), save :: particle_exchange_metrics

   public :: lcm_get_particle_exchange_metrics
   public :: lcm_particle_exchange_metrics_are_enabled
   public :: lcm_record_particle_count_exchange
   public :: lcm_record_particle_outgoing_preparation
   public :: lcm_record_particle_ownership_update
   public :: lcm_record_particle_payload_exchange
   public :: lcm_record_received_particle_insertion
   public :: lcm_reset_particle_exchange_metrics
   public :: lcm_set_particle_exchange_metrics_enabled
   public :: lcm_start_particle_exchange_metrics_timer
   public :: lcm_stop_particle_exchange_metrics_timer

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

   subroutine lcm_start_particle_exchange_metrics_timer(timer)
      type(lcm_particle_exchange_metrics_timer_t), intent(out) :: timer

      call system_clock(timer%count_start, timer%count_rate, timer%count_max)
   end subroutine lcm_start_particle_exchange_metrics_timer

   real(real64) function lcm_stop_particle_exchange_metrics_timer(timer) &
      result(elapsed_seconds)
      type(lcm_particle_exchange_metrics_timer_t), intent(in) :: timer
      integer(int64) :: count_stop
      integer(int64) :: elapsed_counts

      call system_clock(count_stop)
      elapsed_seconds = 0.0_real64
      if (timer%count_rate <= 0_int64) return

      if (count_stop >= timer%count_start) then
         elapsed_counts = count_stop - timer%count_start
      else
         elapsed_counts = (timer%count_max - timer%count_start) + &
                          count_stop + 1_int64
      endif
      elapsed_seconds = real(elapsed_counts, real64) / &
                        real(timer%count_rate, real64)
   end function lcm_stop_particle_exchange_metrics_timer

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

   subroutine lcm_record_particle_count_exchange(phase, elapsed_seconds)
      integer, intent(in) :: phase
      real(real64), intent(in) :: elapsed_seconds

      call lcm_require_valid_exchange_phase(phase)
      if (.not. metrics_enabled) return

      select case (phase)
      case (lcm_exchange_x_phase)
         particle_exchange_metrics%x_count_exchange_seconds = &
            particle_exchange_metrics%x_count_exchange_seconds + elapsed_seconds
         particle_exchange_metrics%x_count_exchange_calls = &
            particle_exchange_metrics%x_count_exchange_calls + 1_int64
      case (lcm_exchange_y_phase)
         particle_exchange_metrics%y_count_exchange_seconds = &
            particle_exchange_metrics%y_count_exchange_seconds + elapsed_seconds
         particle_exchange_metrics%y_count_exchange_calls = &
            particle_exchange_metrics%y_count_exchange_calls + 1_int64
      end select
   end subroutine lcm_record_particle_count_exchange

   subroutine lcm_record_particle_payload_exchange(phase, elapsed_seconds)
      integer, intent(in) :: phase
      real(real64), intent(in) :: elapsed_seconds

      call lcm_require_valid_exchange_phase(phase)
      if (.not. metrics_enabled) return

      select case (phase)
      case (lcm_exchange_x_phase)
         particle_exchange_metrics%x_payload_exchange_seconds = &
            particle_exchange_metrics%x_payload_exchange_seconds + elapsed_seconds
         particle_exchange_metrics%x_payload_exchange_calls = &
            particle_exchange_metrics%x_payload_exchange_calls + 1_int64
      case (lcm_exchange_y_phase)
         particle_exchange_metrics%y_payload_exchange_seconds = &
            particle_exchange_metrics%y_payload_exchange_seconds + elapsed_seconds
         particle_exchange_metrics%y_payload_exchange_calls = &
            particle_exchange_metrics%y_payload_exchange_calls + 1_int64
      end select
   end subroutine lcm_record_particle_payload_exchange

   subroutine lcm_record_received_particle_insertion(phase, elapsed_seconds, &
                                                      particles_inserted)
      integer, intent(in) :: phase
      real(real64), intent(in) :: elapsed_seconds
      integer(int64), intent(in) :: particles_inserted

      call lcm_require_valid_exchange_phase(phase)
      if (.not. metrics_enabled) return

      select case (phase)
      case (lcm_exchange_x_phase)
         particle_exchange_metrics%x_received_insertion_seconds = &
            particle_exchange_metrics%x_received_insertion_seconds + elapsed_seconds
         particle_exchange_metrics%x_received_insertion_calls = &
            particle_exchange_metrics%x_received_insertion_calls + 1_int64
         particle_exchange_metrics%x_particles_inserted = &
            particle_exchange_metrics%x_particles_inserted + particles_inserted
      case (lcm_exchange_y_phase)
         particle_exchange_metrics%y_received_insertion_seconds = &
            particle_exchange_metrics%y_received_insertion_seconds + elapsed_seconds
         particle_exchange_metrics%y_received_insertion_calls = &
            particle_exchange_metrics%y_received_insertion_calls + 1_int64
         particle_exchange_metrics%y_particles_inserted = &
            particle_exchange_metrics%y_particles_inserted + particles_inserted
      end select
   end subroutine lcm_record_received_particle_insertion

   subroutine lcm_record_particle_ownership_update(elapsed_seconds, &
                                                   particle_examinations, &
                                                   particles_moved)
      real(real64), intent(in) :: elapsed_seconds
      integer(int64), intent(in) :: particle_examinations
      integer(int64), intent(in) :: particles_moved

      if (.not. metrics_enabled) return

      particle_exchange_metrics%ownership_update_seconds = &
         particle_exchange_metrics%ownership_update_seconds + elapsed_seconds
      particle_exchange_metrics%ownership_update_calls = &
         particle_exchange_metrics%ownership_update_calls + 1_int64
      particle_exchange_metrics%ownership_particle_examinations = &
         particle_exchange_metrics%ownership_particle_examinations + &
         particle_examinations
      particle_exchange_metrics%ownership_particles_moved = &
         particle_exchange_metrics%ownership_particles_moved + particles_moved
   end subroutine lcm_record_particle_ownership_update

end module lcm_particle_exchange_metrics
