!> High-level orchestration of one cardinal particle-exchange phase.
module lcm_particle_exchange
   use lcm_exchange_directions, only: lcm_require_valid_exchange_phase
   use lcm_particle_count_exchange, only: lcm_exchange_particle_counts
   use lcm_particle_mpi_transport, only: lcm_exchange_particle_payloads
   use lcm_particle_state, only: lcm_extract_outgoing_particles, &
                                 lcm_insert_received_particles

   implicit none

   private

   public :: lcm_exchange_particle_phase

contains

   subroutine lcm_exchange_particle_phase(phase)
      integer, intent(in) :: phase

      call lcm_require_valid_exchange_phase(phase)
      call lcm_extract_outgoing_particles(phase)
      call lcm_exchange_particle_counts(phase)
      call lcm_exchange_particle_payloads(phase)
      call lcm_insert_received_particles(phase)
   end subroutine lcm_exchange_particle_phase

end module lcm_particle_exchange
