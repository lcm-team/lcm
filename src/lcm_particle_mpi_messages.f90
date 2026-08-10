!> Typed temporary messages used by standalone particle MPI exchange.
module lcm_particle_mpi_messages
   use lcm_types, only: lcm_particle_mpi_message_t, lcm_particle_t

   implicit none

   private

   public :: lcm_pack_particle_mpi_message
   public :: lcm_unpack_particle_mpi_message

contains

   !> Assemble particle data and temporary destination metadata for MPI transfer.
   !> This is a typed copy, not the later byte-level MPI serialization step.
   subroutine lcm_pack_particle_mpi_message(global_i, global_j, k, particle, message)
      integer, intent(in) :: global_i
      integer, intent(in) :: global_j
      integer, intent(in) :: k
      type(lcm_particle_t), intent(in) :: particle
      type(lcm_particle_mpi_message_t), intent(out) :: message

      message%global_i = global_i
      message%global_j = global_j
      message%k = k
      message%particle = particle
   end subroutine lcm_pack_particle_mpi_message

   !> Split a received typed MPI message into cell metadata and particle data.
   subroutine lcm_unpack_particle_mpi_message(message, global_i, global_j, k, particle)
      type(lcm_particle_mpi_message_t), intent(in) :: message
      integer, intent(out) :: global_i
      integer, intent(out) :: global_j
      integer, intent(out) :: k
      type(lcm_particle_t), intent(out) :: particle

      global_i = message%global_i
      global_j = message%global_j
      k = message%k
      particle = message%particle
   end subroutine lcm_unpack_particle_mpi_message

end module lcm_particle_mpi_messages
