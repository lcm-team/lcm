!> MPI count exchange for directional particle-message buffers.
module lcm_particle_count_exchange
   use mpi_f08, only: MPI_INTEGER, MPI_PROC_NULL, MPI_SENDRECV, MPI_Status, &
                      MPI_SUCCESS
   use lcm_exchange_directions, only: lcm_exchange_direction_count, &
                                      lcm_exchange_neighbor_rank, &
                                      lcm_exchange_opposite_direction
   use lcm_parallel_state, only: parallel_state
   use lcm_particle_exchange_buffers, only: lcm_reserve_receive_particle_messages, &
                                            lcm_send_particle_message_count

   implicit none

   private

   public :: lcm_exchange_particle_counts

contains

   !> Exchange one count per outgoing direction.
   !>
   !> A count sent toward direction D is received by the neighbor on its
   !> opposite side.  Therefore, the count received with tag D belongs in the
   !> local receive buffer for opposite(D), which denotes the source direction.
   subroutine lcm_exchange_particle_counts()
      integer :: destination_rank
      integer :: direction
      integer :: ierr
      integer :: receive_count
      integer :: receive_direction
      integer :: send_count
      integer :: source_rank
      type(MPI_Status) :: status

      if (.not. parallel_state%is_set) return

      do direction = 1, lcm_exchange_direction_count
         receive_direction = lcm_exchange_opposite_direction(direction)
         destination_rank = lcm_exchange_neighbor_rank(direction)
         source_rank = lcm_exchange_neighbor_rank(receive_direction)
         send_count = lcm_send_particle_message_count(direction)
         receive_count = 0

         if (destination_rank < 0) destination_rank = MPI_PROC_NULL
         if (source_rank < 0) source_rank = MPI_PROC_NULL

         call MPI_SENDRECV(send_count, 1, MPI_INTEGER, destination_rank, direction, &
                           receive_count, 1, MPI_INTEGER, source_rank, direction, &
                           parallel_state%parallel%communicator, status, ierr)
         if (ierr /= MPI_SUCCESS) then
            error stop 'LCM failed to exchange particle-message counts'
         endif
         call lcm_reserve_receive_particle_messages(receive_direction, receive_count)
      enddo
   end subroutine lcm_exchange_particle_counts

end module lcm_particle_count_exchange
