!> MPI count exchange for directional particle-message buffers.
module lcm_particle_count_exchange
   use mpi_f08, only: MPI_INTEGER, MPI_IRECV, MPI_ISEND, MPI_PROC_NULL, &
                      MPI_Request, MPI_Status, MPI_SUCCESS, MPI_WAITALL
   use lcm_exchange_directions, only: lcm_exchange_direction_for_phase, &
                                      lcm_exchange_directions_per_phase, &
                                      lcm_exchange_neighbor_rank, &
                                      lcm_exchange_opposite_direction, &
                                      lcm_require_valid_exchange_phase
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
   !> Both receives are posted before both sends, and the complete cardinal
   !> phase is synchronized with one wait.
   subroutine lcm_exchange_particle_counts(phase)
      integer, intent(in) :: phase
      integer :: destination_ranks(lcm_exchange_directions_per_phase)
      integer :: directions(lcm_exchange_directions_per_phase)
      integer :: ierr
      integer :: phase_direction
      integer, asynchronous :: receive_counts(lcm_exchange_directions_per_phase)
      integer :: receive_directions(lcm_exchange_directions_per_phase)
      integer, parameter :: request_count = 2 * lcm_exchange_directions_per_phase
      type(MPI_Request) :: requests(request_count)
      integer, asynchronous :: send_counts(lcm_exchange_directions_per_phase)
      integer :: source_ranks(lcm_exchange_directions_per_phase)
      type(MPI_Status) :: statuses(request_count)

      if (.not. parallel_state%is_set) return
      call lcm_require_valid_exchange_phase(phase)

      receive_counts = 0
      do phase_direction = 1, lcm_exchange_directions_per_phase
         directions(phase_direction) = &
            lcm_exchange_direction_for_phase(phase, phase_direction)
         receive_directions(phase_direction) = &
            lcm_exchange_opposite_direction(directions(phase_direction))
         destination_ranks(phase_direction) = &
            lcm_exchange_neighbor_rank(directions(phase_direction))
         source_ranks(phase_direction) = &
            lcm_exchange_neighbor_rank(receive_directions(phase_direction))
         send_counts(phase_direction) = &
            lcm_send_particle_message_count(directions(phase_direction))

         if (destination_ranks(phase_direction) < 0) then
            destination_ranks(phase_direction) = MPI_PROC_NULL
         endif
         if (source_ranks(phase_direction) < 0) then
            source_ranks(phase_direction) = MPI_PROC_NULL
         endif
      enddo

      do phase_direction = 1, lcm_exchange_directions_per_phase
         call MPI_IRECV(receive_counts(phase_direction), 1, MPI_INTEGER, &
                        source_ranks(phase_direction), directions(phase_direction), &
                        parallel_state%parallel%communicator, &
                        requests(phase_direction), ierr)
         call require_mpi_success(ierr, 'post a particle-count receive')
      enddo

      do phase_direction = 1, lcm_exchange_directions_per_phase
         call MPI_ISEND(send_counts(phase_direction), 1, MPI_INTEGER, &
                        destination_ranks(phase_direction), directions(phase_direction), &
                        parallel_state%parallel%communicator, &
                        requests(lcm_exchange_directions_per_phase + phase_direction), &
                        ierr)
         call require_mpi_success(ierr, 'post a particle-count send')
      enddo

      call MPI_WAITALL(request_count, requests, statuses, ierr)
      call require_mpi_success(ierr, 'complete particle-count exchange')

      do phase_direction = 1, lcm_exchange_directions_per_phase
         call lcm_reserve_receive_particle_messages( &
            receive_directions(phase_direction), receive_counts(phase_direction))
      enddo
   end subroutine lcm_exchange_particle_counts

   subroutine require_mpi_success(ierr, operation)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: operation

      if (ierr /= MPI_SUCCESS) error stop 'LCM failed to ' // operation
   end subroutine require_mpi_success

end module lcm_particle_count_exchange
