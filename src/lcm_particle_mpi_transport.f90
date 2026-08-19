!> Low-level typed MPI transport for directional particle-message payloads.
module lcm_particle_mpi_transport
   use iso_fortran_env, only: int64
   use mpi_f08, only: MPI_IRECV, MPI_ISEND, MPI_PROC_NULL, MPI_Request, &
                      MPI_Status, MPI_SUCCESS, MPI_TYPE_SIZE, MPI_WAITALL
   use lcm_exchange_directions, only: lcm_exchange_direction_for_phase, &
                                      lcm_exchange_directions_per_phase, &
                                      lcm_exchange_neighbor_rank, &
                                      lcm_exchange_opposite_direction, &
                                      lcm_require_valid_exchange_phase
   use lcm_parallel_state, only: parallel_state
   use lcm_particle_exchange_buffers, only: &
      lcm_expected_receive_particle_message_count, &
      lcm_get_receive_particle_message_view, &
      lcm_get_send_particle_message_view, &
      lcm_mark_receive_particle_messages_received, &
      lcm_send_particle_message_count
   use lcm_particle_mpi_datatype, only: lcm_get_particle_mpi_datatype
   use lcm_particle_exchange_metrics, only: &
      lcm_particle_exchange_metrics_are_enabled, &
      lcm_particle_exchange_metrics_timer_t, &
      lcm_record_particle_payload_exchange, &
      lcm_record_particle_payload_sent, &
      lcm_start_particle_exchange_metrics_timer, &
      lcm_stop_particle_exchange_metrics_timer
   use lcm_types, only: lcm_particle_mpi_message_t

   implicit none

   private

   public :: lcm_exchange_particle_payloads

   integer, parameter :: particle_payload_tag_base = 1000

contains

   !> Exchange both directional payload buffers belonging to one phase.
   !>
   !> Receives are posted before sends and all operations complete together.
   !> MPI reads and writes the reusable directional buffers directly; no
   !> transport-only copy of the particle payload is created.
   subroutine lcm_exchange_particle_payloads(phase)
      integer, intent(in) :: phase
      integer :: destination_rank
      integer :: direction
      integer :: ierr
      integer :: phase_direction
      integer :: receive_count
      integer :: receive_direction
      integer :: request_count
      integer :: send_counts(lcm_exchange_directions_per_phase)
      integer :: particle_payload_bytes
      integer :: source_rank
      type(lcm_particle_mpi_message_t), contiguous, pointer :: messages(:)
      type(MPI_Request) :: requests( &
         2 * lcm_exchange_directions_per_phase)
      type(MPI_Status) :: statuses( &
         2 * lcm_exchange_directions_per_phase)
      logical :: measure_payload_exchange
      type(lcm_particle_exchange_metrics_timer_t) :: metrics_timer

      if (.not. parallel_state%is_set) return
      call lcm_require_valid_exchange_phase(phase)

      request_count = 0
      send_counts = 0
      particle_payload_bytes = 0
      ! Measure payload request preparation, posting, and completion. Marking
      ! completed receive buffers below remains outside this metric.
      measure_payload_exchange = lcm_particle_exchange_metrics_are_enabled()
      if (measure_payload_exchange) then
         call lcm_start_particle_exchange_metrics_timer(metrics_timer)
      endif
      do phase_direction = 1, lcm_exchange_directions_per_phase
         direction = lcm_exchange_direction_for_phase(phase, phase_direction)
         receive_direction = lcm_exchange_opposite_direction(direction)
         source_rank = lcm_exchange_neighbor_rank(receive_direction)
         receive_count = &
            lcm_expected_receive_particle_message_count(receive_direction)

         if (receive_count > 0) then
            if (source_rank == MPI_PROC_NULL .or. source_rank < 0) then
               error stop 'LCM particle payload has no valid source rank'
            endif
            call lcm_get_receive_particle_message_view(receive_direction, messages)
            request_count = request_count + 1
            call MPI_IRECV(messages, receive_count, &
                           lcm_get_particle_mpi_datatype(), source_rank, &
                           particle_payload_tag_base + direction, &
                           parallel_state%parallel%communicator, &
                           requests(request_count), ierr)
            call require_mpi_success(ierr, 'post a particle-payload receive')
         endif
      enddo

      do phase_direction = 1, lcm_exchange_directions_per_phase
         direction = lcm_exchange_direction_for_phase(phase, phase_direction)
         destination_rank = lcm_exchange_neighbor_rank(direction)
         send_counts(phase_direction) = &
            lcm_send_particle_message_count(direction)

         if (send_counts(phase_direction) > 0) then
            if (destination_rank == MPI_PROC_NULL .or. destination_rank < 0) then
               error stop 'LCM particle payload has no valid destination rank'
            endif
            call lcm_get_send_particle_message_view(direction, messages)
            request_count = request_count + 1
            call MPI_ISEND(messages, send_counts(phase_direction), &
                           lcm_get_particle_mpi_datatype(), destination_rank, &
                           particle_payload_tag_base + direction, &
                           parallel_state%parallel%communicator, &
                           requests(request_count), ierr)
            call require_mpi_success(ierr, 'post a particle-payload send')
         endif
      enddo

      if (request_count > 0) then
         call MPI_WAITALL(request_count, requests(1:request_count), &
                          statuses(1:request_count), ierr)
         call require_mpi_success(ierr, 'complete particle-payload transport')
      endif

      if (measure_payload_exchange) then
         call lcm_record_particle_payload_exchange( &
            phase, lcm_stop_particle_exchange_metrics_timer(metrics_timer))
         if (any(send_counts > 0)) then
            call MPI_TYPE_SIZE(lcm_get_particle_mpi_datatype(), &
                               particle_payload_bytes, ierr)
            call require_mpi_success(ierr, 'query the particle-payload size')
         endif
         do phase_direction = 1, lcm_exchange_directions_per_phase
            direction = lcm_exchange_direction_for_phase(phase, phase_direction)
            call lcm_record_particle_payload_sent( &
               direction, int(send_counts(phase_direction), int64), &
               int(send_counts(phase_direction), int64) * &
                  int(particle_payload_bytes, int64))
         enddo
      endif

      do phase_direction = 1, lcm_exchange_directions_per_phase
         direction = lcm_exchange_direction_for_phase(phase, phase_direction)
         receive_direction = lcm_exchange_opposite_direction(direction)
         call lcm_mark_receive_particle_messages_received(receive_direction)
      enddo
   end subroutine lcm_exchange_particle_payloads

   subroutine require_mpi_success(ierr, operation)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: operation

      if (ierr /= MPI_SUCCESS) error stop 'LCM failed to ' // operation
   end subroutine require_mpi_success

end module lcm_particle_mpi_transport
