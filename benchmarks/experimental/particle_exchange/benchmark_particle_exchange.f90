program benchmark_particle_exchange
   use iso_fortran_env, only: int64, real64
   use mpi_f08, only: MPI_ADDRESS_KIND, MPI_BARRIER, MPI_CART_CREATE, &
                      MPI_COMM_FREE, MPI_COMM_RANK, MPI_COMM_SIZE, &
                      MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_DIMS_CREATE, &
                      MPI_FINALIZE, MPI_INIT, MPI_MAX, MPI_REDUCE, MPI_SUCCESS, &
                      MPI_TYPE_GET_EXTENT, MPI_TYPE_SIZE, MPI_WTIME, MPI_Comm
   use lcm_exchange_directions, only: lcm_exchange_direction_for_phase, &
                                      lcm_exchange_directions_per_phase, &
                                      lcm_exchange_phase_count, &
                                      lcm_exchange_x_phase
   use lcm_parallel_state, only: lcm_init_parallel, lcm_parallel_t, &
                                 lcm_set_parallel
   use lcm_particle_count_exchange, only: lcm_exchange_particle_counts
   use lcm_particle_exchange_buffers, only: &
      lcm_append_send_particle_message, &
      lcm_clear_particle_exchange_phase_buffers, &
      lcm_receive_particle_message_count
   use lcm_particle_mpi_datatype, only: lcm_free_particle_mpi_datatype, &
                                        lcm_get_particle_mpi_datatype
   use lcm_particle_mpi_transport, only: lcm_exchange_particle_payloads
   use lcm_types, only: lcm_particle_mpi_message_t

   implicit none

   integer, parameter :: default_iterations = 20
   integer, parameter :: default_maximum_message_count = 10000
   integer, parameter :: warmup_iterations = 3
   integer, parameter :: candidate_message_counts(6) = &
      [1, 10, 100, 1000, 10000, 100000]

   type(MPI_Comm) :: communicator
   integer :: dimensions(2)
   integer :: ierr
   integer :: iterations
   integer :: maximum_message_count
   integer :: message_count_index
   integer :: mpi_size
   integer :: phase
   type(lcm_parallel_t) :: parallel
   logical :: periodicity(2)
   integer :: rank

   call MPI_INIT(ierr)
   call require_mpi_success(ierr, 'initialize MPI')
   call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
   call require_mpi_success(ierr, 'query the MPI size')

   if (mpi_size < 2) then
      error stop 'particle-exchange benchmark requires at least two MPI ranks'
   endif

   dimensions = 0
   call MPI_DIMS_CREATE(mpi_size, 2, dimensions, ierr)
   call require_mpi_success(ierr, 'create Cartesian dimensions')
   periodicity = [.true., .true.]
   call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dimensions, periodicity, .false., &
                        communicator, ierr)
   call require_mpi_success(ierr, 'create the Cartesian communicator')
   call MPI_COMM_RANK(communicator, rank, ierr)
   call require_mpi_success(ierr, 'query the Cartesian rank')

   call lcm_init_parallel(parallel, communicator)
   call lcm_set_parallel(parallel)
   call read_positive_argument(1, default_iterations, iterations)
   call read_positive_argument(2, default_maximum_message_count, &
                               maximum_message_count)

   if (rank == 0) then
      write(*,'(a)') 'ranks,dimensions_x,dimensions_y,phase,' // &
         'messages_per_direction,iterations,message_payload_bytes,' // &
         'message_extent_bytes,count_seconds,payload_seconds,total_seconds,' // &
         'payload_messages_per_second,payload_mib_per_second'
   endif

   do message_count_index = 1, size(candidate_message_counts)
      if (candidate_message_counts(message_count_index) > &
          maximum_message_count) cycle

      do phase = 1, lcm_exchange_phase_count
         call benchmark_phase(phase, &
                              candidate_message_counts(message_count_index), &
                              iterations)
      enddo
   enddo

   call lcm_free_particle_mpi_datatype()
   call MPI_COMM_FREE(communicator, ierr)
   call require_mpi_success(ierr, 'free the Cartesian communicator')
   call MPI_FINALIZE(ierr)

contains

   subroutine benchmark_phase(phase, messages_per_direction, iteration_count)
      integer, intent(in) :: phase
      integer, intent(in) :: messages_per_direction
      integer, intent(in) :: iteration_count
      real(real64) :: count_elapsed
      real(real64) :: count_elapsed_max
      real(real64) :: count_elapsed_sum
      integer(MPI_ADDRESS_KIND) :: datatype_extent
      integer(MPI_ADDRESS_KIND) :: datatype_lower_bound
      integer :: datatype_size
      integer :: iteration
      real(real64) :: message_rate
      real(real64) :: payload_elapsed
      real(real64) :: payload_elapsed_max
      real(real64) :: payload_elapsed_sum
      real(real64) :: payload_mib_rate
      real(real64) :: start_time
      real(real64) :: total_elapsed
      real(real64) :: total_elapsed_max
      real(real64) :: total_elapsed_sum

      call MPI_TYPE_GET_EXTENT(lcm_get_particle_mpi_datatype(), &
                               datatype_lower_bound, datatype_extent, ierr)
      call require_mpi_success(ierr, 'query the particle-message extent')
      call MPI_TYPE_SIZE(lcm_get_particle_mpi_datatype(), datatype_size, ierr)
      call require_mpi_success(ierr, 'query the particle-message size')

      do iteration = 1, warmup_iterations
         call prepare_directional_messages(phase, messages_per_direction)
         call exchange_once(phase, count_elapsed, payload_elapsed)
         call verify_receive_counts(phase, messages_per_direction)
      enddo

      count_elapsed_sum = 0.0_real64
      payload_elapsed_sum = 0.0_real64
      total_elapsed_sum = 0.0_real64
      do iteration = 1, iteration_count
         call prepare_directional_messages(phase, messages_per_direction)
         call MPI_BARRIER(communicator, ierr)
         call require_mpi_success(ierr, 'synchronize a measured iteration')

         start_time = MPI_WTIME()
         call lcm_exchange_particle_counts(phase)
         count_elapsed = MPI_WTIME() - start_time

         start_time = MPI_WTIME()
         call lcm_exchange_particle_payloads(phase)
         payload_elapsed = MPI_WTIME() - start_time
         total_elapsed = count_elapsed + payload_elapsed

         call MPI_REDUCE(count_elapsed, count_elapsed_max, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MAX, 0, communicator, ierr)
         call require_mpi_success(ierr, 'reduce the count-exchange time')
         call MPI_REDUCE(payload_elapsed, payload_elapsed_max, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MAX, 0, communicator, ierr)
         call require_mpi_success(ierr, 'reduce the payload-exchange time')
         call MPI_REDUCE(total_elapsed, total_elapsed_max, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MAX, 0, communicator, ierr)
         call require_mpi_success(ierr, 'reduce the total exchange time')
         call verify_receive_counts(phase, messages_per_direction)

         if (rank == 0) then
            count_elapsed_sum = count_elapsed_sum + count_elapsed_max
            payload_elapsed_sum = payload_elapsed_sum + payload_elapsed_max
            total_elapsed_sum = total_elapsed_sum + total_elapsed_max
         endif
      enddo

      if (rank /= 0) return

      count_elapsed = count_elapsed_sum / real(iteration_count, real64)
      payload_elapsed = payload_elapsed_sum / real(iteration_count, real64)
      total_elapsed = total_elapsed_sum / real(iteration_count, real64)
      if (payload_elapsed > 0.0_real64) then
         message_rate = real(2_int64 * int(messages_per_direction, int64) * &
                             int(mpi_size, int64), real64) / payload_elapsed
         payload_mib_rate = message_rate * real(datatype_size, real64) / &
                            (1024.0_real64 * 1024.0_real64)
      else
         message_rate = 0.0_real64
         payload_mib_rate = 0.0_real64
      endif

      write(*,'(i0,a,i0,a,i0,a,a,a,i0,a,i0,a,i0,a,i0,a,es16.8,a,' // &
              'es16.8,a,es16.8,a,es16.8,a,es16.8)') &
         mpi_size, ',', dimensions(1), ',', dimensions(2), ',', &
         phase_name(phase), ',', messages_per_direction, ',', iteration_count, &
         ',', datatype_size, ',', datatype_extent, ',', count_elapsed, ',', &
         payload_elapsed, ',', total_elapsed, ',', message_rate, ',', &
         payload_mib_rate
   end subroutine benchmark_phase

   subroutine exchange_once(phase, count_elapsed, payload_elapsed)
      integer, intent(in) :: phase
      real(real64), intent(out) :: count_elapsed
      real(real64), intent(out) :: payload_elapsed
      real(real64) :: start_time

      call MPI_BARRIER(communicator, ierr)
      call require_mpi_success(ierr, 'synchronize a warm-up iteration')
      start_time = MPI_WTIME()
      call lcm_exchange_particle_counts(phase)
      count_elapsed = MPI_WTIME() - start_time
      start_time = MPI_WTIME()
      call lcm_exchange_particle_payloads(phase)
      payload_elapsed = MPI_WTIME() - start_time
   end subroutine exchange_once

   subroutine prepare_directional_messages(phase, messages_per_direction)
      integer, intent(in) :: phase
      integer, intent(in) :: messages_per_direction
      integer :: direction
      integer :: message_index
      type(lcm_particle_mpi_message_t) :: message
      integer :: phase_direction

      call lcm_clear_particle_exchange_phase_buffers(phase)
      do phase_direction = 1, lcm_exchange_directions_per_phase
         direction = lcm_exchange_direction_for_phase(phase, phase_direction)
         do message_index = 1, messages_per_direction
            message = benchmark_message(direction, message_index)
            call lcm_append_send_particle_message(direction, message)
         enddo
      enddo
   end subroutine prepare_directional_messages

   function benchmark_message(direction, message_index) result(message)
      integer, intent(in) :: direction
      integer, intent(in) :: message_index
      type(lcm_particle_mpi_message_t) :: message

      message = lcm_particle_mpi_message_t()
      message%global_i = direction
      message%global_j = rank
      message%k = 1
      message%particle%id = int(message_index, int64)
      message%particle%weight = int(rank + 1, int64)
      message%particle%x = real(message_index, real64)
      message%particle%y = real(direction, real64)
      message%particle%z = real(rank, real64)
      message%particle%temperature = 273.15_real64
      message%particle%mass_water = 1.0e-12_real64
   end function benchmark_message

   subroutine verify_receive_counts(phase, expected_count)
      integer, intent(in) :: phase
      integer, intent(in) :: expected_count
      integer :: direction
      integer :: phase_direction

      do phase_direction = 1, lcm_exchange_directions_per_phase
         direction = lcm_exchange_direction_for_phase(phase, phase_direction)
         if (lcm_receive_particle_message_count(direction) /= expected_count) then
            error stop 'particle-exchange benchmark received an unexpected count'
         endif
      enddo
   end subroutine verify_receive_counts

   function phase_name(phase) result(name)
      integer, intent(in) :: phase
      character(len=1) :: name

      if (phase == lcm_exchange_x_phase) then
         name = 'x'
      else
         name = 'y'
      endif
   end function phase_name

   subroutine read_positive_argument(argument_index, default_value, value)
      integer, intent(in) :: argument_index
      integer, intent(in) :: default_value
      integer, intent(out) :: value
      character(len=32) :: argument
      integer :: read_status

      value = default_value
      if (command_argument_count() < argument_index) return

      call get_command_argument(argument_index, argument)
      read(argument, *, iostat=read_status) value
      if (read_status /= 0 .or. value < 1) then
         error stop 'benchmark arguments must be positive integers'
      endif
   end subroutine read_positive_argument

   subroutine require_mpi_success(mpi_error, operation)
      integer, intent(in) :: mpi_error
      character(len=*), intent(in) :: operation

      if (mpi_error /= MPI_SUCCESS) error stop 'benchmark failed to ' // operation
   end subroutine require_mpi_success

end program benchmark_particle_exchange
