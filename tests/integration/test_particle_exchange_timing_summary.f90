program test_particle_exchange_timing_summary
   use iso_fortran_env, only: int64, real64
   use mpi_f08, only: MPI_CART_CREATE, MPI_COMM_FREE, MPI_COMM_RANK, &
                      MPI_COMM_SIZE, MPI_COMM_WORLD, MPI_Comm, MPI_FINALIZE, &
                      MPI_INIT
   use lcm_exchange_directions, only: lcm_exchange_x_phase, &
                                      lcm_exchange_y_phase
   use lcm_parallel_state, only: lcm_init_parallel, lcm_parallel_t, &
                                 lcm_set_parallel
   use lcm_particle_exchange_metrics, only: &
      lcm_record_particle_count_exchange, &
      lcm_record_particle_outgoing_preparation, &
      lcm_record_particle_ownership_update, &
      lcm_record_particle_payload_exchange, &
      lcm_record_received_particle_insertion, &
      lcm_reset_particle_exchange_metrics, &
      lcm_set_particle_exchange_metrics_enabled
   use lcm_particle_exchange_timing_summary, only: &
      lcm_get_particle_exchange_timing_summary, &
      lcm_particle_exchange_timing_summary_t, &
      lcm_rank_timing_summary_t

   implicit none

   type(MPI_Comm) :: communicator
   integer :: dimensions(2)
   real(real64) :: factor
   integer :: ierr
   integer :: rank
   integer :: size
   logical :: periodicity(2)
   type(lcm_parallel_t) :: parallel
   logical :: serial_mode
   type(lcm_particle_exchange_timing_summary_t) :: summary
   character(len=16) :: test_mode

   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

   test_mode = ''
   if (command_argument_count() > 0) call get_command_argument(1, test_mode)
   serial_mode = trim(test_mode) == 'serial'

   if (serial_mode) then
      if (size /= 1) error stop 'serial timing summary test requires one rank'
      factor = 2.0_real64
   else
      if (size /= 4) error stop 'MPI timing summary test requires four ranks'
      dimensions = [2, 2]
      periodicity = [.false., .false.]
      call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dimensions, periodicity, .false., &
                           communicator, ierr)
      call MPI_COMM_RANK(communicator, rank, ierr)
      call lcm_init_parallel(parallel, communicator)
      call lcm_set_parallel(parallel)
      factor = real(rank + 1, real64)
   endif

   call lcm_reset_particle_exchange_metrics()
   call lcm_set_particle_exchange_metrics_enabled(.true.)
   call record_synthetic_timings(factor)
   call lcm_get_particle_exchange_timing_summary(summary)

   if (serial_mode) then
      call assert_equal(summary%rank_count, 1, 'serial summary rank count')
      call assert_all_timings(summary, 2.0_real64, 2.0_real64, 2.0_real64)
   else
      call assert_equal(summary%rank_count, 4, 'MPI summary rank count')
      call assert_all_timings(summary, 1.0_real64, 4.0_real64, 2.5_real64)
      call MPI_COMM_FREE(communicator, ierr)
   endif

   call MPI_FINALIZE(ierr)

contains

   subroutine record_synthetic_timings(rank_factor)
      real(real64), intent(in) :: rank_factor

      call lcm_record_particle_outgoing_preparation( &
         lcm_exchange_x_phase, 1.0_real64 * rank_factor, 0_int64, 0_int64)
      call lcm_record_particle_outgoing_preparation( &
         lcm_exchange_y_phase, 2.0_real64 * rank_factor, 0_int64, 0_int64)
      call lcm_record_particle_count_exchange( &
         lcm_exchange_x_phase, 3.0_real64 * rank_factor)
      call lcm_record_particle_count_exchange( &
         lcm_exchange_y_phase, 4.0_real64 * rank_factor)
      call lcm_record_particle_payload_exchange( &
         lcm_exchange_x_phase, 5.0_real64 * rank_factor)
      call lcm_record_particle_payload_exchange( &
         lcm_exchange_y_phase, 6.0_real64 * rank_factor)
      call lcm_record_received_particle_insertion( &
         lcm_exchange_x_phase, 7.0_real64 * rank_factor, 0_int64)
      call lcm_record_received_particle_insertion( &
         lcm_exchange_y_phase, 8.0_real64 * rank_factor, 0_int64)
      call lcm_record_particle_ownership_update( &
         9.0_real64 * rank_factor, 0_int64, 0_int64)
   end subroutine record_synthetic_timings

   subroutine assert_all_timings(actual, minimum_factor, maximum_factor, &
                                 mean_factor)
      type(lcm_particle_exchange_timing_summary_t), intent(in) :: actual
      real(real64), intent(in) :: minimum_factor
      real(real64), intent(in) :: maximum_factor
      real(real64), intent(in) :: mean_factor

      call assert_timing(actual%x_outgoing_preparation, 1.0_real64, &
                         minimum_factor, maximum_factor, mean_factor, &
                         'x outgoing preparation')
      call assert_timing(actual%y_outgoing_preparation, 2.0_real64, &
                         minimum_factor, maximum_factor, mean_factor, &
                         'y outgoing preparation')
      call assert_timing(actual%x_count_exchange, 3.0_real64, minimum_factor, &
                         maximum_factor, mean_factor, 'x count exchange')
      call assert_timing(actual%y_count_exchange, 4.0_real64, minimum_factor, &
                         maximum_factor, mean_factor, 'y count exchange')
      call assert_timing(actual%x_payload_exchange, 5.0_real64, minimum_factor, &
                         maximum_factor, mean_factor, 'x payload exchange')
      call assert_timing(actual%y_payload_exchange, 6.0_real64, minimum_factor, &
                         maximum_factor, mean_factor, 'y payload exchange')
      call assert_timing(actual%x_received_insertion, 7.0_real64, &
                         minimum_factor, maximum_factor, mean_factor, &
                         'x received insertion')
      call assert_timing(actual%y_received_insertion, 8.0_real64, &
                         minimum_factor, maximum_factor, mean_factor, &
                         'y received insertion')
      call assert_timing(actual%ownership_update, 9.0_real64, minimum_factor, &
                         maximum_factor, mean_factor, 'ownership update')
   end subroutine assert_all_timings

   subroutine assert_timing(actual, coefficient, minimum_factor, &
                            maximum_factor, mean_factor, description)
      type(lcm_rank_timing_summary_t), intent(in) :: actual
      real(real64), intent(in) :: coefficient
      real(real64), intent(in) :: minimum_factor
      real(real64), intent(in) :: maximum_factor
      real(real64), intent(in) :: mean_factor
      character(len=*), intent(in) :: description

      call assert_close(actual%minimum_seconds, coefficient * minimum_factor, &
                        description // ' minimum')
      call assert_close(actual%maximum_seconds, coefficient * maximum_factor, &
                        description // ' maximum')
      call assert_close(actual%mean_seconds, coefficient * mean_factor, &
                        description // ' mean')
   end subroutine assert_timing

   subroutine assert_equal(actual, expected, description)
      integer, intent(in) :: actual
      integer, intent(in) :: expected
      character(len=*), intent(in) :: description

      if (actual /= expected) error stop 'assertion failed: ' // description
   end subroutine assert_equal

   subroutine assert_close(actual, expected, description)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: description

      if (abs(actual - expected) > 1.0e-12_real64) then
         error stop 'assertion failed: ' // description
      endif
   end subroutine assert_close

end program test_particle_exchange_timing_summary
