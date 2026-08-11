program test_parallel_state
   use mpi_f08, only: MPI_CART_CREATE, MPI_COMM_COMPARE, MPI_COMM_FREE, &
                      MPI_COMM_WORLD, MPI_Comm, MPI_FINALIZE, MPI_IDENT, &
                      MPI_INIT
   use lcm_parallel_state, only: lcm_init_parallel, lcm_parallel_t, &
                                 lcm_set_parallel, parallel_state

   implicit none

   type(lcm_parallel_t) :: parallel
   type(MPI_Comm) :: communicator
   integer :: comparison
   integer :: dimensions(2)
   integer :: ierr
   logical :: periodicity(2)

   call MPI_INIT(ierr)
   dimensions = [1, 1]
   periodicity = [.true., .true.]
   call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dimensions, periodicity, .true., &
                        communicator, ierr)

   call lcm_init_parallel(parallel, communicator)
   call lcm_set_parallel(parallel)

   call assert_true(parallel_state%is_set, 'parallel state is set')
   call MPI_COMM_COMPARE(parallel_state%parallel%communicator, communicator, &
                         comparison, ierr)
   call assert_true(comparison == MPI_IDENT, 'communicator')
   call assert_true(parallel_state%southwest_rank == 0, 'south-west rank')
   call assert_true(parallel_state%southeast_rank == 0, 'south-east rank')
   call assert_true(parallel_state%northwest_rank == 0, 'north-west rank')
   call assert_true(parallel_state%northeast_rank == 0, 'north-east rank')

   call MPI_COMM_FREE(communicator, ierr)
   call MPI_FINALIZE(ierr)

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         error stop 'assertion failed: ' // message
      endif
   end subroutine assert_true

end program test_parallel_state
