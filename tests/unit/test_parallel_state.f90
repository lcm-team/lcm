program test_parallel_state
   use lcm_parallel_state, only: lcm_set_parallel, parallel_state
   use lcm_types, only: lcm_init_parallel, lcm_parallel_t

   implicit none

   include 'mpif.h'

   type(lcm_parallel_t) :: parallel
   integer :: communicator
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
   call assert_true(parallel_state%parallel%communicator == communicator, 'communicator')
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
