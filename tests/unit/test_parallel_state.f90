program test_parallel_state
   use lcm_parallel_state, only: lcm_set_parallel, parallel_state
   use lcm_types, only: lcm_init_parallel, lcm_parallel_t

   implicit none

   type(lcm_parallel_t) :: parallel

   call lcm_init_parallel(parallel, communicator=17)
   call lcm_set_parallel(parallel)

   call assert_true(parallel_state%is_set, 'parallel state is set')
   call assert_true(parallel_state%parallel%communicator == 17, 'communicator')

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         error stop 'assertion failed: ' // message
      endif
   end subroutine assert_true

end program test_parallel_state
