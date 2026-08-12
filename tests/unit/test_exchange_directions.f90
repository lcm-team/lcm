program test_exchange_directions
   use lcm_exchange_directions, only: lcm_exchange_direction_for_phase, &
                                      lcm_exchange_direction_from_offset, &
                                      lcm_exchange_opposite_direction, &
                                      lcm_exchange_east, lcm_exchange_west, &
                                      lcm_exchange_south, lcm_exchange_north, &
                                      lcm_exchange_x_phase, lcm_exchange_y_phase

   implicit none

   call assert_direction(lcm_exchange_x_phase, 1, lcm_exchange_east, lcm_exchange_west)
   call assert_direction(lcm_exchange_x_phase, -1, lcm_exchange_west, lcm_exchange_east)
   call assert_direction(lcm_exchange_y_phase, -1, lcm_exchange_south, lcm_exchange_north)
   call assert_direction(lcm_exchange_y_phase, 1, lcm_exchange_north, lcm_exchange_south)

   call assert_true(lcm_exchange_direction_for_phase(lcm_exchange_x_phase, 1) == &
                    lcm_exchange_east, 'first x-phase direction')
   call assert_true(lcm_exchange_direction_for_phase(lcm_exchange_x_phase, 2) == &
                    lcm_exchange_west, 'second x-phase direction')
   call assert_true(lcm_exchange_direction_for_phase(lcm_exchange_y_phase, 1) == &
                    lcm_exchange_south, 'first y-phase direction')
   call assert_true(lcm_exchange_direction_for_phase(lcm_exchange_y_phase, 2) == &
                    lcm_exchange_north, 'second y-phase direction')

contains

   subroutine assert_direction(phase, offset, expected, opposite)
      integer, intent(in) :: phase
      integer, intent(in) :: offset
      integer, intent(in) :: expected
      integer, intent(in) :: opposite

      call assert_true(lcm_exchange_direction_from_offset(phase, offset) == expected, &
                       'direction from phase offset')
      call assert_true(lcm_exchange_opposite_direction(expected) == opposite, &
                       'opposite direction')
   end subroutine assert_direction

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) error stop 'assertion failed: ' // message
   end subroutine assert_true

end program test_exchange_directions
