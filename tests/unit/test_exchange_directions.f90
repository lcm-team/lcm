program test_exchange_directions
   use lcm_exchange_directions, only: lcm_exchange_direction_from_offsets, &
                                      lcm_exchange_opposite_direction, &
                                      lcm_exchange_east, lcm_exchange_west, &
                                      lcm_exchange_southeast, lcm_exchange_northwest, &
                                      lcm_exchange_south, lcm_exchange_north, &
                                      lcm_exchange_southwest, lcm_exchange_northeast

   implicit none

   call assert_direction(1, 0, lcm_exchange_east, lcm_exchange_west)
   call assert_direction(-1, 0, lcm_exchange_west, lcm_exchange_east)
   call assert_direction(1, -1, lcm_exchange_southeast, lcm_exchange_northwest)
   call assert_direction(-1, 1, lcm_exchange_northwest, lcm_exchange_southeast)
   call assert_direction(0, -1, lcm_exchange_south, lcm_exchange_north)
   call assert_direction(0, 1, lcm_exchange_north, lcm_exchange_south)
   call assert_direction(-1, -1, lcm_exchange_southwest, lcm_exchange_northeast)
   call assert_direction(1, 1, lcm_exchange_northeast, lcm_exchange_southwest)

contains

   subroutine assert_direction(x_offset, y_offset, expected, opposite)
      integer, intent(in) :: x_offset
      integer, intent(in) :: y_offset
      integer, intent(in) :: expected
      integer, intent(in) :: opposite

      call assert_true(lcm_exchange_direction_from_offsets(x_offset, y_offset) == expected, &
                       'direction from offsets')
      call assert_true(lcm_exchange_opposite_direction(expected) == opposite, &
                       'opposite direction')
   end subroutine assert_direction

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) error stop 'assertion failed: ' // message
   end subroutine assert_true

end program test_exchange_directions
