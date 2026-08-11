!> Direction identifiers and neighbor lookup for particle MPI exchange.
module lcm_exchange_directions
   use lcm_parallel_state, only: parallel_state

   implicit none

   private

   integer, parameter, public :: lcm_exchange_east = 1
   integer, parameter, public :: lcm_exchange_west = 2
   integer, parameter, public :: lcm_exchange_southeast = 3
   integer, parameter, public :: lcm_exchange_northwest = 4
   integer, parameter, public :: lcm_exchange_south = 5
   integer, parameter, public :: lcm_exchange_north = 6
   integer, parameter, public :: lcm_exchange_southwest = 7
   integer, parameter, public :: lcm_exchange_northeast = 8
   integer, parameter, public :: lcm_exchange_direction_count = 8

   public :: lcm_exchange_direction_from_offsets
   public :: lcm_exchange_neighbor_rank
   public :: lcm_exchange_opposite_direction

contains

   integer function lcm_exchange_direction_from_offsets(x_offset, y_offset) result(direction)
      integer, intent(in) :: x_offset
      integer, intent(in) :: y_offset

      select case (x_offset)
      case (1)
         select case (y_offset)
         case (0)
            direction = lcm_exchange_east
         case (-1)
            direction = lcm_exchange_southeast
         case (1)
            direction = lcm_exchange_northeast
         case default
            error stop 'LCM exchange x offset must be -1, 0, or 1'
         end select
      case (-1)
         select case (y_offset)
         case (0)
            direction = lcm_exchange_west
         case (-1)
            direction = lcm_exchange_southwest
         case (1)
            direction = lcm_exchange_northwest
         case default
            error stop 'LCM exchange x offset must be -1, 0, or 1'
         end select
      case (0)
         select case (y_offset)
         case (-1)
            direction = lcm_exchange_south
         case (1)
            direction = lcm_exchange_north
         case default
            error stop 'LCM exchange direction cannot be zero or exceed one cell'
         end select
      case default
         error stop 'LCM exchange x offset must be -1, 0, or 1'
      end select
   end function lcm_exchange_direction_from_offsets

   integer function lcm_exchange_neighbor_rank(direction) result(rank)
      integer, intent(in) :: direction

      select case (direction)
      case (lcm_exchange_east)
         rank = parallel_state%parallel%east_rank
      case (lcm_exchange_west)
         rank = parallel_state%parallel%west_rank
      case (lcm_exchange_south)
         rank = parallel_state%parallel%south_rank
      case (lcm_exchange_north)
         rank = parallel_state%parallel%north_rank
      case (lcm_exchange_southeast)
         rank = parallel_state%parallel%southeast_rank
      case (lcm_exchange_northwest)
         rank = parallel_state%parallel%northwest_rank
      case (lcm_exchange_southwest)
         rank = parallel_state%parallel%southwest_rank
      case (lcm_exchange_northeast)
         rank = parallel_state%parallel%northeast_rank
      case default
         error stop 'LCM exchange direction is invalid'
      end select
   end function lcm_exchange_neighbor_rank

   integer function lcm_exchange_opposite_direction(direction) result(opposite)
      integer, intent(in) :: direction

      select case (direction)
      case (lcm_exchange_east)
         opposite = lcm_exchange_west
      case (lcm_exchange_west)
         opposite = lcm_exchange_east
      case (lcm_exchange_southeast)
         opposite = lcm_exchange_northwest
      case (lcm_exchange_northwest)
         opposite = lcm_exchange_southeast
      case (lcm_exchange_south)
         opposite = lcm_exchange_north
      case (lcm_exchange_north)
         opposite = lcm_exchange_south
      case (lcm_exchange_southwest)
         opposite = lcm_exchange_northeast
      case (lcm_exchange_northeast)
         opposite = lcm_exchange_southwest
      case default
         error stop 'LCM exchange direction is invalid'
      end select
   end function lcm_exchange_opposite_direction

end module lcm_exchange_directions
