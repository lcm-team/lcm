!> Cardinal direction and phase identifiers for particle MPI exchange.
module lcm_exchange_directions
   use lcm_parallel_state, only: parallel_state

   implicit none

   private

   integer, parameter, public :: lcm_exchange_east = 1
   integer, parameter, public :: lcm_exchange_west = 2
   integer, parameter, public :: lcm_exchange_south = 3
   integer, parameter, public :: lcm_exchange_north = 4
   integer, parameter, public :: lcm_exchange_direction_count = 4

   integer, parameter, public :: lcm_exchange_x_phase = 1
   integer, parameter, public :: lcm_exchange_y_phase = 2
   integer, parameter, public :: lcm_exchange_phase_count = 2
   integer, parameter, public :: lcm_exchange_directions_per_phase = 2

   public :: lcm_exchange_direction_for_phase
   public :: lcm_exchange_direction_from_offset
   public :: lcm_exchange_neighbor_rank
   public :: lcm_exchange_opposite_direction
   public :: lcm_require_valid_exchange_phase

contains

   integer function lcm_exchange_direction_for_phase(phase, phase_direction) &
      result(direction)
      integer, intent(in) :: phase
      integer, intent(in) :: phase_direction

      call lcm_require_valid_exchange_phase(phase)
      if (phase_direction < 1 .or. &
          phase_direction > lcm_exchange_directions_per_phase) then
         error stop 'LCM exchange phase direction is invalid'
      endif

      select case (phase)
      case (lcm_exchange_x_phase)
         if (phase_direction == 1) then
            direction = lcm_exchange_east
         else
            direction = lcm_exchange_west
         endif
      case (lcm_exchange_y_phase)
         if (phase_direction == 1) then
            direction = lcm_exchange_south
         else
            direction = lcm_exchange_north
         endif
      end select
   end function lcm_exchange_direction_for_phase

   integer function lcm_exchange_direction_from_offset(phase, offset) result(direction)
      integer, intent(in) :: phase
      integer, intent(in) :: offset

      call lcm_require_valid_exchange_phase(phase)
      if (offset /= -1 .and. offset /= 1) then
         error stop 'LCM exchange offset must be -1 or 1'
      endif

      select case (phase)
      case (lcm_exchange_x_phase)
         if (offset < 0) then
            direction = lcm_exchange_west
         else
            direction = lcm_exchange_east
         endif
      case (lcm_exchange_y_phase)
         if (offset < 0) then
            direction = lcm_exchange_south
         else
            direction = lcm_exchange_north
         endif
      end select
   end function lcm_exchange_direction_from_offset

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
      case (lcm_exchange_south)
         opposite = lcm_exchange_north
      case (lcm_exchange_north)
         opposite = lcm_exchange_south
      case default
         error stop 'LCM exchange direction is invalid'
      end select
   end function lcm_exchange_opposite_direction

   subroutine lcm_require_valid_exchange_phase(phase)
      integer, intent(in) :: phase

      if (phase < 1 .or. phase > lcm_exchange_phase_count) then
         error stop 'LCM exchange phase is invalid'
      endif
   end subroutine lcm_require_valid_exchange_phase

end module lcm_exchange_directions
