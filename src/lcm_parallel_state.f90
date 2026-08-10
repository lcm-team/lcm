!> LCM-owned communication context supplied by the host model.
module lcm_parallel_state
   use lcm_types, only: lcm_parallel_t

   implicit none

   private

   public :: t_lcm_parallel_state
   public :: parallel_state
   public :: lcm_set_parallel

   type :: t_lcm_parallel_state
      type(lcm_parallel_t) :: parallel
      logical :: is_set = .false.
      integer :: southwest_rank = -1
      integer :: southeast_rank = -1
      integer :: northwest_rank = -1
      integer :: northeast_rank = -1
   end type t_lcm_parallel_state

   type(t_lcm_parallel_state), protected, save :: parallel_state

contains

   subroutine lcm_set_parallel(host_parallel)
      type(lcm_parallel_t), intent(in) :: host_parallel

      if (parallel_state%is_set) then
         error stop 'LCM parallel state has already been initialized'
      endif

      parallel_state%parallel = host_parallel
      call derive_diagonal_neighbor_ranks()
      parallel_state%is_set = .true.
   end subroutine lcm_set_parallel

   subroutine derive_diagonal_neighbor_ranks()
      include 'mpif.h'

      integer :: coordinates(2)
      integer :: dimensions(2)
      integer :: ierr
      logical :: periodicity(2)

      call MPI_CART_GET(parallel_state%parallel%communicator, 2, dimensions, &
                        periodicity, coordinates, ierr)

      parallel_state%southwest_rank = cartesian_neighbor_rank(-1, -1, &
         dimensions, periodicity, coordinates)
      parallel_state%southeast_rank = cartesian_neighbor_rank(1, -1, &
         dimensions, periodicity, coordinates)
      parallel_state%northwest_rank = cartesian_neighbor_rank(-1, 1, &
         dimensions, periodicity, coordinates)
      parallel_state%northeast_rank = cartesian_neighbor_rank(1, 1, &
         dimensions, periodicity, coordinates)
   end subroutine derive_diagonal_neighbor_ranks

   integer function cartesian_neighbor_rank(x_offset, y_offset, dimensions, &
                                            periodicity, coordinates) result(rank)
      include 'mpif.h'

      integer, intent(in) :: x_offset
      integer, intent(in) :: y_offset
      integer, intent(in) :: dimensions(2)
      logical, intent(in) :: periodicity(2)
      integer, intent(in) :: coordinates(2)
      integer :: ierr
      integer :: neighbor_coordinates(2)

      neighbor_coordinates = coordinates + [x_offset, y_offset]
      call wrap_or_reject_coordinate(neighbor_coordinates(1), dimensions(1), &
                                     periodicity(1), rank)
      if (rank < 0) return
      call wrap_or_reject_coordinate(neighbor_coordinates(2), dimensions(2), &
                                     periodicity(2), rank)
      if (rank < 0) return

      call MPI_CART_RANK(parallel_state%parallel%communicator, &
                         neighbor_coordinates, rank, ierr)
   end function cartesian_neighbor_rank

   subroutine wrap_or_reject_coordinate(coordinate, dimension, periodic, status)
      integer, intent(inout) :: coordinate
      integer, intent(in) :: dimension
      logical, intent(in) :: periodic
      integer, intent(out) :: status

      status = 0
      if (coordinate >= 0 .and. coordinate < dimension) return
      if (.not. periodic) then
         status = -1
         return
      endif

      coordinate = modulo(coordinate, dimension)
   end subroutine wrap_or_reject_coordinate

end module lcm_parallel_state
