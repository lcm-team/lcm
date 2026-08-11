!> LCM-owned communication context supplied by the host model.
module lcm_parallel_state
   use mpi_f08, only: MPI_CART_GET, MPI_CART_RANK, MPI_Comm, MPI_COMM_NULL, &
                      MPI_PROC_NULL, MPI_SUCCESS

   implicit none

   private

   public :: t_lcm_parallel_state
   public :: lcm_parallel_t
   public :: lcm_init_parallel
   public :: parallel_state
   public :: lcm_set_parallel

   !> Host-neutral MPI execution context supplied by an LES adapter.
   !>
   !> Grid geometry remains in lcm_grid_t. The host supplies its existing
   !> Cartesian communicator; LCM does not create a separate decomposition.
   type :: lcm_parallel_t
      type(MPI_Comm) :: communicator = MPI_COMM_NULL
   end type lcm_parallel_t

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

   subroutine lcm_init_parallel(parallel, communicator)
      type(lcm_parallel_t), intent(out) :: parallel
      type(MPI_Comm), intent(in) :: communicator

      parallel%communicator = communicator
   end subroutine lcm_init_parallel

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
      integer :: coordinates(2)
      integer :: dimensions(2)
      integer :: ierr
      logical :: periodicity(2)

      call MPI_CART_GET(parallel_state%parallel%communicator, 2, dimensions, &
                        periodicity, coordinates, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the host Cartesian communicator'
      endif

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
      integer, intent(in) :: x_offset
      integer, intent(in) :: y_offset
      integer, intent(in) :: dimensions(2)
      logical, intent(in) :: periodicity(2)
      integer, intent(in) :: coordinates(2)
      integer :: ierr
      integer :: neighbor_coordinates(2)
      logical :: is_valid

      neighbor_coordinates = coordinates + [x_offset, y_offset]
      call wrap_or_reject_coordinate(neighbor_coordinates(1), dimensions(1), &
                                     periodicity(1), is_valid)
      if (.not. is_valid) then
         rank = MPI_PROC_NULL
         return
      endif
      call wrap_or_reject_coordinate(neighbor_coordinates(2), dimensions(2), &
                                     periodicity(2), is_valid)
      if (.not. is_valid) then
         rank = MPI_PROC_NULL
         return
      endif

      call MPI_CART_RANK(parallel_state%parallel%communicator, &
                         neighbor_coordinates, rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to derive a Cartesian neighbor rank'
      endif
   end function cartesian_neighbor_rank

   subroutine wrap_or_reject_coordinate(coordinate, dimension, periodic, is_valid)
      integer, intent(inout) :: coordinate
      integer, intent(in) :: dimension
      logical, intent(in) :: periodic
      logical, intent(out) :: is_valid

      is_valid = coordinate >= 0 .and. coordinate < dimension
      if (is_valid) return
      if (.not. periodic) then
         return
      endif

      coordinate = modulo(coordinate, dimension)
      is_valid = .true.
   end subroutine wrap_or_reject_coordinate

end module lcm_parallel_state
