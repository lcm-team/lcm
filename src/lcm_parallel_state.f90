!> LCM-owned communication context supplied by the host model.
module lcm_parallel_state
   use mpi_f08, only: MPI_CART_GET, MPI_CART_SHIFT, MPI_Comm, &
                      MPI_COMM_NULL, MPI_COMM_RANK, MPI_COMM_SIZE, &
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
   !> All members after communicator are derived by lcm_set_parallel and must
   !> not be populated independently by a host adapter.
   type :: lcm_parallel_t
      type(MPI_Comm) :: communicator = MPI_COMM_NULL
      integer :: rank = 0
      integer :: size = 1
      integer :: coordinates(2) = [0, 0]
      integer :: dimensions(2) = [1, 1]
      logical :: periodicity(2) = [.false., .false.]
      integer :: west_rank = MPI_PROC_NULL
      integer :: east_rank = MPI_PROC_NULL
      integer :: south_rank = MPI_PROC_NULL
      integer :: north_rank = MPI_PROC_NULL
   end type lcm_parallel_t

   type :: t_lcm_parallel_state
      type(lcm_parallel_t) :: parallel
      logical :: is_set = .false.
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
      call initialize_cartesian_topology(parallel_state%parallel)
      parallel_state%is_set = .true.
   end subroutine lcm_set_parallel

   subroutine initialize_cartesian_topology(parallel)
      type(lcm_parallel_t), intent(inout) :: parallel
      integer :: ierr

      call MPI_COMM_RANK(parallel%communicator, parallel%rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the host communicator rank'
      endif
      call MPI_COMM_SIZE(parallel%communicator, parallel%size, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the host communicator size'
      endif

      call MPI_CART_GET(parallel%communicator, 2, parallel%dimensions, &
                        parallel%periodicity, parallel%coordinates, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the host Cartesian communicator'
      endif

      call MPI_CART_SHIFT(parallel%communicator, 0, 1, parallel%west_rank, &
                          parallel%east_rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to derive west/east Cartesian neighbors'
      endif
      call MPI_CART_SHIFT(parallel%communicator, 1, 1, parallel%south_rank, &
                          parallel%north_rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to derive south/north Cartesian neighbors'
      endif
   end subroutine initialize_cartesian_topology

end module lcm_parallel_state
