!> LCM-owned communication context derived from supplied communicator from the host model.
module lcm_parallel_state
   use mpi_f08, only: MPI_Comm, MPI_Comm_rank, MPI_Comm_size, MPI_Cart_get, MPI_Cart_shift,  &
                      MPI_Comm_split_type, MPI_Comm_split, MPI_Barrier,                      &
                      MPI_SUCCESS, MPI_COMM_NULL, MPI_PROC_NULL, MPI_INFO_NULL,              &
                      MPI_COMM_TYPE_SHARED

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
      type(lcm_parallel_t) :: global      !> based on global 2d cartesian topology
      type(lcm_parallel_t) :: intranode   !> based on intranode 2d cartesian topology
      type(lcm_parallel_t) :: internode   !> based on no-topology pool of all node heads
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

      parallel_state%global = host_parallel
      call initialize_global_topology(parallel_state%global)
      call initialize_intranode_topology(parallel_state%global, parallel_state%intranode)
      call initialize_internode_communicator(parallel_state%global, parallel_state%intranode, &
                                             parallel_state%internode)
      !> synchronize all processes here once
      call MPI_Barrier(parallel_state%global%communicator)

      parallel_state%is_set = .true.
   end subroutine lcm_set_parallel

   !> fills the parallel state info of the global topology
   subroutine initialize_global_topology(parallel)
      type(lcm_parallel_t), intent(inout) :: parallel
      integer :: ierr

      call MPI_Comm_rank(parallel%communicator, parallel%rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the host communicator rank'
      endif
      call MPI_Comm_size(parallel%communicator, parallel%size, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the host communicator size'
      endif

      call MPI_Cart_get(parallel%communicator, 2, parallel%dimensions, &
                        parallel%periodicity, parallel%coordinates, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the host communicator information'
      endif

      call MPI_Cart_shift(parallel%communicator, 0, 1, parallel%west_rank, &
                          parallel%east_rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to derive west/east Cartesian neighbors'
      endif
      call MPI_Cart_shift(parallel%communicator, 1, 1, parallel%south_rank, &
                          parallel%north_rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to derive south/north Cartesian neighbors'
      endif
   end subroutine initialize_global_topology

   !> splits the global topology nodewise into subcommunicators,
   !> each containing all processes running on a node (=intranode).
   !> the returned parallel state is again a 2d cartesian topology,
   !> describing the layout of the processes inside a node.
   subroutine initialize_intranode_topology(global_topo, intranode_topo)
      type(lcm_parallel_t), intent(inout) :: global_topo
      type(lcm_parallel_t), intent(inout) :: intranode_topo
      type(lcm_parallel_t)                :: intranode_notopo
      integer :: ierr, nprocy
      integer :: dims(2)
      logical :: periods(2)

      call MPI_Comm_split_type(global_topo%communicator, MPI_COMM_TYPE_SHARED, global_topo%rank, &
                               MPI_INFO_NULL, intranode_notopo%communicator, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to split the global host communicator' // &
                    ' into intranode no-topology communicators'
      endif

      call MPI_Comm_rank(intranode_notopo%communicator, intranode_notopo%rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the intranode no-topology communicator rank'
      endif
      call MPI_Comm_size(intranode_notopo%communicator, intranode_notopo%size, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the intranode no-topology communicator size'
      endif

      nprocy = global_topo%dimensions(2)
      if (mod(intranode_notopo%size, nprocy) /= 0) then
         if (intranode_notopo%rank == 0) then
            error stop 'LCM checked the runtime partition of the MPI processes onto nodes' // &
                       ' and found that processes-per-node is not a multiple of nprocy'    // &
                       ' leading to ragged lines that cannot be handled in the output'
         endif
      endif

      dims = [intranode_notopo%size / nprocy, nprocy]

      !> take over periodicity info provided by host model communicator.
      !> check if running on more than one node and if yes, set to
      !> non-periodic boundaries in x-direction (=between nodewise stripes of subdomains).
      periods = global_topo%periodicity
      if (global_topo%size /= intranode_notopo%size) then
         periods(1) = .false.
      endif

      call MPI_Cart_create(intranode_notopo%communicator, 2, dims, periods, .false., &
                           intranode_topo%communicator, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to create intranode cartesian topology'
      endif

      call MPI_Comm_rank(intranode_topo%communicator, intranode_topo%rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the intranode topo communicator rank'
      endif
      call MPI_Comm_size(intranode_topo%communicator, intranode_topo%size, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the intranode topo communicator size'
      endif

      call MPI_Cart_get(intranode_topo%communicator, 2, intranode_topo%dimensions, &
                        intranode_topo%periodicity, intranode_topo%coordinates, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the intranode topo communicator information'
      endif
   end subroutine initialize_intranode_topology

   !> creates a new communicator as a no-topology pool of all node head processes (=internode)
   !> needed for (nodewise) parallel file I/O, internode rank 0 = global rank 0
   !> the rank in this communicator tells the position of each node head process
   subroutine initialize_internode_communicator(global_topo, intranode_topo, internode_notopo)
      type(lcm_parallel_t), intent(in)    :: global_topo
      type(lcm_parallel_t), intent(in)    :: intranode_topo
      type(lcm_parallel_t), intent(inout) :: internode_notopo
      integer :: ierr

      call MPI_Comm_split(global_topo%communicator, intranode_topo%rank, global_topo%rank,   &
                          internode_notopo%communicator, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to split the global host communicator' // &
                    ' into an internode no-topology communicator'
      endif

      call MPI_Comm_rank(internode_notopo%communicator, internode_notopo%rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the internode topo communicator rank'
      endif
      call MPI_Comm_size(internode_notopo%communicator, internode_notopo%size, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to query the internode topo communicator size'
      endif
   end subroutine initialize_internode_communicator

end module lcm_parallel_state
