program test_particle_exchange_mpi_edge_cases
   use iso_fortran_env, only: int64, real64
   use mpi_f08, only: MPI_ALLREDUCE, MPI_CART_COORDS, MPI_CART_CREATE, &
                      MPI_CART_RANK, MPI_COMM_FREE, MPI_COMM_RANK, MPI_COMM_SIZE, &
                      MPI_COMM_WORLD, MPI_Comm, MPI_FINALIZE, MPI_INIT, MPI_INTEGER, &
                      MPI_SUM
   use lcm_exchange_directions, only: lcm_exchange_x_phase, lcm_exchange_y_phase
   use lcm_grid_state, only: lcm_set_grid
   use lcm_parallel_state, only: lcm_init_parallel, lcm_parallel_t, lcm_set_parallel
   use lcm_particle_exchange, only: lcm_exchange_particle_phase
   use lcm_particle_state, only: cell_particles, lcm_add_particle, &
                                 lcm_allocate_particle_cells
   use lcm_types, only: lcm_grid_t, lcm_init_grid, lcm_particle_t

   implicit none

   type(MPI_Comm) :: communicator
   integer :: coordinates(2)
   integer :: dimensions(2)
   real(real64) :: dz(1)
   type(lcm_grid_t) :: grid
   integer :: ierr
   integer :: mpi_size
   integer :: rank
   integer :: source_rank
   integer :: destination_rank
   integer :: source_coordinates(2)
   integer :: destination_coordinates(2)
   integer :: local_count
   integer :: global_count
   logical :: periodicity(2)
   type(lcm_parallel_t) :: parallel
   type(lcm_particle_t) :: particle
   real(real64) :: z_center(1)
   real(real64) :: z_face(2)

   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
   if (mpi_size /= 4) error stop 'edge-case exchange test requires four ranks'

   dimensions = [2, 2]
   periodicity = [.false., .false.]
   call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dimensions, periodicity, .false., &
                        communicator, ierr)
   call MPI_COMM_RANK(communicator, rank, ierr)
   call MPI_CART_COORDS(communicator, rank, 2, coordinates, ierr)

   call lcm_init_parallel(parallel, communicator)
   call lcm_set_parallel(parallel)

   dz = [10.0_real64]
   z_center = [5.0_real64]
   z_face = [0.0_real64, 10.0_real64]
   call lcm_init_grid(grid, 2, 2, 1, 4, 4, &
                      2, 3, 2, 3, 10.0_real64, 10.0_real64, &
                      z_center, z_face, dz, 1, 1, 0, &
                      2 * coordinates(1) + 1, 2 * coordinates(2) + 1, 1)
   call lcm_set_grid(grid)
   call lcm_allocate_particle_cells()

   ! A completely empty exchange must be a no-op on every rank.
   call lcm_exchange_particle_phase(lcm_exchange_x_phase)
   call lcm_exchange_particle_phase(lcm_exchange_y_phase)
   call assert_equal(total_particle_count(), 0, 'zero-particle exchange is a no-op')
   local_count = total_particle_count()
   call MPI_ALLREDUCE(local_count, global_count, 1, MPI_INTEGER, MPI_SUM, &
                      communicator, ierr)
   call assert_equal(global_count, 0, 'zero-particle global conservation')

   ! Only the lower-left rank sends one particle east. All other ranks have
   ! zero outgoing traffic, which exercises one-sided communication.
   source_coordinates = [0, 0]
   destination_coordinates = [1, 0]
   call MPI_CART_RANK(communicator, source_coordinates, source_rank, ierr)
   call MPI_CART_RANK(communicator, destination_coordinates, destination_rank, ierr)
   if (rank == source_rank) then
      particle = lcm_particle_t()
      particle%id = 9001_int64
      particle%weight = 17_int64
      particle%x = 25.0_real64
      particle%y = 5.0_real64
      particle%z = 5.0_real64
      call lcm_add_particle(2, 1, 1, particle)
   endif

   call lcm_exchange_particle_phase(lcm_exchange_x_phase)
   call assert_equal(total_particle_count(), merge(1, 0, rank == destination_rank), &
                     'one-sided x migration ownership')
   call assert_equal(cell_particles(1, 1, 1)%count, &
                     merge(1, 0, rank == destination_rank), &
                     'one-sided x migration destination cell')
   if (rank == destination_rank) then
      call assert_true(cell_particles(1, 1, 1)%particles(1)%id == 9001_int64, &
                       'one-sided migration preserves particle id')
      call assert_true(cell_particles(1, 1, 1)%particles(1)%weight == 17_int64, &
                       'one-sided migration preserves particle weight')
   endif

   call lcm_exchange_particle_phase(lcm_exchange_y_phase)
   call assert_equal(total_particle_count(), merge(1, 0, rank == destination_rank), &
                     'one-sided y phase preserves ownership')
   local_count = total_particle_count()
   call MPI_ALLREDUCE(local_count, global_count, 1, MPI_INTEGER, MPI_SUM, &
                      communicator, ierr)
   call assert_equal(global_count, 1, 'one-sided global conservation')

   call MPI_COMM_FREE(communicator, ierr)
   call MPI_FINALIZE(ierr)

contains

   integer function total_particle_count() result(count)
      integer :: i
      integer :: j
      integer :: k

      count = 0
      do k = 1, size(cell_particles, 3)
         do j = 1, size(cell_particles, 2)
            do i = 1, size(cell_particles, 1)
               count = count + cell_particles(i,j,k)%count
            enddo
         enddo
      enddo
   end function total_particle_count

   subroutine assert_true(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description

      if (.not. condition) error stop 'assertion failed: ' // description
   end subroutine assert_true

   subroutine assert_equal(actual, expected, description)
      integer, intent(in) :: actual
      integer, intent(in) :: expected
      character(len=*), intent(in) :: description

      if (actual /= expected) error stop 'assertion failed: ' // description
   end subroutine assert_equal

end program test_particle_exchange_mpi_edge_cases
