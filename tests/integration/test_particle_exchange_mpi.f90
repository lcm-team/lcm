program test_particle_exchange_mpi
   use iso_fortran_env, only: int64, real64
   use mpi_f08, only: MPI_CART_COORDS, MPI_CART_CREATE, MPI_CART_RANK, &
                      MPI_COMM_FREE, MPI_COMM_RANK, MPI_COMM_SIZE, &
                      MPI_COMM_WORLD, MPI_Comm, MPI_FINALIZE, MPI_INIT
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
   integer :: initial_i
   integer :: initial_j
   integer :: expected_source_coordinates(2)
   integer :: expected_source_rank
   integer :: final_i
   integer :: final_j
   integer :: mpi_size
   integer, parameter :: particles_per_route = 3
   type(lcm_parallel_t) :: parallel
   type(lcm_particle_t) :: particle
   logical :: periodicity(2)
   integer :: rank
   integer :: sequence
   real(real64) :: z_center(1)
   real(real64) :: z_face(2)

   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
   if (mpi_size /= 4) error stop 'corner-route payload test requires four ranks'

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

   ! Every rank sends one particle toward the diagonally opposite rank:
   ! NE, NW, SE, and SW are therefore exercised in the same 2x2 run.
   if (coordinates(1) == 0) then
      initial_i = 2
      final_i = 2
   else
      initial_i = 1
      final_i = 1
   endif
   if (coordinates(2) == 0) then
      initial_j = 2
      final_j = 2
   else
      initial_j = 1
      final_j = 1
   endif

   do sequence = 1, particles_per_route
      particle = lcm_particle_t()
      particle%id = int(1000 + 10 * rank + sequence, int64)
      particle%weight = int(2000 + 10 * rank + sequence, int64)
      particle%mass_water = real(10 * rank + sequence, real64) + 0.25_real64
      particle%temperature = 270.0_real64 + real(10 * rank + sequence, real64)
      particle%x = merge(25.0_real64, 15.0_real64, coordinates(1) == 0)
      particle%y = merge(25.0_real64, 15.0_real64, coordinates(2) == 0)
      particle%z = 5.0_real64
      call lcm_add_particle(initial_i, initial_j, 1, particle)
   enddo

   call lcm_exchange_particle_phase(lcm_exchange_x_phase)
   call assert_equal(total_particle_count(), particles_per_route, &
                     'x phase preserves every route particle')
   call lcm_exchange_particle_phase(lcm_exchange_y_phase)

   expected_source_coordinates = 1 - coordinates
   call MPI_CART_RANK(communicator, expected_source_coordinates, &
                      expected_source_rank, ierr)

   call assert_equal(total_particle_count(), particles_per_route, &
                     'final local particle count')
   call assert_equal(cell_particles(final_i,final_j,1)%count, &
                     particles_per_route, &
                     'particles inserted in final diagonal destination cell')
   do sequence = 1, particles_per_route
      particle = cell_particles(final_i,final_j,1)%particles(sequence)
      call assert_equal(int(particle%id), &
                        1000 + 10 * expected_source_rank + sequence, &
                        'ordered diagonal source particle id')
      call assert_equal(int(particle%weight), &
                        2000 + 10 * expected_source_rank + sequence, &
                        'diagonal source particle weight')
      call assert_close(particle%mass_water, &
                        real(10 * expected_source_rank + sequence, real64) + &
                           0.25_real64, &
                        'diagonal source water mass')
      call assert_close(particle%temperature, &
                        270.0_real64 + &
                           real(10 * expected_source_rank + sequence, real64), &
                        'diagonal source temperature')
   enddo

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

   subroutine assert_equal(actual, expected, description)
      integer, intent(in) :: actual
      integer, intent(in) :: expected
      character(len=*), intent(in) :: description

      if (actual /= expected) error stop 'assertion failed: ' // description
   end subroutine assert_equal

   subroutine assert_close(actual, expected, description)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: description

      if (actual /= expected) error stop 'assertion failed: ' // description
   end subroutine assert_close

end program test_particle_exchange_mpi
