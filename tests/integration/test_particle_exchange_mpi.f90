program test_particle_exchange_mpi
   use iso_fortran_env, only: int64, real64
   use mpi_f08, only: MPI_CART_COORDS, MPI_CART_CREATE, MPI_CART_RANK, &
                      MPI_COMM_FREE, MPI_COMM_RANK, MPI_COMM_SIZE, &
                      MPI_COMM_WORLD, MPI_Comm, MPI_FINALIZE, MPI_INIT
   use lcm_exchange_directions, only: lcm_exchange_north, lcm_exchange_south, &
                                      lcm_exchange_x_phase, lcm_exchange_y_phase
   use lcm_grid_state, only: lcm_set_grid
   use lcm_parallel_state, only: lcm_init_parallel, lcm_parallel_t, lcm_set_parallel
   use lcm_particle_exchange, only: lcm_exchange_particle_phase
   use lcm_particle_exchange_buffers, only: lcm_send_particle_message_count
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
   type(lcm_particle_t) :: expected_particle
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
      particle = populated_particle(rank, sequence, coordinates)
      call lcm_add_particle(initial_i, initial_j, 1, particle)
   enddo

   call lcm_exchange_particle_phase(lcm_exchange_x_phase)
   call assert_equal(total_particle_count(), 0, &
                     'x-received corner particles remain staged for y')
   call assert_equal( &
      lcm_send_particle_message_count(lcm_exchange_south) + &
      lcm_send_particle_message_count(lcm_exchange_north), &
      particles_per_route, 'x-received corner particles populate y buffers')
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
      expected_particle = populated_particle(expected_source_rank, sequence, &
                                             expected_source_coordinates)
      call assert_particle_equal(particle, expected_particle, &
                                 'diagonal particle state')
   enddo

   call MPI_COMM_FREE(communicator, ierr)
   call MPI_FINALIZE(ierr)

contains

   function populated_particle(source_rank, source_sequence, source_coordinates) &
      result(result_particle)
      integer, intent(in) :: source_rank
      integer, intent(in) :: source_sequence
      integer, intent(in) :: source_coordinates(2)
      type(lcm_particle_t) :: result_particle
      integer :: integer_seed
      real(real64) :: seed

      integer_seed = 100 * source_rank + source_sequence
      seed = real(integer_seed, real64)
      result_particle%weight = int(2000 + 10 * source_rank + source_sequence, &
                                   int64)
      result_particle%id = int(1000 + 10 * source_rank + source_sequence, int64)
      result_particle%particle_type = integer_seed + 2
      result_particle%x = merge(25.0_real64, 15.0_real64, &
                                source_coordinates(1) == 0)
      result_particle%y = merge(25.0_real64, 15.0_real64, &
                                source_coordinates(2) == 0)
      result_particle%z = 5.0_real64
      result_particle%u = seed + 1.0_real64
      result_particle%v = seed + 2.0_real64
      result_particle%w = seed + 3.0_real64
      result_particle%u_sgs = seed + 4.0_real64
      result_particle%v_sgs = seed + 5.0_real64
      result_particle%w_sgs = seed + 6.0_real64
      result_particle%mass_water = seed + 7.25_real64
      result_particle%mass_aero = seed + 8.25_real64
      result_particle%w_sedi = seed + 9.0_real64
      result_particle%temperature = seed + 270.0_real64
      result_particle%qv = seed + 10.0_real64
      result_particle%eta = seed + 11.0_real64
      result_particle%axis_a = seed + 12.0_real64
      result_particle%axis_c = seed + 13.0_real64
      result_particle%density = seed + 14.0_real64
      result_particle%rimed_mass = seed + 15.0_real64
      result_particle%freezing_temperature = seed + 16.0_real64
      result_particle%in_surface_area = seed + &
         [17.0_real64, 18.0_real64, 19.0_real64, 20.0_real64]
      result_particle%aerosol_portions = seed + &
         [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
      result_particle%cloudy_age = seed + 21.0_real64
      result_particle%cell_particle_index = integer_seed + 22
      result_particle%last_i = integer_seed + 23
      result_particle%last_j = integer_seed + 24
      result_particle%last_k = integer_seed + 25
      result_particle%n_monomer = int(integer_seed + 2000, int64)
      result_particle%spray_marker = integer_seed + 26
      result_particle%track_marker = integer_seed + 27
   end function populated_particle

   subroutine assert_particle_equal(actual, expected, description)
      type(lcm_particle_t), intent(in) :: actual
      type(lcm_particle_t), intent(in) :: expected
      character(len=*), intent(in) :: description

      call assert_true(actual%weight == expected%weight, description // ' weight')
      call assert_true(actual%id == expected%id, description // ' id')
      call assert_true(actual%particle_type == expected%particle_type, &
                       description // ' type')
      call assert_true(all([actual%x, actual%y, actual%z, actual%u, actual%v, &
                            actual%w, actual%u_sgs, actual%v_sgs, actual%w_sgs] == &
                           [expected%x, expected%y, expected%z, expected%u, &
                            expected%v, expected%w, expected%u_sgs, &
                            expected%v_sgs, expected%w_sgs]), &
                       description // ' motion')
      call assert_true(all([actual%mass_water, actual%mass_aero, actual%w_sedi, &
                            actual%temperature, actual%qv, actual%eta] == &
                           [expected%mass_water, expected%mass_aero, &
                            expected%w_sedi, expected%temperature, expected%qv, &
                            expected%eta]), description // ' microphysics')
      call assert_true(all([actual%axis_a, actual%axis_c, actual%density, &
                            actual%rimed_mass, actual%freezing_temperature] == &
                           [expected%axis_a, expected%axis_c, expected%density, &
                            expected%rimed_mass, expected%freezing_temperature]), &
                       description // ' ice state')
      call assert_true(all(actual%in_surface_area == expected%in_surface_area), &
                       description // ' surface area')
      call assert_true(all(actual%aerosol_portions == expected%aerosol_portions), &
                       description // ' aerosol portions')
      call assert_true(actual%cloudy_age == expected%cloudy_age, &
                       description // ' cloudy age')
      call assert_true(all([actual%cell_particle_index, actual%last_i, &
                            actual%last_j, actual%last_k, actual%spray_marker, &
                            actual%track_marker] == &
                           [expected%cell_particle_index, expected%last_i, &
                            expected%last_j, expected%last_k, &
                            expected%spray_marker, expected%track_marker]), &
                       description // ' integer metadata')
      call assert_true(actual%n_monomer == expected%n_monomer, &
                       description // ' monomers')
   end subroutine assert_particle_equal

   subroutine assert_true(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description

      if (.not. condition) error stop 'assertion failed: ' // description
   end subroutine assert_true

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

end program test_particle_exchange_mpi
