program test_outgoing_particle_message_packing
   use iso_fortran_env, only: int64, real64
   use mpi_f08, only: MPI_CART_CREATE, MPI_COMM_FREE, MPI_COMM_WORLD, MPI_Comm, &
                      MPI_FINALIZE, MPI_INIT
   use lcm_exchange_directions, only: lcm_exchange_east, lcm_exchange_north, &
                                      lcm_exchange_south, lcm_exchange_west, &
                                      lcm_exchange_x_phase
   use lcm_parallel_state, only: lcm_init_parallel, lcm_parallel_t, &
                                 lcm_set_parallel
   use lcm_particle_exchange_buffers, only: lcm_get_send_particle_message, &
                                            lcm_send_particle_message_count
   use lcm_particle_state, only: cell_particles, lcm_add_particle, &
                                 lcm_allocate_particle_cells, &
                                 lcm_extract_outgoing_particles
   use lcm_types, only: lcm_grid_t, lcm_init_grid, lcm_particle_mpi_message_t, &
                        lcm_particle_t
   use lcm_grid_state, only: lcm_set_grid

   implicit none

   type(lcm_grid_t) :: grid
   type(lcm_parallel_t) :: parallel
   type(MPI_Comm) :: communicator
   integer :: dimensions(2)
   integer :: ierr
   logical :: periodicity(2)
   type(lcm_particle_mpi_message_t) :: message
   type(lcm_particle_t) :: expected_local_particle
   type(lcm_particle_t) :: particle
   real(real64) :: dz(1)
   real(real64) :: z_center(1)
   real(real64) :: z_face(2)

   dz = [10.0_real64]
   z_center = [5.0_real64]
   z_face = [0.0_real64, 10.0_real64]
   call MPI_INIT(ierr)
   dimensions = [1, 1]
   periodicity = [.true., .true.]
   call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dimensions, periodicity, .false., &
                        communicator, ierr)
   call lcm_init_parallel(parallel, communicator)
   call lcm_set_parallel(parallel)
   call lcm_init_grid(grid, 3, 2, 1, 20, 20, &
                      2, 4, 2, 3, 10.0_real64, 10.0_real64, &
                      z_center, z_face, dz, 1, 1, 0, 11, 7, 1)
   call lcm_set_grid(grid)
   call lcm_allocate_particle_cells()

   expected_local_particle = populated_particle(503_int64, 115.0_real64, &
                                                75.0_real64)
   call add_particle(501_int64, 105.0_real64, 65.0_real64)
   call add_particle(502_int64, 130.0_real64, 65.0_real64)
   call add_particle(503_int64, 115.0_real64, 75.0_real64)
   call add_particle(504_int64, 130.0_real64, 80.0_real64)
   call add_particle(505_int64, 125.0_real64, 65.0_real64)
   call add_particle(506_int64, 95.0_real64, 65.0_real64)
   call add_particle(507_int64, 95.0_real64, 80.0_real64)
   call add_particle(508_int64, 105.0_real64, 90.0_real64)
   call add_particle(509_int64, 105.0_real64, 55.0_real64)

   call lcm_extract_outgoing_particles(lcm_exchange_x_phase)

   call assert_equal(cell_particles(1,1,1)%count, 1, 'same-cell particle count')
   call assert_equal(int(cell_particles(1,1,1)%particles(1)%id), 501, &
                     'same-cell particle order')
   call assert_equal(cell_particles(2,2,1)%count, 1, &
                     'local two-dimensional move count')
   call assert_particle_equal(cell_particles(2,2,1)%particles(1), &
                              expected_local_particle, &
                              'local two-dimensional move')
   call assert_equal(cell_particles(3,1,1)%count, 1, 'local x move count')
   call assert_equal(int(cell_particles(3,1,1)%particles(1)%id), 505, &
                     'local x move id')
   call assert_equal(int(cell_particles(1,1,1)%particles(2)%id), 0, &
                     'first vacated slot cleared')
   call assert_equal(int(cell_particles(1,1,1)%particles(9)%id), 0, &
                     'last vacated slot cleared')
   call assert_equal(lcm_send_particle_message_count(lcm_exchange_east), 2, &
                     'east send count')
   call assert_equal(lcm_send_particle_message_count(lcm_exchange_west), 2, &
                     'west send count')
   call assert_equal(lcm_send_particle_message_count(lcm_exchange_north), 1, &
                     'staged north send count')
   call assert_equal(lcm_send_particle_message_count(lcm_exchange_south), 1, &
                     'staged south send count')
   call lcm_get_send_particle_message(lcm_exchange_east, 1, message)
   call assert_equal(int(message%particle%id), 502, 'east particle id')
   call lcm_get_send_particle_message(lcm_exchange_east, 2, message)
   call assert_equal(message%global_i, 14, 'packed global i')
   call assert_equal(message%global_j, 9, 'packed global j')
   call assert_equal(message%k, 1, 'packed k')
   call assert_equal(int(message%particle%id), 504, 'packed particle id')
   call assert_close(message%particle%x, 130.0_real64, 'internal crossing x')
   call assert_close(message%particle%y, 80.0_real64, 'internal crossing y')
   call lcm_get_send_particle_message(lcm_exchange_west, 1, message)
   call assert_equal(int(message%particle%id), 506, 'first west particle id')
   call lcm_get_send_particle_message(lcm_exchange_west, 2, message)
   call assert_equal(int(message%particle%id), 507, 'second west particle id')
   call lcm_get_send_particle_message(lcm_exchange_north, 1, message)
   call assert_equal(int(message%particle%id), 508, 'north particle id')
   call lcm_get_send_particle_message(lcm_exchange_south, 1, message)
   call assert_equal(int(message%particle%id), 509, 'south particle id')

   call MPI_COMM_FREE(communicator, ierr)
   call MPI_FINALIZE(ierr)

contains

   subroutine add_particle(id, x, y)
      integer(int64), intent(in) :: id
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y

      particle = populated_particle(id, x, y)
      call lcm_add_particle(1, 1, 1, particle)
   end subroutine add_particle

   function populated_particle(id, x, y) result(result_particle)
      integer(int64), intent(in) :: id
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      type(lcm_particle_t) :: result_particle
      real(real64) :: seed

      seed = real(id - 500_int64, real64)
      result_particle%weight = id + 1000_int64
      result_particle%id = id
      result_particle%particle_type = int(seed) + 2
      result_particle%x = x
      result_particle%y = y
      result_particle%z = 5.0_real64
      result_particle%u = seed + 1.0_real64
      result_particle%v = seed + 2.0_real64
      result_particle%w = seed + 3.0_real64
      result_particle%u_sgs = seed + 4.0_real64
      result_particle%v_sgs = seed + 5.0_real64
      result_particle%w_sgs = seed + 6.0_real64
      result_particle%mass_water = seed + 7.0_real64
      result_particle%mass_aero = seed + 8.0_real64
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
      result_particle%cell_particle_index = int(seed) + 22
      result_particle%last_i = int(seed) + 23
      result_particle%last_j = int(seed) + 24
      result_particle%last_k = int(seed) + 25
      result_particle%n_monomer = id + 2000_int64
      result_particle%spray_marker = int(seed) + 26
      result_particle%track_marker = int(seed) + 27
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

      if (abs(actual - expected) > 1.0e-12_real64) then
         error stop 'assertion failed: ' // description
      endif
   end subroutine assert_close

end program test_outgoing_particle_message_packing
