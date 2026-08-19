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
   call assert_equal(int(cell_particles(2,2,1)%particles(1)%id), 503, &
                     'local two-dimensional move id')
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

      particle = lcm_particle_t()
      particle%id = id
      particle%x = x
      particle%y = y
      particle%z = 5.0_real64
      call lcm_add_particle(1, 1, 1, particle)
   end subroutine add_particle

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
