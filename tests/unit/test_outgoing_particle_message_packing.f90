program test_outgoing_particle_message_packing
   use iso_fortran_env, only: int64, real64
   use lcm_exchange_directions, only: lcm_exchange_northeast
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
   type(lcm_particle_mpi_message_t) :: message
   type(lcm_particle_t) :: particle
   real(real64) :: dz(1)
   real(real64) :: z_center(1)
   real(real64) :: z_face(2)

   dz = [10.0_real64]
   z_center = [5.0_real64]
   z_face = [0.0_real64, 10.0_real64]
   call lcm_init_grid(grid, 3, 2, 1, 20, 20, &
                      2, 4, 2, 3, 10.0_real64, 10.0_real64, &
                      z_center, z_face, dz, 1, 1, 0, 0, 0, &
                      2, 2, 11, 7, 0, 0, 0, 0, &
                      .false., .false., 1)
   call lcm_set_grid(grid)
   call lcm_allocate_particle_cells()

   particle = lcm_particle_t()
   particle%id = 501_int64
   particle%x = 130.0_real64
   particle%y = 80.0_real64
   particle%z = 5.0_real64
   call lcm_add_particle(1, 1, 1, particle)

   call lcm_extract_outgoing_particles()

   call assert_equal(cell_particles(1,1,1)%count, 0, 'outgoing particle removed')
   call assert_equal(lcm_send_particle_message_count(lcm_exchange_northeast), 1, &
                     'north-east send count')
   call lcm_get_send_particle_message(lcm_exchange_northeast, 1, message)
   call assert_equal(message%global_i, 14, 'packed global i')
   call assert_equal(message%global_j, 9, 'packed global j')
   call assert_equal(message%k, 1, 'packed k')
   call assert_equal(int(message%particle%id), 501, 'packed particle id')

contains

   subroutine assert_equal(actual, expected, description)
      integer, intent(in) :: actual
      integer, intent(in) :: expected
      character(len=*), intent(in) :: description

      if (actual /= expected) error stop 'assertion failed: ' // description
   end subroutine assert_equal

end program test_outgoing_particle_message_packing
