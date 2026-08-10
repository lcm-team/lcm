program test_local_subdomain_bounds
   use iso_fortran_env, only: real64
   use lcm_grid_state, only: lcm_local_subdomain_bounds, lcm_set_grid
   use lcm_particle_location, only: lcm_locate_particle_cell
   use lcm_types, only: lcm_grid_t, lcm_init_grid

   implicit none

   type(lcm_grid_t) :: grid
   real(real64) :: dz(1)
   real(real64) :: z_center(1)
   real(real64) :: z_face(2)
   integer :: i
   integer :: j
   integer :: k
   integer :: global_i_first
   integer :: global_i_last
   integer :: global_j_first
   integer :: global_j_last
   logical :: is_local

   dz = [10.0_real64]
   z_center = [5.0_real64]
   z_face = [0.0_real64, 10.0_real64]

   ! Host-array interior indices deliberately differ from global LCM indices.
   call lcm_init_grid(grid, 3, 2, 1, 20, 20, &
                      2, 4, 2, 3, 10.0_real64, 10.0_real64, &
                      z_center, z_face, dz, 1, 1, 0, 0, 0, &
                      2, 2, 11, 7, -1, -1, -1, -1, &
                      .false., .false., 1)
   call lcm_set_grid(grid)

   call lcm_local_subdomain_bounds(global_i_first, global_i_last, &
                                   global_j_first, global_j_last)
   call assert_equal(global_i_first, 11, 'global i first')
   call assert_equal(global_i_last, 13, 'global i last')
   call assert_equal(global_j_first, 7, 'global j first')
   call assert_equal(global_j_last, 8, 'global j last')

   call lcm_locate_particle_cell(100.1_real64, 60.1_real64, 5.0_real64, &
                                 i, j, k, is_local)
   call assert_equal(i, 1, 'first local i')
   call assert_equal(j, 1, 'first local j')
   call assert_equal(k, 1, 'local k')
   call assert_true(is_local, 'first global cell is local')

   call lcm_locate_particle_cell(130.0_real64, 80.1_real64, 5.0_real64, &
                                 i, j, k, is_local)
   call assert_equal(i, 4, 'east outgoing i')
   call assert_equal(j, 3, 'north outgoing j')
   call assert_true(.not. is_local, 'north-east cell is outside local subdomain')

contains

   subroutine assert_equal(actual, expected, message)
      integer, intent(in) :: actual
      integer, intent(in) :: expected
      character(len=*), intent(in) :: message

      if (actual /= expected) error stop 'assertion failed: ' // message
   end subroutine assert_equal

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) error stop 'assertion failed: ' // message
   end subroutine assert_true

end program test_local_subdomain_bounds
