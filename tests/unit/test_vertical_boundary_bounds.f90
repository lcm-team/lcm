!> Verify that LCM physical vertical limits come from host-grid face heights.
program test_vertical_boundary_bounds
   use iso_fortran_env, only: real64
   use lcm_grid_state, only: grid_state, lcm_set_grid, &
                             lcm_vertical_boundary_inset
   use lcm_types, only: lcm_grid_t, lcm_init_grid

   implicit none

   type(lcm_grid_t) :: grid
   real(real64) :: dz(3)
   real(real64) :: z_center(3)
   real(real64) :: z_face(4)

   dz = [15.0_real64, 40.0_real64, 55.0_real64]
   z_face = [10.0_real64, 25.0_real64, 65.0_real64, 120.0_real64]
   z_center = [17.5_real64, 45.0_real64, 92.5_real64]

   call lcm_init_grid(grid, 1, 1, 3, 1, 1, &
                      1, 1, 1, 1, 10.0_real64, 10.0_real64, &
                      z_center, z_face, dz, 0, 0, 0, 0, 0, &
                      1, 1, 1, 1, -1, -1, -1, -1, &
                      .false., .false., 1)
   call lcm_set_grid(grid)

   call assert_close(grid_state%zbottom, z_face(1) + &
                     lcm_vertical_boundary_inset, "bottom from first z face")
   call assert_close(grid_state%ztop, z_face(4) - &
                     lcm_vertical_boundary_inset, "top from final z face")

contains

   subroutine assert_close(actual, expected, message)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: message

      if (abs(actual - expected) > 1.0e-12_real64) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_close

end program test_vertical_boundary_bounds
