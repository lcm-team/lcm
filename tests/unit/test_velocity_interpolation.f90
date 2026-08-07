program test_velocity_interpolation
   use iso_fortran_env, only: real64
   use lcm_field_state, only: lcm_attach_field
   use lcm_grid_state, only: lcm_set_grid
   use lcm_types, only: lcm_grid_t, lcm_init_grid
   use lcm_velocity_interpolation, only: lcm_interpolate_particle_velocity

   implicit none

   type(lcm_grid_t) :: grid
   real(real64), target :: u(3,2,2)
   real(real64), target :: v(2,3,2)
   real(real64), target :: w(2,2,3)

   call configure_grid()
   call set_uniform_velocity_fields()
   call lcm_attach_field("u", u)
   call lcm_attach_field("v", v)
   call lcm_attach_field("w", w)

   call assert_velocity(2.5_real64, 5.0_real64, 50.0_real64, &
                        10.0_real64, 20.0_real64, 30.0_real64, &
                        "uniform first cell")
   call set_linear_velocity_fields()
   call assert_velocity(15.0_real64, 25.0_real64, 150.0_real64, &
                        115.0_real64, 262.5_real64, 450.0_real64, &
                        "linear interpolation")
   call assert_velocity(15.0_real64, 25.0_real64, 175.0_real64, &
                        115.0_real64, 262.5_real64, 475.0_real64, &
                        "top vertical face")

contains

   subroutine configure_grid()
      real(real64) :: dz(2)
      real(real64) :: z_center(2)
      real(real64) :: z_face(3)

      dz = [100.0_real64, 100.0_real64]
      z_center = [50.0_real64, 150.0_real64]
      z_face = [0.0_real64, 100.0_real64, 200.0_real64]

      call lcm_init_grid(grid, 2, 2, 2, 2, 2, &
                         1, 2, 1, 2, 10.0_real64, 20.0_real64, &
                         z_center, z_face, dz, 0, 0, 0, 0, 0, &
                         1, 1, 1, 1, -1, -1, -1, -1, &
                         .false., .false., 1)
      call lcm_set_grid(grid)
   end subroutine configure_grid

   subroutine set_uniform_velocity_fields()
      u = 10.0_real64
      v = 20.0_real64
      w = 30.0_real64
   end subroutine set_uniform_velocity_fields

   subroutine set_linear_velocity_fields()
      integer :: i
      integer :: j
      integer :: k

      do k = 1, 2
         do j = 1, 2
            do i = 1, 3
               u(i, j, k) = 100.0_real64 + real(i - 1, real64) * 10.0_real64
            enddo
         enddo
      enddo
      do k = 1, 2
         do j = 1, 3
            do i = 1, 2
               v(i, j, k) = 200.0_real64 + real(j - 1, real64) * 50.0_real64
            enddo
         enddo
      enddo
      do k = 1, 3
         do j = 1, 2
            do i = 1, 2
               w(i, j, k) = 300.0_real64 + real(k - 1, real64) * 100.0_real64
            enddo
         enddo
      enddo
   end subroutine set_linear_velocity_fields

   subroutine assert_velocity(x, y, z, expected_u, expected_v, expected_w, message)
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      real(real64), intent(in) :: z
      real(real64), intent(in) :: expected_u
      real(real64), intent(in) :: expected_v
      real(real64), intent(in) :: expected_w
      character(len=*), intent(in) :: message

      real(real64) :: actual_u
      real(real64) :: actual_v
      real(real64) :: actual_w
      logical :: is_local

      call lcm_interpolate_particle_velocity(x, y, z, actual_u, actual_v, &
                                             actual_w, is_local)

      call assert_true(is_local, message // " local")
      call assert_close(actual_u, expected_u, message // " u")
      call assert_close(actual_v, expected_v, message // " v")
      call assert_close(actual_w, expected_w, message // " w")
   end subroutine assert_velocity

   subroutine assert_close(actual, expected, message)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: message

      if (abs(actual - expected) > 1.0e-12_real64) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_close

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_true

end program test_velocity_interpolation
