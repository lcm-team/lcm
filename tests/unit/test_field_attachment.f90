program test_field_attachment
   use iso_fortran_env, only: real64
   use lcm_field_state, only: field_state, lcm_attach_field, lcm_attach_fields
   use lcm_grid_state, only: lcm_set_grid
   use lcm_types, only: lcm_grid_t, lcm_init_grid

   implicit none

   type(lcm_grid_t) :: grid
   real(real64), target :: u_face(3,2,2)
   real(real64), target :: v_face(2,3,2)
   real(real64), target :: w_face(2,2,3)
   real(real64), target :: temperature(2,2,2)
   real(real64), target :: qv(2,2,2)

   call configure_grid()

   u_face = 1.0_real64
   v_face = 2.0_real64
   w_face = 3.0_real64
   temperature = 280.0_real64
   qv = 0.008_real64

   call lcm_attach_field("u", u_face)
   call lcm_attach_field("v", v_face)
   call lcm_attach_field("w", w_face)
   call lcm_attach_field("temperature", temperature)
   call lcm_attach_fields(qv=qv)

   call assert_true(field_state%u_is_attached, "u attached")
   call assert_true(field_state%v_is_attached, "v attached")
   call assert_true(field_state%w_is_attached, "w attached")
   call assert_true(field_state%temperature_is_attached, "temperature attached")
   call assert_true(field_state%qv_is_attached, "qv attached")

   call assert_true(all(shape(field_state%u) == [3, 2, 2]), "u shape")
   call assert_true(all(shape(field_state%v) == [2, 3, 2]), "v shape")
   call assert_true(all(shape(field_state%w) == [2, 2, 3]), "w shape")
   call assert_true(all(shape(field_state%temperature) == [2, 2, 2]), &
                    "temperature shape")
   call assert_true(all(shape(field_state%qv) == [2, 2, 2]), "qv shape")

contains

   subroutine configure_grid()
      real(real64) :: dz(2)
      real(real64) :: z_center(2)
      real(real64) :: z_face(3)

      dz = [100.0_real64, 100.0_real64]
      z_center = [50.0_real64, 150.0_real64]
      z_face = [0.0_real64, 100.0_real64, 200.0_real64]

      call lcm_init_grid(grid, 2, 2, 2, 4, 4, &
                         2, 3, 2, 3, 10.0_real64, 20.0_real64, &
                         z_center, z_face, dz, 0, 0, 0, 1, 1, 1)
      call lcm_set_grid(grid)
   end subroutine configure_grid

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_true

end program test_field_attachment
