program test_field_copy
   use iso_fortran_env, only: real64
   use lcm_field_state, only: field_state, lcm_set_field
   use lcm_grid_state, only: lcm_set_grid
   use lcm_types, only: lcm_grid_t, lcm_init_grid

   implicit none

   type(lcm_grid_t) :: grid
   real(real64) :: u_face(3,2,2)
   real(real64) :: v_face(2,3,2)
   real(real64) :: w_face(2,2,3)
   real(real64) :: temperature(2,2,2)
   real(real64) :: qv(2,2,2)
   real(real64) :: pressure(2), density(2), dry_density(2,2,2)
   real(real64) :: static_energy(2,2,2), liquid(2,2,2)
   real(real64) :: strided_temperature(4,2,2)
   real(real64), allocatable :: temporary_temperature(:,:,:)

   ! Metadata exists before storage is supplied; allocation defines readiness.
   call assert_true(.not. allocated(field_state%temperature%values), "temperature initially absent")
   call assert_true(field_state%temperature%var_name == "temperature", "temperature name")
   call assert_true(field_state%temperature%units == "K", "temperature units")
   call assert_true(field_state%temperature%long_name == "Air temperature", "temperature description")
   call assert_true(all(field_state%temperature%dimensions == [character(len=8) :: "x","y","z"]), &
                    "volume physical dimensions")
   call assert_true(all(field_state%pressure%dimensions == [character(len=8) :: "z","",""]), &
                    "profile physical dimension")
   call assert_true(all(field_state%u%dimensions == [character(len=8) :: "x_face","y","z"]), &
                    "u staggering")
   call assert_true(all(field_state%v%dimensions == [character(len=8) :: "x","y_face","z"]), &
                    "v staggering")
   call assert_true(all(field_state%w%dimensions == [character(len=8) :: "x","y","z_face"]), &
                    "w staggering")
   call configure_grid()

   u_face = 1.0_real64
   v_face = 2.0_real64
   w_face = 3.0_real64
   temperature = 280.0_real64
   qv = 0.008_real64

   ! Each field keeps the shape required by the particle velocity interpolator.
   call lcm_set_field("u", u_face)
   call lcm_set_field("v", v_face)
   call lcm_set_field("w", w_face)
   call lcm_set_field("qv", qv)
   call lcm_set_field("temperature", temperature)

   call assert_true(allocated(field_state%u%values), "u supplied")
   call assert_true(allocated(field_state%v%values), "v supplied")
   call assert_true(allocated(field_state%w%values), "w supplied")
   call assert_true(allocated(field_state%temperature%values), "temperature supplied")
   call assert_true(allocated(field_state%qv%values), "qv supplied")

   call assert_true(all(shape(field_state%u%values) == [3, 2, 2]), "u shape")
   call assert_true(all(shape(field_state%v%values) == [2, 3, 2]), "v shape")
   call assert_true(all(shape(field_state%w%values) == [2, 2, 3]), "w shape")
   call assert_true(all(shape(field_state%temperature%values) == [2, 2, 2]), &
                    "temperature shape")
   call assert_true(all(shape(field_state%qv%values) == [2, 2, 2]), "qv shape")

   pressure = [100000.0_real64, 98000.0_real64]
   density = [1.2_real64, 1.1_real64]
   dry_density = 1.0_real64
   static_energy = 300.0_real64
   liquid = 0.001_real64
   call lcm_set_field("pressure", reshape(pressure, [size(pressure),1,1]))
   call lcm_set_field("density", reshape(density, [size(density),1,1]))
   call assert_true(all(shape(field_state%pressure%values) == [2,1,1]), "pressure profile storage")
   call assert_true(all(shape(field_state%density%values) == [2,1,1]), "density profile storage")
   call lcm_set_field("dry_air_density", dry_density)
   call lcm_set_field("static_energy", static_energy)
   call lcm_set_field("cloud_liquid", liquid)

   ! Changes to host memory are invisible until another explicit copy.
   u_face = -1.0_real64
   v_face = -2.0_real64
   w_face = -3.0_real64
   temperature = -280.0_real64
   qv = -0.008_real64
   pressure = 0.0_real64
   density = 0.0_real64
   dry_density = 0.0_real64
   static_energy = 0.0_real64
   liquid = 0.0_real64
   call assert_true(all(field_state%u%values == 1.0_real64), "u independent of host")
   call assert_true(all(field_state%v%values == 2.0_real64), "v independent of host")
   call assert_true(all(field_state%w%values == 3.0_real64), "w independent of host")
   call assert_true(all(field_state%temperature%values == 280.0_real64), "temperature independent of host")
   call assert_true(all(field_state%qv%values == 0.008_real64), "qv independent of host")
   call assert_true(all(field_state%pressure%values(:,1,1) == [100000.0_real64, 98000.0_real64]), "pressure independent of host")
   call assert_true(all(field_state%density%values(:,1,1) == [1.2_real64, 1.1_real64]), "density independent of host")
   call assert_true(all(field_state%dry_air_density%values == 1.0_real64), "dry density independent of host")
   call assert_true(all(field_state%static_energy%values == 300.0_real64), "static energy independent of host")
   call assert_true(all(field_state%cloud_liquid%values == 0.001_real64), "liquid independent of host")

   call lcm_set_field("u", u_face)
   call lcm_set_field("v", v_face)
   call lcm_set_field("w", w_face)
   call lcm_set_field("qv", qv)
   call lcm_set_field("pressure", reshape(pressure, [size(pressure),1,1]))
   call lcm_set_field("density", reshape(density, [size(density),1,1]))
   call assert_true(all(field_state%u%values == -1.0_real64), "u explicitly updated")
   call assert_true(all(field_state%v%values == -2.0_real64), "v explicitly updated")
   call assert_true(all(field_state%w%values == -3.0_real64), "w explicitly updated")
   call assert_true(all(field_state%qv%values == -0.008_real64), "qv explicitly updated")
   call assert_true(all(field_state%pressure%values == 0.0_real64), "pressure explicitly updated")
   call assert_true(all(field_state%density%values == 0.0_real64), "density explicitly updated")
   call assert_true(all(field_state%temperature%values == 280.0_real64), "other field unchanged")

   ! The host can release its arrays and pass noncontiguous sections or expressions.
   allocate(temporary_temperature(2,2,2))
   temporary_temperature = 281.0_real64
   call lcm_set_field("temperature", temporary_temperature)
   deallocate(temporary_temperature)
   call assert_true(all(field_state%temperature%values == 281.0_real64), "input lifetime independent")
   strided_temperature = 282.0_real64
   call lcm_set_field("temperature", strided_temperature(1:4:2,:,:))
   strided_temperature = 0.0_real64
   call assert_true(all(field_state%temperature%values == 282.0_real64), "noncontiguous input copied")
   call lcm_set_field("temperature", field_state%temperature%values + 1.0_real64)
   call assert_true(all(field_state%temperature%values == 283.0_real64), "temporary expression copied")

   ! Hosts explicitly pack profiles into the same rank-3 storage convention.
   pressure = [91000.0_real64, 90000.0_real64]
   call lcm_set_field("pressure", reshape(pressure, [2,1,1]))
   pressure = 0.0_real64
   call assert_true(all(field_state%pressure%values(:,1,1) == [91000.0_real64, 90000.0_real64]), &
                    "packed profile independent of host")
   call lcm_set_field("pressure", reshape(pressure, [size(pressure),1,1]))
   call assert_true(all(field_state%pressure%values == 0.0_real64), "reshaped profile updates same storage")

   call assert_true(field_state%temperature%units == "K", "updates preserve units")
   call assert_true(field_state%temperature%long_name == "Air temperature", "updates preserve description")
   call assert_true(all(field_state%pressure%dimensions == [character(len=8) :: "z","",""]), &
                    "updates preserve physical dimensions")
   call assert_true(count(len_trim(field_state%pressure%dimensions) > 0) == 1, &
                    "writer sees one physical profile dimension")

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

end program test_field_copy
