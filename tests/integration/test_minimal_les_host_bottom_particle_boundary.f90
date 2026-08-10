!> Verify SAM-compatible bottom removal and reflection through lcm_advance.
program test_minimal_les_host_bottom_particle_boundary
   use iso_fortran_env, only: int64, real64
   use lcm_grid_state, only: grid_state
   use lcm_host_interface, only: lcm_advance
   use lcm_parameters, only: pi, rho_liq
   use lcm_particle_properties, only: lcm_particle_radius
   use lcm_particle_state, only: cell_particles, lcm_add_particle, &
                                 lcm_reset_particle_cells
   use lcm_types, only: lcm_particle_t
   use minimal_les_host, only: minimal_host_dx, minimal_host_dy, &
                               minimal_host_u, minimal_host_v, minimal_host_w, &
                               setup_minimal_les_host

   implicit none

   integer(int64), parameter :: reflected_particle_id = 9301_int64
   integer(int64), parameter :: deleted_particle_id = 9302_int64
   integer(int64), parameter :: threshold_particle_id = 9303_int64
   real(real64) :: threshold_radius

   threshold_radius = water_particle_radius(1.0e-6_real64)
   call setup_minimal_les_host("bin", bottom_removal_radius=threshold_radius)
   call lcm_reset_particle_cells()
   call configure_downward_host_velocity()
   call add_reflected_particle()
   call add_deleted_particle()
   call add_threshold_particle()

   call lcm_advance(1.0_real64)

   call assert_true(cell_particles(1,1,1)%count == 2, &
                    "small and threshold particles remain after bottom handling")
   call assert_particle(reflected_particle_id, 1.251_real64, 0.25_real64, &
                        "small particle is reflected at the bottom")
   call assert_particle(threshold_particle_id, 1.501_real64, 0.50_real64, &
                        "threshold-radius particle is reflected at the bottom")
   call assert_particle_absent(deleted_particle_id, &
                               "large particle is removed at the bottom")

contains

   subroutine configure_downward_host_velocity()
      minimal_host_u = 0.0_real64
      minimal_host_v = 0.0_real64
      minimal_host_w = -2.0_real64
   end subroutine configure_downward_host_velocity

   subroutine add_reflected_particle()
      type(lcm_particle_t) :: particle

      particle = lcm_particle_t()
      particle%id = reflected_particle_id
      particle%x = 0.5_real64 * minimal_host_dx
      particle%y = 0.5_real64 * minimal_host_dy
      particle%z = grid_state%zbottom + 1.0_real64
      particle%w_sgs = -0.25_real64
      call set_water_radius(particle, 0.5e-6_real64)
      call lcm_add_particle(1, 1, 1, particle)
   end subroutine add_reflected_particle

   subroutine add_deleted_particle()
      type(lcm_particle_t) :: particle

      particle = lcm_particle_t()
      particle%id = deleted_particle_id
      particle%x = 0.5_real64 * minimal_host_dx
      particle%y = 0.5_real64 * minimal_host_dy
      particle%z = grid_state%zbottom + 1.0_real64
      call set_water_radius(particle, 2.0e-6_real64)
      call lcm_add_particle(1, 1, 1, particle)
   end subroutine add_deleted_particle

   subroutine add_threshold_particle()
      type(lcm_particle_t) :: particle

      particle = lcm_particle_t()
      particle%id = threshold_particle_id
      particle%x = 0.5_real64 * minimal_host_dx
      particle%y = 0.5_real64 * minimal_host_dy
      particle%z = grid_state%zbottom + 1.0_real64
      particle%w_sgs = -0.50_real64
      call set_water_radius(particle, 1.0e-6_real64)
      call lcm_add_particle(1, 1, 1, particle)
   end subroutine add_threshold_particle

   real(real64) function water_particle_radius(nominal_radius) result(radius)
      real(real64), intent(in) :: nominal_radius
      type(lcm_particle_t) :: particle

      particle = lcm_particle_t()
      call set_water_radius(particle, nominal_radius)
      radius = lcm_particle_radius(particle, [0.0_real64, 0.0_real64, &
                                              0.0_real64, 0.0_real64])
   end function water_particle_radius

   subroutine set_water_radius(particle, radius)
      type(lcm_particle_t), intent(inout) :: particle
      real(real64), intent(in) :: radius

      particle%weight = 10_int64
      particle%mass_water = real(particle%weight, real64) * 4.0_real64 / 3.0_real64 * &
                            pi * rho_liq * radius**3
   end subroutine set_water_radius

   subroutine assert_particle(id, expected_z, expected_w_sgs, message)
      integer(int64), intent(in) :: id
      real(real64), intent(in) :: expected_z
      real(real64), intent(in) :: expected_w_sgs
      character(len=*), intent(in) :: message

      integer :: n

      do n = 1, cell_particles(1,1,1)%count
         if (cell_particles(1,1,1)%particles(n)%id == id) then
            call assert_close(cell_particles(1,1,1)%particles(n)%z, expected_z, &
                              message // " z")
            call assert_close(cell_particles(1,1,1)%particles(n)%w_sgs, &
                              expected_w_sgs, message // " w_sgs")
            return
         endif
      enddo

      error stop "assertion failed: " // message // " present"
   end subroutine assert_particle

   subroutine assert_particle_absent(id, message)
      integer(int64), intent(in) :: id
      character(len=*), intent(in) :: message

      integer :: n

      do n = 1, cell_particles(1,1,1)%count
         if (cell_particles(1,1,1)%particles(n)%id == id) then
            error stop "assertion failed: " // message
         endif
      enddo
   end subroutine assert_particle_absent

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

end program test_minimal_les_host_bottom_particle_boundary
