!> Verify top reflection through the minimal LES host and public LCM advance.
program test_minimal_les_host_top_particle_boundary
   use iso_fortran_env, only: int64, real64
   use lcm_grid_state, only: grid_state
   use lcm_particle_state, only: cell_particles, lcm_add_particle
   use lcm_types, only: lcm_particle_t
   use minimal_les_host, only: minimal_host_dx, minimal_host_dy, &
                               minimal_host_u, minimal_host_v, minimal_host_w, &
                               setup_minimal_les_host, &
                               advance_minimal_les_transport

   implicit none

   integer(int64), parameter :: top_crossing_particle_id = 9201_int64
   integer(int64), parameter :: in_domain_particle_id = 9202_int64

   call setup_minimal_les_host("bin")
   call configure_upward_host_velocity()
   call add_top_crossing_particle()
   call add_in_domain_particle()

   call advance_minimal_les_transport(1.0_real64)

   call assert_true(cell_particles(1,1,2)%count == 2, &
                    "both integration-test particles remain in the upper cell")
   call assert_particle(top_crossing_particle_id, grid_state%ztop - 1.25_real64, &
                        -0.25_real64, "particle crossing the top is reflected")
   call assert_particle(in_domain_particle_id, 152.0_real64, 0.0_real64, &
                        "in-domain particle only follows host motion")

contains

   subroutine configure_upward_host_velocity()
      minimal_host_u = 0.0_real64
      minimal_host_v = 0.0_real64
      minimal_host_w = 2.0_real64
   end subroutine configure_upward_host_velocity

   subroutine add_top_crossing_particle()
      type(lcm_particle_t) :: particle

      particle = lcm_particle_t()
      particle%id = top_crossing_particle_id
      particle%x = 0.5_real64 * minimal_host_dx
      particle%y = 0.5_real64 * minimal_host_dy
      particle%z = grid_state%ztop - 1.0_real64
      particle%w_sgs = 0.25_real64
      call lcm_add_particle(1, 1, 2, particle)
   end subroutine add_top_crossing_particle

   subroutine add_in_domain_particle()
      type(lcm_particle_t) :: particle

      particle = lcm_particle_t()
      particle%id = in_domain_particle_id
      particle%x = 0.5_real64 * minimal_host_dx
      particle%y = 0.5_real64 * minimal_host_dy
      particle%z = 150.0_real64
      call lcm_add_particle(1, 1, 2, particle)
   end subroutine add_in_domain_particle

   subroutine assert_particle(id, expected_z, expected_w_sgs, message)
      integer(int64), intent(in) :: id
      real(real64), intent(in) :: expected_z
      real(real64), intent(in) :: expected_w_sgs
      character(len=*), intent(in) :: message

      integer :: n
      logical :: found

      found = .false.
      do n = 1, cell_particles(1,1,2)%count
         if (cell_particles(1,1,2)%particles(n)%id == id) then
            found = .true.
            call assert_close(cell_particles(1,1,2)%particles(n)%z, expected_z, &
                              message // " z")
            call assert_close(cell_particles(1,1,2)%particles(n)%w_sgs, &
                              expected_w_sgs, message // " w_sgs")
         endif
      enddo

      call assert_true(found, message // " present")
   end subroutine assert_particle

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

end program test_minimal_les_host_top_particle_boundary
