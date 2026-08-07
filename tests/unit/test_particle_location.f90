program test_particle_location
   use iso_fortran_env, only: real64
   use lcm_particle_location, only: lcm_locate_particle_cell
   use minimal_les_host, only: minimal_host_dx, minimal_host_dy, &
                               minimal_host_z_face, setup_minimal_les_host

   implicit none

   call setup_minimal_les_host()

   call assert_location(0.0_real64, 0.0_real64, minimal_host_z_face(1), &
                        1, 1, 1, .true., "lower local corner")
   call assert_location(minimal_host_dx - 0.001_real64, &
                        minimal_host_dy - 0.001_real64, &
                        minimal_host_z_face(2) - 0.001_real64, &
                        1, 1, 1, .true., "inside first cell")
   call assert_location(minimal_host_dx, minimal_host_dy, &
                        minimal_host_z_face(2), &
                        2, 2, 2, .true., "exact inner faces")
   call assert_location(2.0_real64 * minimal_host_dx - 0.001_real64, &
                        2.0_real64 * minimal_host_dy - 0.001_real64, &
                        minimal_host_z_face(3) - 0.001_real64, &
                        2, 2, 2, .true., "inside upper local cell")
   call assert_location(-0.001_real64, 0.0_real64, minimal_host_z_face(1), &
                        0, 1, 1, .false., "west of local subdomain")
   call assert_location(2.0_real64 * minimal_host_dx, 0.0_real64, &
                        minimal_host_z_face(1), &
                        3, 1, 1, .false., "east of local subdomain")
   call assert_location(0.0_real64, -0.001_real64, minimal_host_z_face(1), &
                        1, 0, 1, .false., "south of local subdomain")
   call assert_location(0.0_real64, 2.0_real64 * minimal_host_dy, &
                        minimal_host_z_face(1), &
                        1, 3, 1, .false., "north of local subdomain")
   call assert_location(0.0_real64, 0.0_real64, &
                        minimal_host_z_face(1) - 0.001_real64, &
                        1, 1, 0, .false., "below vertical domain")
   call assert_location(0.0_real64, 0.0_real64, minimal_host_z_face(3), &
                        1, 1, 3, .false., "above vertical domain")

contains

   subroutine assert_location(x, y, z, expected_i, expected_j, expected_k, &
                              expected_is_local, message)
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      real(real64), intent(in) :: z
      integer, intent(in) :: expected_i
      integer, intent(in) :: expected_j
      integer, intent(in) :: expected_k
      logical, intent(in) :: expected_is_local
      character(len=*), intent(in) :: message

      integer :: i
      integer :: j
      integer :: k
      logical :: is_local

      call lcm_locate_particle_cell(x, y, z, i, j, k, is_local)

      call assert_true(i == expected_i, message // " i")
      call assert_true(j == expected_j, message // " j")
      call assert_true(k == expected_k, message // " k")
      call assert_true(is_local .eqv. expected_is_local, &
                       message // " local status")
   end subroutine assert_location

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_true

end program test_particle_location
