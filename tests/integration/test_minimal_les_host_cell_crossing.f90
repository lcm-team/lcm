program test_minimal_les_host_cell_crossing
   use iso_fortran_env, only: int64, real64
   use lcm_host_interface, only: lcm_advance
   use lcm_particle_state, only: cell_particles, lcm_add_particle
   use lcm_types, only: lcm_particle_t
   use minimal_les_host, only: minimal_host_dx, minimal_host_dy, &
                               minimal_host_np_xy, minimal_host_u, &
                               minimal_host_v, minimal_host_w, &
                               minimal_host_z_face, setup_minimal_les_host

   implicit none

   integer(int64), parameter :: crossing_particle_id = 9001_int64
   character(len=8) :: direction

   call get_command_argument(1, direction)
   if (len_trim(direction) == 0) direction = "x"
   call setup_minimal_les_host("bin")
   call configure_crossing_velocity_field(trim(direction))
   call add_particle_near_cell_face(trim(direction))

   call assert_initial_cell_counts(trim(direction))

   call lcm_advance(1.0_real64)

   call assert_final_cell_counts(trim(direction))
   call assert_crossing_particle_in_new_cell(trim(direction))

contains

   subroutine configure_crossing_velocity_field(direction)
      character(len=*), intent(in) :: direction

      minimal_host_u = 0.0_real64
      minimal_host_v = 0.0_real64
      minimal_host_w = 0.0_real64

      select case (direction)
      case ("x")
         minimal_host_u = 2.0_real64
      case ("y")
         minimal_host_v = 2.0_real64
      case ("z")
         minimal_host_w = 2.0_real64
      case default
         error stop "unknown minimal LES host crossing direction"
      end select
   end subroutine configure_crossing_velocity_field

   subroutine add_particle_near_cell_face(direction)
      character(len=*), intent(in) :: direction
      type(lcm_particle_t) :: particle

      particle = lcm_particle_t()
      particle%id = crossing_particle_id
      particle%y = 0.5_real64 * minimal_host_dy
      particle%x = 0.5_real64 * minimal_host_dx
      particle%z = 0.5_real64 * (minimal_host_z_face(1) + minimal_host_z_face(2))

      select case (direction)
      case ("x")
         particle%x = minimal_host_dx - 1.0_real64
      case ("y")
         particle%y = minimal_host_dy - 1.0_real64
      case ("z")
         particle%z = minimal_host_z_face(2) - 1.0_real64
      case default
         error stop "unknown minimal LES host crossing direction"
      end select

      call lcm_add_particle(1, 1, 1, particle)
   end subroutine add_particle_near_cell_face

   subroutine assert_initial_cell_counts(direction)
      character(len=*), intent(in) :: direction

      call assert_true(cell_particles(1,1,1)%count == minimal_host_np_xy + 1, &
                       "crossing particle starts in source cell")

      select case (direction)
      case ("x")
         call assert_true(cell_particles(2,1,1)%count == minimal_host_np_xy, &
                          "x target starts with initialized particles")
      case ("y")
         call assert_true(cell_particles(1,2,1)%count == minimal_host_np_xy, &
                          "y target starts with initialized particles")
      case ("z")
         call assert_true(cell_particles(1,1,2)%count == 0, &
                          "z target starts empty above initialized levels")
      end select
   end subroutine assert_initial_cell_counts

   subroutine assert_final_cell_counts(direction)
      character(len=*), intent(in) :: direction

      call assert_true(cell_particles(1,1,1)%count == minimal_host_np_xy, &
                       "crossing particle removed from old cell")

      select case (direction)
      case ("x")
         call assert_true(cell_particles(2,1,1)%count == minimal_host_np_xy + 1, &
                          "x crossing particle appended to new local cell")
      case ("y")
         call assert_true(cell_particles(1,2,1)%count == minimal_host_np_xy + 1, &
                          "y crossing particle appended to new local cell")
      case ("z")
         call assert_true(cell_particles(1,1,2)%count == 1, &
                          "z crossing particle appended to new local cell")
      end select
   end subroutine assert_final_cell_counts

   subroutine assert_crossing_particle_in_new_cell(direction)
      character(len=*), intent(in) :: direction
      integer :: n
      integer :: target_i
      integer :: target_j
      integer :: target_k
      logical :: found
      real(real64) :: expected_x
      real(real64) :: expected_y
      real(real64) :: expected_z

      call expected_crossing_state(direction, target_i, target_j, target_k, &
                                   expected_x, expected_y, expected_z)
      found = .false.
      do n = 1, cell_particles(target_i,target_j,target_k)%count
         if (cell_particles(target_i,target_j,target_k)%particles(n)%id == &
             crossing_particle_id) then
            found = .true.
            call assert_near(cell_particles(target_i,target_j,target_k)%particles(n)%x, &
                             expected_x, &
                             "crossing particle x after advance")
            call assert_near(cell_particles(target_i,target_j,target_k)%particles(n)%y, &
                             expected_y, &
                             "crossing particle y after advance")
            call assert_near(cell_particles(target_i,target_j,target_k)%particles(n)%z, &
                             expected_z, &
                             "crossing particle z after advance")
         endif
      enddo

      call assert_true(found, "crossing particle found in new cell")
   end subroutine assert_crossing_particle_in_new_cell

   subroutine expected_crossing_state(direction, target_i, target_j, target_k, &
                                     expected_x, expected_y, expected_z)
      character(len=*), intent(in) :: direction
      integer, intent(out) :: target_i
      integer, intent(out) :: target_j
      integer, intent(out) :: target_k
      real(real64), intent(out) :: expected_x
      real(real64), intent(out) :: expected_y
      real(real64), intent(out) :: expected_z

      target_i = 1
      target_j = 1
      target_k = 1
      expected_x = 0.5_real64 * minimal_host_dx
      expected_y = 0.5_real64 * minimal_host_dy
      expected_z = 0.5_real64 * (minimal_host_z_face(1) + minimal_host_z_face(2))

      select case (direction)
      case ("x")
         target_i = 2
         expected_x = minimal_host_dx + 1.0_real64
      case ("y")
         target_j = 2
         expected_y = minimal_host_dy + 1.0_real64
      case ("z")
         target_k = 2
         expected_z = minimal_host_z_face(2) + 1.0_real64
      case default
         error stop "unknown minimal LES host crossing direction"
      end select
   end subroutine expected_crossing_state

   subroutine assert_near(actual, expected, message)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: message

      call assert_true(abs(actual - expected) <= 1.0e-12_real64, message)
   end subroutine assert_near

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_true

end program test_minimal_les_host_cell_crossing
