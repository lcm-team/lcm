program test_minimal_les_host_outgoing_particles
   use iso_fortran_env, only: int64, real64
   use lcm_host_interface, only: lcm_advance
   use lcm_exchange_directions, only: lcm_exchange_east, lcm_exchange_north, &
                                      lcm_exchange_northeast
   use lcm_particle_state, only: cell_particles, lcm_add_particle, &
                                 lcm_reset_particle_cells, &
                                 n_outgoing_particles, outgoing_particles
   use lcm_types, only: lcm_particle_t
   use minimal_les_host, only: minimal_host_u, minimal_host_v, &
                               minimal_host_w, setup_minimal_les_host

   implicit none

   call setup_minimal_les_host("bin", use_internal_subdomain=.true.)
   call lcm_reset_particle_cells()

   minimal_host_u = 2.0_real64
   minimal_host_v = 2.0_real64
   minimal_host_w = 0.0_real64

   call add_outgoing_particle(9101_int64, 2, 1, 199.0_real64, 125.0_real64)
   call add_outgoing_particle(9102_int64, 1, 2, 125.0_real64, 199.0_real64)
   call add_outgoing_particle(9103_int64, 2, 2, 199.0_real64, 199.0_real64)

   call lcm_advance(1.0_real64)

   call assert_true(cell_particles(2,1,1)%count == 0, "east particle removed")
   call assert_true(cell_particles(1,2,1)%count == 0, "north particle removed")
   call assert_true(cell_particles(2,2,1)%count == 0, "north-east particle removed")
   call assert_true(n_outgoing_particles == 3, "outgoing particle count")

   call assert_outgoing_particle(9101_int64, lcm_exchange_east, "east outgoing particle")
   call assert_outgoing_particle(9102_int64, lcm_exchange_north, "north outgoing particle")
   call assert_outgoing_particle(9103_int64, lcm_exchange_northeast, "north-east outgoing particle")

contains

   subroutine add_outgoing_particle(id, i, j, x, y)
      integer(int64), intent(in) :: id
      integer, intent(in) :: i
      integer, intent(in) :: j
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      type(lcm_particle_t) :: particle

      particle = lcm_particle_t()
      particle%id = id
      particle%x = x
      particle%y = y
      particle%z = 50.0_real64
      call lcm_add_particle(i, j, 1, particle)
   end subroutine add_outgoing_particle

   subroutine assert_outgoing_particle(id, expected_direction, message)
      integer(int64), intent(in) :: id
      integer, intent(in) :: expected_direction
      character(len=*), intent(in) :: message
      integer :: n
      logical :: found

      found = .false.
      do n = 1, n_outgoing_particles
         if (outgoing_particles(n)%particle%id == id) then
            found = .true.
            call assert_true(outgoing_particles(n)%direction == expected_direction, &
                             message // " direction")
         endif
      enddo

      call assert_true(found, message // " present")
   end subroutine assert_outgoing_particle

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_true

end program test_minimal_les_host_outgoing_particles
