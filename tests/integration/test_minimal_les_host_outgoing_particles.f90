program test_minimal_les_host_outgoing_particles
   use iso_fortran_env, only: int64, real64
   use lcm_host_interface, only: lcm_advance
   use lcm_exchange_directions, only: lcm_exchange_direction_count
   use lcm_particle_exchange_buffers, only: lcm_send_particle_message_count
   use lcm_particle_state, only: cell_particles, lcm_add_particle, &
                                 lcm_reset_particle_cells, &
                                 n_outgoing_particles
   use lcm_types, only: lcm_particle_t
   use minimal_les_host, only: minimal_host_u, minimal_host_v, &
                               minimal_host_w, setup_minimal_les_host

   implicit none

   integer :: direction

   call setup_minimal_les_host("bin", use_internal_subdomain=.true.)
   call lcm_reset_particle_cells()

   minimal_host_u = 2.0_real64
   minimal_host_v = 2.0_real64
   minimal_host_w = 0.0_real64

   call add_outgoing_particle(9101_int64, 2, 1, 199.0_real64, 125.0_real64)
   call add_outgoing_particle(9102_int64, 1, 2, 125.0_real64, 199.0_real64)
   call add_outgoing_particle(9103_int64, 2, 2, 199.0_real64, 199.0_real64)

   call lcm_advance(1.0_real64)

   call assert_true(total_particle_count() == 3, &
                    "serial exchange preserves all particles")
   call assert_true(n_outgoing_particles == 0, &
                    "serial exchange does not classify outgoing particles")
   do direction = 1, lcm_exchange_direction_count
      call assert_true(lcm_send_particle_message_count(direction) == 0, &
                       "serial exchange leaves send buffers empty")
   enddo

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

   integer function total_particle_count() result(count)
      integer :: i
      integer :: j
      integer :: k

      count = 0
      do k = 1, size(cell_particles, 3)
         do j = 1, size(cell_particles, 2)
            do i = 1, size(cell_particles, 1)
               count = count + cell_particles(i,j,k)%count
            enddo
         enddo
      enddo
   end function total_particle_count

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_true

end program test_minimal_les_host_outgoing_particles
