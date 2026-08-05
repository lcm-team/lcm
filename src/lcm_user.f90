!> Standalone counterpart of SAM-LCM micro_user.f90.
!!
!! These hooks mirror SAM-LCM's compile-time editable user routines while
!! exposing standalone LCM data structures instead of SAM's flat particle(:)
!! storage and n_start/n_end cell slices.
module lcm_user
   use lcm_types, only: lcm_particle_cell_t

   implicit none

   private

   public :: lcm_user_init_particles

contains

   subroutine lcm_user_init_particles(i, j, k, particle_cell)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      type(lcm_particle_cell_t), intent(inout) :: particle_cell

      integer :: n

      ! Default hook is intentionally empty, matching SAM-LCM
      ! micro_user_init_particles. User code may loop over:
      ! do n = 1, particle_cell%count
      !    particle_cell%particles(n)%... = ...
      ! enddo
      do n = 1, particle_cell%count
      enddo
   end subroutine lcm_user_init_particles

end module lcm_user
