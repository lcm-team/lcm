!> Host-neutral particle location helpers for the standalone cell-owned state.
module lcm_particle_location
   use iso_fortran_env, only: real64
   use lcm_grid_state, only: grid_state, lcm_local_subdomain_bounds

   implicit none

   private

   public :: lcm_locate_particle_cell

contains

   subroutine lcm_locate_particle_cell(x, y, z, i, j, k, is_local, &
                                       global_i, global_j)
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      real(real64), intent(in) :: z
      integer, intent(out) :: i
      integer, intent(out) :: j
      integer, intent(out) :: k
      logical, intent(out) :: is_local
      integer, intent(out), optional :: global_i
      integer, intent(out), optional :: global_j
      integer :: global_i_first
      integer :: global_i_last
      integer :: global_j_first
      integer :: global_j_last

      call require_grid_set()

      call lcm_local_subdomain_bounds(global_i_first, global_i_last, &
                                      global_j_first, global_j_last)

      i = locate_horizontal_cell(x, grid_state%grid%dx, global_i_first, &
                                 global_i_last - global_i_first + 1)
      if (present(global_i)) global_i = global_horizontal_cell(x, grid_state%grid%dx)
      if (grid_state%grid%y_dimension == 1) then
         j = locate_horizontal_cell(y, grid_state%grid%dy, global_j_first, &
                                    global_j_last - global_j_first + 1)
         if (present(global_j)) global_j = global_horizontal_cell(y, grid_state%grid%dy)
      else
         j = 1
         if (present(global_j)) global_j = global_j_first
      endif
      k = locate_vertical_cell(z)

      is_local = i >= 1 .and. i <= grid_state%grid%nx_local .and. &
                 j >= 1 .and. j <= grid_state%grid%ny_local .and. &
                 k >= 1 .and. k <= grid_state%grid%nz
   end subroutine lcm_locate_particle_cell

   integer function locate_horizontal_cell(position, spacing, global_start, &
                                           local_count) result(index)
      real(real64), intent(in) :: position
      real(real64), intent(in) :: spacing
      integer, intent(in) :: global_start
      integer, intent(in) :: local_count

      integer :: global_index

      global_index = global_horizontal_cell(position, spacing)
      index = global_index - global_start + 1

      if (index < 1) then
         index = 0
      elseif (index > local_count) then
         index = local_count + 1
      endif
   end function locate_horizontal_cell

   integer function global_horizontal_cell(position, spacing) result(global_index)
      real(real64), intent(in) :: position
      real(real64), intent(in) :: spacing

      if (spacing <= 0.0_real64) then
         error stop "LCM grid spacing must be positive to locate particles"
      endif

      global_index = floor(position / spacing) + 1
   end function global_horizontal_cell

   integer function locate_vertical_cell(z) result(index)
      real(real64), intent(in) :: z

      integer :: finish
      integer :: mid
      integer :: start

      if (z < grid_state%grid%z_face(1)) then
         index = 0
         return
      endif
      if (z >= grid_state%grid%z_face(grid_state%grid%nz + 1)) then
         index = grid_state%grid%nz + 1
         return
      endif

      start = 1
      finish = grid_state%grid%nz + 1

      do while (finish - start > 1)
         mid = (start + finish) / 2
         if (z >= grid_state%grid%z_face(mid)) then
            start = mid
         else
            finish = mid
         endif
      enddo

      index = start
   end function locate_vertical_cell

   subroutine require_grid_set()
      if (.not. grid_state%is_set) then
         error stop "LCM grid must be initialized before locating particles"
      endif
   end subroutine require_grid_set

end module lcm_particle_location
