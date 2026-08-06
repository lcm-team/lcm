!> SAM-LCM-compatible interpolation from attached host velocity fields.
module lcm_velocity_interpolation
   use iso_fortran_env, only: real64
   use lcm_field_state, only: field_state
   use lcm_grid_state, only: grid_state
   use lcm_particle_location, only: lcm_locate_particle_cell

   implicit none

   private

   public :: lcm_interpolate_particle_velocity

contains

   subroutine lcm_interpolate_particle_velocity(x, y, z, u, v, w, is_local, &
                                                clamp_to_local)
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      real(real64), intent(in) :: z
      real(real64), intent(out) :: u
      real(real64), intent(out) :: v
      real(real64), intent(out) :: w
      logical, intent(out) :: is_local
      logical, intent(in), optional :: clamp_to_local

      integer :: i
      integer :: j
      integer :: k
      logical :: use_clamped_cell
      real(real64) :: alpha
      real(real64) :: beta
      real(real64) :: gamma

      call require_velocity_fields_attached()

      call lcm_locate_particle_cell(x, y, z, i, j, k, is_local)
      if (.not. is_local) then
         use_clamped_cell = .false.
         if (present(clamp_to_local)) use_clamped_cell = clamp_to_local
         if (use_clamped_cell) then
            i = max(min(i, grid_state%grid%nx_local), 1)
            j = max(min(j, grid_state%grid%ny_local), 1)
            k = max(min(k, grid_state%grid%nz), 1)
         else
            u = 0.0_real64
            v = 0.0_real64
            w = 0.0_real64
            return
         endif
      endif

      alpha = clamp_unit_interval((x - real(grid_state%gpxs(i), real64)) / &
                                  grid_state%grid%dx)
      u = field_state%u(i, j, k) * (1.0_real64 - alpha) + &
          field_state%u(i + 1, j, k) * alpha

      if (grid_state%grid%y_dimension == 1) then
         beta = clamp_unit_interval((y - real(grid_state%gpys(j), real64)) / &
                                    grid_state%grid%dy)
         v = field_state%v(i, j, k) * (1.0_real64 - beta) + &
             field_state%v(i, j + 1, k) * beta
      else
         v = field_state%v(i, j, k)
      endif

      gamma = clamp_unit_interval((z - grid_state%grid%z_face(k)) / &
                                  (grid_state%grid%z_face(k + 1) - &
                                   grid_state%grid%z_face(k)))
      w = field_state%w(i, j, k) * (1.0_real64 - gamma) + &
          field_state%w(i, j, k + 1) * gamma
   end subroutine lcm_interpolate_particle_velocity

   real(real64) function clamp_unit_interval(value) result(clamped)
      real(real64), intent(in) :: value

      clamped = max(min(value, 1.0_real64), 0.0_real64)
   end function clamp_unit_interval

   subroutine require_velocity_fields_attached()
      if (.not. grid_state%is_set) then
         error stop "LCM grid must be initialized before interpolating velocity"
      endif
      if (.not. field_state%u_is_attached) then
         error stop "LCM u field must be attached before interpolating velocity"
      endif
      if (.not. field_state%v_is_attached) then
         error stop "LCM v field must be attached before interpolating velocity"
      endif
      if (.not. field_state%w_is_attached) then
         error stop "LCM w field must be attached before interpolating velocity"
      endif
   end subroutine require_velocity_fields_attached

end module lcm_velocity_interpolation
