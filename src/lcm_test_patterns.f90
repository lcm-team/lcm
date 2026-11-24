!> Contains routines for creating (eulerian) test patterns
!!
module lcm_test_patterns
   use iso_fortran_env, only: real32, real64

   use lcm_types,          only: lcm_grid_t
   use lcm_grid_state,     only: grid_state
   use lcm_field_state,    only: field_state
   use lcm_parallel_state, only: parallel_state, lcm_parallel_t

   implicit none

   private

   public :: lcm_set_qv_test_pattern

contains

   !> adds to qv fields a global test pattern
   !> the option 'smiley' creates a horizontal smiley that rotates with height
   subroutine lcm_set_qv_test_pattern(pattern_name)
      character(*), intent(in) :: pattern_name
      integer :: i, j, ig, jg, k
      real(real64), parameter :: pi = 3.14159265358d0
      real(real64) :: s, u, v, phi, up, vp, r, rm, w, f
      real(real64) :: outline, eyes, mouth

      if (trim(pattern_name) == 'smiley') then
         do k = 1, grid_state%grid%nz
            do j = 1, grid_state%grid%ny_local
               jg = j + grid_state%grid%global_j_start - 1
               do i = 1, grid_state%grid%nx_local
                  ig = i + grid_state%grid%global_i_start - 1

                  s   = real(max(grid_state%grid%nx_global, grid_state%grid%ny_global), real64)
                  w   = max(0.015d0, 2.0d0/s)
                  u   = (ig - 0.5d0 - 0.5d0*grid_state%grid%nx_global) / s
                  v   = (jg - 0.5d0 - 0.5d0*grid_state%grid%ny_global) / s

                  phi =  2.0d0*pi * real(k-1, real64) / real(grid_state%grid%nz, real64)
                  up  =  cos(phi)*u + sin(phi)*v
                  vp  = -sin(phi)*u + cos(phi)*v
                  r   = sqrt(up**2 + vp**2)

                  outline = exp(-((r - 0.40d0)/w)**2)

                  eyes = exp(-((up + 0.14d0)**2 + (vp - 0.12d0)**2) / 0.04d0**2) &
                       + exp(-((up - 0.14d0)/0.07d0)**2 - ((vp - 0.12d0)/w)**2)

                  rm    = sqrt(up**2 + (vp - 0.02d0)**2)
                  mouth = exp(-((rm - 0.24d0)/w)**2) * 0.5d0*(1.0d0 - tanh((vp - 0.02d0)/w))

                  f = outline + eyes + mouth + 0.2d0*u + 0.1d0*v + 0.1d0*real(k-1, real64)/real(grid_state%grid%nz, real64)
                  field_state%qv(i,j,k) = field_state%qv(i,j,k) + f * 1.0e-4
               enddo
            enddo
         enddo
      endif
   end subroutine lcm_set_qv_test_pattern

end module lcm_test_patterns
