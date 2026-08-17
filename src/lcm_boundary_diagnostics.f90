!> Per-call physical particle-boundary loss diagnostics.
!!
!! Values are local to the current rank and are reset for every call to
!! lcm_apply_particle_boundaries. Host feedback decides later whether and how
!! to reduce, normalize, or convert these losses to a precipitation flux.
module lcm_boundary_diagnostics
   use iso_fortran_env, only: int64, real64

   implicit none

   private

   public :: boundary_diagnostics
   public :: lcm_record_boundary_removals
   public :: lcm_reset_boundary_diagnostics

   type :: lcm_boundary_diagnostics_t
      integer(int64) :: deleted_particle_count = 0_int64
      real(real64) :: deleted_water_mass = 0.0_real64
      real(real64) :: deleted_aerosol_mass = 0.0_real64
   end type lcm_boundary_diagnostics_t

   type(lcm_boundary_diagnostics_t), protected, save :: boundary_diagnostics

contains

   subroutine lcm_reset_boundary_diagnostics()
      boundary_diagnostics = lcm_boundary_diagnostics_t()
   end subroutine lcm_reset_boundary_diagnostics

   subroutine lcm_record_boundary_removals(deleted_particle_count, &
                                           deleted_water_mass, &
                                           deleted_aerosol_mass)
      integer(int64), intent(in) :: deleted_particle_count
      real(real64), intent(in) :: deleted_water_mass
      real(real64), intent(in) :: deleted_aerosol_mass

      boundary_diagnostics%deleted_particle_count = &
         boundary_diagnostics%deleted_particle_count + deleted_particle_count
      boundary_diagnostics%deleted_water_mass = &
         boundary_diagnostics%deleted_water_mass + deleted_water_mass
      boundary_diagnostics%deleted_aerosol_mass = &
         boundary_diagnostics%deleted_aerosol_mass + deleted_aerosol_mass
   end subroutine lcm_record_boundary_removals

end module lcm_boundary_diagnostics
