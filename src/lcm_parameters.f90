!> Shared physical constants for standalone LCM.
module lcm_parameters
   use iso_fortran_env, only: real64

   implicit none

   private

   public :: pi
   public :: epsv
   public :: latent_heat_condensation
   public :: specific_heat_dry_air
   public :: fac_cond
   public :: molecular_weight_water
   public :: rho_liq
   public :: rv
   public :: surface_tension_water

   real(real64), parameter :: pi = &
      3.141592653589793238462643383279502884_real64
   real(real64), parameter :: epsv = 0.622_real64
   real(real64), parameter :: latent_heat_condensation = 2.5e6_real64
   real(real64), parameter :: specific_heat_dry_air = 1004.0_real64
   real(real64), parameter :: fac_cond = &
      latent_heat_condensation / specific_heat_dry_air
   real(real64), parameter :: molecular_weight_water = 18.0_real64
   real(real64), parameter :: rho_liq = 1000.0_real64
   real(real64), parameter :: rv = 461.5_real64
   real(real64), parameter :: surface_tension_water = 0.0728_real64

end module lcm_parameters
