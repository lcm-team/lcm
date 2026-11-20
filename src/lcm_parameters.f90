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
   public :: gravity
   public :: molecular_viscosity_air
   public :: thermal_conductivity_air
   public :: water_vapor_diffusivity_air
   public :: latent_heat_sublimation
   public :: latent_heat_fusion
   public :: fac_sub
   public :: fac_fus
   public :: r_act
   public :: alpha_mordy

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
   real(real64), parameter :: gravity = 9.81_real64
   real(real64), parameter :: molecular_viscosity_air = 1.8e-5_real64
   real(real64), parameter :: thermal_conductivity_air = 2.40e-2_real64
   real(real64), parameter :: water_vapor_diffusivity_air = 2.21e-5_real64
   real(real64), parameter :: latent_heat_sublimation = 2.834e6_real64
   real(real64), parameter :: latent_heat_fusion = 3.34e5_real64
   real(real64), parameter :: fac_sub = &
      latent_heat_sublimation / specific_heat_dry_air
   real(real64), parameter :: fac_fus = &
      latent_heat_fusion / specific_heat_dry_air
   real(real64), parameter :: r_act = 1.0e-6_real64
   real(real64), parameter :: alpha_mordy = 0.036_real64

end module lcm_parameters
