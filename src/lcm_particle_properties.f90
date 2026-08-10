!> Named-field physical properties of standalone LCM particles.
module lcm_particle_properties
   use iso_fortran_env, only: real64
   use lcm_parameters, only: pi, rho_liq
   use lcm_types, only: lcm_particle_t

   implicit none

   private

   public :: lcm_particle_radius

contains

   !> Return the physical radius of one computational particle.
   !>
   !> This is the named-field counterpart of SAM-LCM's get_r_particle(n).
   !> The water mass and aerosol mass are total super-droplet masses, so both
   !> are divided by weight before being converted to volume.  In
   !> micro_vars.f90, SAM adds aerosol volume only when iN / iM > 0.1; this
   !> routine preserves that approximation as mass_aero / mass_water > 0.1.
   pure real(real64) function lcm_particle_radius(particle, aerosol_density) &
      result(radius)
      type(lcm_particle_t), intent(in) :: particle
      real(real64), intent(in) :: aerosol_density(:)

      real(real64) :: particle_volume
      real(real64) :: aerosol_volume
      real(real64) :: weight
      integer :: mode

      if (particle%weight <= 0) then
         radius = 0.0_real64
         return
      endif

      weight = real(particle%weight, real64)
      particle_volume = particle%mass_water / weight / rho_liq

      if (particle%mass_water > 0.0_real64 .and. &
          particle%mass_aero / particle%mass_water > 0.1_real64) then
         aerosol_volume = 0.0_real64
         do mode = 1, min(size(particle%aerosol_portions), size(aerosol_density))
            if (aerosol_density(mode) > 0.0_real64) then
               aerosol_volume = aerosol_volume + particle%aerosol_portions(mode) * &
                                particle%mass_aero / weight / aerosol_density(mode)
            endif
         enddo
         particle_volume = particle_volume + aerosol_volume
      endif

      radius = (3.0_real64 * max(particle_volume, 0.0_real64) / &
                (4.0_real64 * pi))**(1.0_real64 / 3.0_real64)
   end function lcm_particle_radius

end module lcm_particle_properties
