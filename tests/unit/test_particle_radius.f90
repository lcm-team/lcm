!> Verify the standalone named-field replacement for SAM get_r_particle.
program test_particle_radius
   use iso_fortran_env, only: int64, real64
   use lcm_parameters, only: pi, rho_liq
   use lcm_particle_properties, only: lcm_particle_radius
   use lcm_types, only: lcm_particle_t

   implicit none

   call test_aqueous_radius_ignores_small_aerosol_contribution()
   call test_aqueous_and_aerosol_volume_define_radius()
   call test_zero_weight_has_zero_radius()

contains

   subroutine test_aqueous_radius_ignores_small_aerosol_contribution()
      type(lcm_particle_t) :: particle
      real(real64) :: aerosol_density(4)
      real(real64), parameter :: aqueous_radius = 2.0e-6_real64

      aerosol_density = [2000.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
      call set_aqueous_particle(particle, aqueous_radius)
      particle%mass_aero = 0.05_real64 * particle%mass_water
      particle%aerosol_portions(1) = 1.0_real64

      call assert_close(lcm_particle_radius(particle, aerosol_density), &
                        aqueous_radius, "small aerosol mass follows SAM approximation")
   end subroutine test_aqueous_radius_ignores_small_aerosol_contribution

   subroutine test_aqueous_and_aerosol_volume_define_radius()
      type(lcm_particle_t) :: particle
      real(real64) :: aerosol_density(4)
      real(real64), parameter :: aqueous_radius = 2.0e-6_real64
      real(real64), parameter :: combined_radius = 3.0e-6_real64

      aerosol_density = [2000.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
      call set_aqueous_particle(particle, aqueous_radius)
      particle%mass_aero = real(particle%weight, real64) * 4.0_real64 / 3.0_real64 * &
                            pi * aerosol_density(1) * &
                            (combined_radius**3 - aqueous_radius**3)
      particle%aerosol_portions(1) = 1.0_real64

      call assert_close(lcm_particle_radius(particle, aerosol_density), &
                        combined_radius, "large aerosol mass adds aerosol volume")
   end subroutine test_aqueous_and_aerosol_volume_define_radius

   subroutine test_zero_weight_has_zero_radius()
      type(lcm_particle_t) :: particle
      real(real64) :: aerosol_density(4)

      aerosol_density = [2000.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
      call set_aqueous_particle(particle, 2.0e-6_real64)
      particle%weight = 0_int64

      call assert_close(lcm_particle_radius(particle, aerosol_density), 0.0_real64, &
                        "zero-weight particle")
   end subroutine test_zero_weight_has_zero_radius

   subroutine set_aqueous_particle(particle, radius)
      type(lcm_particle_t), intent(out) :: particle
      real(real64), intent(in) :: radius

      particle%weight = 10_int64
      particle%mass_water = real(particle%weight, real64) * 4.0_real64 / 3.0_real64 * &
                            pi * rho_liq * radius**3
   end subroutine set_aqueous_particle

   subroutine assert_close(actual, expected, message)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: message

      if (abs(actual - expected) > 1.0e-18_real64) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_close

end program test_particle_radius
