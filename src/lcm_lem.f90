!> LEM-specific supersaturation and particle feedback for diffusional growth.
!!
!! The coordinator and shared droplet-growth solver remain in
!! lcm_condensation. This module owns the LEM policy on each side of that
!! shared solver: it computes particle supersaturation, then applies the
!! solver's weighted mass change to particle vapor, temperature, and eta.
module lcm_lem
   use iso_fortran_env, only: real64
   use lcm_types, only: lcm_particle_t
   use lcm_parameters, only: pi, fac_cond, gravity, specific_heat_dry_air
   use lcm_particle_state, only: qsatw_liquid

   implicit none
   private

   public :: lcm_lem_growth_context_t
   public :: lcm_lem_particle_supersaturation
   public :: lcm_lem_apply_growth_feedback

   !> Explicit process inputs needed by the LEM growth policy.
   type :: lcm_lem_growth_context_t
      real(real64) :: time_step
      real(real64) :: pressure
      real(real64) :: cell_volume
      real(real64) :: air_density
      real(real64) :: qsatw_pre
      real(real64) :: c_liq_pre
      real(real64) :: a_liq_pre
      real(real64) :: radius
      real(real64) :: r0
      real(real64) :: ventilation
      integer :: particle_count
   end type lcm_lem_growth_context_t

   !> Callback for the shared TFL mean-supersaturation equations.
   abstract interface
      real(real64) function tfl_mean_supersaturation(eta_last, forcing, tau, &
                                                      time_step) result(eta_mean)
         import :: real64
         real(real64), intent(in) :: eta_last, forcing, tau, time_step
      end function tfl_mean_supersaturation
   end interface

contains

   !> Calculate one particle's mean supersaturation over the timestep.
   real(real64) function lcm_lem_particle_supersaturation(particle, context, &
      mean_supersaturation) result(supersaturation)
      type(lcm_particle_t), intent(in) :: particle
      type(lcm_lem_growth_context_t), intent(in) :: context
      procedure(tfl_mean_supersaturation) :: mean_supersaturation

      real(real64) :: particle_temperature, pre_liq, tau
      real(real64) :: eta_dyn, forcing, eta_mean

      particle_temperature = particle%temperature - &
         particle%z * gravity / specific_heat_dry_air
      pre_liq = context%a_liq_pre * 4.0_real64 * pi * context%c_liq_pre / &
                (context%qsatw_pre * &
                 (context%cell_volume / real(context%particle_count, real64)))
      tau = 1.0_real64 / (pre_liq * &
         real(particle%weight, real64) * context%radius**2 / &
         (context%radius + context%r0) * context%ventilation)
      eta_dyn = particle%qv - qsatw_liquid(particle_temperature, context%pressure)
      forcing = (eta_dyn - particle%eta) / context%time_step
      eta_mean = mean_supersaturation(particle%eta, forcing, tau, context%time_step)
      supersaturation = eta_mean / &
                        qsatw_liquid(particle_temperature, context%pressure)
   end function lcm_lem_particle_supersaturation

   !> Apply growth feedback to particle-local vapor, heat, and eta.
   subroutine lcm_lem_apply_growth_feedback(particle, context, weighted_mass_change)
      type(lcm_particle_t), intent(inout) :: particle
      type(lcm_lem_growth_context_t), intent(in) :: context
      real(real64), intent(in) :: weighted_mass_change

      real(real64) :: mixing_ratio_change, particle_temperature

      mixing_ratio_change = weighted_mass_change / &
         ((context%cell_volume / real(context%particle_count, real64)) * &
          context%air_density)
      particle%temperature = particle%temperature + mixing_ratio_change * fac_cond
      particle%qv = particle%qv - mixing_ratio_change
      particle_temperature = particle%temperature - &
         particle%z * gravity / specific_heat_dry_air
      particle%eta = particle%qv - &
         qsatw_liquid(particle_temperature, context%pressure)
   end subroutine lcm_lem_apply_growth_feedback

end module lcm_lem
