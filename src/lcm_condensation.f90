!> Standalone port of the liquid-growth portion of SAM-LCM micro_condensation.
!!
!! The implementation keeps the SAM supersaturation and radius-growth
!! equations. Ice growth and phase transitions remain explicit blockers until
!! their complete SAM branches and dependencies are ported; they are never
!! silently skipped when enabled or when ice particles are present.
module lcm_condensation
   use iso_fortran_env, only: int64, real64
   use lcm_types, only: lcm_particle_cell_t, lcm_particle_t, &
                        liquid_no_in, &
                        liquid_with_in, &
                        ice_no_in, &
                        ice_with_in
   use lcm_parameters, only: pi, rho_liq, rv, specific_heat_dry_air, &
                             molecular_viscosity_air, thermal_conductivity_air, &
                             water_vapor_diffusivity_air, &
                             surface_tension_water, latent_heat_condensation, &
                             alpha_mordy, r_act
   use lcm_config_state, only: config_state
   use lcm_field_state, only: field_state
   use lcm_grid_state, only: grid_state
   use lcm_microphysics_state, only: lcm_get_cell_mean_history, &
      lcm_save_cell_mean_history, lcm_reset_condensation_diagnostics, &
      lcm_add_condensation_diagnostics
   use lcm_particle_state, only: cell_particles, &
                                 particle_aerosol_radius, particle_bfactor, &
                                 cell_tabs_lcm, qsatw_liquid, esatw_liquid
   use lcm_lem, only: lcm_lem_growth_context_t, &
                      lcm_lem_particle_supersaturation, &
                      lcm_lem_apply_growth_feedback

   implicit none
   private

   public :: lcm_apply_condensation

contains

   !> Apply liquid condensational growth and evaporation to local particles.
   !!
   !! `current_step` is SAM's `nstep` counterpart and is used only for its
   !! configured initial supersaturation cap. The host fields are inputs and
   !! remain unchanged.
   subroutine lcm_apply_condensation(time_step, current_step)
      real(real64), intent(in) :: time_step
      integer, intent(in) :: current_step

      integer :: i, j, k
      integer :: nx, ny, nz

      call lcm_reset_condensation_diagnostics()

      ! SAM's micro_condensation always executes ice growth and phase changes
      ! when ice particles exist. Those branches are not yet ported here, so
      ! fail before changing any particle instead of omitting their physics.
      if (config_state%switch_ice) then
         error stop 'LCM condensation: SAM ice-growth/freezing branches are not yet ported'
      endif

      nx = grid_state%grid%nx_local
      ny = grid_state%grid%ny_local
      nz = grid_state%grid%nz

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               call micro_condensation_cell(i, j, k, time_step, current_step)
            enddo
         enddo
      enddo
   end subroutine lcm_apply_condensation

   subroutine micro_condensation_cell(i, j, k, time_step, current_step)
      integer, intent(in) :: i, j, k, current_step
      real(real64), intent(in) :: time_step

      integer :: n, count
      real(real64) :: temperature, pressure, air_density, cell_volume
      real(real64) :: qsatw_pre, esatw_pre, c_liq_pre, a_liq_pre
      real(real64) :: thermal_conductivity, vapor_diffusivity
      real(real64) :: alpha_factor
      real(real64) :: weighted_mass_change

      count = cell_particles(i,j,k)%count
      if (count == 0) then
         call lcm_save_cell_mean_history(i, j, k)
         return
      endif

      temperature = cell_tabs_lcm(i, j, k)
      pressure = field_state%pressure(k)
      air_density = field_state%density(k)
      cell_volume = real(grid_state%dV(k), real64)
      call get_growth_coefficients(temperature, pressure, thermal_conductivity, &
                                   vapor_diffusivity)

      qsatw_pre = qsatw_liquid(temperature, pressure)
      esatw_pre = esatw_liquid(temperature)
      c_liq_pre = 1.0_real64 / &
         (rv * temperature / (100.0_real64 * esatw_pre * vapor_diffusivity) + &
          (latent_heat_condensation / (rv * temperature) - 1.0_real64) * &
          latent_heat_condensation / (thermal_conductivity * temperature))
      a_liq_pre = 1.0_real64 + &
         dtqsatw_sam(temperature, pressure) * &
         latent_heat_condensation / specific_heat_dry_air
      if (config_state%switch_koehler_micro) then
         alpha_factor = 2.0_real64 * surface_tension_water / &
                        (rho_liq * rv * temperature)
      else
         alpha_factor = 0.0_real64
      endif

      do n = 1, count
         call require_liquid_particle(cell_particles(i,j,k)%particles(n))
      enddo

      if (config_state%switch_lem_micro) then
         call apply_lem_condensation_cell(i, j, k, count, time_step, &
            temperature, pressure, air_density, cell_volume, qsatw_pre, &
            c_liq_pre, a_liq_pre, alpha_factor, esatw_pre, &
            thermal_conductivity, vapor_diffusivity)
      else
         call apply_non_lem_condensation_cell(i, j, k, count, time_step, &
            current_step, temperature, pressure, air_density, cell_volume, &
            qsatw_pre, c_liq_pre, a_liq_pre, alpha_factor, esatw_pre, &
            thermal_conductivity, vapor_diffusivity)
      endif

      call lcm_save_cell_mean_history(i, j, k)
   end subroutine micro_condensation_cell

   !> Apply cell-mean supersaturation depletion and growth for non-LEM mode.
   subroutine apply_non_lem_condensation_cell(i, j, k, count, time_step, &
      current_step, temperature, pressure, air_density, cell_volume, qsatw_pre, &
      c_liq_pre, a_liq_pre, alpha_factor, esatw_pre, thermal_conductivity, &
      vapor_diffusivity)
      integer, intent(in) :: i, j, k, count, current_step
      real(real64), intent(in) :: time_step, temperature, pressure, air_density
      real(real64), intent(in) :: cell_volume, qsatw_pre, c_liq_pre, a_liq_pre
      real(real64), intent(in) :: alpha_factor, esatw_pre
      real(real64), intent(in) :: thermal_conductivity, vapor_diffusivity

      integer :: n
      real(real64) :: previous_temperature, previous_pressure
      real(real64) :: radius, r0, ventilation, depletion_sum
      real(real64) :: pre_liq, tau, eta_dyn, eta_last, forcing
      real(real64) :: mean_absolute_supersaturation, supersaturation
      real(real64) :: weighted_mass_change

      call lcm_get_cell_mean_history(i, j, k, previous_temperature, previous_pressure)

      depletion_sum = 0.0_real64
      do n = 1, count
         radius = liquid_radius(cell_particles(i,j,k)%particles(n))
         r0 = growth_length_scale(temperature, esatw_pre, &
            thermal_conductivity, vapor_diffusivity)
         ventilation = liquid_ventilation(radius, &
            cell_particles(i,j,k)%particles(n)%w_sedi, air_density)
         depletion_sum = depletion_sum + &
            real(cell_particles(i,j,k)%particles(n)%weight, real64) * &
            radius**2 / (radius + r0) * ventilation
      enddo

      if (depletion_sum > 0.0_real64) then
         pre_liq = a_liq_pre * 4.0_real64 * pi * c_liq_pre / &
                   (qsatw_pre * cell_volume)
         tau = 1.0_real64 / (pre_liq * depletion_sum)
         eta_dyn = field_state%qv(i,j,k) - qsatw_liquid(temperature, pressure)
         eta_last = field_state%qv(i,j,k) - &
                    qsatw_liquid(previous_temperature, previous_pressure)
         forcing = (eta_dyn - eta_last) / time_step
         mean_absolute_supersaturation = mean_supersaturation( &
            eta_last, forcing, tau, time_step)
         supersaturation = mean_absolute_supersaturation / qsatw_pre
      else
         supersaturation = 0.0_real64
      endif

      if (current_step <= config_state%nstep_cut_supersat .and. &
          supersaturation > real(config_state%supersat_threshold, real64)) then
         supersaturation = real(config_state%supersat_threshold, real64)
      endif

      do n = 1, count
         call integrate_liquid_particle(i, j, k, n, time_step, temperature, &
            air_density, c_liq_pre, alpha_factor, supersaturation, &
            esatw_pre, thermal_conductivity, vapor_diffusivity, &
            cell_particles(i,j,k), weighted_mass_change)
      enddo
   end subroutine apply_non_lem_condensation_cell

   !> Apply particle-local supersaturation depletion, growth, and LEM updates.
   subroutine apply_lem_condensation_cell(i, j, k, count, time_step, &
      temperature, pressure, air_density, cell_volume, qsatw_pre, c_liq_pre, &
      a_liq_pre, alpha_factor, esatw_pre, thermal_conductivity, vapor_diffusivity)
      integer, intent(in) :: i, j, k, count
      real(real64), intent(in) :: time_step, temperature, pressure, air_density
      real(real64), intent(in) :: cell_volume, qsatw_pre, c_liq_pre, a_liq_pre
      real(real64), intent(in) :: alpha_factor, esatw_pre
      real(real64), intent(in) :: thermal_conductivity, vapor_diffusivity

      integer :: n
      real(real64) :: radius, r0
      real(real64) :: bfactor, ventilation, supersaturation
      real(real64) :: weighted_mass_change
      type(lcm_lem_growth_context_t) :: lem_context

      do n = 1, count
         radius = liquid_radius(cell_particles(i,j,k)%particles(n))
         r0 = growth_length_scale(temperature, esatw_pre, &
            thermal_conductivity, vapor_diffusivity)
         if (config_state%switch_koehler_micro) then
            bfactor = particle_bfactor(i, j, k, n)
         else
            bfactor = 0.0_real64
         endif
         ventilation = liquid_ventilation(radius, &
            cell_particles(i,j,k)%particles(n)%w_sedi, air_density)
         lem_context%time_step = time_step
         lem_context%pressure = pressure
         lem_context%cell_volume = cell_volume
         lem_context%air_density = air_density
         lem_context%qsatw_pre = qsatw_pre
         lem_context%c_liq_pre = c_liq_pre
         lem_context%a_liq_pre = a_liq_pre
         lem_context%radius = radius
         lem_context%r0 = r0
         lem_context%ventilation = ventilation
         lem_context%particle_count = count
         supersaturation = lcm_lem_particle_supersaturation( &
            cell_particles(i,j,k)%particles(n), lem_context, mean_supersaturation)

         call integrate_liquid_particle(i, j, k, n, time_step, temperature, &
            air_density, c_liq_pre, alpha_factor, supersaturation, &
            esatw_pre, thermal_conductivity, vapor_diffusivity, &
            cell_particles(i,j,k), weighted_mass_change, &
            local_r0=r0, local_bfactor=bfactor, local_ventilation=ventilation)
         call lcm_lem_apply_growth_feedback(cell_particles(i,j,k)%particles(n), &
                                            lem_context, weighted_mass_change)
      enddo
   end subroutine apply_lem_condensation_cell

   !> Mean absolute supersaturation over one timestep using SAM's TFL formulas.
   real(real64) function mean_supersaturation(eta_last, forcing, tau, time_step) &
      result(eta_mean)
      real(real64), intent(in) :: eta_last, forcing, tau, time_step

      real(real64) :: eta_end, eta_integral, fsw

      eta_end = eta_new(eta_last, forcing, tau, time_step)
      eta_integral = eta_int(eta_last, forcing, tau, time_step)
      fsw = supersaturation_weight(eta_last, eta_end, eta_integral, time_step)
      eta_mean = (eta_integral * (1.0_real64 - fsw) + &
         fsw * (eta_last + eta_end) * time_step * 0.5_real64) / time_step
   end function mean_supersaturation

   real(real64) function eta_int(eta_last, forcing, tau, time_step) result(value)
      real(real64), intent(in) :: eta_last, forcing, tau, time_step

      value = tau * (1.0_real64 - exp(-time_step / tau)) * &
         (eta_last - tau * forcing) + tau * forcing * time_step
   end function eta_int

   real(real64) function eta_new(eta_last, forcing, tau, time_step) result(value)
      real(real64), intent(in) :: eta_last, forcing, tau, time_step

      value = tau * forcing + (eta_last - tau * forcing) * &
         exp(-time_step / tau)
   end function eta_new

   real(real64) function supersaturation_weight(eta_last, eta_end, eta_integral, &
                                                time_step) result(fsw)
      real(real64), intent(in) :: eta_last, eta_end, eta_integral, time_step

      fsw = max(min(1.0_real64, max(&
         (min(eta_last, eta_end) * time_step - eta_integral) * 1.0e35_real64, &
         0.0_real64)), min(1.0_real64, max(&
         (eta_integral - max(eta_last, eta_end) * time_step) * 1.0e35_real64, &
         0.0_real64)))
   end function supersaturation_weight

   subroutine integrate_liquid_particle(i, j, k, n, time_step, temperature, &
      air_density, c_liq, alpha_factor, supersaturation, &
      esatw_pre, thermal_conductivity, &
      vapor_diffusivity, particle_cell, weighted_mass_change, &
      local_r0, local_bfactor, &
      local_ventilation)
      integer, intent(in) :: i, j, k, n
      real(real64), intent(in) :: time_step, temperature
      real(real64), intent(in) :: air_density, c_liq
      real(real64), intent(in) :: alpha_factor, supersaturation
      real(real64), intent(in) :: esatw_pre, thermal_conductivity
      real(real64), intent(in) :: vapor_diffusivity
      type(lcm_particle_cell_t), intent(inout) :: particle_cell
      real(real64), intent(out) :: weighted_mass_change
      real(real64), intent(in), optional :: local_r0, local_bfactor
      real(real64), intent(in), optional :: local_ventilation

      type(lcm_particle_t) :: particle
      real(real64) :: radius, new_radius, aerosol_radius, r0, bfactor
      real(real64) :: ventilation, dm

      particle = particle_cell%particles(n)
      radius = liquid_radius(particle)
      aerosol_radius = particle_aerosol_radius(i, j, k, n)
      r0 = growth_length_scale(temperature, esatw_pre, &
         thermal_conductivity, vapor_diffusivity)
      bfactor = 0.0_real64
      if (config_state%switch_koehler_micro) bfactor = &
         particle_bfactor(i, j, k, n)
      ventilation = liquid_ventilation(radius, particle%w_sedi, air_density)
      if (present(local_r0)) r0 = local_r0
      if (present(local_bfactor)) bfactor = local_bfactor
      if (present(local_ventilation)) ventilation = local_ventilation

      if (config_state%switch_analytical_micro) then
         new_radius = sqrt(max((radius + r0)**2 + &
            2.0_real64 * time_step * c_liq / rho_liq * &
            supersaturation * ventilation, (aerosol_radius + r0)**2)) - r0
      else
         new_radius = radius
         call radius_liquid_euler(new_radius, time_step, r0, c_liq, &
            supersaturation, ventilation, alpha_factor, bfactor, aerosol_radius)
      endif

      particle%mass_water = (4.0_real64 / 3.0_real64) * pi * rho_liq * &
         new_radius**3 * real(particle%weight, real64)
      dm = (4.0_real64 / 3.0_real64) * pi * rho_liq * &
           (new_radius**3 - radius**3)
      weighted_mass_change = dm * real(particle%weight, real64)

      if (new_radius > radius) then
         call lcm_add_condensation_diagnostics(i, j, k, &
            dm * real(particle%weight, real64), 0.0_real64, &
            merge(real(particle%weight, real64), 0.0_real64, &
               radius <= r_act .and. new_radius > r_act), 0.0_real64)
      else
         call lcm_add_condensation_diagnostics(i, j, k, 0.0_real64, &
            dm * real(particle%weight, real64), 0.0_real64, &
            merge(real(particle%weight, real64), 0.0_real64, &
               radius > r_act .and. new_radius <= r_act))
      endif

      particle_cell%particles(n) = particle
   end subroutine integrate_liquid_particle

   subroutine radius_liquid_euler(radius, time_step, r0, c_pre, supersaturation, &
      ventilation, alpha_factor, beta_factor, aerosol_radius)
      real(real64), intent(inout) :: radius
      real(real64), intent(in) :: time_step, r0, c_pre, supersaturation
      real(real64), intent(in) :: ventilation, alpha_factor, beta_factor
      real(real64), intent(in) :: aerosol_radius

      integer :: iteration
      real(real64) :: r_euler, r_euler_old, relative_change
      real(real64) :: dr2dt, d2r2dtdr2, f, dfdr2, dt_euler, t_euler, g_pre

      g_pre = c_pre / rho_liq
      r_euler = radius
      r_euler_old = radius
      dt_euler = time_step
      t_euler = 0.0_real64

      do while (t_euler < time_step - 1.0e-20_real64)
         do iteration = 1, 500
            dr2dt = 2.0_real64 * g_pre * ventilation * &
               (supersaturation - alpha_factor / r_euler + &
                beta_factor / r_euler**3) * &
               r_euler / (r_euler + r0)
            d2r2dtdr2 = g_pre * ventilation * &
               (alpha_factor * r_euler**3 - &
                beta_factor * (3.0_real64 * r_euler + 2.0_real64 * r0) - &
                r_euler**3 * r0 * supersaturation) / &
               (r_euler**4 * (r_euler + r0)**2)

            dt_euler = min(time_step - t_euler, time_step)
            f = r_euler**2 - radius**2 - dt_euler * dr2dt
            dfdr2 = 1.0_real64 - dt_euler * d2r2dtdr2
            r_euler = sqrt(max(r_euler**2 - f / dfdr2, aerosol_radius**2))

            relative_change = abs(r_euler - r_euler_old) / r_euler_old
            r_euler_old = r_euler
            if (relative_change < 1.0e-12_real64) exit
         enddo
         t_euler = t_euler + dt_euler
         radius = r_euler
      enddo
   end subroutine radius_liquid_euler

   real(real64) function liquid_radius(particle) result(radius)
      type(lcm_particle_t), intent(in) :: particle

      radius = (particle%mass_water / &
         ((4.0_real64 / 3.0_real64) * pi * &
          real(particle%weight, real64) * rho_liq))**(1.0_real64 / 3.0_real64)
   end function liquid_radius

   real(real64) function growth_length_scale(temperature, esatw, &
      thermal_conductivity, vapor_diffusivity) result(r0)
      real(real64), intent(in) :: temperature, esatw
      real(real64), intent(in) :: thermal_conductivity, vapor_diffusivity

      r0 = vapor_diffusivity / alpha_mordy * &
         sqrt(2.0_real64 * pi / (rv * temperature)) / &
         (1.0_real64 + vapor_diffusivity * &
          latent_heat_condensation**2 * 100.0_real64 * esatw / &
          (thermal_conductivity * rv**2 * temperature**3))
   end function growth_length_scale

   real(real64) function liquid_ventilation(radius, settling_velocity, &
                                             air_density) result(value)
      real(real64), intent(in) :: radius, settling_velocity, air_density
      real(real64) :: reynolds

      if (radius > 40.0e-6_real64) then
         reynolds = 2.0_real64 * radius * settling_velocity / &
                    (molecular_viscosity_air / air_density)
         if (reynolds > 2.5_real64) then
            value = 0.78_real64 + 0.28_real64 * sqrt(reynolds)
         else
            value = 1.0_real64 + 0.09_real64 * reynolds
         endif
      else
         value = 1.0_real64
      endif
   end function liquid_ventilation

   subroutine get_growth_coefficients(temperature, pressure, thermal_conductivity, &
                                      vapor_diffusivity)
      real(real64), intent(in) :: temperature, pressure
      real(real64), intent(out) :: thermal_conductivity, vapor_diffusivity

      thermal_conductivity = thermal_conductivity_air
      vapor_diffusivity = water_vapor_diffusivity_air
      if (config_state%switch_variable_therco) then
         thermal_conductivity = 2.3823e-2_real64 + &
            7.1177e-5_real64 * (temperature - 273.15_real64)
      endif
      if (config_state%switch_variable_diffelq) then
         vapor_diffusivity = 2.11e-5_real64 * &
            (temperature / 273.15_real64)**1.94_real64 * &
            101325.0_real64 / pressure
      endif
   end subroutine get_growth_coefficients

   real(real64) function dtqsatw_sam(temperature, pressure) result(value)
      real(real64), intent(in) :: temperature, pressure

      value = (qsatw_liquid(temperature + 0.01_real64, pressure) - &
               qsatw_liquid(temperature - 0.01_real64, pressure)) / 0.02_real64
   end function dtqsatw_sam


   subroutine require_liquid_particle(particle)
      type(lcm_particle_t), intent(in) :: particle

      if (particle%particle_type == ice_no_in .or. &
          particle%particle_type == ice_with_in) then
         error stop 'LCM condensation: SAM ice-growth/phase-change branches are not yet ported'
      endif
      if (particle%particle_type /= liquid_no_in .and. &
          particle%particle_type /= liquid_with_in) then
         error stop 'LCM condensation: invalid SAM-compatible particle_type'
      endif
      if (particle%weight <= 0_int64) then
         error stop 'LCM condensation requires positive particle multiplicity'
      endif
      if (particle%mass_water < 0.0_real64) then
         error stop 'LCM condensation requires nonnegative particle water mass'
      endif
   end subroutine require_liquid_particle

end module lcm_condensation
