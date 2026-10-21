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
   use lcm_parameters, only: pi, rho_liq, rv, specific_heat_dry_air, fac_cond, gravity, &
                             molecular_viscosity_air, thermal_conductivity_air, &
                             water_vapor_diffusivity_air, &
                             surface_tension_water, latent_heat_condensation, &
                             alpha_mordy, molecular_weight_ffc, f_ffc, &
                             n_avogadro, r_act
   use lcm_config_state, only: config_state
   use lcm_field_state, only: field_state
   use lcm_grid_state, only: grid_state
   use lcm_microphysics_state, only: lcm_get_cell_mean_history, &
      lcm_save_cell_mean_history, lcm_reset_condensation_diagnostics, &
      lcm_add_condensation_diagnostics
   use lcm_particle_state, only: cell_particles, &
                                 particle_aerosol_radius, particle_bfactor, &
                                 cell_tabs_lcm, qsatw_liquid, esatw_liquid

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
      real(real64) :: qsatw_pre, esatw_pre, c_liq_pre, a_liq_pre, d_pre
      real(real64) :: thermal_conductivity, vapor_diffusivity
      real(real64) :: ql_local, radiation, radiation_integral, qrad_liq
      real(real64) :: alpha_factor, beta_factor
      real(real64) :: radius, new_radius, aerosol_radius, r0
      real(real64) :: ventilation, q_absorb, dsd_liq, pre_liq, tau
      real(real64) :: eta_dyn, eta_last, forcing, eta_end, eta_integral
      real(real64) :: eta_mean, fsw, supersaturation
      real(real64) :: particle_temperature, r0_local
      real(real64) :: bfactor, previous_temperature_cell, previous_pressure_cell

      count = cell_particles(i,j,k)%count
      if (count == 0) then
         call lcm_save_cell_mean_history(i, j, k)
         return
      endif

      temperature = cell_tabs_lcm(i, j, k)
      pressure = field_state%pressure(k)
      air_density = field_state%density(k)
      cell_volume = real(grid_state%dV(k), real64)
      call lcm_get_cell_mean_history(i, j, k, previous_temperature_cell, &
                                     previous_pressure_cell)
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
      d_pre = (latent_heat_condensation / (rv * temperature) - 1.0_real64) / &
              (thermal_conductivity * temperature)

      if (config_state%switch_koehler_micro) then
         alpha_factor = 2.0_real64 * surface_tension_water / &
                        (rho_liq * rv * temperature)
      else
         alpha_factor = 0.0_real64
      endif

      ql_local = 0.0_real64
      radiation_integral = 0.0_real64
      do n = 1, count
         call require_liquid_particle(cell_particles(i,j,k)%particles(n))
         radius = liquid_radius(cell_particles(i,j,k)%particles(n))
         ql_local = ql_local + cell_particles(i,j,k)%particles(n)%mass_water / &
                    (cell_volume * air_density)
         if (config_state%switch_radiation_micro) then
            q_absorb = 1.18_real64 * (1.0_real64 - exp(-0.28e6_real64 * radius))
            radiation_integral = radiation_integral + &
               4.0_real64 * pi * radius**2 * q_absorb * &
               real(cell_particles(i,j,k)%particles(n)%weight, real64) / cell_volume
         endif
      enddo

      radiation = 0.0_real64
      if (config_state%switch_radiation_micro .and. ql_local >= 1.0e-5_real64 .and. &
          radiation_integral > 0.0_real64) then
         ! SAM sets k_liq/k_tot to 1.0 in this routine.
         qrad_liq = field_state%qrad(i,j,k)
         radiation = qrad_liq * specific_heat_dry_air * air_density / radiation_integral
      endif

      dsd_liq = 0.0_real64
      if (.not. config_state%switch_lem_micro) then
         do n = 1, count
            radius = liquid_radius(cell_particles(i,j,k)%particles(n))
            r0 = growth_length_scale(temperature, esatw_pre, &
                 thermal_conductivity, vapor_diffusivity, &
                 cell_particles(i,j,k)%particles(n), alpha_mordy)
            ventilation = liquid_ventilation(radius, &
               cell_particles(i,j,k)%particles(n)%w_sedi, air_density)
            dsd_liq = dsd_liq + &
               real(cell_particles(i,j,k)%particles(n)%weight, real64) * &
               radius**2 / (radius + r0) * ventilation
         enddo

         if (dsd_liq > 0.0_real64) then
            pre_liq = a_liq_pre * 4.0_real64 * pi * c_liq_pre / &
                      (qsatw_pre * cell_volume)
            tau = 1.0_real64 / (pre_liq * dsd_liq)
            eta_dyn = field_state%qv(i,j,k) - &
                      qsatw_liquid(temperature, pressure)
            eta_last = field_state%qv(i,j,k) - &
                       qsatw_liquid(previous_temperature_cell, previous_pressure_cell)
            forcing = (eta_dyn - eta_last) / time_step
            eta_end = eta_new(eta_last, forcing, tau, time_step)
            eta_integral = eta_int(eta_last, forcing, tau, time_step)
            fsw = supersaturation_weight(eta_last, eta_end, eta_integral, time_step)
            eta_mean = (eta_integral * (1.0_real64 - fsw) + &
               fsw * (eta_last + eta_end) * time_step * 0.5_real64) / time_step
            supersaturation = eta_mean / qsatw_pre
         else
            supersaturation = 0.0_real64
         endif

         if (current_step <= config_state%nstep_cut_supersat .and. &
             supersaturation > real(config_state%supersat_threshold, real64)) then
            supersaturation = real(config_state%supersat_threshold, real64)
         endif

         do n = 1, count
            call integrate_liquid_particle(i, j, k, n, time_step, temperature, &
               pressure, air_density, cell_volume, c_liq_pre, d_pre, &
               alpha_factor, supersaturation, radiation, esatw_pre, &
               thermal_conductivity, vapor_diffusivity, cell_particles(i,j,k))
         enddo
      else
         do n = 1, count
            radius = liquid_radius(cell_particles(i,j,k)%particles(n))
            aerosol_radius = particle_aerosol_radius(i, j, k, n)
            particle_temperature = cell_particles(i,j,k)%particles(n)%temperature - &
               cell_particles(i,j,k)%particles(n)%z * gravity / specific_heat_dry_air
            r0 = growth_length_scale(temperature, esatw_pre, thermal_conductivity, &
               vapor_diffusivity, cell_particles(i,j,k)%particles(n), alpha_mordy)
            if (config_state%switch_koehler_micro) then
               bfactor = particle_bfactor(i, j, k, n)
            else
               bfactor = 0.0_real64
            endif
            ventilation = liquid_ventilation(radius, &
               cell_particles(i,j,k)%particles(n)%w_sedi, air_density)
            pre_liq = a_liq_pre * 4.0_real64 * pi * c_liq_pre / &
                      (qsatw_pre * (cell_volume / real(count, real64)))
            tau = 1.0_real64 / (pre_liq * &
               real(cell_particles(i,j,k)%particles(n)%weight, real64) * &
               radius**2 / (radius + r0) * ventilation)
            eta_dyn = cell_particles(i,j,k)%particles(n)%qv - &
                      qsatw_liquid(particle_temperature, pressure)
            eta_last = cell_particles(i,j,k)%particles(n)%eta
            forcing = (eta_dyn - eta_last) / time_step
            eta_end = eta_new(eta_last, forcing, tau, time_step)
            eta_integral = eta_int(eta_last, forcing, tau, time_step)
            fsw = supersaturation_weight(eta_last, eta_end, eta_integral, time_step)
            eta_mean = (eta_integral * (1.0_real64 - fsw) + &
               fsw * (eta_last + eta_end) * time_step * 0.5_real64) / time_step
            supersaturation = eta_mean / &
                              qsatw_liquid(particle_temperature, pressure)

            r0_local = growth_length_scale(temperature, esatw_pre, &
               thermal_conductivity, &
               vapor_diffusivity, cell_particles(i,j,k)%particles(n), alpha_mordy)
            call integrate_liquid_particle(i, j, k, n, time_step, temperature, &
               pressure, air_density, cell_volume, c_liq_pre, d_pre, &
               alpha_factor, supersaturation, radiation, esatw_pre, &
               thermal_conductivity, vapor_diffusivity, cell_particles(i,j,k), &
               local_r0=r0_local, local_bfactor=bfactor, &
               local_ventilation=ventilation)
         enddo

         ! SAM saves post-condensation LEM supersaturation for the next step.
         do n = 1, count
            particle_temperature = cell_particles(i,j,k)%particles(n)%temperature - &
               cell_particles(i,j,k)%particles(n)%z * gravity / specific_heat_dry_air
            cell_particles(i,j,k)%particles(n)%eta = cell_particles(i,j,k)%particles(n)%qv - &
               qsatw_liquid(particle_temperature, pressure)
         enddo
      endif

      call lcm_save_cell_mean_history(i, j, k)
   end subroutine micro_condensation_cell

   subroutine integrate_liquid_particle(i, j, k, n, time_step, temperature, &
      pressure, air_density, cell_volume, c_liq, d_pre, alpha_factor, &
      supersaturation, radiation, esatw_pre, thermal_conductivity, &
      vapor_diffusivity, particle_cell, local_r0, local_bfactor, &
      local_ventilation)
      integer, intent(in) :: i, j, k, n
      real(real64), intent(in) :: time_step, temperature, pressure
      real(real64), intent(in) :: air_density, cell_volume, c_liq, d_pre
      real(real64), intent(in) :: alpha_factor, supersaturation, radiation
      real(real64), intent(in) :: esatw_pre, thermal_conductivity
      real(real64), intent(in) :: vapor_diffusivity
      type(lcm_particle_cell_t), intent(inout) :: particle_cell
      real(real64), intent(in), optional :: local_r0, local_bfactor
      real(real64), intent(in), optional :: local_ventilation

      type(lcm_particle_t) :: particle
      real(real64) :: radius, new_radius, aerosol_radius, r0, bfactor
      real(real64) :: ventilation, q_absorb, dm, dmass_cell
      real(real64) :: particle_temperature
      integer :: count

      particle = particle_cell%particles(n)
      count = particle_cell%count
      radius = liquid_radius(particle)
      aerosol_radius = particle_aerosol_radius(i, j, k, n)
      r0 = growth_length_scale(temperature, esatw_pre, &
         thermal_conductivity, vapor_diffusivity, particle, alpha_mordy)
      bfactor = 0.0_real64
      if (config_state%switch_koehler_micro) bfactor = &
         particle_bfactor(i, j, k, n)
      ventilation = liquid_ventilation(radius, particle%w_sedi, air_density)
      if (present(local_r0)) r0 = local_r0
      if (present(local_bfactor)) bfactor = local_bfactor
      if (present(local_ventilation)) ventilation = local_ventilation

      q_absorb = 0.0_real64
      if (config_state%switch_radiation_micro) then
         q_absorb = 1.18_real64 * (1.0_real64 - exp(-0.28e6_real64 * radius))
      endif

      if (config_state%switch_analytical_micro) then
         new_radius = sqrt(max((radius + r0)**2 + &
            2.0_real64 * time_step * c_liq / rho_liq * &
            supersaturation * ventilation, (aerosol_radius + r0)**2)) - r0
      else
         new_radius = radius
         call radius_liquid_euler(new_radius, time_step, r0, c_liq, &
            supersaturation, ventilation, alpha_factor, bfactor, aerosol_radius, &
            d_pre, q_absorb * radiation)
      endif

      particle%mass_water = (4.0_real64 / 3.0_real64) * pi * rho_liq * &
         new_radius**3 * real(particle%weight, real64)
      dm = (4.0_real64 / 3.0_real64) * pi * rho_liq * &
           (new_radius**3 - radius**3)
      dmass_cell = dm * real(particle%weight, real64)

      if (config_state%switch_lem_micro) then
         dmass_cell = dmass_cell / &
            ((cell_volume / real(count, real64)) * air_density)
         particle%temperature = particle%temperature + dmass_cell * fac_cond
         particle%qv = particle%qv - dmass_cell
      endif

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
      ventilation, alpha_factor, beta_factor, aerosol_radius, d_pre, radiation)
      real(real64), intent(inout) :: radius
      real(real64), intent(in) :: time_step, r0, c_pre, supersaturation
      real(real64), intent(in) :: ventilation, alpha_factor, beta_factor
      real(real64), intent(in) :: aerosol_radius, d_pre, radiation

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
                beta_factor / r_euler**3 - d_pre * radiation * r_euler) * &
               r_euler / (r_euler + r0)
            d2r2dtdr2 = g_pre * ventilation * &
               (alpha_factor * r_euler**3 - &
                beta_factor * (3.0_real64 * r_euler + 2.0_real64 * r0) - &
                r_euler**3 * (d_pre * radiation * r_euler * &
                   (r_euler + 2.0_real64 * r0) - r0 * supersaturation)) / &
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

   real(real64) function liquid_radius(particle) result(radius)
      type(lcm_particle_t), intent(in) :: particle

      radius = (particle%mass_water / &
         ((4.0_real64 / 3.0_real64) * pi * &
          real(particle%weight, real64) * rho_liq))**(1.0_real64 / 3.0_real64)
   end function liquid_radius

   real(real64) function growth_length_scale(temperature, esatw, &
      thermal_conductivity, vapor_diffusivity, particle, alpha) result(r0)
      real(real64), intent(in) :: temperature, esatw
      real(real64), intent(in) :: thermal_conductivity, vapor_diffusivity, alpha
      type(lcm_particle_t), intent(in) :: particle
      real(real64) :: accommodation

      accommodation = alpha
      if (config_state%switch_ffc_micro) then
         accommodation = ffc_accommodation(particle, temperature)
      endif
      r0 = vapor_diffusivity / accommodation * &
         sqrt(2.0_real64 * pi / (rv * temperature)) / &
         (1.0_real64 + vapor_diffusivity * &
          latent_heat_condensation**2 * 100.0_real64 * esatw / &
          (thermal_conductivity * rv**2 * temperature**3))
   end function growth_length_scale

   real(real64) function ffc_accommodation(particle, temperature) result(alpha)
      type(lcm_particle_t), intent(in) :: particle
      real(real64), intent(in) :: temperature
      real(real64) :: droplet_area, n_ffc, radius

      radius = liquid_radius(particle)
      droplet_area = 4.0_real64 * pi * radius**2
      n_ffc = f_ffc * particle%mass_aero / &
         real(particle%weight, real64) / molecular_weight_ffc * n_avogadro
      if (n_ffc <= 0.0_real64) then
         error stop 'LCM condensation: SAM FFC path requires positive aerosol mass'
      endif
      if (droplet_area / n_ffc < 22.0e-20_real64) then
         alpha = 1.0e-5_real64
      else
         alpha = alpha_mordy
      endif
      ! Keep the temperature argument explicit because SAM evaluates this
      ! accommodation length at the cell/particle temperature for each path.
      if (temperature <= 0.0_real64) error stop 'LCM condensation: invalid temperature'
   end function ffc_accommodation

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
