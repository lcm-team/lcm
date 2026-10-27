program test_condensation_growth_options
   use iso_fortran_env, only: int64, real64
   use lcm_condensation, only: lcm_apply_condensation
   use lcm_config_state, only: config_state, lcm_set_config
   use lcm_parameters, only: pi, rho_liq, gravity, specific_heat_dry_air
   use lcm_particle_state, only: cell_particles, lcm_add_particle, &
                                 lcm_reset_particle_cells, qsatw_liquid
   use lcm_types, only: lcm_particle_t, liquid_no_in
   use minimal_les_host, only: minimal_host_config, minimal_host_static_energy, &
                               minimal_host_qv, minimal_host_pressure, &
                               minimal_host_z_center, setup_minimal_les_host

   implicit none

   real(real64), parameter :: initial_radius = 0.5e-6_real64
   real(real64), parameter :: aerosol_radius = 0.05e-6_real64
   real(real64), parameter :: test_temperature = 280.0_real64
   real(real64), parameter :: relative_tolerance = 1.0e-10_real64
   real(real64) :: baseline_mass, no_koehler_mass
   real(real64) :: variable_thermal_mass, variable_diffusivity_mass
   real(real64) :: higher_supersaturation_mass, evaporation_mass
   real(real64) :: evaporation_radius, lem_growth_mass
   real(real64) :: initial_lem_qv, initial_lem_temperature

   call setup_minimal_les_host()
   call prepare_supersaturated_cell()

   minimal_host_config%switch_lem_micro = .false.
   minimal_host_config%switch_analytical_micro = .false.
   minimal_host_config%switch_koehler_micro = .true.
   minimal_host_config%switch_variable_therco = .false.
   minimal_host_config%switch_variable_diffelq = .false.
   call lcm_set_config(minimal_host_config)

   if (.not. config_state%switch_koehler_micro) then
      error stop 'Kohler effects must be enabled by default'
   endif

   baseline_mass = run_growth_case(0.04_real64)
   call assert_true(baseline_mass > initial_particle_mass(), &
                    'supersaturated liquid droplet should grow')
   call assert_true(cell_particles(1,1,1)%particles(1)%qv == &
                    minimal_host_qv(1,1,1), &
                    'non-LEM growth should not change particle-local vapor')
   call assert_true(cell_particles(1,1,1)%particles(1)%temperature == &
                    test_temperature, &
                    'non-LEM growth should not change particle-local temperature')

   minimal_host_config%switch_koehler_micro = .false.
   call lcm_set_config(minimal_host_config)
   no_koehler_mass = run_growth_case(0.04_real64)
   call assert_greater(no_koehler_mass, baseline_mass, &
                       'Kohler curvature and solute should restrain growth here')

   minimal_host_config%switch_koehler_micro = .true.
   call lcm_set_config(minimal_host_config)
   higher_supersaturation_mass = run_growth_case(0.08_real64)
   call assert_greater(higher_supersaturation_mass, baseline_mass, &
                       'higher supersaturation should produce more growth')

   minimal_host_config%switch_variable_therco = .true.
   call lcm_set_config(minimal_host_config)
   variable_thermal_mass = run_growth_case(0.04_real64)
   call assert_greater(variable_thermal_mass, baseline_mass, &
                       'higher variable thermal conductivity should increase growth')

   minimal_host_config%switch_variable_therco = .false.
   minimal_host_config%switch_variable_diffelq = .true.
   call lcm_set_config(minimal_host_config)
   variable_diffusivity_mass = run_growth_case(0.04_real64)
   call assert_greater(variable_diffusivity_mass, baseline_mass, &
                       'higher variable vapor diffusivity should increase growth')

   minimal_host_config%switch_variable_diffelq = .false.
   call lcm_set_config(minimal_host_config)
   evaporation_mass = run_growth_case(-0.20_real64)
   evaporation_radius = radius_from_mass(evaporation_mass)
   call assert_true(evaporation_mass < initial_particle_mass(), &
                    'subsaturated liquid droplet should evaporate')
   call assert_true(evaporation_radius >= aerosol_radius * (1.0_real64 - 1.0e-12_real64), &
                    'evaporation should not shrink below dry aerosol radius')

   minimal_host_config%switch_lem_micro = .true.
   call lcm_set_config(minimal_host_config)
   initial_lem_qv = 1.04_real64 * &
      qsatw_liquid(test_temperature, minimal_host_pressure(1))
   initial_lem_temperature = test_temperature
   lem_growth_mass = run_growth_case(0.04_real64, 1000000_int64)
   call assert_true(lem_growth_mass > initial_particle_mass(1000000_int64), &
                    'LEM supersaturated droplet should grow')
   call assert_true(cell_particles(1,1,1)%particles(1)%qv < initial_lem_qv, &
                    'LEM condensation should reduce particle-local vapor')
   call assert_true(cell_particles(1,1,1)%particles(1)%temperature > &
                    initial_lem_temperature, &
                    'LEM condensation should warm the particle')
   call assert_true(abs(cell_particles(1,1,1)%particles(1)%eta - &
                    (cell_particles(1,1,1)%particles(1)%qv - &
                     qsatw_liquid(cell_particles(1,1,1)%particles(1)%temperature, &
                                  minimal_host_pressure(1)))) < 1.0e-14_real64, &
                    'LEM eta should reflect post-growth vapor and temperature')

contains

   subroutine prepare_supersaturated_cell()
      real(real64) :: gamaz

      gamaz = gravity / specific_heat_dry_air * minimal_host_z_center(1)
      minimal_host_static_energy(1,1,1) = test_temperature + gamaz
   end subroutine prepare_supersaturated_cell

   real(real64) function run_growth_case(supersaturation_ratio, multiplicity) &
      result(mass_after)
      real(real64), intent(in) :: supersaturation_ratio
      integer(int64), intent(in), optional :: multiplicity
      type(lcm_particle_t) :: particle

      call lcm_reset_particle_cells()
      minimal_host_qv(1,1,1) = (1.0_real64 + supersaturation_ratio) * &
         qsatw_liquid(test_temperature, minimal_host_pressure(1))
      particle%weight = 1_int64
      if (present(multiplicity)) particle%weight = multiplicity
      particle%particle_type = liquid_no_in
      particle%mass_water = (4.0_real64 / 3.0_real64) * pi * rho_liq * &
                            initial_radius**3 * real(particle%weight, real64)
      particle%mass_aero = minimal_host_config%rho_aero(1) * &
         (4.0_real64 / 3.0_real64) * pi * aerosol_radius**3 * &
         real(particle%weight, real64)
      particle%aerosol_portions(1) = 1.0_real64
      particle%temperature = test_temperature
      particle%qv = minimal_host_qv(1,1,1)
      particle%eta = particle%qv - &
         qsatw_liquid(test_temperature, minimal_host_pressure(1))
      particle%w_sedi = 0.0_real64
      call lcm_add_particle(1, 1, 1, particle)

      call lcm_apply_condensation(1.0_real64, minimal_host_config%nstep_lcm_init + 1)
      mass_after = cell_particles(1,1,1)%particles(1)%mass_water
   end function run_growth_case

   real(real64) function initial_particle_mass(multiplicity) result(mass)
      integer(int64), intent(in), optional :: multiplicity
      real(real64) :: weight

      weight = 1.0_real64
      if (present(multiplicity)) weight = real(multiplicity, real64)
      mass = (4.0_real64 / 3.0_real64) * pi * rho_liq * initial_radius**3
      mass = mass * weight
   end function initial_particle_mass

   real(real64) function radius_from_mass(mass) result(radius)
      real(real64), intent(in) :: mass

      radius = (mass / ((4.0_real64 / 3.0_real64) * pi * rho_liq))** &
               (1.0_real64 / 3.0_real64)
   end function radius_from_mass

   subroutine assert_greater(value, reference, message)
      real(real64), intent(in) :: value, reference
      character(len=*), intent(in) :: message

      call assert_true(value > reference + &
         relative_tolerance * max(abs(reference), initial_particle_mass()), message)
   end subroutine assert_greater

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) error stop message
   end subroutine assert_true

end program test_condensation_growth_options
