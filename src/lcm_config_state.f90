!> LCM configuration supplied through the standalone host interface.
module lcm_config_state
   use lcm_types, only: lcm_config_t
   use micro_vars, only: np_xy, nzl, nstep_LCM_init, nstep_cut_supersat, &
                         supersat_threshold, n_uvw_iteration,           &
                         switch_sgs_velocities, switch_sedimentation,   &
                         start_sedimentation, switch_steady_aerosol,    &
                         inject_bottom, inject_top, switch_collection,  &
                         start_collection,                              &
                         switch_linear_sampling_collection,             &
                         switch_breakup, breakup_conserve_mode,         &
                         switch_analytical_micro,                       &
                         switch_radiation_micro, switch_LEM_micro,      &
                         switch_koehler_micro, switch_FFC_micro,        &
                         switch_steady_pichamber,                       &
                         switch_reflective_walls_pichamber,             &
                         switch_walls_pichamber, switch_cloud_tracking, &
                         switch_phase_sources, switch_particle_source,  &
                         switch_moving_source, n_aero, rm_aero,         &
                         sigma_aero, shares_aero, n_ice, rm_in,         &
                         sigma_in, rho_aero, molecular_weight_aero,     &
                         vanthoff_aero, init_aero, unit_aero,           &
                         r_separate_micro, nz_b_source_aero,            &
                         nz_t_source_aero, ny_s_source_aero,            &
                         ny_n_source_aero, nx_w_source_aero,            &
                         nx_e_source_aero, nz_b_source_in,              &
                         nz_t_source_in, ny_s_source_in, ny_n_source_in,&
                         nx_w_source_in, nx_e_source_in,                &
                         n_start_source_aero, n_end_source_aero,        &
                         n_source_aero, n_start_source_in,              &
                         n_end_source_in, n_source_in, r_source_aero,   &
                         r_source_in, sigma_source_aero,                &
                         sigma_source_in, type_source_in,               &
                         p_source_aero, p_source_in,                    &
                         frac_in_source_particles,                      &
                         frac_in_source_weight, unit_source,            &
                         type_source, salinity_source, vx_source_aero,  &
                         vy_source_aero, vz_source_aero, vx_source_in,  &
                         vy_source_in, vz_source_in, alpha_spec,        &
                         r_start_spec, r_end_spec, spec_nsave,          &
                         spec_nsavestart, switch_spec_aero,             &
                         horiz_average_spec, track_some_particles,      &
                         n_track_particles, particle_ntrack,            &
                         particle_ntrackstart, particle_ntrackstop,     &
                         frac_ice, switch_ice, switch_habit_ice,        &
                         habitparam, start_ice, switch_icebreak,        &
                         switch_rs, switch_ds, ds_param, icebreak_param,&
                         switch_diskice_sizethresh, sizethresh_chen,    &
                         sizethresh_diskice

   implicit none

   private

   public :: lcm_set_config
   public :: config_is_set
   public :: config_state

   logical, protected, save :: config_is_set = .false.
   type(lcm_config_t), protected, save :: config_state

contains

   subroutine lcm_set_config(config)
      type(lcm_config_t), intent(in) :: config

      call validate_config(config)
      config_state = config

      np_xy = config%np_xy
      nzl = config%nzl
      nstep_LCM_init = config%nstep_lcm_init
      nstep_cut_supersat = config%nstep_cut_supersat
      supersat_threshold = config%supersat_threshold
      n_uvw_iteration = config%n_uvw_iteration
      switch_sgs_velocities = config%switch_sgs_velocities
      switch_sedimentation = config%switch_sedimentation
      start_sedimentation = config%start_sedimentation
      switch_steady_aerosol = config%switch_steady_aerosol
      inject_bottom = config%inject_bottom
      inject_top = config%inject_top
      switch_collection = config%switch_collection
      start_collection = config%start_collection
      switch_linear_sampling_collection = config%switch_linear_sampling_collection
      switch_breakup = config%switch_breakup
      breakup_conserve_mode = config%breakup_conserve_mode
      switch_analytical_micro = config%switch_analytical_micro
      switch_radiation_micro = config%switch_radiation_micro
      switch_LEM_micro = config%switch_lem_micro
      switch_koehler_micro = config%switch_koehler_micro
      switch_FFC_micro = config%switch_ffc_micro
      switch_steady_pichamber = config%switch_steady_pichamber
      switch_reflective_walls_pichamber = config%switch_reflective_walls_pichamber
      switch_walls_pichamber = config%switch_walls_pichamber
      switch_cloud_tracking = config%switch_cloud_tracking
      switch_phase_sources = config%switch_phase_sources
      switch_particle_source = config%switch_particle_source
      switch_moving_source = config%switch_moving_source
      switch_ice = config%switch_ice
      switch_habit_ice = config%switch_habit_ice
      habitparam = config%habitparam
      start_ice = config%start_ice
      switch_icebreak = config%switch_icebreak
      switch_rs = config%switch_rs
      switch_ds = config%switch_ds
      ds_param = config%ds_param
      icebreak_param = config%icebreak_param
      switch_diskice_sizethresh = config%switch_diskice_sizethresh
      sizethresh_chen = config%sizethresh_chen
      sizethresh_diskice = config%sizethresh_diskice
      n_aero = config%n_aero
      rm_aero = config%rm_aero
      sigma_aero = config%sigma_aero
      shares_aero = config%shares_aero
      n_ice = config%n_ice
      rm_in = config%rm_in
      sigma_in = config%sigma_in
      rho_aero = config%rho_aero
      molecular_weight_aero = config%molecular_weight_aero
      vanthoff_aero = config%vanthoff_aero
      init_aero = config%init_aero
      unit_aero = config%unit_aero
      r_separate_micro = config%r_separate_micro
      nz_b_source_aero = config%nz_b_source_aero
      nz_t_source_aero = config%nz_t_source_aero
      ny_s_source_aero = config%ny_s_source_aero
      ny_n_source_aero = config%ny_n_source_aero
      nx_w_source_aero = config%nx_w_source_aero
      nx_e_source_aero = config%nx_e_source_aero
      nz_b_source_in = config%nz_b_source_in
      nz_t_source_in = config%nz_t_source_in
      ny_s_source_in = config%ny_s_source_in
      ny_n_source_in = config%ny_n_source_in
      nx_w_source_in = config%nx_w_source_in
      nx_e_source_in = config%nx_e_source_in
      n_start_source_aero = config%n_start_source_aero
      n_end_source_aero = config%n_end_source_aero
      n_source_aero = config%n_source_aero
      n_start_source_in = config%n_start_source_in
      n_end_source_in = config%n_end_source_in
      n_source_in = config%n_source_in
      r_source_aero = config%r_source_aero
      r_source_in = config%r_source_in
      sigma_source_aero = config%sigma_source_aero
      sigma_source_in = config%sigma_source_in
      type_source_in = config%type_source_in
      p_source_aero = config%p_source_aero
      p_source_in = config%p_source_in
      frac_in_source_particles = config%frac_in_source_particles
      frac_in_source_weight = config%frac_in_source_weight
      unit_source = config%unit_source
      type_source = config%type_source
      salinity_source = config%salinity_source
      vx_source_aero = config%vx_source_aero
      vy_source_aero = config%vy_source_aero
      vz_source_aero = config%vz_source_aero
      vx_source_in = config%vx_source_in
      vy_source_in = config%vy_source_in
      vz_source_in = config%vz_source_in
      alpha_spec = config%alpha_spec
      r_start_spec = config%r_start_spec
      r_end_spec = config%r_end_spec
      spec_nsave = config%spec_nsave
      spec_nsavestart = config%spec_nsavestart
      switch_spec_aero = config%switch_spec_aero
      horiz_average_spec = config%horiz_average_spec
      track_some_particles = config%track_some_particles
      n_track_particles = config%n_track_particles
      particle_ntrack = config%particle_ntrack
      particle_ntrackstart = config%particle_ntrackstart
      particle_ntrackstop = config%particle_ntrackstop
      config_is_set = .true.
   end subroutine lcm_set_config

   subroutine validate_config(config)
      type(lcm_config_t), intent(in) :: config

      if (config%np_xy <= 0) then
         error stop 'LCM config np_xy must be positive'
      endif

      if (config%nzl <= 0) then
         error stop 'LCM config nzl must be positive'
      endif

      if (config%nstep_lcm_init <= 0) then
         error stop 'LCM config nstep_lcm_init must be positive'
      endif

      if (config%n_uvw_iteration < 0) then
         error stop 'LCM config n_uvw_iteration must not be negative'
      endif

      if (config%bottom_removal_radius <= 0.0) then
         error stop 'LCM config bottom_removal_radius must be positive'
      endif

      if (trim(config%init_aero) /= 'bin' .and. &
          trim(config%init_aero) /= 'random' .and. &
          trim(config%init_aero) /= 'USM') then
         error stop 'LCM config init_aero must be bin, random, or USM'
      endif

      if (trim(config%unit_aero) /= '1/kg' .and. &
          trim(config%unit_aero) /= '1/m3') then
         error stop 'LCM config unit_aero must be 1/kg or 1/m3'
      endif
   end subroutine validate_config

end module lcm_config_state
