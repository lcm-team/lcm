!> LCM configuration supplied through the standalone host interface.
module lcm_config_state
   use lcm_types, only: lcm_config_t
   use micro_vars, only: np_xy, nzl, nstep_LCM_init, n_uvw_iteration, &
                         switch_sgs_velocities, switch_sedimentation, &
                         switch_collection, switch_analytical_micro, &
                         switch_radiation_micro, switch_LEM_micro, &
                         switch_koehler_micro, switch_FFC_micro, &
                         switch_ice, n_aero, rm_aero, sigma_aero, &
                         shares_aero, rho_aero, molecular_weight_aero, &
                         vanthoff_aero, init_aero, unit_aero

   implicit none

   private

   public :: lcm_set_config
   public :: config_is_set

   logical, protected, save :: config_is_set = .false.

contains

   subroutine lcm_set_config(config)
      type(lcm_config_t), intent(in) :: config

      call validate_config(config)

      np_xy = config%np_xy
      nzl = config%nzl
      nstep_LCM_init = config%nstep_lcm_init
      n_uvw_iteration = config%n_uvw_iteration
      switch_sgs_velocities = config%switch_sgs_velocities
      switch_sedimentation = config%switch_sedimentation
      switch_collection = config%switch_collection
      switch_analytical_micro = config%switch_analytical_micro
      switch_radiation_micro = config%switch_radiation_micro
      switch_LEM_micro = config%switch_lem_micro
      switch_koehler_micro = config%switch_koehler_micro
      switch_FFC_micro = config%switch_ffc_micro
      switch_ice = config%switch_ice
      n_aero = config%n_aero
      rm_aero = config%rm_aero
      sigma_aero = config%sigma_aero
      shares_aero = config%shares_aero
      rho_aero = config%rho_aero
      molecular_weight_aero = config%molecular_weight_aero
      vanthoff_aero = config%vanthoff_aero
      init_aero = config%init_aero
      unit_aero = config%unit_aero
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
