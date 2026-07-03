!> LCM configuration supplied through the standalone host interface.
module lcm_config_state
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

   public :: lcm_config_t
   public :: lcm_init_config
   public :: lcm_set_config

   type :: lcm_config_t
      integer :: np_xy = 50
      integer :: nzl = 50
      integer :: nstep_lcm_init = 3
      integer :: n_uvw_iteration = 1
      logical :: switch_sgs_velocities = .true.
      logical :: switch_sedimentation = .true.
      logical :: switch_collection = .true.
      logical :: switch_analytical_micro = .false.
      logical :: switch_radiation_micro = .false.
      logical :: switch_lem_micro = .true.
      logical :: switch_koehler_micro = .true.
      logical :: switch_ffc_micro = .false.
      logical :: switch_ice = .false.
      real :: n_aero(4) = [100.0e6, 0.0, 0.0, 0.0]
      real :: rm_aero(4) = [0.04e-6, 0.0, 0.0, 0.0]
      real :: sigma_aero(4) = [1.4, 0.0, 0.0, 0.0]
      real :: shares_aero(4) = [1.0, 1.0, 1.0, 1.0]
      real :: rho_aero(4) = [2170.0, 0.0, 0.0, 0.0]
      real :: molecular_weight_aero(4) = [58.443, 0.0, 0.0, 0.0]
      real :: vanthoff_aero(4) = [2.0, 0.0, 0.0, 0.0]
      character(len=6) :: init_aero = 'bin'
      character(len=12) :: unit_aero = '1/kg'
   end type lcm_config_t

contains

   subroutine lcm_init_config(config)
      type(lcm_config_t), intent(out) :: config

      config = lcm_config_t()
   end subroutine lcm_init_config

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
