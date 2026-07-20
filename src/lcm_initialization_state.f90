!> Standalone LCM initialization lifecycle.
!!
!! Hosts call lcm_initialize after config, grid, and atmospheric fields have
!! been supplied. The timestep path may then assume LCM-owned storage exists.
module lcm_initialization_state
   use iso_fortran_env, only: int64
   use lcm_config_state, only: config_is_set
   use lcm_field_state, only: field_state
   use lcm_grid_state, only: grid_state
   use lcm_particle_state, only: lcm_allocate_particle_cells
   use micro_vars, only: frac_ice, n_source_aero, np_xy, nzl, &
                         p_source_aero, p_source_in,         &
                         switch_analytical_micro, switch_ice,&
                         switch_LEM_micro, switch_particle_source, &
                         switch_radiation_micro, frac_in_source_particles

   implicit none

   private

   public :: lcm_initialize
   public :: lcm_is_initialized
   public :: lcm_initial_particle_count_local
   public :: lcm_initial_particle_count_global
   public :: lcm_initialization_mode

   logical, protected, save :: lcm_is_initialized = .false.
   integer(int64), protected, save :: lcm_initial_particle_count_local = 0_int64
   integer(int64), protected, save :: lcm_initial_particle_count_global = 0_int64
   character(len=16), protected, save :: lcm_initialization_mode = "not_initialized"

contains

   subroutine lcm_initialize(from_restart)
      logical, intent(in), optional :: from_restart

      if (lcm_is_initialized) then
         error stop "LCM has already been initialized"
      endif

      call validate_initialization_mode(from_restart)
      call validate_initialization_inputs()
      call initialize_particle_count_state()
      call lcm_allocate_particle_cells(initial_capacity_per_cell=np_xy)

      ! SAM-LCM micro_set also initialized grid-helper geometry. In standalone
      ! LCM that responsibility belongs to lcm_set_grid.
      lcm_initialization_mode = "fresh"
      lcm_is_initialized = .true.
   end subroutine lcm_initialize

   subroutine validate_initialization_mode(from_restart)
      logical, intent(in), optional :: from_restart

      if (present(from_restart)) then
         if (from_restart) then
            error stop "LCM restart initialization is not implemented"
         endif
      endif
   end subroutine validate_initialization_mode

   subroutine validate_initialization_inputs()
      if (.not. config_is_set) then
         error stop "LCM config must be set before lcm_initialize"
      endif

      if (.not. grid_state%is_set) then
         error stop "LCM grid must be set before lcm_initialize"
      endif

      if (nzl > grid_state%grid%nz) then
         error stop "LCM nzl must not exceed the vertical grid size"
      endif

      if (switch_radiation_micro .and. switch_analytical_micro) then
         error stop "LCM radiation microphysics requires analytical microphysics to be disabled"
      endif

      if (switch_LEM_micro .and. has_stretched_vertical_grid()) then
         error stop "LCM LEM SGS model requires an unstretched vertical grid"
      endif

      if (switch_ice .and. switch_LEM_micro) then
         error stop "LCM ice microphysics is not compatible with the LEM SGS model"
      endif

      if (switch_ice .and. switch_radiation_micro) then
         error stop "LCM radiation effects on ice are not implemented"
      endif

      if (switch_ice .and. (frac_ice < 0.0 .or. frac_ice > 1.0)) then
         error stop "LCM frac_ice must be between 0 and 1"
      endif

      if ((.not. switch_ice) .and. switch_particle_source .and. &
          (maxval(p_source_in) > 0.0 .or. &
           maxval(p_source_aero * n_source_aero * frac_in_source_particles) /= 0.0)) then
         error stop "LCM ice particle source is present but ice microphysics is disabled"
      endif

      call require_attached_field(field_state%u_is_attached, "u")
      call require_attached_field(field_state%v_is_attached, "v")
      call require_attached_field(field_state%w_is_attached, "w")
      call require_attached_field(field_state%static_energy_is_attached, &
                                  "static_energy")
      call require_attached_field(field_state%temperature_is_attached, &
                                  "temperature")
      call require_attached_field(field_state%qv_is_attached, "qv")
      call require_attached_field(field_state%pressure_is_attached, "pressure")
      call require_attached_field(field_state%density_is_attached, "density")
   end subroutine validate_initialization_inputs

   subroutine require_attached_field(is_attached, field_name)
      logical, intent(in) :: is_attached
      character(len=*), intent(in) :: field_name

      if (.not. is_attached) then
         error stop "LCM required field " // field_name // &
                    " must be attached before lcm_initialize"
      endif
   end subroutine require_attached_field

   logical function has_stretched_vertical_grid()
      real :: dz_reference
      real :: tolerance

      dz_reference = real(grid_state%grid%dz(1))
      tolerance = 100.0 * epsilon(dz_reference) * max(1.0, abs(dz_reference))

      has_stretched_vertical_grid = any(abs(real(grid_state%grid%dz) - &
                                           dz_reference) > tolerance)
   end function has_stretched_vertical_grid

   subroutine initialize_particle_count_state()
      integer :: initialization_levels

      initialization_levels = min(nzl, grid_state%grid%nz)

      lcm_initial_particle_count_local = &
         int(grid_state%grid%nx_local, int64) * &
         int(grid_state%grid%ny_local, int64) * &
         int(initialization_levels, int64) * &
         int(np_xy, int64)

      lcm_initial_particle_count_global = &
         int(grid_state%grid%nx_global, int64) * &
         int(grid_state%grid%ny_global, int64) * &
         int(initialization_levels, int64) * &
         int(np_xy, int64)
   end subroutine initialize_particle_count_state

end module lcm_initialization_state
