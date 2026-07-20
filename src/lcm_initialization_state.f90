!> Standalone LCM initialization lifecycle.
!!
!! Hosts call lcm_initialize after config, grid, and atmospheric fields have
!! been supplied. The timestep path may then assume LCM-owned storage exists.
module lcm_initialization_state
   use iso_fortran_env, only: int64
   use lcm_config_state, only: config_is_set, config_state
   use lcm_field_state, only: field_state
   use lcm_grid_state, only: grid_state
   use lcm_particle_state, only: lcm_allocate_particle_cells

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

   subroutine lcm_initialize()
      if (lcm_is_initialized) then
         error stop "LCM has already been initialized"
      endif

      call validate_initialization_inputs()

      if (config_state%start_from_restart) then
         call initialize_lcm_from_restart()
      else
         call initialize_lcm_from_config()
      endif

      lcm_is_initialized = .true.
   end subroutine lcm_initialize

   subroutine initialize_lcm_from_config()
      call initialize_particle_count_state()
      call lcm_allocate_particle_cells(initial_capacity_per_cell=config_state%np_xy)

      ! SAM-LCM micro_set also initialized grid-helper geometry. In standalone
      ! LCM that responsibility belongs to lcm_set_grid.
      lcm_initialization_mode = "fresh"
   end subroutine initialize_lcm_from_config

   subroutine initialize_lcm_from_restart()
      lcm_initialization_mode = "restart"
      error stop "LCM restart initialization is not implemented"
   end subroutine initialize_lcm_from_restart

   subroutine validate_initialization_inputs()
      if (.not. config_is_set) then
         error stop "LCM config must be set before lcm_initialize"
      endif

      if (.not. grid_state%is_set) then
         error stop "LCM grid must be set before lcm_initialize"
      endif

      if (config_state%nzl > grid_state%grid%nz) then
         error stop "LCM nzl must not exceed the vertical grid size"
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

   subroutine initialize_particle_count_state()
      integer :: initialization_levels

      initialization_levels = min(config_state%nzl, grid_state%grid%nz)

      lcm_initial_particle_count_local = &
         int(grid_state%grid%nx_local, int64) * &
         int(grid_state%grid%ny_local, int64) * &
         int(initialization_levels, int64) * &
         int(config_state%np_xy, int64)

      lcm_initial_particle_count_global = &
         int(grid_state%grid%nx_global, int64) * &
         int(grid_state%grid%ny_global, int64) * &
         int(initialization_levels, int64) * &
         int(config_state%np_xy, int64)
   end subroutine initialize_particle_count_state

end module lcm_initialization_state
