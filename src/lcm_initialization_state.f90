!> Standalone LCM initialization lifecycle.
!!
!! Hosts call lcm_initialize after config, grid, and atmospheric fields have
!! been supplied. The timestep path may then assume LCM-owned storage exists.
module lcm_initialization_state
   use lcm_config_state, only: config_is_set
   use lcm_field_state, only: field_state
   use lcm_grid_state, only: grid_state
   use lcm_particle_state, only: lcm_allocate_particle_cells
   use micro_vars, only: np_xy

   implicit none

   private

   public :: lcm_initialize
   public :: lcm_is_initialized

   logical, protected, save :: lcm_is_initialized = .false.

contains

   subroutine lcm_initialize()
      if (lcm_is_initialized) then
         error stop "LCM has already been initialized"
      endif

      call validate_initialization_inputs()
      call lcm_allocate_particle_cells(initial_capacity_per_cell=np_xy)

      ! SAM-LCM micro_set also initialized grid-helper geometry. In standalone
      ! LCM that responsibility belongs to lcm_set_grid.
      lcm_is_initialized = .true.
   end subroutine lcm_initialize

   subroutine validate_initialization_inputs()
      if (.not. config_is_set) then
         error stop "LCM config must be set before lcm_initialize"
      endif

      if (.not. grid_state%is_set) then
         error stop "LCM grid must be set before lcm_initialize"
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

end module lcm_initialization_state
