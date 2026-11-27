!> LCM-owned atmospheric field state supplied by the host model.
module lcm_field_state
   use iso_fortran_env, only: real64
   use lcm_grid_state, only: grid_state
   use lcm_types, only: lcm_field_state_t

   implicit none

   private

   public :: field_state
   public :: lcm_set_field
   public :: lcm_set_fields

   type(lcm_field_state_t), protected, save :: field_state

contains

   !> Copy supplied fields into LCM-owned storage. Omitted fields are unchanged.
   !! Supply updated values before initialization and each advance; LCM never
   !! retains a reference to the input arrays. Existing allocations are reused.
   !! Vapor qv is kg vapor / kg dry air. Optional dry_air_density supplies
   !! the matching mass basis for condensation; density remains the air profile
   !! used by particle initialization and ventilation. Without the optional
   !! field, condensation retains the density profile as its mass basis.
   subroutine lcm_set_fields(u, v, w, static_energy, temperature, qv, &
                                cloud_liquid, pressure, density, dry_air_density)
      real(real64), intent(in), optional :: u(:,:,:)
      real(real64), intent(in), optional :: v(:,:,:)
      real(real64), intent(in), optional :: w(:,:,:)
      real(real64), intent(in), optional :: static_energy(:,:,:)
      real(real64), intent(in), optional :: temperature(:,:,:)
      real(real64), intent(in), optional :: qv(:,:,:)
      real(real64), intent(in), optional :: cloud_liquid(:,:,:)
      real(real64), intent(in), optional :: pressure(:)
      real(real64), intent(in), optional :: density(:)
      real(real64), intent(in), optional :: dry_air_density(:,:,:)
      logical :: has_field

      if (.not. grid_state%is_set) then
         error stop "LCM grid must be initialized before setting fields"
      endif

      has_field = present(u) .or. present(v) .or. present(w) .or. &
                  present(static_energy) .or. present(temperature) .or. &
                  present(qv) .or. present(cloud_liquid) .or. &
                  present(pressure) .or. present(density) .or. present(dry_air_density)

      if (.not. has_field) then
         error stop "LCM set fields requires at least one field"
      endif

      if (present(u)) then
         call copy_3d_field(u, field_state%u, field_state%u_is_set, &
                              "u", grid_state%grid%nx_local + 1, &
                              grid_state%grid%ny_local, grid_state%grid%nz)
      endif
      if (present(v)) then
         call copy_3d_field(v, field_state%v, field_state%v_is_set, &
                              "v", grid_state%grid%nx_local, &
                              grid_state%grid%ny_local + 1, grid_state%grid%nz)
      endif
      if (present(w)) then
         call copy_3d_field(w, field_state%w, field_state%w_is_set, &
                              "w", grid_state%grid%nx_local, &
                              grid_state%grid%ny_local, grid_state%grid%nz + 1)
      endif
      if (present(static_energy)) then
         call copy_cell_center_field("static_energy", static_energy, &
                                       field_state%static_energy, &
                                       field_state%static_energy_is_set)
      endif
      if (present(temperature)) then
         call copy_cell_center_field("temperature", temperature, &
                                       field_state%temperature, &
                                       field_state%temperature_is_set)
      endif
      if (present(qv)) then
         call copy_cell_center_field("qv", qv, field_state%qv, &
                                       field_state%qv_is_set)
      endif
      if (present(cloud_liquid)) then
         call copy_cell_center_field("cloud_liquid", cloud_liquid, &
                                       field_state%cloud_liquid, &
                                       field_state%cloud_liquid_is_set)
      endif
      if (present(pressure)) then
         call copy_profile_field(pressure, field_state%pressure, &
                                   field_state%pressure_is_set, "pressure")
      endif
      if (present(density)) then
         call copy_profile_field(density, field_state%density, &
                                   field_state%density_is_set, "density")
      endif
      if (present(dry_air_density)) then
         call copy_cell_center_field("dry_air_density", dry_air_density, &
                                       field_state%dry_air_density, &
                                       field_state%dry_air_density_is_set)
      endif
   end subroutine lcm_set_fields

   subroutine lcm_set_field(name, field)
      character(len=*), intent(in) :: name
      real(real64), intent(in) :: field(:,:,:)

      if (.not. grid_state%is_set) then
         error stop "LCM grid must be initialized before setting fields"
      endif

      select case (trim(name))
      case ("u")
         call copy_3d_field(field, field_state%u, field_state%u_is_set, &
                              "u", grid_state%grid%nx_local + 1, &
                              grid_state%grid%ny_local, grid_state%grid%nz)
      case ("v")
         call copy_3d_field(field, field_state%v, field_state%v_is_set, &
                              "v", grid_state%grid%nx_local, &
                              grid_state%grid%ny_local + 1, grid_state%grid%nz)
      case ("w")
         call copy_3d_field(field, field_state%w, field_state%w_is_set, &
                              "w", grid_state%grid%nx_local, &
                              grid_state%grid%ny_local, grid_state%grid%nz + 1)
      case ("static_energy")
         call copy_cell_center_field("static_energy", field, &
                                       field_state%static_energy, &
                                       field_state%static_energy_is_set)
      case ("temperature")
         call copy_cell_center_field("temperature", field, &
                                       field_state%temperature, &
                                       field_state%temperature_is_set)
      case ("qv")
         call copy_cell_center_field("qv", field, field_state%qv, &
                                       field_state%qv_is_set)
      case ("cloud_liquid")
         call copy_cell_center_field("cloud_liquid", field, &
                                       field_state%cloud_liquid, &
                                       field_state%cloud_liquid_is_set)
      case ("dry_air_density")
         call copy_cell_center_field("dry_air_density", field, &
                                       field_state%dry_air_density, &
                                       field_state%dry_air_density_is_set)
      case default
         error stop "LCM cannot set unknown field: " // trim(name)
      end select
   end subroutine lcm_set_field

   subroutine copy_cell_center_field(field_name, field, field_values, is_set)
      character(len=*), intent(in) :: field_name
      real(real64), intent(in) :: field(:,:,:)
      real(real64), allocatable, intent(inout) :: field_values(:,:,:)
      logical, intent(inout) :: is_set

      call copy_3d_field(field, field_values, is_set, field_name, &
                           grid_state%grid%nx_local, grid_state%grid%ny_local, &
                           grid_state%grid%nz)
   end subroutine copy_cell_center_field

   subroutine copy_3d_field(field, field_values, is_set, field_name, &
                              expected_nx, expected_ny, expected_nz)
      real(real64), intent(in) :: field(:,:,:)
      real(real64), allocatable, intent(inout) :: field_values(:,:,:)
      logical, intent(inout) :: is_set
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: expected_nx
      integer, intent(in) :: expected_ny
      integer, intent(in) :: expected_nz

      call validate_3d_field_shape(field, field_name, expected_nx, expected_ny, &
                                   expected_nz)
      if (.not. allocated(field_values)) allocate(field_values, mold=field)
      field_values = field
      is_set = .true.
   end subroutine copy_3d_field

   subroutine copy_profile_field(field, field_values, is_set, field_name)
      real(real64), intent(in) :: field(:)
      real(real64), allocatable, intent(inout) :: field_values(:)
      logical, intent(inout) :: is_set
      character(len=*), intent(in) :: field_name

      call validate_profile_field_shape(field, field_name)
      if (.not. allocated(field_values)) allocate(field_values, mold=field)
      field_values = field
      is_set = .true.
   end subroutine copy_profile_field

   subroutine validate_3d_field_shape(field, field_name, expected_nx, expected_ny, &
                                      expected_nz)
      real(real64), intent(in) :: field(:,:,:)
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: expected_nx
      integer, intent(in) :: expected_ny
      integer, intent(in) :: expected_nz

      if (any(shape(field) /= [expected_nx, expected_ny, expected_nz])) then
         error stop "LCM " // field_name // " field shape does not match the grid"
      endif
   end subroutine validate_3d_field_shape

   subroutine validate_profile_field_shape(field, field_name)
      real(real64), intent(in) :: field(:)
      character(len=*), intent(in) :: field_name

      if (size(field) /= grid_state%grid%nz) then
         error stop "LCM " // field_name // " profile shape does not match the grid"
      endif
   end subroutine validate_profile_field_shape

end module lcm_field_state
