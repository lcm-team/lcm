!> LCM-owned atmospheric field state supplied by the host model.
module lcm_field_state
   use iso_fortran_env, only: real64
   use lcm_grid_state, only: grid_state
   use lcm_types, only: lcm_field_state_t

   implicit none

   private

   public :: field_state
   public :: lcm_attach_field
   public :: lcm_attach_fields

   type(lcm_field_state_t), protected, save :: field_state

contains

   subroutine lcm_attach_fields(u, v, w, static_energy, temperature, qv, &
                                cloud_liquid, pressure, density)
      real(real64), intent(in), target, optional :: u(:,:,:)
      real(real64), intent(in), target, optional :: v(:,:,:)
      real(real64), intent(in), target, optional :: w(:,:,:)
      real(real64), intent(in), target, optional :: static_energy(:,:,:)
      real(real64), intent(in), target, optional :: temperature(:,:,:)
      real(real64), intent(in), target, optional :: qv(:,:,:)
      real(real64), intent(in), target, optional :: cloud_liquid(:,:,:)
      real(real64), intent(in), target, optional :: pressure(:)
      real(real64), intent(in), target, optional :: density(:)
      logical :: has_field

      if (.not. grid_state%is_set) then
         error stop "LCM grid must be initialized before attaching fields"
      endif

      has_field = present(u) .or. present(v) .or. present(w) .or. &
                  present(static_energy) .or. present(temperature) .or. &
                  present(qv) .or. present(cloud_liquid) .or. &
                  present(pressure) .or. present(density)

      if (.not. has_field) then
         error stop "LCM attach fields requires at least one field"
      endif

      if (present(u)) then
         call attach_3d_field(u, field_state%u, field_state%u_is_attached, &
                              "u", grid_state%grid%nx_local, &
                              grid_state%grid%ny_local, grid_state%grid%nz)
      endif
      if (present(v)) then
         call attach_3d_field(v, field_state%v, field_state%v_is_attached, &
                              "v", grid_state%grid%nx_local, &
                              grid_state%grid%ny_local, grid_state%grid%nz)
      endif
      if (present(w)) then
         call attach_3d_field(w, field_state%w, field_state%w_is_attached, &
                              "w", grid_state%grid%nx_local, &
                              grid_state%grid%ny_local, grid_state%grid%nz)
      endif
      if (present(static_energy)) then
         call attach_cell_center_field("static_energy", static_energy, &
                                       field_state%static_energy, &
                                       field_state%static_energy_is_attached)
      endif
      if (present(temperature)) then
         call attach_cell_center_field("temperature", temperature, &
                                       field_state%temperature, &
                                       field_state%temperature_is_attached)
      endif
      if (present(qv)) then
         call attach_cell_center_field("qv", qv, field_state%qv, &
                                       field_state%qv_is_attached)
      endif
      if (present(cloud_liquid)) then
         call attach_cell_center_field("cloud_liquid", cloud_liquid, &
                                       field_state%cloud_liquid, &
                                       field_state%cloud_liquid_is_attached)
      endif
      if (present(pressure)) then
         call attach_profile_field(pressure, field_state%pressure, &
                                   field_state%pressure_is_attached, "pressure")
      endif
      if (present(density)) then
         call attach_profile_field(density, field_state%density, &
                                   field_state%density_is_attached, "density")
      endif
   end subroutine lcm_attach_fields

   subroutine lcm_attach_field(name, field)
      character(len=*), intent(in) :: name
      real(real64), intent(in), target :: field(:,:,:)

      if (.not. grid_state%is_set) then
         error stop "LCM grid must be initialized before attaching fields"
      endif

      select case (trim(name))
      case ("u")
         call attach_3d_field(field, field_state%u, field_state%u_is_attached, &
                              "u", grid_state%grid%nx_local + 1, &
                              grid_state%grid%ny_local, grid_state%grid%nz)
      case ("v")
         call attach_3d_field(field, field_state%v, field_state%v_is_attached, &
                              "v", grid_state%grid%nx_local, &
                              grid_state%grid%ny_local + 1, grid_state%grid%nz)
      case ("w")
         call attach_3d_field(field, field_state%w, field_state%w_is_attached, &
                              "w", grid_state%grid%nx_local, &
                              grid_state%grid%ny_local, grid_state%grid%nz + 1)
      case ("static_energy")
         call attach_cell_center_field("static_energy", field, &
                                       field_state%static_energy, &
                                       field_state%static_energy_is_attached)
      case ("temperature")
         call attach_cell_center_field("temperature", field, &
                                       field_state%temperature, &
                                       field_state%temperature_is_attached)
      case ("qv")
         call attach_cell_center_field("qv", field, field_state%qv, &
                                       field_state%qv_is_attached)
      case ("cloud_liquid")
         call attach_cell_center_field("cloud_liquid", field, &
                                       field_state%cloud_liquid, &
                                       field_state%cloud_liquid_is_attached)
      case default
         error stop "LCM cannot attach unknown field: " // trim(name)
      end select
   end subroutine lcm_attach_field

   subroutine attach_cell_center_field(field_name, field, field_pointer, is_attached)
      character(len=*), intent(in) :: field_name
      real(real64), intent(in), target :: field(:,:,:)
      real(real64), pointer, intent(inout) :: field_pointer(:,:,:)
      logical, intent(inout) :: is_attached

      call attach_3d_field(field, field_pointer, is_attached, field_name, &
                           grid_state%grid%nx_local, grid_state%grid%ny_local, &
                           grid_state%grid%nz)
   end subroutine attach_cell_center_field

   subroutine attach_3d_field(field, field_pointer, is_attached, field_name, &
                              expected_nx, expected_ny, expected_nz)
      real(real64), intent(in), target :: field(:,:,:)
      real(real64), pointer, intent(inout) :: field_pointer(:,:,:)
      logical, intent(inout) :: is_attached
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: expected_nx
      integer, intent(in) :: expected_ny
      integer, intent(in) :: expected_nz

      call validate_3d_field_shape(field, field_name, expected_nx, expected_ny, &
                                   expected_nz)
      field_pointer => field
      is_attached = .true.
   end subroutine attach_3d_field

   subroutine attach_profile_field(field, field_pointer, is_attached, field_name)
      real(real64), intent(in), target :: field(:)
      real(real64), pointer, intent(inout) :: field_pointer(:)
      logical, intent(inout) :: is_attached
      character(len=*), intent(in) :: field_name

      call validate_profile_field_shape(field, field_name)
      field_pointer => field
      is_attached = .true.
   end subroutine attach_profile_field

   subroutine validate_3d_field_shape(field, field_name, expected_nx, expected_ny, &
                                      expected_nz)
      real(real64), intent(in), target :: field(:,:,:)
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: expected_nx
      integer, intent(in) :: expected_ny
      integer, intent(in) :: expected_nz

      if (any(shape(field) /= [expected_nx, expected_ny, expected_nz])) then
         error stop "LCM " // field_name // " field shape does not match the grid"
      endif
   end subroutine validate_3d_field_shape

   subroutine validate_profile_field_shape(field, field_name)
      real(real64), intent(in), target :: field(:)
      character(len=*), intent(in) :: field_name

      if (size(field) /= grid_state%grid%nz) then
         error stop "LCM " // field_name // " profile shape does not match the grid"
      endif
   end subroutine validate_profile_field_shape

end module lcm_field_state
