!> LCM-owned atmospheric field state supplied by the host model.
module lcm_field_state
   use iso_fortran_env, only: real64
   use lcm_grid_state, only: grid_state
   use lcm_types, only: lcm_field_state_t

   implicit none

   private

   public :: field_state
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
         call attach_3d_field(u, field_state%u, field_state%u_is_attached, "u")
      endif
      if (present(v)) then
         call attach_3d_field(v, field_state%v, field_state%v_is_attached, "v")
      endif
      if (present(w)) then
         call attach_3d_field(w, field_state%w, field_state%w_is_attached, "w")
      endif
      if (present(static_energy)) then
         call attach_3d_field(static_energy, field_state%static_energy, &
                              field_state%static_energy_is_attached, "static_energy")
      endif
      if (present(temperature)) then
         call attach_3d_field(temperature, field_state%temperature, &
                              field_state%temperature_is_attached, "temperature")
      endif
      if (present(qv)) then
         call attach_3d_field(qv, field_state%qv, field_state%qv_is_attached, "qv")
      endif
      if (present(cloud_liquid)) then
         call attach_3d_field(cloud_liquid, field_state%cloud_liquid, &
                              field_state%cloud_liquid_is_attached, "cloud_liquid")
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

   subroutine attach_3d_field(field, field_pointer, is_attached, field_name)
      real(real64), intent(in), target :: field(:,:,:)
      real(real64), pointer, intent(inout) :: field_pointer(:,:,:)
      logical, intent(inout) :: is_attached
      character(len=*), intent(in) :: field_name

      call validate_3d_field_shape(field, field_name)
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

   subroutine validate_3d_field_shape(field, field_name)
      real(real64), intent(in), target :: field(:,:,:)
      character(len=*), intent(in) :: field_name

      if (any(shape(field) /= [grid_state%grid%nx_local, &
                               grid_state%grid%ny_local, &
                               grid_state%grid%nz])) then
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
