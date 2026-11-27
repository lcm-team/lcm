!> LCM-owned atmospheric field state supplied by the host model.
module lcm_field_state
   use iso_fortran_env, only: real64
   use lcm_grid_state, only: grid_state
   use lcm_types, only: lcm_variable_t, lcm_fields_t

   implicit none
   private

   public :: field_state
   public :: lcm_set_field

   type(lcm_fields_t), protected, save :: field_state

contains

   !> Copy one field into LCM-owned storage, reusing its existing allocation.
   !! The host supplies updated fields before initialization and each advance.
   !! Volume fields retain (x,y,z) ordering and the existing velocity staggering.
   !! Pressure and density are vertical profiles stored as (nz,1,1).
   !! Hosts must reshape 1D profiles before calling this 3D-only interface.
   !! Vapor qv is kg vapor / kg dry air. Supplying dry_air_density gives
   !! condensation the matching mass basis; otherwise it uses density.
   subroutine lcm_set_field(name, field)
      character(len=*), intent(in) :: name
      real(real64), intent(in) :: field(:,:,:)
      integer :: nx, ny, nz

      if (.not. grid_state%is_set) then
         error stop "LCM grid must be initialized before setting fields"
      endif
      nx = grid_state%grid%nx_local
      ny = grid_state%grid%ny_local
      nz = grid_state%grid%nz

      select case (trim(name))
      case ("u")
         call copy_3d_field(field, field_state%u, [nx+1,ny,nz])
      case ("v")
         call copy_3d_field(field, field_state%v, [nx,ny+1,nz])
      case ("w")
         call copy_3d_field(field, field_state%w, [nx,ny,nz+1])
      case ("static_energy")
         call copy_3d_field(field, field_state%static_energy, [nx,ny,nz])
      case ("temperature")
         call copy_3d_field(field, field_state%temperature, [nx,ny,nz])
      case ("qv")
         call copy_3d_field(field, field_state%qv, [nx,ny,nz])
      case ("cloud_liquid")
         call copy_3d_field(field, field_state%cloud_liquid, [nx,ny,nz])
      case ("dry_air_density")
         call copy_3d_field(field, field_state%dry_air_density, [nx,ny,nz])
      case ("pressure")
         call copy_3d_field(field, field_state%pressure, [nz,1,1])
      case ("density")
         call copy_3d_field(field, field_state%density, [nz,1,1])
      case default
         error stop "LCM cannot set unknown field: " // trim(name)
      end select
   end subroutine lcm_set_field

   subroutine copy_3d_field(field, variable, expected_shape)
      real(real64), intent(in) :: field(:,:,:)
      type(lcm_variable_t), intent(inout) :: variable
      integer, intent(in) :: expected_shape(3)

      if (any(shape(field) /= expected_shape)) then
         error stop "LCM " // trim(variable%var_name) // " field shape does not match the grid"
      endif
      if (.not. allocated(variable%values)) allocate(variable%values, mold=field)
      variable%values = field
   end subroutine copy_3d_field

end module lcm_field_state
