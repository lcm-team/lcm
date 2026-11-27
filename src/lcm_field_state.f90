!> LCM-owned atmospheric field state supplied by the host model.
module lcm_field_state
   use iso_fortran_env, only: real64
   use lcm_grid_state, only: grid_state
   use lcm_types, only: lcm_field_state_t

   implicit none
   private

   public :: field_state
   public :: lcm_set_field

   type(lcm_field_state_t), protected, save :: field_state

   interface lcm_set_field
      module procedure set_field_3d
      module procedure set_field_profile
   end interface lcm_set_field

contains

   !> Copy one field into LCM-owned storage, reusing its existing allocation.
   !! The host supplies updated fields before initialization and each advance.
   !! Volume fields retain (x,y,z) ordering and the existing velocity staggering.
   !! Pressure and density are vertical profiles stored as (nz,1,1).
   !! Vapor qv is kg vapor / kg dry air. Supplying dry_air_density gives
   !! condensation the matching mass basis; otherwise it uses density.
   subroutine set_field_3d(name, field)
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
         call copy_3d_field(field, field_state%u, &
                            field_state%u_is_set, name, [nx+1,ny,nz])
      case ("v")
         call copy_3d_field(field, field_state%v, &
                            field_state%v_is_set, name, [nx,ny+1,nz])
      case ("w")
         call copy_3d_field(field, field_state%w, &
                            field_state%w_is_set, name, [nx,ny,nz+1])
      case ("static_energy")
         call copy_3d_field(field, field_state%static_energy, &
                            field_state%static_energy_is_set, name, [nx,ny,nz])
      case ("temperature")
         call copy_3d_field(field, field_state%temperature, &
                            field_state%temperature_is_set, name, [nx,ny,nz])
      case ("qv")
         call copy_3d_field(field, field_state%qv, &
                            field_state%qv_is_set, name, [nx,ny,nz])
      case ("cloud_liquid")
         call copy_3d_field(field, field_state%cloud_liquid, &
                            field_state%cloud_liquid_is_set, name, [nx,ny,nz])
      case ("dry_air_density")
         call copy_3d_field(field, field_state%dry_air_density, &
                            field_state%dry_air_density_is_set, name, [nx,ny,nz])
      case ("pressure")
         call copy_3d_field(field, field_state%pressure, &
                            field_state%pressure_is_set, name, [nz,1,1])
      case ("density")
         call copy_3d_field(field, field_state%density, &
                            field_state%density_is_set, name, [nz,1,1])
      case default
         error stop "LCM cannot set unknown field: " // trim(name)
      end select
   end subroutine set_field_3d

   !> Accept a host's 1D profile without changing its physical dimension.
   !! Packing into (nz,1,1) only standardizes LCM storage; it adds no values.
   subroutine set_field_profile(name, field)
      character(len=*), intent(in) :: name
      real(real64), intent(in) :: field(:)

      select case (trim(name))
      case ("pressure", "density")
         call set_field_3d(name, reshape(field, [size(field),1,1]))
      case default
         error stop "LCM 1D input requires a vertical profile field: " // trim(name)
      end select
   end subroutine set_field_profile

   subroutine copy_3d_field(field, values, is_set, name, expected_shape)
      real(real64), intent(in) :: field(:,:,:)
      real(real64), allocatable, intent(inout) :: values(:,:,:)
      logical, intent(inout) :: is_set
      character(len=*), intent(in) :: name
      integer, intent(in) :: expected_shape(3)

      if (any(shape(field) /= expected_shape)) then
         error stop "LCM " // trim(name) // " field shape does not match the grid"
      endif
      if (.not. allocated(values)) allocate(values, mold=field)
      values = field
      is_set = .true.
   end subroutine copy_3d_field

end module lcm_field_state
