!> LCM-owned process diagnostics prepared for later output and analysis.
!!
!! These arrays hold raw per-step increments before normalization. Their future
!! output interface must define reduction, normalization, and accumulation scope.
module lcm_diagnostics_state
   use iso_fortran_env, only: real64
   use lcm_grid_state, only: grid_state

   implicit none
   private

   public :: lcm_ensure_condensation_diagnostics
   public :: lcm_reset_condensation_diagnostics
   public :: lcm_add_condensation_diagnostics

   real(real64), allocatable :: condensation_mass(:,:,:)
   real(real64), allocatable :: evaporation_mass(:,:,:)
   real(real64), allocatable :: activation_weight(:,:,:)
   real(real64), allocatable :: deactivation_weight(:,:,:)

contains

   subroutine lcm_ensure_condensation_diagnostics()
      integer :: nx, ny, nz

      nx = grid_state%grid%nx_local
      ny = grid_state%grid%ny_local
      nz = grid_state%grid%nz
      if (allocated(condensation_mass)) then
         if (all(shape(condensation_mass) == [nx,ny,nz])) return
         deallocate(condensation_mass, evaporation_mass, &
                    activation_weight, deactivation_weight)
      endif
      allocate(condensation_mass(nx,ny,nz), evaporation_mass(nx,ny,nz), &
               activation_weight(nx,ny,nz), deactivation_weight(nx,ny,nz))
   end subroutine lcm_ensure_condensation_diagnostics

   subroutine lcm_reset_condensation_diagnostics()
      if (.not. allocated(condensation_mass)) then
         error stop 'LCM condensation diagnostics have not been allocated'
      endif
      condensation_mass = 0.0_real64
      evaporation_mass = 0.0_real64
      activation_weight = 0.0_real64
      deactivation_weight = 0.0_real64
   end subroutine lcm_reset_condensation_diagnostics

   subroutine lcm_add_condensation_diagnostics(i, j, k, condensation, &
                                                evaporation, activation, deactivation)
      integer, intent(in) :: i, j, k
      real(real64), intent(in) :: condensation, evaporation
      real(real64), intent(in) :: activation, deactivation

      if (.not. allocated(condensation_mass)) then
         error stop 'LCM condensation diagnostics have not been allocated'
      endif
      if (i < 1 .or. i > size(condensation_mass,1) .or. &
          j < 1 .or. j > size(condensation_mass,2) .or. &
          k < 1 .or. k > size(condensation_mass,3)) then
         error stop 'LCM condensation diagnostic index is outside the local grid'
      endif
      condensation_mass(i,j,k) = condensation_mass(i,j,k) + condensation
      evaporation_mass(i,j,k) = evaporation_mass(i,j,k) + evaporation
      activation_weight(i,j,k) = activation_weight(i,j,k) + activation
      deactivation_weight(i,j,k) = deactivation_weight(i,j,k) + deactivation
   end subroutine lcm_add_condensation_diagnostics

end module lcm_diagnostics_state
