!> Shared LCM microphysics state for physical-process history and diagnostics.
!!
!! Persistent thermodynamic history is working state used between process
!! calls. Raw process increments are diagnostic state with a separate reset
!! lifecycle. Host atmospheric fields remain owned by lcm_field_state.
module lcm_microphysics_state
   use iso_fortran_env, only: real64
   use lcm_field_state, only: field_state
   use lcm_grid_state, only: grid_state

   implicit none
   private

   public :: lcm_ensure_cell_mean_history
   public :: lcm_get_cell_mean_history
   public :: lcm_save_cell_mean_history
   public :: lcm_ensure_condensation_diagnostics
   public :: lcm_reset_condensation_diagnostics
   public :: lcm_add_condensation_diagnostics

   real(real64), allocatable :: previous_temperature(:,:,:)
   real(real64), allocatable :: previous_pressure(:,:,:)

   ! Raw per-step increments, retained for a future diagnostics/output reader.
   real(real64), allocatable :: condensation_mass(:,:,:)
   real(real64), allocatable :: evaporation_mass(:,:,:)
   real(real64), allocatable :: activation_weight(:,:,:)
   real(real64), allocatable :: deactivation_weight(:,:,:)

contains

   !> Allocate history for the current local grid and seed it from host fields.
   subroutine lcm_ensure_cell_mean_history()
      integer :: nx, ny, nz, i, j, k

      nx = grid_state%grid%nx_local
      ny = grid_state%grid%ny_local
      nz = grid_state%grid%nz
      if (allocated(previous_temperature)) then
         if (all(shape(previous_temperature) == [nx,ny,nz])) return
         deallocate(previous_temperature, previous_pressure)
      endif

      allocate(previous_temperature(nx,ny,nz), previous_pressure(nx,ny,nz))
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               previous_temperature(i,j,k) = field_state%temperature(i,j,k)
               previous_pressure(i,j,k) = field_state%pressure(k)
            enddo
         enddo
      enddo
   end subroutine lcm_ensure_cell_mean_history

   !> Return the saved cell-mean temperature and pressure for process forcing.
   subroutine lcm_get_cell_mean_history(i, j, k, temperature, pressure)
      integer, intent(in) :: i, j, k
      real(real64), intent(out) :: temperature, pressure

      if (.not. allocated(previous_temperature)) then
         error stop 'LCM physical-process history has not been initialized'
      endif
      if (i < 1 .or. i > size(previous_temperature,1) .or. &
          j < 1 .or. j > size(previous_temperature,2) .or. &
          k < 1 .or. k > size(previous_temperature,3)) then
         error stop 'LCM physical-process history index is outside the local grid'
      endif
      temperature = previous_temperature(i,j,k)
      pressure = previous_pressure(i,j,k)
   end subroutine lcm_get_cell_mean_history

   !> Save the current host thermodynamic state for the next process call.
   subroutine lcm_save_cell_mean_history(i, j, k)
      integer, intent(in) :: i, j, k

      if (.not. allocated(previous_temperature)) then
         error stop 'LCM physical-process history has not been initialized'
      endif
      if (i < 1 .or. i > size(previous_temperature,1) .or. &
          j < 1 .or. j > size(previous_temperature,2) .or. &
          k < 1 .or. k > size(previous_temperature,3)) then
         error stop 'LCM physical-process history index is outside the local grid'
      endif
      previous_temperature(i,j,k) = field_state%temperature(i,j,k)
      previous_pressure(i,j,k) = field_state%pressure(k)
   end subroutine lcm_save_cell_mean_history

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
         error stop 'LCM microphysics diagnostics have not been allocated'
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
         error stop 'LCM microphysics diagnostics have not been allocated'
      endif
      if (i < 1 .or. i > size(condensation_mass,1) .or. &
          j < 1 .or. j > size(condensation_mass,2) .or. &
          k < 1 .or. k > size(condensation_mass,3)) then
         error stop 'LCM microphysics diagnostic index is outside the local grid'
      endif
      condensation_mass(i,j,k) = condensation_mass(i,j,k) + condensation
      evaporation_mass(i,j,k) = evaporation_mass(i,j,k) + evaporation
      activation_weight(i,j,k) = activation_weight(i,j,k) + activation
      deactivation_weight(i,j,k) = deactivation_weight(i,j,k) + deactivation
   end subroutine lcm_add_condensation_diagnostics

end module lcm_microphysics_state
