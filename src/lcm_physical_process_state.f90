!> Persistent state shared by LCM physical-process calculations.
!!
!! This module owns cell-mean thermodynamic history used between process calls.
!! Host atmospheric fields remain owned by lcm_field_state.
module lcm_physical_process_state
   use iso_fortran_env, only: real64
   use lcm_field_state, only: field_state
   use lcm_grid_state, only: grid_state

   implicit none
   private

   public :: lcm_ensure_cell_mean_history
   public :: lcm_get_cell_mean_history
   public :: lcm_save_cell_mean_history

   real(real64), allocatable :: previous_temperature(:,:,:)
   real(real64), allocatable :: previous_pressure(:,:,:)

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

end module lcm_physical_process_state
