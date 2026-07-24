module micro_task_util_MPI
   implicit none
contains
   subroutine micro_task_abort()
      error stop "SAM-LCM probe micro_task_abort"
   end subroutine micro_task_abort
end module micro_task_util_MPI

module microphysics
   implicit none
   integer, parameter :: imv = 1
   integer, parameter :: nmicro_fields = 6
   integer, parameter :: nmicro_fields_diag = 76
   real, allocatable, save :: micro_field(:,:,:,:)
end module microphysics

subroutine micro_user_init_particles(i, j, k)
   implicit none
   integer, intent(in) :: i
   integer, intent(in) :: j
   integer, intent(in) :: k
end subroutine micro_user_init_particles
