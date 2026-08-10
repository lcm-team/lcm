!> LCM-owned communication context supplied by the host model.
module lcm_parallel_state
   use lcm_types, only: lcm_parallel_t

   implicit none

   private

   public :: t_lcm_parallel_state
   public :: parallel_state
   public :: lcm_set_parallel

   type :: t_lcm_parallel_state
      type(lcm_parallel_t) :: parallel
      logical :: is_set = .false.
   end type t_lcm_parallel_state

   type(t_lcm_parallel_state), protected, save :: parallel_state

contains

   subroutine lcm_set_parallel(host_parallel)
      type(lcm_parallel_t), intent(in) :: host_parallel

      if (parallel_state%is_set) then
         error stop 'LCM parallel state has already been initialized'
      endif

      parallel_state%parallel = host_parallel
      parallel_state%is_set = .true.
   end subroutine lcm_set_parallel

end module lcm_parallel_state
