program test_boundary_config
   use iso_fortran_env, only: real64
   use lcm_config_state, only: config_state, lcm_set_config
   use lcm_types, only: lcm_config_t, lcm_init_config

   implicit none

   type(lcm_config_t) :: config

   call lcm_init_config(config)
   call assert_close(config%bottom_removal_radius, 1.0e-6_real64, &
                     "default bottom-removal radius")

   config%bottom_removal_radius = 2.5e-6_real64
   call lcm_set_config(config)
   call assert_close(config_state%bottom_removal_radius, 2.5e-6_real64, &
                     "configured bottom-removal radius")

contains

   subroutine assert_close(actual, expected, message)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: message

      if (abs(actual - expected) > 1.0e-18_real64) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_close

end program test_boundary_config
