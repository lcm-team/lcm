program test_particle_message
   use iso_fortran_env, only: int64, real64
   use lcm_particle_mpi_messages, only: lcm_pack_particle_mpi_message, &
                                        lcm_unpack_particle_mpi_message
   use lcm_types, only: lcm_particle_mpi_message_t, lcm_particle_t

   implicit none

   type(lcm_particle_mpi_message_t) :: message
   type(lcm_particle_t) :: particle
   type(lcm_particle_t) :: unpacked_particle
   integer :: global_i
   integer :: global_j
   integer :: k

   call populate_particle(particle)
   call lcm_pack_particle_mpi_message(21, 34, 5, particle, message)
   call lcm_unpack_particle_mpi_message(message, global_i, global_j, k, unpacked_particle)

   call assert_true(global_i == 21, 'global i')
   call assert_true(global_j == 34, 'global j')
   call assert_true(k == 5, 'vertical cell')
   call assert_particle_equal(unpacked_particle, particle)

contains

   subroutine populate_particle(particle)
      type(lcm_particle_t), intent(out) :: particle

      particle%weight = 17_int64
      particle%id = 42_int64
      particle%particle_type = 2
      particle%x = 1.0_real64
      particle%y = 2.0_real64
      particle%z = 3.0_real64
      particle%u = 4.0_real64
      particle%v = 5.0_real64
      particle%w = 6.0_real64
      particle%u_sgs = 7.0_real64
      particle%v_sgs = 8.0_real64
      particle%w_sgs = 9.0_real64
      particle%mass_water = 10.0_real64
      particle%mass_aero = 11.0_real64
      particle%w_sedi = 12.0_real64
      particle%temperature = 13.0_real64
      particle%qv = 14.0_real64
      particle%eta = 15.0_real64
      particle%axis_a = 16.0_real64
      particle%axis_c = 17.0_real64
      particle%density = 18.0_real64
      particle%rimed_mass = 19.0_real64
      particle%freezing_temperature = 20.0_real64
      particle%in_surface_area = [21.0_real64, 22.0_real64, 23.0_real64, 24.0_real64]
      particle%aerosol_portions = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
      particle%cloudy_age = 25.0_real64
      particle%cell_particle_index = 26
      particle%last_i = 27
      particle%last_j = 28
      particle%last_k = 29
      particle%n_monomer = 30_int64
      particle%spray_marker = 31
      particle%track_marker = 32
   end subroutine populate_particle

   subroutine assert_particle_equal(actual, expected)
      type(lcm_particle_t), intent(in) :: actual
      type(lcm_particle_t), intent(in) :: expected

      call assert_true(actual%weight == expected%weight, 'weight')
      call assert_true(actual%id == expected%id, 'id')
      call assert_true(actual%particle_type == expected%particle_type, 'type')
      call assert_true(all([actual%x, actual%y, actual%z, actual%u, actual%v, actual%w, &
                            actual%u_sgs, actual%v_sgs, actual%w_sgs] == &
                           [expected%x, expected%y, expected%z, expected%u, expected%v, expected%w, &
                            expected%u_sgs, expected%v_sgs, expected%w_sgs]), 'motion')
      call assert_true(all([actual%mass_water, actual%mass_aero, actual%w_sedi, &
                            actual%temperature, actual%qv, actual%eta] == &
                           [expected%mass_water, expected%mass_aero, expected%w_sedi, &
                            expected%temperature, expected%qv, expected%eta]), 'microphysics')
      call assert_true(all([actual%axis_a, actual%axis_c, actual%density, actual%rimed_mass, &
                            actual%freezing_temperature] == &
                           [expected%axis_a, expected%axis_c, expected%density, expected%rimed_mass, &
                            expected%freezing_temperature]), 'ice state')
      call assert_true(all(actual%in_surface_area == expected%in_surface_area), 'surface area')
      call assert_true(all(actual%aerosol_portions == expected%aerosol_portions), 'aerosol portions')
      call assert_true(actual%cloudy_age == expected%cloudy_age, 'cloudy age')
      call assert_true(all([actual%cell_particle_index, actual%last_i, actual%last_j, actual%last_k, &
                            actual%spray_marker, actual%track_marker] == &
                           [expected%cell_particle_index, expected%last_i, expected%last_j, expected%last_k, &
                            expected%spray_marker, expected%track_marker]), 'integer metadata')
      call assert_true(actual%n_monomer == expected%n_monomer, 'monomers')
   end subroutine assert_particle_equal

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) error stop 'assertion failed: ' // message
   end subroutine assert_true

end program test_particle_message
