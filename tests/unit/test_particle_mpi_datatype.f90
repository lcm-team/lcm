program test_particle_mpi_datatype
   use iso_fortran_env, only: int64, real64
   use mpi_f08, only: MPI_COMM_RANK, MPI_COMM_SIZE, MPI_COMM_WORLD, MPI_Datatype, &
                      MPI_FINALIZE, MPI_INIT, MPI_SENDRECV, MPI_Status, MPI_SUCCESS
   use lcm_particle_mpi_datatype, only: lcm_free_particle_mpi_datatype, &
                                        lcm_get_particle_mpi_datatype
   use lcm_types, only: lcm_particle_mpi_message_t

   implicit none

   type(lcm_particle_mpi_message_t) :: expected_message
   type(lcm_particle_mpi_message_t) :: receive_messages(2)
   type(lcm_particle_mpi_message_t) :: send_messages(2)
   type(MPI_Datatype) :: particle_message_datatype
   type(MPI_Status) :: status
   integer :: ierr
   integer :: message_index
   integer :: partner_rank
   integer :: rank
   integer :: mpi_size

   call MPI_INIT(ierr)
   call require_mpi_success(ierr, 'initialize MPI')
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   call require_mpi_success(ierr, 'query the test rank')
   call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
   call require_mpi_success(ierr, 'query the test size')
   if (mpi_size /= 2) error stop 'particle MPI datatype test requires two ranks'

   partner_rank = 1 - rank
   do message_index = 1, size(send_messages)
      call fill_message(send_messages(message_index), rank, message_index)
   enddo

   particle_message_datatype = lcm_get_particle_mpi_datatype()
   call MPI_SENDRECV(send_messages, size(send_messages), &
                     particle_message_datatype, partner_rank, 701, &
                     receive_messages, size(receive_messages), &
                     particle_message_datatype, partner_rank, 701, &
                     MPI_COMM_WORLD, status, ierr)
   call require_mpi_success(ierr, 'exchange typed particle messages')

   do message_index = 1, size(receive_messages)
      call fill_message(expected_message, partner_rank, message_index)
      call assert_message_equal(receive_messages(message_index), expected_message)
   enddo

   call lcm_free_particle_mpi_datatype()
   call MPI_FINALIZE(ierr)
   call require_mpi_success(ierr, 'finalize MPI')

contains

   subroutine fill_message(message, owner_rank, sequence)
      type(lcm_particle_mpi_message_t), intent(out) :: message
      integer, intent(in) :: owner_rank
      integer, intent(in) :: sequence
      integer :: base
      real(real64) :: real_base

      base = 1000 * owner_rank + 100 * sequence
      real_base = real(base, real64)

      message%global_i = base + 1
      message%global_j = base + 2
      message%k = base + 3
      message%particle%weight = int(base + 4, int64)
      message%particle%id = int(base + 5, int64)
      message%particle%particle_type = base + 6
      message%particle%x = real_base + 0.25_real64
      message%particle%y = real_base + 0.50_real64
      message%particle%z = real_base + 0.75_real64
      message%particle%u = real_base + 1.00_real64
      message%particle%v = real_base + 1.25_real64
      message%particle%w = real_base + 1.50_real64
      message%particle%u_sgs = real_base + 1.75_real64
      message%particle%v_sgs = real_base + 2.00_real64
      message%particle%w_sgs = real_base + 2.25_real64
      message%particle%mass_water = real_base + 2.50_real64
      message%particle%mass_aero = real_base + 2.75_real64
      message%particle%w_sedi = real_base + 3.00_real64
      message%particle%temperature = real_base + 3.25_real64
      message%particle%qv = real_base + 3.50_real64
      message%particle%eta = real_base + 3.75_real64
      message%particle%axis_a = real_base + 4.00_real64
      message%particle%axis_c = real_base + 4.25_real64
      message%particle%density = real_base + 4.50_real64
      message%particle%rimed_mass = real_base + 4.75_real64
      message%particle%freezing_temperature = real_base + 5.00_real64
      message%particle%in_surface_area = real_base + &
         [5.25_real64, 5.50_real64, 5.75_real64, 6.00_real64]
      message%particle%aerosol_portions = real_base + &
         [6.25_real64, 6.50_real64, 6.75_real64, 7.00_real64]
      message%particle%cloudy_age = real_base + 7.25_real64
      message%particle%cell_particle_index = base + 7
      message%particle%last_i = base + 8
      message%particle%last_j = base + 9
      message%particle%last_k = base + 10
      message%particle%n_monomer = int(base + 11, int64)
      message%particle%spray_marker = base + 12
      message%particle%track_marker = base + 13
   end subroutine fill_message

   subroutine assert_message_equal(actual, expected)
      type(lcm_particle_mpi_message_t), intent(in) :: actual
      type(lcm_particle_mpi_message_t), intent(in) :: expected

      call assert_integer(actual%global_i, expected%global_i, 'global_i')
      call assert_integer(actual%global_j, expected%global_j, 'global_j')
      call assert_integer(actual%k, expected%k, 'k')
      call assert_int64(actual%particle%weight, expected%particle%weight, 'weight')
      call assert_int64(actual%particle%id, expected%particle%id, 'id')
      call assert_integer(actual%particle%particle_type, &
                          expected%particle%particle_type, 'particle_type')
      call assert_real64(actual%particle%x, expected%particle%x, 'x')
      call assert_real64(actual%particle%y, expected%particle%y, 'y')
      call assert_real64(actual%particle%z, expected%particle%z, 'z')
      call assert_real64(actual%particle%u, expected%particle%u, 'u')
      call assert_real64(actual%particle%v, expected%particle%v, 'v')
      call assert_real64(actual%particle%w, expected%particle%w, 'w')
      call assert_real64(actual%particle%u_sgs, expected%particle%u_sgs, 'u_sgs')
      call assert_real64(actual%particle%v_sgs, expected%particle%v_sgs, 'v_sgs')
      call assert_real64(actual%particle%w_sgs, expected%particle%w_sgs, 'w_sgs')
      call assert_real64(actual%particle%mass_water, expected%particle%mass_water, &
                         'mass_water')
      call assert_real64(actual%particle%mass_aero, expected%particle%mass_aero, &
                         'mass_aero')
      call assert_real64(actual%particle%w_sedi, expected%particle%w_sedi, 'w_sedi')
      call assert_real64(actual%particle%temperature, &
                         expected%particle%temperature, 'temperature')
      call assert_real64(actual%particle%qv, expected%particle%qv, 'qv')
      call assert_real64(actual%particle%eta, expected%particle%eta, 'eta')
      call assert_real64(actual%particle%axis_a, expected%particle%axis_a, 'axis_a')
      call assert_real64(actual%particle%axis_c, expected%particle%axis_c, 'axis_c')
      call assert_real64(actual%particle%density, expected%particle%density, 'density')
      call assert_real64(actual%particle%rimed_mass, expected%particle%rimed_mass, &
                         'rimed_mass')
      call assert_real64(actual%particle%freezing_temperature, &
                         expected%particle%freezing_temperature, &
                         'freezing_temperature')
      call assert_real64_array(actual%particle%in_surface_area, &
                               expected%particle%in_surface_area, &
                               'in_surface_area')
      call assert_real64_array(actual%particle%aerosol_portions, &
                               expected%particle%aerosol_portions, &
                               'aerosol_portions')
      call assert_real64(actual%particle%cloudy_age, &
                         expected%particle%cloudy_age, 'cloudy_age')
      call assert_integer(actual%particle%cell_particle_index, &
                          expected%particle%cell_particle_index, &
                          'cell_particle_index')
      call assert_integer(actual%particle%last_i, expected%particle%last_i, 'last_i')
      call assert_integer(actual%particle%last_j, expected%particle%last_j, 'last_j')
      call assert_integer(actual%particle%last_k, expected%particle%last_k, 'last_k')
      call assert_int64(actual%particle%n_monomer, expected%particle%n_monomer, &
                        'n_monomer')
      call assert_integer(actual%particle%spray_marker, &
                          expected%particle%spray_marker, 'spray_marker')
      call assert_integer(actual%particle%track_marker, &
                          expected%particle%track_marker, 'track_marker')
   end subroutine assert_message_equal

   subroutine assert_integer(actual, expected, component)
      integer, intent(in) :: actual
      integer, intent(in) :: expected
      character(len=*), intent(in) :: component

      if (actual /= expected) error stop 'integer mismatch: ' // component
   end subroutine assert_integer

   subroutine assert_int64(actual, expected, component)
      integer(int64), intent(in) :: actual
      integer(int64), intent(in) :: expected
      character(len=*), intent(in) :: component

      if (actual /= expected) error stop 'int64 mismatch: ' // component
   end subroutine assert_int64

   subroutine assert_real64(actual, expected, component)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: component

      if (actual /= expected) error stop 'real64 mismatch: ' // component
   end subroutine assert_real64

   subroutine assert_real64_array(actual, expected, component)
      real(real64), intent(in) :: actual(:)
      real(real64), intent(in) :: expected(:)
      character(len=*), intent(in) :: component

      if (any(actual /= expected)) error stop 'real64-array mismatch: ' // component
   end subroutine assert_real64_array

   subroutine require_mpi_success(error_code, operation)
      integer, intent(in) :: error_code
      character(len=*), intent(in) :: operation

      if (error_code /= MPI_SUCCESS) error stop 'failed to ' // operation
   end subroutine require_mpi_success

end program test_particle_mpi_datatype
