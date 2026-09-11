!> MPI datatype for typed particle-exchange messages.
module lcm_particle_mpi_datatype
   use iso_fortran_env, only: int64, real64
   use mpi_f08, only: MPI_ADDRESS_KIND, MPI_Datatype, MPI_DATATYPE_NULL, &
                      MPI_GET_ADDRESS, MPI_INTEGER, MPI_SUCCESS, &
                      MPI_TYPE_COMMIT, MPI_TYPE_CREATE_RESIZED, &
                      MPI_TYPE_CREATE_STRUCT, MPI_TYPE_FREE, &
                      MPI_TYPE_MATCH_SIZE, MPI_TYPECLASS_INTEGER, &
                      MPI_TYPECLASS_REAL
   use lcm_types, only: lcm_particle_mpi_message_t

   implicit none

   private

   public :: lcm_free_particle_mpi_datatype
   public :: lcm_get_particle_mpi_datatype

   integer, parameter :: particle_message_block_count = 36

   type(MPI_Datatype), save :: particle_message_datatype = MPI_DATATYPE_NULL
   logical, save :: particle_message_datatype_is_initialized = .false.

   interface component_address
      module procedure default_integer_address
      module procedure int64_address
      module procedure real64_address
      module procedure real64_array_address
   end interface component_address

contains

   function lcm_get_particle_mpi_datatype() result(datatype)
      type(MPI_Datatype) :: datatype

      if (.not. particle_message_datatype_is_initialized) then
         call create_particle_message_datatype()
      endif
      datatype = particle_message_datatype
   end function lcm_get_particle_mpi_datatype

   subroutine lcm_free_particle_mpi_datatype()
      integer :: ierr

      if (.not. particle_message_datatype_is_initialized) return

      call MPI_TYPE_FREE(particle_message_datatype, ierr)
      call require_mpi_success(ierr, 'free the particle-message datatype')
      particle_message_datatype = MPI_DATATYPE_NULL
      particle_message_datatype_is_initialized = .false.
   end subroutine lcm_free_particle_mpi_datatype

   subroutine create_particle_message_datatype()
      type(lcm_particle_mpi_message_t), target :: sample(2)
      integer :: block_lengths(particle_message_block_count)
      integer(MPI_ADDRESS_KIND) :: displacements(particle_message_block_count)
      type(MPI_Datatype) :: block_types(particle_message_block_count)
      type(MPI_Datatype) :: int64_datatype
      type(MPI_Datatype) :: real64_datatype
      type(MPI_Datatype) :: unresized_datatype
      integer(MPI_ADDRESS_KIND) :: base_address
      integer(MPI_ADDRESS_KIND) :: next_message_address
      integer(MPI_ADDRESS_KIND) :: message_extent
      integer :: ierr

      call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, &
                               storage_size(sample(1)%particle%weight) / 8, &
                               int64_datatype, ierr)
      call require_mpi_success(ierr, 'match the particle int64 datatype')
      call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, &
                               storage_size(sample(1)%particle%x) / 8, &
                               real64_datatype, ierr)
      call require_mpi_success(ierr, 'match the particle real64 datatype')

      block_lengths = 1
      block_lengths(27:28) = 4
      block_types(1:3) = MPI_INTEGER
      block_types(4:5) = int64_datatype
      block_types(6) = MPI_INTEGER
      block_types(7:29) = real64_datatype
      block_types(30:33) = MPI_INTEGER
      block_types(34) = int64_datatype
      block_types(35:36) = MPI_INTEGER

      call component_address(sample(1)%global_i, displacements(1))
      call component_address(sample(1)%global_j, displacements(2))
      call component_address(sample(1)%k, displacements(3))
      call component_address(sample(1)%particle%weight, displacements(4))
      call component_address(sample(1)%particle%id, displacements(5))
      call component_address(sample(1)%particle%particle_type, displacements(6))
      call component_address(sample(1)%particle%x, displacements(7))
      call component_address(sample(1)%particle%y, displacements(8))
      call component_address(sample(1)%particle%z, displacements(9))
      call component_address(sample(1)%particle%u, displacements(10))
      call component_address(sample(1)%particle%v, displacements(11))
      call component_address(sample(1)%particle%w, displacements(12))
      call component_address(sample(1)%particle%u_sgs, displacements(13))
      call component_address(sample(1)%particle%v_sgs, displacements(14))
      call component_address(sample(1)%particle%w_sgs, displacements(15))
      call component_address(sample(1)%particle%mass_water, displacements(16))
      call component_address(sample(1)%particle%mass_aero, displacements(17))
      call component_address(sample(1)%particle%w_sedi, displacements(18))
      call component_address(sample(1)%particle%temperature, displacements(19))
      call component_address(sample(1)%particle%qv, displacements(20))
      call component_address(sample(1)%particle%eta, displacements(21))
      call component_address(sample(1)%particle%axis_a, displacements(22))
      call component_address(sample(1)%particle%axis_c, displacements(23))
      call component_address(sample(1)%particle%density, displacements(24))
      call component_address(sample(1)%particle%rimed_mass, displacements(25))
      call component_address(sample(1)%particle%freezing_temperature, displacements(26))
      call component_address(sample(1)%particle%in_surface_area, displacements(27))
      call component_address(sample(1)%particle%aerosol_portions, displacements(28))
      call component_address(sample(1)%particle%cloudy_age, displacements(29))
      call component_address(sample(1)%particle%cell_particle_index, displacements(30))
      call component_address(sample(1)%particle%last_i, displacements(31))
      call component_address(sample(1)%particle%last_j, displacements(32))
      call component_address(sample(1)%particle%last_k, displacements(33))
      call component_address(sample(1)%particle%n_monomer, displacements(34))
      call component_address(sample(1)%particle%spray_marker, displacements(35))
      call component_address(sample(1)%particle%track_marker, displacements(36))

      base_address = displacements(1)
      displacements = displacements - base_address
      call component_address(sample(2)%global_i, next_message_address)
      message_extent = next_message_address - base_address

      call MPI_TYPE_CREATE_STRUCT(particle_message_block_count, block_lengths, &
                                  displacements, block_types, &
                                  unresized_datatype, ierr)
      call require_mpi_success(ierr, 'create the particle-message datatype')
      call MPI_TYPE_CREATE_RESIZED(unresized_datatype, 0_MPI_ADDRESS_KIND, &
                                   message_extent, particle_message_datatype, ierr)
      call require_mpi_success(ierr, 'resize the particle-message datatype')
      call MPI_TYPE_COMMIT(particle_message_datatype, ierr)
      call require_mpi_success(ierr, 'commit the particle-message datatype')
      call MPI_TYPE_FREE(unresized_datatype, ierr)
      call require_mpi_success(ierr, 'free the unresized particle-message datatype')

      particle_message_datatype_is_initialized = .true.
   end subroutine create_particle_message_datatype

   subroutine default_integer_address(value, address)
      integer, target, intent(in) :: value
      integer(MPI_ADDRESS_KIND), intent(out) :: address
      integer :: ierr

      call MPI_GET_ADDRESS(value, address, ierr)
      call require_mpi_success(ierr, 'query a default-integer component address')
   end subroutine default_integer_address

   subroutine int64_address(value, address)
      integer(int64), target, intent(in) :: value
      integer(MPI_ADDRESS_KIND), intent(out) :: address
      integer :: ierr

      call MPI_GET_ADDRESS(value, address, ierr)
      call require_mpi_success(ierr, 'query an int64 component address')
   end subroutine int64_address

   subroutine real64_address(value, address)
      real(real64), target, intent(in) :: value
      integer(MPI_ADDRESS_KIND), intent(out) :: address
      integer :: ierr

      call MPI_GET_ADDRESS(value, address, ierr)
      call require_mpi_success(ierr, 'query a real64 component address')
   end subroutine real64_address

   subroutine real64_array_address(values, address)
      real(real64), target, intent(in) :: values(:)
      integer(MPI_ADDRESS_KIND), intent(out) :: address
      integer :: ierr

      call MPI_GET_ADDRESS(values(1), address, ierr)
      call require_mpi_success(ierr, 'query a real64-array component address')
   end subroutine real64_array_address

   subroutine require_mpi_success(ierr, operation)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: operation

      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to ' // operation
      endif
   end subroutine require_mpi_success

end module lcm_particle_mpi_datatype
