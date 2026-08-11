program test_particle_count_exchange
   use iso_fortran_env, only: int64, real64
   use lcm_exchange_directions, only: lcm_exchange_direction_count, &
                                      lcm_exchange_neighbor_rank, &
                                      lcm_exchange_opposite_direction
   use lcm_grid_state, only: lcm_set_grid
   use lcm_parallel_state, only: lcm_set_parallel
   use lcm_particle_count_exchange, only: lcm_exchange_particle_counts
   use lcm_particle_exchange_buffers, only: lcm_append_send_particle_message, &
                                            lcm_clear_particle_exchange_buffers, &
                                            lcm_expected_receive_particle_message_count
   use lcm_types, only: lcm_grid_t, lcm_init_grid, lcm_init_parallel, &
                        lcm_parallel_t, lcm_particle_mpi_message_t

   implicit none

   include 'mpif.h'

   integer :: communicator
   integer :: dimensions(2)
   integer :: direction
   integer :: east_rank
   integer :: ierr
   integer :: north_rank
   integer :: opposite_direction
   integer :: rank
   integer :: size
   integer :: south_rank
   integer :: source_rank
   character(len=32) :: topology_argument
   logical :: periodic_x
   logical :: periodicity(2)
   integer :: west_rank
   type(lcm_grid_t) :: grid
   type(lcm_parallel_t) :: parallel
   type(lcm_particle_mpi_message_t) :: message
   real(real64) :: dz(1)
   real(real64) :: z_center(1)
   real(real64) :: z_face(2)

   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
   periodic_x = .true.
   if (command_argument_count() > 0) then
      call get_command_argument(1, topology_argument)
      if (trim(topology_argument) == 'nonperiodic') periodic_x = .false.
   endif
   dimensions = [size, 1]
   periodicity = [periodic_x, .true.]
   call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dimensions, periodicity, .false., &
                        communicator, ierr)
   call MPI_COMM_RANK(communicator, rank, ierr)
   call MPI_CART_SHIFT(communicator, 0, 1, west_rank, east_rank, ierr)
   call MPI_CART_SHIFT(communicator, 1, 1, south_rank, north_rank, ierr)

   dz = [1.0_real64]
   z_center = [0.5_real64]
   z_face = [0.0_real64, 1.0_real64]
   call lcm_init_grid(grid, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                      1.0_real64, 1.0_real64, &
                      z_center, z_face, dz, 0, 0, 0, 0, 0, 1, 1, 1, 1, &
                      west_rank, east_rank, south_rank, north_rank, &
                      periodic_x, .true., 1)
   call lcm_set_grid(grid)
   call lcm_init_parallel(parallel, communicator)
   call lcm_set_parallel(parallel)

   call lcm_clear_particle_exchange_buffers()
   do direction = 1, lcm_exchange_direction_count
      message = lcm_particle_mpi_message_t()
      message%particle%id = int(10 * rank + direction, int64)
      call append_messages(direction, message, 10 * rank + direction)
   enddo

   call lcm_exchange_particle_counts()

   do direction = 1, lcm_exchange_direction_count
      opposite_direction = lcm_exchange_opposite_direction(direction)
      source_rank = lcm_exchange_neighbor_rank(opposite_direction)
      if (source_rank < 0) then
         call assert_equal(lcm_expected_receive_particle_message_count(opposite_direction), &
                           0, 'missing neighbor receive count')
      else
         call assert_equal(lcm_expected_receive_particle_message_count(opposite_direction), &
                           10 * source_rank + direction, &
                        'received count from opposite direction')
      endif
   enddo

   call MPI_COMM_FREE(communicator, ierr)
   call MPI_FINALIZE(ierr)

contains

   subroutine append_messages(direction, message, count)
      integer, intent(in) :: direction
      integer, intent(in) :: count
      type(lcm_particle_mpi_message_t), intent(in) :: message
      integer :: n

      do n = 1, count
         call lcm_append_send_particle_message(direction, message)
      enddo
   end subroutine append_messages

   subroutine assert_equal(actual, expected, description)
      integer, intent(in) :: actual
      integer, intent(in) :: expected
      character(len=*), intent(in) :: description

      if (actual /= expected) error stop 'assertion failed: ' // description
   end subroutine assert_equal

end program test_particle_count_exchange
