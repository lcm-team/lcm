program test_particle_count_exchange
   use iso_fortran_env, only: int64, real64
   use mpi_f08, only: MPI_CART_COORDS, MPI_CART_CREATE, MPI_CART_RANK, &
                      MPI_CART_SHIFT, MPI_COMM_FREE, MPI_COMM_RANK, &
                      MPI_COMM_SIZE, MPI_COMM_WORLD, MPI_Comm, MPI_FINALIZE, &
                      MPI_INIT, MPI_PROC_NULL
   use lcm_exchange_directions, only: lcm_exchange_direction_count, &
                                      lcm_exchange_neighbor_rank, &
                                      lcm_exchange_opposite_direction
   use lcm_grid_state, only: lcm_set_grid
   use lcm_parallel_state, only: lcm_init_parallel, lcm_parallel_t, &
                                 lcm_set_parallel
   use lcm_particle_count_exchange, only: lcm_exchange_particle_counts
   use lcm_particle_exchange_buffers, only: lcm_append_send_particle_message, &
                                            lcm_clear_particle_exchange_buffers, &
                                            lcm_expected_receive_particle_message_count
   use lcm_types, only: lcm_grid_t, lcm_init_grid, lcm_particle_mpi_message_t

   implicit none

   type(MPI_Comm) :: communicator
   integer :: coordinates(2)
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
   logical :: square_topology
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
   square_topology = .false.
   if (command_argument_count() > 0) then
      call get_command_argument(1, topology_argument)
      if (trim(topology_argument) == 'nonperiodic') periodic_x = .false.
      if (trim(topology_argument) == 'square') square_topology = .true.
   endif
   if (square_topology) then
      if (size /= 4) error stop 'square count-exchange test requires four ranks'
      dimensions = [2, 2]
      periodicity = [.true., .true.]
   else
      dimensions = [size, 1]
      periodicity = [periodic_x, .true.]
   endif
   call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dimensions, periodicity, .false., &
                        communicator, ierr)
   call MPI_COMM_RANK(communicator, rank, ierr)
   call MPI_CART_COORDS(communicator, rank, 2, coordinates, ierr)
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

   do direction = 1, lcm_exchange_direction_count
      call assert_equal(normalized_neighbor_rank(lcm_exchange_neighbor_rank(direction)), &
                        normalized_neighbor_rank(expected_neighbor_rank(direction)), &
                        'Cartesian neighbor rank')
   enddo

   call lcm_clear_particle_exchange_buffers()
   do direction = 1, lcm_exchange_direction_count
      message = lcm_particle_mpi_message_t()
      message%particle%id = int(10 * rank + direction, int64)
      call append_messages(direction, message, 10 * rank + direction)
   enddo

   call lcm_exchange_particle_counts()

   do direction = 1, lcm_exchange_direction_count
      opposite_direction = lcm_exchange_opposite_direction(direction)
      source_rank = expected_neighbor_rank(opposite_direction)
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

   integer function expected_neighbor_rank(direction) result(neighbor_rank)
      integer, intent(in) :: direction
      integer :: neighbor_coordinates(2)
      integer :: x_offset
      integer :: y_offset

      call direction_offsets(direction, x_offset, y_offset)
      neighbor_coordinates = coordinates + [x_offset, y_offset]

      if (.not. periodicity(1) .and. &
          (neighbor_coordinates(1) < 0 .or. neighbor_coordinates(1) >= dimensions(1))) then
         neighbor_rank = MPI_PROC_NULL
         return
      endif
      if (.not. periodicity(2) .and. &
          (neighbor_coordinates(2) < 0 .or. neighbor_coordinates(2) >= dimensions(2))) then
         neighbor_rank = MPI_PROC_NULL
         return
      endif

      neighbor_coordinates(1) = modulo(neighbor_coordinates(1), dimensions(1))
      neighbor_coordinates(2) = modulo(neighbor_coordinates(2), dimensions(2))
      call MPI_CART_RANK(communicator, neighbor_coordinates, neighbor_rank, ierr)
   end function expected_neighbor_rank

   subroutine direction_offsets(direction, x_offset, y_offset)
      use lcm_exchange_directions, only: lcm_exchange_east, lcm_exchange_west, &
                                         lcm_exchange_southeast, &
                                         lcm_exchange_northwest, &
                                         lcm_exchange_south, lcm_exchange_north, &
                                         lcm_exchange_southwest, &
                                         lcm_exchange_northeast

      integer, intent(in) :: direction
      integer, intent(out) :: x_offset
      integer, intent(out) :: y_offset

      select case (direction)
      case (lcm_exchange_east)
         x_offset = 1
         y_offset = 0
      case (lcm_exchange_west)
         x_offset = -1
         y_offset = 0
      case (lcm_exchange_southeast)
         x_offset = 1
         y_offset = -1
      case (lcm_exchange_northwest)
         x_offset = -1
         y_offset = 1
      case (lcm_exchange_south)
         x_offset = 0
         y_offset = -1
      case (lcm_exchange_north)
         x_offset = 0
         y_offset = 1
      case (lcm_exchange_southwest)
         x_offset = -1
         y_offset = -1
      case (lcm_exchange_northeast)
         x_offset = 1
         y_offset = 1
      case default
         error stop 'invalid direction in count-exchange test'
      end select
   end subroutine direction_offsets

   integer function normalized_neighbor_rank(neighbor_rank) result(normalized_rank)
      integer, intent(in) :: neighbor_rank

      if (neighbor_rank < 0) then
         normalized_rank = -1
      else
         normalized_rank = neighbor_rank
      endif
   end function normalized_neighbor_rank

   subroutine assert_equal(actual, expected, description)
      integer, intent(in) :: actual
      integer, intent(in) :: expected
      character(len=*), intent(in) :: description

      if (actual /= expected) error stop 'assertion failed: ' // description
   end subroutine assert_equal

end program test_particle_count_exchange
