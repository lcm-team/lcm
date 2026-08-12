program test_periodic_particle_destination_wrapping
   use iso_fortran_env, only: int64, real64
   use mpi_f08, only: MPI_CART_CREATE, MPI_COMM_FREE, MPI_COMM_WORLD, MPI_Comm, &
                      MPI_FINALIZE, MPI_INIT
   use lcm_exchange_directions, only: lcm_exchange_east, lcm_exchange_north, &
                                      lcm_exchange_south, lcm_exchange_west, &
                                      lcm_exchange_x_phase, lcm_exchange_y_phase
   use lcm_grid_state, only: lcm_set_grid
   use lcm_parallel_state, only: lcm_init_parallel, lcm_parallel_t, &
                                 lcm_set_parallel
   use lcm_particle_exchange, only: lcm_exchange_particle_phase
   use lcm_particle_exchange_buffers, only: lcm_get_send_particle_message, &
                                            lcm_send_particle_message_count
   use lcm_particle_state, only: lcm_add_particle, lcm_allocate_particle_cells, &
                                 cell_particles
   use lcm_types, only: lcm_grid_t, lcm_init_grid, lcm_particle_mpi_message_t, &
                        lcm_particle_t

   implicit none

   type(lcm_grid_t) :: grid
   type(lcm_parallel_t) :: parallel
   type(MPI_Comm) :: communicator
   integer :: dimensions(2)
   integer :: ierr
   logical :: periodicity(2)
   type(lcm_particle_mpi_message_t) :: message
   real(real64) :: dz(1)
   real(real64) :: z_center(1)
   real(real64) :: z_face(2)

   dz = [10.0_real64]
   z_center = [5.0_real64]
   z_face = [0.0_real64, 10.0_real64]
   call MPI_INIT(ierr)
   dimensions = [1, 1]
   periodicity = [.true., .true.]
   call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dimensions, periodicity, .false., &
                        communicator, ierr)
   call lcm_init_parallel(parallel, communicator)
   call lcm_set_parallel(parallel)
   call lcm_init_grid(grid, 3, 2, 1, 3, 2, &
                      2, 4, 2, 3, 10.0_real64, 10.0_real64, &
                      z_center, z_face, dz, 1, 1, 0, 1, 1, 1)
   call lcm_set_grid(grid)
   call lcm_allocate_particle_cells()

   call add_particle(601_int64, -1.0_real64, 5.0_real64)
   call add_particle(602_int64, 31.0_real64, 5.0_real64)
   call add_particle(603_int64, 5.0_real64, -1.0_real64)
   call add_particle(604_int64, 5.0_real64, 21.0_real64)
   call add_particle(605_int64, 31.0_real64, 21.0_real64)
   call add_particle(606_int64, -1.0_real64, 21.0_real64)
   call add_particle(607_int64, 31.0_real64, -1.0_real64)
   call add_particle(608_int64, -1.0_real64, -1.0_real64)

   call lcm_exchange_particle_phase(lcm_exchange_x_phase)

   call assert_wrapped_message(lcm_exchange_west, 601, 3, 1, &
                               29.0_real64, 5.0_real64, 'west-to-east')
   call assert_wrapped_message(lcm_exchange_east, 602, 1, 1, &
                               1.0_real64, 5.0_real64, 'east-to-west')

   ! The x phase wraps x only. Corner particles retain their out-of-range y
   ! destination until the receiving rank classifies them in the y phase.
   call assert_wrapped_message(lcm_exchange_east, 605, 1, 3, &
                               1.0_real64, 21.0_real64, 'north-east x phase')
   call assert_wrapped_message(lcm_exchange_west, 606, 3, 3, &
                               29.0_real64, 21.0_real64, 'north-west x phase')
   call assert_wrapped_message(lcm_exchange_east, 607, 1, 0, &
                               1.0_real64, -1.0_real64, 'south-east x phase')
   call assert_wrapped_message(lcm_exchange_west, 608, 3, 0, &
                               29.0_real64, -1.0_real64, 'south-west x phase')

   call assert_equal(total_particle_count(), 8, &
                     'x payload transport reinserts every periodic particle')
   call lcm_exchange_particle_phase(lcm_exchange_y_phase)

   call assert_wrapped_message(lcm_exchange_south, 603, 1, 2, &
                               5.0_real64, 19.0_real64, 'south-to-north')
   call assert_wrapped_message(lcm_exchange_north, 604, 1, 1, &
                               5.0_real64, 1.0_real64, 'north-to-south')
   call assert_wrapped_message(lcm_exchange_north, 605, 1, 1, &
                               1.0_real64, 1.0_real64, 'north-east corner')
   call assert_wrapped_message(lcm_exchange_north, 606, 3, 1, &
                               29.0_real64, 1.0_real64, 'north-west corner')
   call assert_wrapped_message(lcm_exchange_south, 607, 1, 2, &
                               1.0_real64, 19.0_real64, 'south-east corner')
   call assert_wrapped_message(lcm_exchange_south, 608, 3, 2, &
                               29.0_real64, 19.0_real64, 'south-west corner')
   call assert_equal(total_particle_count(), 8, &
                     'y payload transport preserves every periodic particle')

   call MPI_COMM_FREE(communicator, ierr)
   call MPI_FINALIZE(ierr)

contains

   subroutine add_particle(id, x, y)
      integer(int64), intent(in) :: id
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      type(lcm_particle_t) :: particle

      particle = lcm_particle_t()
      particle%id = id
      particle%x = x
      particle%y = y
      particle%z = 5.0_real64
      call lcm_add_particle(1, 1, 1, particle)
   end subroutine add_particle

   integer function total_particle_count() result(count)
      integer :: i
      integer :: j
      integer :: k

      count = 0
      do k = 1, size(cell_particles, 3)
         do j = 1, size(cell_particles, 2)
            do i = 1, size(cell_particles, 1)
               count = count + cell_particles(i,j,k)%count
            enddo
         enddo
      enddo
   end function total_particle_count

   subroutine assert_wrapped_message(direction, expected_id, expected_i, expected_j, &
                                     expected_x, expected_y, description)
      integer, intent(in) :: direction
      integer, intent(in) :: expected_id
      integer, intent(in) :: expected_i
      integer, intent(in) :: expected_j
      real(real64), intent(in) :: expected_x
      real(real64), intent(in) :: expected_y
      character(len=*), intent(in) :: description

      integer :: index
      logical :: found

      found = .false.
      do index = 1, lcm_send_particle_message_count(direction)
         call lcm_get_send_particle_message(direction, index, message)
         if (message%particle%id == int(expected_id, int64)) then
            found = .true.
            exit
         endif
      enddo
      call assert_true(found, trim(description) // ' message found')
      call assert_equal(int(message%particle%id), expected_id, &
                        trim(description) // ' particle id')
      call assert_equal(message%global_i, expected_i, trim(description) // ' global i')
      call assert_equal(message%global_j, expected_j, trim(description) // ' global j')
      call assert_close(message%particle%x, expected_x, trim(description) // ' x')
      call assert_close(message%particle%y, expected_y, trim(description) // ' y')
   end subroutine assert_wrapped_message

   subroutine assert_true(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description

      if (.not. condition) error stop 'assertion failed: ' // description
   end subroutine assert_true

   subroutine assert_equal(actual, expected, description)
      integer, intent(in) :: actual
      integer, intent(in) :: expected
      character(len=*), intent(in) :: description

      if (actual /= expected) error stop 'assertion failed: ' // description
   end subroutine assert_equal

   subroutine assert_close(actual, expected, description)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: description

      if (abs(actual - expected) > 1.0e-12_real64) then
         error stop 'assertion failed: ' // description
      endif
   end subroutine assert_close

end program test_periodic_particle_destination_wrapping
