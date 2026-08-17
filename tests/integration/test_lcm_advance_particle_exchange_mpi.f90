program test_lcm_advance_particle_exchange_mpi
   use iso_fortran_env, only: int64, real64
   use mpi_f08, only: MPI_ALLGATHER, MPI_ALLREDUCE, MPI_CART_COORDS, &
                      MPI_CART_CREATE, MPI_CART_RANK, &
                      MPI_COMM_FREE, MPI_COMM_RANK, MPI_COMM_SIZE, &
                      MPI_COMM_WORLD, MPI_INTEGER, MPI_INTEGER8, &
                      MPI_SUM, MPI_Comm, MPI_FINALIZE, MPI_INIT
   use lcm_host_interface, only: lcm_advance, lcm_attach_field, &
                                 lcm_attach_fields, lcm_config_t, &
                                 lcm_get_particle_exchange_metrics, &
                                 lcm_grid_t, lcm_init_config, &
                                 lcm_init_grid, lcm_init_parallel, &
                                 lcm_initialize, lcm_parallel_t, &
                                 lcm_particle_exchange_metrics_t, &
                                 lcm_reset_particle_exchange_metrics, &
                                 lcm_set_config, lcm_set_grid, &
                                 lcm_set_parallel, &
                                 lcm_set_particle_exchange_metrics_enabled
   use lcm_initialization_state, only: lcm_initial_particle_count_global, &
                                       lcm_initial_particle_count_local
   use lcm_particle_location, only: lcm_locate_particle_cell
   use lcm_particle_state, only: cell_particles

   implicit none

   integer, parameter :: nx_local = 2
   integer, parameter :: ny_local = 2
   integer, parameter :: nz = 1
   integer, parameter :: np_xy = 1
   integer, parameter :: expected_local_count = nx_local * ny_local * np_xy
   integer, parameter :: expected_global_count = 4 * expected_local_count
   real(real64), parameter :: dx = 10.0_real64
   real(real64), parameter :: dy = 10.0_real64
   real(real64), parameter :: time_step = 1.0_real64
   real(real64), parameter :: subdomain_x_width = real(nx_local, real64) * dx
   real(real64), parameter :: subdomain_y_width = real(ny_local, real64) * dy

   type(MPI_Comm) :: communicator
   type(lcm_config_t) :: config
   integer :: coordinates(2)
   integer :: dimensions(2)
   real(real64) :: dz(nz)
   type(lcm_grid_t) :: grid
   integer :: ierr
   integer :: mpi_size
   integer(int64) :: initial_ids_by_rank(expected_local_count, 4)
   type(lcm_parallel_t) :: parallel
   type(lcm_particle_exchange_metrics_t) :: metrics
   logical :: periodicity(2)
   integer :: rank
   real(real64), target :: density(nz)
   real(real64), target :: pressure(nz)
   real(real64), target :: qv(nx_local, ny_local, nz)
   real(real64), target :: static_energy(nx_local, ny_local, nz)
   real(real64), target :: temperature(nx_local, ny_local, nz)
   real(real64), target :: u(nx_local + 1, ny_local, nz)
   real(real64), target :: v(nx_local, ny_local + 1, nz)
   real(real64), target :: w(nx_local, ny_local, nz + 1)
   real(real64) :: z_center(nz)
   real(real64) :: z_face(nz + 1)

   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
   if (mpi_size /= 4) then
      error stop 'lcm_advance particle-exchange test requires four ranks'
   endif

   dimensions = [2, 2]
   periodicity = [.true., .true.]
   call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dimensions, periodicity, .false., &
                        communicator, ierr)
   call MPI_COMM_RANK(communicator, rank, ierr)
   call MPI_CART_COORDS(communicator, rank, 2, coordinates, ierr)

   call configure_lcm()
   call attach_host_fields()
   call lcm_initialize()

   call assert_equal_int64(lcm_initial_particle_count_local, &
                           int(expected_local_count, int64), &
                           'initialized local particle count')
   call assert_equal_int64(lcm_initial_particle_count_global, &
                           int(expected_global_count, int64), &
                           'initialized global particle count')
   call assert_equal(total_particle_count(), expected_local_count, &
                     'cell-owned count before advance')
   call gather_initial_particle_ids()

   call lcm_reset_particle_exchange_metrics()
   call lcm_set_particle_exchange_metrics_enabled(.true.)
   call lcm_advance(time_step)

   call assert_global_particle_conservation()
   call assert_particles_arrived_from_diagonal_rank()
   call assert_final_cell_ownership()
   call lcm_get_particle_exchange_metrics(metrics)
   call assert_equal_int64(metrics%ownership_update_calls, 1_int64, &
                           'one ownership update through lcm_advance')
   call assert_equal_int64(metrics%ownership_particle_examinations, &
                           int(expected_local_count, int64), &
                           'all final local particles examined')
   call assert_equal_int64(metrics%ownership_particles_moved, 0_int64, &
                           'received particles already own final cells')
   call assert_true(metrics%ownership_update_seconds >= 0.0_real64, &
                    'lcm_advance ownership update time is nonnegative')

   call MPI_COMM_FREE(communicator, ierr)
   call MPI_FINALIZE(ierr)

contains

   subroutine configure_lcm()
      call lcm_init_parallel(parallel, communicator)
      call lcm_set_parallel(parallel)

      dz = [10.0_real64]
      z_center = [5.0_real64]
      z_face = [0.0_real64, 10.0_real64]
      call lcm_init_grid(grid, nx_local, ny_local, nz, &
                         dimensions(1) * nx_local, &
                         dimensions(2) * ny_local, &
                         1, nx_local, 1, ny_local, dx, dy, &
                         z_center, z_face, dz, 0, 0, 0, &
                         coordinates(1) * nx_local + 1, &
                         coordinates(2) * ny_local + 1, 1)
      call lcm_set_grid(grid)

      call lcm_init_config(config)
      config%np_xy = np_xy
      config%nzl = 1
      config%n_uvw_iteration = 0
      config%switch_sgs_velocities = .false.
      config%switch_sedimentation = .false.
      config%switch_collection = .false.
      config%switch_breakup = .false.
      config%switch_phase_sources = .false.
      config%switch_particle_source = .false.
      call lcm_set_config(config)
   end subroutine configure_lcm

   subroutine attach_host_fields()
      u = subdomain_x_width / time_step
      v = subdomain_y_width / time_step
      w = 0.0_real64
      static_energy = 300.0_real64
      temperature = 280.0_real64
      qv = 0.01_real64
      pressure = 100000.0_real64
      density = 1.2_real64

      call lcm_attach_field('u', u)
      call lcm_attach_field('v', v)
      call lcm_attach_field('w', w)
      call lcm_attach_fields(static_energy=static_energy, &
                             temperature=temperature, qv=qv, &
                             pressure=pressure, density=density)
   end subroutine attach_host_fields

   subroutine assert_global_particle_conservation()
      integer :: global_count
      integer :: local_count

      local_count = total_particle_count()
      call assert_equal(local_count, expected_local_count, &
                        'local particle count after advance')
      call MPI_ALLREDUCE(local_count, global_count, 1, MPI_INTEGER, MPI_SUM, &
                         communicator, ierr)
      call assert_equal(global_count, expected_global_count, &
                        'global particle count after advance')
   end subroutine assert_global_particle_conservation

   subroutine gather_initial_particle_ids()
      integer :: i
      integer :: j
      integer(int64) :: local_ids(expected_local_count)
      integer :: n

      n = 0
      do j = 1, ny_local
         do i = 1, nx_local
            n = n + 1
            local_ids(n) = cell_particles(i,j,1)%particles(1)%id
         enddo
      enddo

      call MPI_ALLGATHER(local_ids, expected_local_count, MPI_INTEGER8, &
                         initial_ids_by_rank, expected_local_count, MPI_INTEGER8, &
                         communicator, ierr)
   end subroutine gather_initial_particle_ids

   subroutine assert_particles_arrived_from_diagonal_rank()
      integer :: expected_source_coordinates(2)
      integer :: expected_source_rank
      integer :: i
      integer :: j
      integer :: n

      expected_source_coordinates = modulo(coordinates - 1, dimensions)
      call MPI_CART_RANK(communicator, expected_source_coordinates, &
                         expected_source_rank, ierr)

      do j = 1, ny_local
         do i = 1, nx_local
            call assert_equal(cell_particles(i,j,1)%count, np_xy, &
                              'one initialized particle per final cell')
            do n = 1, cell_particles(i,j,1)%count
               call assert_equal_int64(cell_particles(i,j,1)%particles(n)%id, &
                  initial_ids_by_rank((j - 1) * nx_local + i, &
                                      expected_source_rank + 1), &
                  'particle retained its source-cell identity')
            enddo
         enddo
      enddo
   end subroutine assert_particles_arrived_from_diagonal_rank

   subroutine assert_final_cell_ownership()
      integer :: i
      integer :: j
      integer :: k
      integer :: located_i
      integer :: located_j
      integer :: located_k
      integer :: n
      logical :: is_local
      real(real64) :: x
      real(real64) :: y

      do k = 1, nz
         do j = 1, ny_local
            do i = 1, nx_local
               do n = 1, cell_particles(i,j,k)%count
                  x = cell_particles(i,j,k)%particles(n)%x
                  y = cell_particles(i,j,k)%particles(n)%y
                  call assert_true(x >= 0.0_real64 .and. &
                                   x < real(dimensions(1) * nx_local, real64) * dx, &
                                   'particle x is wrapped into global domain')
                  call assert_true(y >= 0.0_real64 .and. &
                                   y < real(dimensions(2) * ny_local, real64) * dy, &
                                   'particle y is wrapped into global domain')
                  call lcm_locate_particle_cell( &
                     x, y, cell_particles(i,j,k)%particles(n)%z, &
                     located_i, located_j, located_k, is_local)
                  call assert_true(is_local, &
                                   'received particle is local after both phases')
                  call assert_equal(located_i, i, &
                                    'received particle x cell ownership')
                  call assert_equal(located_j, j, &
                                    'received particle y cell ownership')
                  call assert_equal(located_k, k, &
                                    'received particle vertical cell ownership')
                  call assert_close(cell_particles(i,j,k)%particles(n)%u, &
                                    subdomain_x_width / time_step, &
                                    'transported particle u velocity')
                  call assert_close(cell_particles(i,j,k)%particles(n)%v, &
                                    subdomain_y_width / time_step, &
                                    'transported particle v velocity')
               enddo
            enddo
         enddo
      enddo
   end subroutine assert_final_cell_ownership

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

   subroutine assert_equal_int64(actual, expected, description)
      integer(int64), intent(in) :: actual
      integer(int64), intent(in) :: expected
      character(len=*), intent(in) :: description

      if (actual /= expected) error stop 'assertion failed: ' // description
   end subroutine assert_equal_int64

   subroutine assert_close(actual, expected, description)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: description

      if (abs(actual - expected) > 1.0e-12_real64) then
         error stop 'assertion failed: ' // description
      endif
   end subroutine assert_close

end program test_lcm_advance_particle_exchange_mpi
