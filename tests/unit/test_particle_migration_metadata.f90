program test_particle_migration_metadata
   use iso_fortran_env, only: int64, real64
   use lcm_exchange_directions, only: lcm_exchange_east, lcm_exchange_north, &
                                      lcm_exchange_south, lcm_exchange_west
   use lcm_grid_state, only: lcm_set_grid
   use lcm_particle_migration, only: lcm_classify_particle_migration, &
                                      lcm_migration_local_cell, &
                                      lcm_migration_outgoing, &
                                      lcm_migration_same_cell, &
                                      lcm_migration_unclassified, &
                                      lcm_particle_migration_is_outgoing, &
                                      lcm_particle_migration_is_valid, &
                                      lcm_particle_migration_t
   use lcm_types, only: lcm_grid_t, lcm_init_grid, lcm_particle_t

   implicit none

   type(lcm_particle_migration_t) :: migration

   call setup_grid()
   call assert_true(migration%disposition == lcm_migration_unclassified, &
                    'migration defaults to unclassified')
   call assert_true(lcm_particle_migration_is_valid(migration), &
                    'default migration metadata is valid')

   migration = lcm_particle_migration_t(7, 8, 9, lcm_migration_same_cell, 0)
   call assert_true(lcm_particle_migration_is_valid(migration), &
                    'same-cell migration has a local destination')
   call assert_true(.not. lcm_particle_migration_is_outgoing(migration), &
                    'same-cell migration is not outgoing')

   migration%disposition = lcm_migration_local_cell
   call assert_true(lcm_particle_migration_is_valid(migration), &
                    'local-cell migration has a local destination')

   migration = lcm_particle_migration_t(0, 8, 9, lcm_migration_outgoing, &
                                        lcm_exchange_west)
   call assert_true(lcm_particle_migration_is_valid(migration), &
                    'outgoing migration may precede periodic wrapping')
   call assert_true(lcm_particle_migration_is_outgoing(migration), &
                    'west migration is outgoing')

   migration%direction = 0
   call assert_true(.not. lcm_particle_migration_is_valid(migration), &
                    'outgoing migration requires a direction')

   migration = lcm_particle_migration_t(7, 8, 9, lcm_migration_local_cell, &
                                        lcm_exchange_east)
   call assert_true(.not. lcm_particle_migration_is_valid(migration), &
                    'local migration cannot have an MPI direction')

   call test_same_cell()
   call test_local_cell_moves()
   call test_cardinal_routes()
   call test_corner_route_precedence()

contains

   subroutine setup_grid()
      integer, parameter :: nx_local = 2
      integer, parameter :: ny_local = 2
      integer, parameter :: nz = 2
      real(real64), parameter :: z_center(nz) = [50.0_real64, 160.0_real64]
      real(real64), parameter :: z_face(nz + 1) = &
         [0.0_real64, 100.0_real64, 220.0_real64]
      real(real64), parameter :: dz(nz) = [100.0_real64, 120.0_real64]
      type(lcm_grid_t) :: grid

      call lcm_init_grid(grid, nx_local, ny_local, nz, 6, 6, &
                         1, nx_local, 1, ny_local, &
                         50.0_real64, 50.0_real64, z_center, z_face, dz, &
                         0, 0, 0, 3, 3, 1)
      call lcm_set_grid(grid)
   end subroutine setup_grid

   subroutine test_same_cell()
      type(lcm_particle_t) :: particle

      particle = particle_at(125.0_real64, 125.0_real64, 50.0_real64)
      particle%id = 41_int64

      call lcm_classify_particle_migration(particle, 1, 1, 1, migration)

      call assert_migration(migration, 3, 3, 1, lcm_migration_same_cell, 0, &
                            'same cell')
      call assert_true(particle%id == 41_int64, 'classifier preserves particle id')
      call assert_true(particle%x == 125.0_real64, &
                       'classifier preserves particle position')
   end subroutine test_same_cell

   subroutine test_local_cell_moves()
      type(lcm_particle_t) :: particle

      particle = particle_at(175.0_real64, 125.0_real64, 50.0_real64)
      call lcm_classify_particle_migration(particle, 1, 1, 1, migration)
      call assert_migration(migration, 4, 3, 1, lcm_migration_local_cell, 0, &
                            'local x move')

      particle = particle_at(125.0_real64, 175.0_real64, 50.0_real64)
      call lcm_classify_particle_migration(particle, 1, 1, 1, migration)
      call assert_migration(migration, 3, 4, 1, lcm_migration_local_cell, 0, &
                            'local y move')

      particle = particle_at(125.0_real64, 125.0_real64, 150.0_real64)
      call lcm_classify_particle_migration(particle, 1, 1, 1, migration)
      call assert_migration(migration, 3, 3, 2, lcm_migration_local_cell, 0, &
                            'local z move')
   end subroutine test_local_cell_moves

   subroutine test_cardinal_routes()
      call assert_particle_route(99.0_real64, 125.0_real64, &
                                 lcm_exchange_west, 2, 3, 'west route')
      call assert_particle_route(200.0_real64, 125.0_real64, &
                                 lcm_exchange_east, 5, 3, 'east route')
      call assert_particle_route(125.0_real64, 99.0_real64, &
                                 lcm_exchange_south, 3, 2, 'south route')
      call assert_particle_route(125.0_real64, 200.0_real64, &
                                 lcm_exchange_north, 3, 5, 'north route')
   end subroutine test_cardinal_routes

   subroutine test_corner_route_precedence()
      call assert_particle_route(99.0_real64, 99.0_real64, &
                                 lcm_exchange_west, 2, 2, 'southwest route')
      call assert_particle_route(99.0_real64, 200.0_real64, &
                                 lcm_exchange_west, 2, 5, 'northwest route')
      call assert_particle_route(200.0_real64, 99.0_real64, &
                                 lcm_exchange_east, 5, 2, 'southeast route')
      call assert_particle_route(200.0_real64, 200.0_real64, &
                                 lcm_exchange_east, 5, 5, 'northeast route')
   end subroutine test_corner_route_precedence

   subroutine assert_particle_route(x, y, direction, global_i, global_j, message)
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      integer, intent(in) :: direction
      integer, intent(in) :: global_i
      integer, intent(in) :: global_j
      character(len=*), intent(in) :: message
      type(lcm_particle_t) :: particle

      particle = particle_at(x, y, 50.0_real64)
      call lcm_classify_particle_migration(particle, 1, 1, 1, migration)
      call assert_migration(migration, global_i, global_j, 1, &
                            lcm_migration_outgoing, direction, message)
   end subroutine assert_particle_route

   subroutine assert_migration(actual, global_i, global_j, k, disposition, &
                               direction, message)
      type(lcm_particle_migration_t), intent(in) :: actual
      integer, intent(in) :: global_i
      integer, intent(in) :: global_j
      integer, intent(in) :: k
      integer, intent(in) :: disposition
      integer, intent(in) :: direction
      character(len=*), intent(in) :: message

      call assert_true(lcm_particle_migration_is_valid(actual), &
                       message // ' validity')
      call assert_true(actual%global_i == global_i, message // ' global i')
      call assert_true(actual%global_j == global_j, message // ' global j')
      call assert_true(actual%k == k, message // ' k')
      call assert_true(actual%disposition == disposition, &
                       message // ' disposition')
      call assert_true(actual%direction == direction, message // ' direction')
   end subroutine assert_migration

   function particle_at(x, y, z) result(particle)
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      real(real64), intent(in) :: z
      type(lcm_particle_t) :: particle

      particle%x = x
      particle%y = y
      particle%z = z
   end function particle_at

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         write(*,'(a)') 'FAIL: ' // trim(message)
         error stop 1
      endif
   end subroutine assert_true

end program test_particle_migration_metadata
