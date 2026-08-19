program test_particle_migration_metadata
   use lcm_exchange_directions, only: lcm_exchange_east, lcm_exchange_west
   use lcm_particle_migration, only: lcm_migration_local_cell, &
                                      lcm_migration_outgoing, &
                                      lcm_migration_same_cell, &
                                      lcm_migration_unclassified, &
                                      lcm_particle_migration_is_outgoing, &
                                      lcm_particle_migration_is_valid, &
                                      lcm_particle_migration_t

   implicit none

   type(lcm_particle_migration_t) :: migration

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

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         write(*,'(a)') 'FAIL: ' // trim(message)
         error stop 1
      endif
   end subroutine assert_true

end program test_particle_migration_metadata
