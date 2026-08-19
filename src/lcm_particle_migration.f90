!> Temporary particle-migration metadata.
!>
!> A migration record is produced after physical particle-boundary handling
!> and consumed while compacting local cells, constructing MPI send buffers,
!> and inserting particles into destination cells.  It is discarded before
!> microphysics and is deliberately separate from lcm_particle_t.
!>
!> The migration planner gives horizontal MPI routing x precedence: a
!> particle outside x is routed west/east first, even when it is also outside
!> y.  Only particles local in x may be routed south/north in the initial
!> plan.  Particles received by the x phase are classified separately for the
!> y phase, preserving the existing two-phase corner-routing semantics.
module lcm_particle_migration
   use lcm_exchange_directions, only: lcm_exchange_direction_count

   implicit none

   private

   integer, parameter, public :: lcm_migration_unclassified = 0
   integer, parameter, public :: lcm_migration_same_cell = 1
   integer, parameter, public :: lcm_migration_local_cell = 2
   integer, parameter, public :: lcm_migration_outgoing = 3

   !> Destination and route for one particle during a single advance step.
   !>
   !> global_i, global_j, and k are temporary ownership metadata.  For an
   !> outgoing particle, direction is one of the cardinal exchange-direction
   !> identifiers.  Non-outgoing records have direction zero.
   type, public :: lcm_particle_migration_t
      integer :: global_i = 0
      integer :: global_j = 0
      integer :: k = 0
      integer :: disposition = lcm_migration_unclassified
      integer :: direction = 0
   end type lcm_particle_migration_t

   public :: lcm_particle_migration_is_outgoing
   public :: lcm_particle_migration_is_valid

contains

   pure logical function lcm_particle_migration_is_outgoing(migration)
      type(lcm_particle_migration_t), intent(in) :: migration

      lcm_particle_migration_is_outgoing = &
         migration%disposition == lcm_migration_outgoing
   end function lcm_particle_migration_is_outgoing

   pure logical function lcm_particle_migration_is_valid(migration)
      type(lcm_particle_migration_t), intent(in) :: migration

      select case (migration%disposition)
      case (lcm_migration_unclassified)
         lcm_particle_migration_is_valid = migration%direction == 0
      case (lcm_migration_same_cell, lcm_migration_local_cell)
         lcm_particle_migration_is_valid = migration%global_i >= 1 .and. &
                                           migration%global_j >= 1 .and. &
                                           migration%k >= 1 .and. &
                                           migration%direction == 0
      case (lcm_migration_outgoing)
         lcm_particle_migration_is_valid = migration%k >= 1 .and. &
                                           migration%direction >= 1 .and. &
                                           migration%direction <= &
                                              lcm_exchange_direction_count
      case default
         lcm_particle_migration_is_valid = .false.
      end select
   end function lcm_particle_migration_is_valid

end module lcm_particle_migration
