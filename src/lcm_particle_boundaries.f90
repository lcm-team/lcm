!> Physical particle-boundary rules for standalone LCM.
module lcm_particle_boundaries
   use lcm_grid_state, only: grid_state
   use lcm_particle_state, only: lcm_reflect_particles_at_top

   implicit none

   private

   public :: lcm_apply_particle_boundaries

contains

   !> Apply physical particle-boundary rules after particle motion.
   !>
   !> This routine deliberately excludes MPI exchange and final cell ownership.
   !> The initial implementation preserves SAM-LCM's top reflection rule;
   !> bottom-boundary behavior is added separately.
   subroutine lcm_apply_particle_boundaries()
      call lcm_reflect_particles_at_top(grid_state%ztop)
   end subroutine lcm_apply_particle_boundaries

end module lcm_particle_boundaries
