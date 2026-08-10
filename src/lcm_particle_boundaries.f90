!> Physical particle-boundary rules for standalone LCM.
module lcm_particle_boundaries
   use iso_fortran_env, only: real64
   use lcm_config_state, only: config_state
   use lcm_grid_state, only: grid_state
   use lcm_particle_state, only: lcm_apply_particles_at_bottom, &
                                 lcm_reflect_particles_at_top

   implicit none

   private

   public :: lcm_apply_particle_boundaries

contains

   !> Apply physical particle-boundary rules after particle motion.
   !>
   !> This routine deliberately excludes MPI exchange and final cell ownership.
   !> The implementation preserves SAM-LCM's top and bottom rules.
   subroutine lcm_apply_particle_boundaries()
      call lcm_reflect_particles_at_top(grid_state%ztop)
      call lcm_apply_particles_at_bottom(grid_state%zbottom, &
                                         config_state%bottom_removal_radius, &
                                         real(config_state%rho_aero, real64))
   end subroutine lcm_apply_particle_boundaries

end module lcm_particle_boundaries
