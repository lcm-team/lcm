!> Public interface for host models coupled to LCM.
module lcm_host_interface
   use lcm_types, only: lcm_grid_t, lcm_config_t, lcm_init_grid, lcm_init_config
   use lcm_config_state, only: lcm_set_config
   use lcm_grid_state, only: lcm_set_grid
   use lcm_parallel_state, only: lcm_parallel_t, lcm_init_parallel, &
                                 lcm_set_parallel
   use lcm_initialization_state, only: lcm_configure_runtime, lcm_initialize
   use lcm_advance_mod, only: lcm_advance
   use lcm_field_state, only: lcm_attach_field, lcm_attach_fields
   use lcm_particle_exchange_metrics, only: lcm_get_particle_exchange_metrics, &
      lcm_particle_exchange_metrics_t, &
      lcm_reset_particle_exchange_metrics, &
      lcm_set_particle_exchange_metrics_enabled
   use lcm_particle_storage_metrics, only: lcm_get_particle_storage_metrics, &
      lcm_particle_storage_metrics_t

   implicit none

   private

   public :: lcm_grid_t
   public :: lcm_parallel_t
   public :: lcm_config_t
   public :: lcm_init_grid
   public :: lcm_init_parallel
   public :: lcm_init_config
   public :: lcm_set_config
   public :: lcm_set_grid
   public :: lcm_set_parallel
   public :: lcm_configure_runtime
   public :: lcm_initialize
   public :: lcm_advance
   public :: lcm_attach_field
   public :: lcm_attach_fields
   public :: lcm_particle_exchange_metrics_t
   public :: lcm_set_particle_exchange_metrics_enabled
   public :: lcm_reset_particle_exchange_metrics
   public :: lcm_get_particle_exchange_metrics
   public :: lcm_particle_storage_metrics_t
   public :: lcm_get_particle_storage_metrics

end module lcm_host_interface
