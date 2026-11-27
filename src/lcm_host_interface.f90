!> Public interface for host models coupled to LCM.
module lcm_host_interface
   use lcm_types, only: lcm_grid_t, lcm_config_t, lcm_variable_t, lcm_fields_t, &
                         lcm_init_grid, lcm_init_config
   use lcm_config_state, only: lcm_set_config
   use lcm_grid_state, only: lcm_set_grid
   use lcm_parallel_state, only: lcm_parallel_t, lcm_init_parallel, &
                                 lcm_set_parallel
   use lcm_initialization_state, only: lcm_configure_runtime, lcm_initialize, &
                                       lcm_pass_timer_routines
   use lcm_advance_mod, only: lcm_advance
   use lcm_field_state, only: lcm_set_field
   use lcm_microphysics_state, only: lcm_set_thermodynamic_history

   use lcm_sgs_mixing, only: lcm_set_lem_environment
   implicit none

   private

   public :: lcm_grid_t
   public :: lcm_parallel_t
   public :: lcm_config_t
   public :: lcm_variable_t
   public :: lcm_fields_t
   public :: lcm_init_grid
   public :: lcm_init_parallel
   public :: lcm_init_config
   public :: lcm_set_config
   public :: lcm_set_grid
   public :: lcm_set_parallel
   public :: lcm_configure_runtime
   public :: lcm_initialize
   public :: lcm_pass_timer_routines
   public :: lcm_set_lem_environment
   public :: lcm_advance
   public :: lcm_set_field
   public :: lcm_set_thermodynamic_history

end module lcm_host_interface
