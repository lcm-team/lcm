!> Explicit cross-rank summaries for optional particle-exchange timings.
!!
!! Calling lcm_get_particle_exchange_timing_summary is collective when an
!! MPI parallel state is configured. It is deliberately separate from
!! lcm_advance so no global collective is added to the timestep path.
module lcm_particle_exchange_timing_summary
   use iso_fortran_env, only: real64
   use mpi_f08, only: MPI_ALLREDUCE, MPI_Datatype, MPI_MAX, MPI_MIN, &
                      MPI_SUCCESS, MPI_SUM, MPI_TYPECLASS_REAL, &
                      MPI_TYPE_MATCH_SIZE
   use lcm_parallel_state, only: parallel_state
   use lcm_particle_exchange_metrics, only: &
      lcm_get_particle_exchange_metrics, lcm_particle_exchange_metrics_t

   implicit none

   private

   integer, parameter :: timing_count = 22
   integer, parameter :: x_outgoing_preparation_index = 1
   integer, parameter :: y_outgoing_preparation_index = 2
   integer, parameter :: x_count_exchange_index = 3
   integer, parameter :: y_count_exchange_index = 4
   integer, parameter :: x_payload_exchange_index = 5
   integer, parameter :: y_payload_exchange_index = 6
   integer, parameter :: x_received_insertion_index = 7
   integer, parameter :: y_received_insertion_index = 8
   integer, parameter :: ownership_update_index = 9
   integer, parameter :: advance_total_index = 10
   integer, parameter :: prepare_particle_velocities_index = 11
   integer, parameter :: refine_particle_velocities_index = 12
   integer, parameter :: update_particle_positions_index = 13
   integer, parameter :: particle_boundaries_index = 14
   integer, parameter :: migration_total_index = 15
   integer, parameter :: migration_buffer_clear_index = 16
   integer, parameter :: migration_plan_allocate_index = 17
   integer, parameter :: migration_plan_classification_index = 18
   integer, parameter :: migration_plan_apply_index = 19
   integer, parameter :: migration_plan_deallocate_index = 20
   integer, parameter :: particle_exchange_total_index = 21
   integer, parameter :: particle_sources_index = 22

   type, public :: lcm_rank_timing_summary_t
      real(real64) :: minimum_seconds = 0.0_real64
      real(real64) :: maximum_seconds = 0.0_real64
      real(real64) :: mean_seconds = 0.0_real64
   end type lcm_rank_timing_summary_t

   type, public :: lcm_particle_exchange_timing_summary_t
      integer :: rank_count = 1
      type(lcm_rank_timing_summary_t) :: x_outgoing_preparation
      type(lcm_rank_timing_summary_t) :: y_outgoing_preparation
      type(lcm_rank_timing_summary_t) :: x_count_exchange
      type(lcm_rank_timing_summary_t) :: y_count_exchange
      type(lcm_rank_timing_summary_t) :: x_payload_exchange
      type(lcm_rank_timing_summary_t) :: y_payload_exchange
      type(lcm_rank_timing_summary_t) :: x_received_insertion
      type(lcm_rank_timing_summary_t) :: y_received_insertion
      type(lcm_rank_timing_summary_t) :: ownership_update
      type(lcm_rank_timing_summary_t) :: advance_total
      type(lcm_rank_timing_summary_t) :: prepare_particle_velocities
      type(lcm_rank_timing_summary_t) :: refine_particle_velocities
      type(lcm_rank_timing_summary_t) :: update_particle_positions
      type(lcm_rank_timing_summary_t) :: particle_boundaries
      type(lcm_rank_timing_summary_t) :: migration_total
      type(lcm_rank_timing_summary_t) :: migration_buffer_clear
      type(lcm_rank_timing_summary_t) :: migration_plan_allocate
      type(lcm_rank_timing_summary_t) :: migration_plan_classification
      type(lcm_rank_timing_summary_t) :: migration_plan_apply
      type(lcm_rank_timing_summary_t) :: migration_plan_deallocate
      type(lcm_rank_timing_summary_t) :: particle_exchange_total
      type(lcm_rank_timing_summary_t) :: particle_sources
   end type lcm_particle_exchange_timing_summary_t

   public :: lcm_get_particle_exchange_timing_summary

contains

   !> Summarize cumulative rank-local timings over the configured communicator.
   !!
   !! Every rank in the communicator must call this routine. The same summary
   !! is returned on every rank. Without configured parallel state, the local
   !! value is returned as the minimum, maximum, and mean for one rank.
   subroutine lcm_get_particle_exchange_timing_summary(summary)
      type(lcm_particle_exchange_timing_summary_t), intent(out) :: summary
      type(lcm_particle_exchange_metrics_t) :: metrics
      real(real64) :: local_values(timing_count)
      real(real64) :: minimum_values(timing_count)
      real(real64) :: maximum_values(timing_count)
      real(real64) :: summed_values(timing_count)
      type(MPI_Datatype) :: real64_datatype
      integer :: ierr

      call lcm_get_particle_exchange_metrics(metrics)
      call pack_timing_values(metrics, local_values)

      summary = lcm_particle_exchange_timing_summary_t()
      if (parallel_state%is_set) then
         summary%rank_count = parallel_state%parallel%size
         call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, &
                                  storage_size(local_values(1)) / 8, &
                                  real64_datatype, ierr)
         call require_mpi_success(ierr, 'match the timing real datatype')
         call MPI_ALLREDUCE(local_values, minimum_values, timing_count, &
                            real64_datatype, MPI_MIN, &
                            parallel_state%parallel%communicator, ierr)
         call require_mpi_success(ierr, 'reduce minimum exchange timings')
         call MPI_ALLREDUCE(local_values, maximum_values, timing_count, &
                            real64_datatype, MPI_MAX, &
                            parallel_state%parallel%communicator, ierr)
         call require_mpi_success(ierr, 'reduce maximum exchange timings')
         call MPI_ALLREDUCE(local_values, summed_values, timing_count, &
                            real64_datatype, MPI_SUM, &
                            parallel_state%parallel%communicator, ierr)
         call require_mpi_success(ierr, 'sum exchange timings')
      else
         minimum_values = local_values
         maximum_values = local_values
         summed_values = local_values
      endif

      call unpack_timing_summary( &
         summary, minimum_values, maximum_values, &
         summed_values / real(summary%rank_count, real64))
   end subroutine lcm_get_particle_exchange_timing_summary

   subroutine pack_timing_values(metrics, values)
      type(lcm_particle_exchange_metrics_t), intent(in) :: metrics
      real(real64), intent(out) :: values(timing_count)

      values(x_outgoing_preparation_index) = &
         metrics%x_outgoing_preparation_seconds
      values(y_outgoing_preparation_index) = &
         metrics%y_outgoing_preparation_seconds
      values(x_count_exchange_index) = metrics%x_count_exchange_seconds
      values(y_count_exchange_index) = metrics%y_count_exchange_seconds
      values(x_payload_exchange_index) = metrics%x_payload_exchange_seconds
      values(y_payload_exchange_index) = metrics%y_payload_exchange_seconds
      values(x_received_insertion_index) = &
         metrics%x_received_insertion_seconds
      values(y_received_insertion_index) = &
         metrics%y_received_insertion_seconds
      values(ownership_update_index) = metrics%ownership_update_seconds
      values(advance_total_index) = metrics%advance_total_seconds
      values(prepare_particle_velocities_index) = &
         metrics%prepare_particle_velocities_seconds
      values(refine_particle_velocities_index) = &
         metrics%refine_particle_velocities_seconds
      values(update_particle_positions_index) = &
         metrics%update_particle_positions_seconds
      values(particle_boundaries_index) = metrics%particle_boundaries_seconds
      values(migration_total_index) = metrics%migration_total_seconds
      values(migration_buffer_clear_index) = &
         metrics%migration_buffer_clear_seconds
      values(migration_plan_allocate_index) = &
         metrics%migration_plan_allocate_seconds
      values(migration_plan_classification_index) = &
         metrics%migration_plan_classification_seconds
      values(migration_plan_apply_index) = &
         metrics%migration_plan_apply_seconds
      values(migration_plan_deallocate_index) = &
         metrics%migration_plan_deallocate_seconds
      values(particle_exchange_total_index) = &
         metrics%particle_exchange_total_seconds
      values(particle_sources_index) = metrics%particle_sources_seconds
   end subroutine pack_timing_values

   subroutine unpack_timing_summary(summary, minimum_values, maximum_values, &
                                    mean_values)
      type(lcm_particle_exchange_timing_summary_t), intent(inout) :: summary
      real(real64), intent(in) :: minimum_values(timing_count)
      real(real64), intent(in) :: maximum_values(timing_count)
      real(real64), intent(in) :: mean_values(timing_count)

      call set_timing_summary(summary%x_outgoing_preparation, &
                              x_outgoing_preparation_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%y_outgoing_preparation, &
                              y_outgoing_preparation_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%x_count_exchange, &
                              x_count_exchange_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%y_count_exchange, &
                              y_count_exchange_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%x_payload_exchange, &
                              x_payload_exchange_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%y_payload_exchange, &
                              y_payload_exchange_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%x_received_insertion, &
                              x_received_insertion_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%y_received_insertion, &
                              y_received_insertion_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%ownership_update, &
                              ownership_update_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%advance_total, advance_total_index, &
                              minimum_values, maximum_values, mean_values)
      call set_timing_summary(summary%prepare_particle_velocities, &
                              prepare_particle_velocities_index, &
                              minimum_values, maximum_values, mean_values)
      call set_timing_summary(summary%refine_particle_velocities, &
                              refine_particle_velocities_index, &
                              minimum_values, maximum_values, mean_values)
      call set_timing_summary(summary%update_particle_positions, &
                              update_particle_positions_index, &
                              minimum_values, maximum_values, mean_values)
      call set_timing_summary(summary%particle_boundaries, &
                              particle_boundaries_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%migration_total, migration_total_index, &
                              minimum_values, maximum_values, mean_values)
      call set_timing_summary(summary%migration_buffer_clear, &
                              migration_buffer_clear_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%migration_plan_allocate, &
                              migration_plan_allocate_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%migration_plan_classification, &
                              migration_plan_classification_index, &
                              minimum_values, maximum_values, mean_values)
      call set_timing_summary(summary%migration_plan_apply, &
                              migration_plan_apply_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%migration_plan_deallocate, &
                              migration_plan_deallocate_index, &
                              minimum_values, maximum_values, mean_values)
      call set_timing_summary(summary%particle_exchange_total, &
                              particle_exchange_total_index, minimum_values, &
                              maximum_values, mean_values)
      call set_timing_summary(summary%particle_sources, particle_sources_index, &
                              minimum_values, maximum_values, mean_values)
   end subroutine unpack_timing_summary

   subroutine set_timing_summary(timing, index, minimum_values, &
                                 maximum_values, mean_values)
      type(lcm_rank_timing_summary_t), intent(out) :: timing
      integer, intent(in) :: index
      real(real64), intent(in) :: minimum_values(timing_count)
      real(real64), intent(in) :: maximum_values(timing_count)
      real(real64), intent(in) :: mean_values(timing_count)

      timing%minimum_seconds = minimum_values(index)
      timing%maximum_seconds = maximum_values(index)
      timing%mean_seconds = mean_values(index)
   end subroutine set_timing_summary

   subroutine require_mpi_success(ierr, operation)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: operation

      if (ierr /= MPI_SUCCESS) error stop 'LCM failed to ' // operation
   end subroutine require_mpi_success

end module lcm_particle_exchange_timing_summary
