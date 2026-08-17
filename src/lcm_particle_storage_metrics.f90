!> On-demand, rank-local metrics for cell-owned particle storage.
!!
!! The byte values describe lcm_particle_t payload storage only. They exclude
!! allocatable-array descriptors, allocator bookkeeping, cell containers,
!! exchange buffers, and other process memory. Total resident memory must be
!! measured separately by the host or system tooling.
module lcm_particle_storage_metrics
   use iso_fortran_env, only: int64
   use lcm_particle_state, only: cell_particles
   use lcm_types, only: lcm_particle_t

   implicit none

   private

   type, public :: lcm_particle_storage_metrics_t
      integer(int64) :: cell_count = 0_int64
      integer(int64) :: allocated_cell_count = 0_int64
      integer(int64) :: empty_allocated_cell_count = 0_int64
      integer(int64) :: active_particle_count = 0_int64
      integer(int64) :: reserved_particle_capacity = 0_int64
      integer(int64) :: unused_reserved_particle_capacity = 0_int64
      integer(int64) :: particle_payload_bytes_per_slot = 0_int64
      integer(int64) :: active_particle_payload_bytes = 0_int64
      integer(int64) :: reserved_particle_payload_bytes = 0_int64
   end type lcm_particle_storage_metrics_t

   public :: lcm_get_particle_storage_metrics

contains

   subroutine lcm_get_particle_storage_metrics(metrics)
      type(lcm_particle_storage_metrics_t), intent(out) :: metrics
      integer :: i
      integer :: j
      integer :: k
      integer(int64) :: particle_payload_bits

      metrics = lcm_particle_storage_metrics_t()
      particle_payload_bits = int(storage_size(lcm_particle_t()), int64)
      metrics%particle_payload_bytes_per_slot = &
         (particle_payload_bits + 7_int64) / 8_int64

      if (.not. allocated(cell_particles)) return

      metrics%cell_count = int(size(cell_particles), int64)
      do k = 1, size(cell_particles, 3)
         do j = 1, size(cell_particles, 2)
            do i = 1, size(cell_particles, 1)
               metrics%active_particle_count = metrics%active_particle_count + &
                  int(cell_particles(i,j,k)%count, int64)
               metrics%reserved_particle_capacity = &
                  metrics%reserved_particle_capacity + &
                  int(cell_particles(i,j,k)%capacity, int64)
               if (allocated(cell_particles(i,j,k)%particles)) then
                  metrics%allocated_cell_count = &
                     metrics%allocated_cell_count + 1_int64
                  if (cell_particles(i,j,k)%count == 0) then
                     metrics%empty_allocated_cell_count = &
                        metrics%empty_allocated_cell_count + 1_int64
                  endif
               endif
            enddo
         enddo
      enddo

      metrics%unused_reserved_particle_capacity = &
         metrics%reserved_particle_capacity - metrics%active_particle_count
      metrics%active_particle_payload_bytes = metrics%active_particle_count * &
         metrics%particle_payload_bytes_per_slot
      metrics%reserved_particle_payload_bytes = &
         metrics%reserved_particle_capacity * &
         metrics%particle_payload_bytes_per_slot
   end subroutine lcm_get_particle_storage_metrics

end module lcm_particle_storage_metrics
