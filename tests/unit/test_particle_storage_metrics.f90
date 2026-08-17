!> Verify on-demand accounting of active and reserved particle storage.
program test_particle_storage_metrics
   use iso_fortran_env, only: int64, real64
   use lcm_grid_state, only: lcm_set_grid
   use lcm_particle_state, only: lcm_add_particle, lcm_allocate_particle_cells
   use lcm_particle_storage_metrics, only: lcm_get_particle_storage_metrics, &
      lcm_particle_storage_metrics_t
   use lcm_types, only: lcm_grid_t, lcm_init_grid, lcm_particle_t

   implicit none

   type(lcm_grid_t) :: grid
   type(lcm_particle_storage_metrics_t) :: metrics
   type(lcm_particle_t) :: particle
   integer :: particle_index

   call lcm_get_particle_storage_metrics(metrics)
   call assert_equal(metrics%cell_count, 0_int64, &
                     'unallocated storage has no cells')
   call assert_equal(metrics%reserved_particle_capacity, 0_int64, &
                     'unallocated storage has no reserved capacity')
   call assert_true(metrics%particle_payload_bytes_per_slot > 0_int64, &
                    'particle payload size is available before allocation')

   call configure_grid()
   call lcm_allocate_particle_cells(4)

   do particle_index = 1, 5
      particle = lcm_particle_t()
      particle%id = int(particle_index, int64)
      call lcm_add_particle(1, 1, 1, particle)
   enddo
   particle = lcm_particle_t()
   particle%id = 6_int64
   call lcm_add_particle(2, 1, 1, particle)

   call lcm_get_particle_storage_metrics(metrics)
   call assert_equal(metrics%cell_count, 4_int64, 'local cell count')
   call assert_equal(metrics%allocated_cell_count, 4_int64, &
                     'cells with particle allocations')
   call assert_equal(metrics%empty_allocated_cell_count, 2_int64, &
                     'empty cells retaining particle allocations')
   call assert_equal(metrics%active_particle_count, 6_int64, &
                     'active particle count')
   call assert_equal(metrics%reserved_particle_capacity, 20_int64, &
                     'reserved particle capacity after growth')
   call assert_equal(metrics%unused_reserved_particle_capacity, 14_int64, &
                     'unused reserved particle slots')
   call assert_equal(metrics%active_particle_payload_bytes, &
                     6_int64 * metrics%particle_payload_bytes_per_slot, &
                     'active particle payload bytes')
   call assert_equal(metrics%reserved_particle_payload_bytes, &
                     20_int64 * metrics%particle_payload_bytes_per_slot, &
                     'reserved particle payload bytes')

contains

   subroutine configure_grid()
      real(real64) :: dz(1)
      real(real64) :: z_center(1)
      real(real64) :: z_face(2)

      dz = [100.0_real64]
      z_center = [50.0_real64]
      z_face = [0.0_real64, 100.0_real64]

      call lcm_init_grid(grid, 2, 2, 1, 2, 2, &
                         1, 2, 1, 2, 10.0_real64, 10.0_real64, &
                         z_center, z_face, dz, 0, 0, 0, 1, 1, 1)
      call lcm_set_grid(grid)
   end subroutine configure_grid

   subroutine assert_equal(actual, expected, description)
      integer(int64), intent(in) :: actual
      integer(int64), intent(in) :: expected
      character(len=*), intent(in) :: description

      if (actual /= expected) error stop 'assertion failed: ' // description
   end subroutine assert_equal

   subroutine assert_true(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description

      if (.not. condition) error stop 'assertion failed: ' // description
   end subroutine assert_true

end program test_particle_storage_metrics
