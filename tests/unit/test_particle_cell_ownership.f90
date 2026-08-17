program test_particle_cell_ownership
   use iso_fortran_env, only: int64, real64
   use lcm_grid_state, only: lcm_set_grid
   use lcm_particle_exchange_metrics, only: lcm_get_particle_exchange_metrics, &
      lcm_particle_exchange_metrics_t, lcm_reset_particle_exchange_metrics, &
      lcm_set_particle_exchange_metrics_enabled
   use lcm_particle_state, only: cell_particles, lcm_add_particle, &
                                 lcm_allocate_particle_cells, &
                                 lcm_update_particle_cell_ownership
   use lcm_types, only: lcm_grid_t, lcm_init_grid, lcm_particle_t

   implicit none

   type(lcm_grid_t) :: grid
   type(lcm_particle_exchange_metrics_t) :: metrics

   call configure_grid()
   call lcm_allocate_particle_cells(1)
   call add_particles()
   call lcm_reset_particle_exchange_metrics()
   call lcm_set_particle_exchange_metrics_enabled(.true.)

   call lcm_update_particle_cell_ownership()

   call lcm_get_particle_exchange_metrics(metrics)
   call assert_true(metrics%ownership_update_calls == 1_int64, &
                    'one ownership update measured')
   call assert_true(metrics%ownership_particle_examinations == 7_int64, &
                    'actual ownership examination operations')
   call assert_true(metrics%ownership_particles_moved == 3_int64, &
                    'particles moved between local cells')
   call assert_true(metrics%ownership_update_seconds >= 0.0_real64, &
                    'ownership update time is nonnegative')

   call assert_true(cell_particles(1,1,1)%count == 1, "original cell count")
   call assert_true(cell_particles(2,1,1)%count == 1, "east cell count")
   call assert_true(cell_particles(1,2,1)%count == 1, "north cell count")
   call assert_true(cell_particles(1,1,2)%count == 1, "upper cell count")

   call assert_particle_id(1, 1, 1, 1, 10_int64, "stationary particle id")
   call assert_particle_id(2, 1, 1, 1, 20_int64, "eastward particle id")
   call assert_particle_id(1, 2, 1, 1, 30_int64, "northward particle id")
   call assert_particle_id(1, 1, 2, 1, 40_int64, "upward particle id")

contains

   subroutine configure_grid()
      real(real64) :: dz(2)
      real(real64) :: z_center(2)
      real(real64) :: z_face(3)

      dz = [100.0_real64, 100.0_real64]
      z_center = [50.0_real64, 150.0_real64]
      z_face = [0.0_real64, 100.0_real64, 200.0_real64]

      call lcm_init_grid(grid, 2, 2, 2, 2, 2, &
                         1, 2, 1, 2, 50.0_real64, 50.0_real64, &
                         z_center, z_face, dz, 0, 0, 0, 1, 1, 1)
      call lcm_set_grid(grid)
   end subroutine configure_grid

   subroutine add_particles()
      call add_particle(10_int64, 25.0_real64, 25.0_real64, 50.0_real64)
      call add_particle(20_int64, 75.0_real64, 25.0_real64, 50.0_real64)
      call add_particle(30_int64, 25.0_real64, 75.0_real64, 50.0_real64)
      call add_particle(40_int64, 25.0_real64, 25.0_real64, 150.0_real64)
   end subroutine add_particles

   subroutine add_particle(id, x, y, z)
      integer(int64), intent(in) :: id
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      real(real64), intent(in) :: z

      type(lcm_particle_t) :: particle

      particle = lcm_particle_t()
      particle%id = id
      particle%x = x
      particle%y = y
      particle%z = z
      call lcm_add_particle(1, 1, 1, particle)
   end subroutine add_particle

   subroutine assert_particle_id(i, j, k, n, expected_id, message)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n
      integer(int64), intent(in) :: expected_id
      character(len=*), intent(in) :: message

      call assert_true(cell_particles(i,j,k)%particles(n)%id == expected_id, &
                       message)
   end subroutine assert_particle_id

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_true

end program test_particle_cell_ownership
