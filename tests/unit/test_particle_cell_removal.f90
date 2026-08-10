!> Verify safe in-place deletion while traversing one cell container.
program test_particle_cell_removal
   use iso_fortran_env, only: int64, real64
   use lcm_grid_state, only: lcm_set_grid
   use lcm_particle_state, only: cell_particles, lcm_add_particle, &
                                 lcm_allocate_particle_cells, lcm_remove_particle
   use lcm_types, only: lcm_grid_t, lcm_init_grid, lcm_particle_t

   implicit none

   type(lcm_grid_t) :: grid

   call configure_grid()
   call lcm_allocate_particle_cells(4)
   call add_particle(1_int64, 50.0_real64)
   call add_particle(2_int64, -1.0_real64)
   call add_particle(3_int64, -2.0_real64)
   call add_particle(4_int64, 75.0_real64)

   call remove_particles_below_zero()

   call assert_true(cell_particles(1,1,1)%count == 2, &
                    "two in-domain particles remain")
   call assert_true(cell_particles(1,1,1)%particles(1)%id == 1_int64, &
                    "first retained particle keeps its order")
   call assert_true(cell_particles(1,1,1)%particles(2)%id == 4_int64, &
                    "second retained particle is not skipped after deletion")

contains

   subroutine configure_grid()
      real(real64) :: dz(1)
      real(real64) :: z_center(1)
      real(real64) :: z_face(2)

      dz = [100.0_real64]
      z_center = [50.0_real64]
      z_face = [0.0_real64, 100.0_real64]

      call lcm_init_grid(grid, 1, 1, 1, 1, 1, &
                         1, 1, 1, 1, 10.0_real64, 10.0_real64, &
                         z_center, z_face, dz, 0, 0, 0, 0, 0, &
                         1, 1, 1, 1, -1, -1, -1, -1, &
                         .false., .false., 1)
      call lcm_set_grid(grid)
   end subroutine configure_grid

   subroutine add_particle(id, z)
      integer(int64), intent(in) :: id
      real(real64), intent(in) :: z

      type(lcm_particle_t) :: particle

      particle%id = id
      particle%z = z
      call lcm_add_particle(1, 1, 1, particle)
   end subroutine add_particle

   subroutine remove_particles_below_zero()
      integer :: particle_index

      particle_index = 1
      do while (particle_index <= cell_particles(1,1,1)%count)
         if (cell_particles(1,1,1)%particles(particle_index)%z < 0.0_real64) then
            call lcm_remove_particle(1, 1, 1, particle_index)
         else
            particle_index = particle_index + 1
         endif
      enddo
   end subroutine remove_particles_below_zero

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_true

end program test_particle_cell_removal
