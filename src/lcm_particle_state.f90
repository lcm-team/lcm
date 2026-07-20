!> LCM-owned computational particle storage.
!!
!! Each local grid cell owns a dynamic list of computational particles.
!! cell_particles(i,j,k)%count is the number of active particles in that cell.
!! cell_particles(i,j,k)%capacity is the allocated storage length. When a cell
!! is full, lcm_add_particle resizes that cell's particle array and preserves
!! the existing particles.
module lcm_particle_state
   use iso_fortran_env, only: int64
   use lcm_grid_state, only: grid_state
   use lcm_types, only: lcm_particle_cell_t, lcm_particle_t

   implicit none

   private

   public :: cell_particles
   public :: lcm_allocate_particle_cells
   public :: lcm_add_particle
   public :: lcm_reset_particle_cells
   public :: lcm_particle_id_for_cell_member

   type(lcm_particle_cell_t), allocatable, protected, save :: cell_particles(:,:,:)

contains

   subroutine lcm_allocate_particle_cells(initial_capacity_per_cell)
      integer, intent(in), optional :: initial_capacity_per_cell
      integer :: capacity
      integer :: i, j, k

      if (.not. grid_state%is_set) then
         error stop "LCM grid must be initialized before allocating particle cells"
      endif
      if (allocated(cell_particles)) then
         error stop "LCM particle cells are already allocated"
      endif

      capacity = 0
      if (present(initial_capacity_per_cell)) then
         if (initial_capacity_per_cell < 0) then
            error stop "LCM initial particle cell capacity must not be negative"
         endif
         capacity = initial_capacity_per_cell
      endif

      ! Allocate one particle container per local grid cell. Particle counts
      ! are stored inside each cell and grow dynamically.
      allocate(cell_particles(grid_state%grid%nx_local, &
                              grid_state%grid%ny_local, &
                              grid_state%grid%nz))

      if (capacity > 0) then
         do k = 1, grid_state%grid%nz
            do j = 1, grid_state%grid%ny_local
               do i = 1, grid_state%grid%nx_local
                  call resize_cell_storage(cell_particles(i,j,k), capacity)
               enddo
            enddo
         enddo
      endif
   end subroutine lcm_allocate_particle_cells

   subroutine lcm_add_particle(i, j, k, particle)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      type(lcm_particle_t), intent(in) :: particle

      call require_particle_cells_allocated()
      call validate_cell_index(i, j, k)

      if (cell_particles(i,j,k)%count == cell_particles(i,j,k)%capacity) then
         call resize_cell_storage(cell_particles(i,j,k), &
                                  next_capacity(cell_particles(i,j,k)%capacity))
      endif

      cell_particles(i,j,k)%count = cell_particles(i,j,k)%count + 1
      cell_particles(i,j,k)%particles(cell_particles(i,j,k)%count) = particle
   end subroutine lcm_add_particle

   subroutine lcm_reset_particle_cells()
      integer :: i, j, k

      call require_particle_cells_allocated()

      do k = 1, size(cell_particles, 3)
         do j = 1, size(cell_particles, 2)
            do i = 1, size(cell_particles, 1)
               cell_particles(i,j,k)%count = 0
            enddo
         enddo
      enddo
   end subroutine lcm_reset_particle_cells

   integer(int64) function lcm_particle_id_for_cell_member( &
      i, j, k, cell_particle_index, particles_per_cell) result(particle_id)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: cell_particle_index
      integer, intent(in) :: particles_per_cell

      integer(int64) :: global_i
      integer(int64) :: global_j
      integer(int64) :: global_cell_id

      call require_grid_set()
      call validate_grid_cell_index(i, j, k)

      if (particles_per_cell <= 0) then
         error stop "LCM particles per cell must be positive for particle IDs"
      endif
      if (cell_particle_index < 1 .or. &
          cell_particle_index > particles_per_cell) then
         error stop "LCM cell particle index is outside the ID range"
      endif

      global_i = int(grid_state%grid%global_i_start + i - 1, int64)
      global_j = int(grid_state%grid%global_j_start + j - 1, int64)

      global_cell_id = &
         (int(k - 1, int64) * int(grid_state%grid%ny_global, int64) + &
          global_j - 1_int64) * int(grid_state%grid%nx_global, int64) + &
         global_i

      particle_id = (global_cell_id - 1_int64) * &
                    int(particles_per_cell, int64) + &
                    int(cell_particle_index, int64)
   end function lcm_particle_id_for_cell_member

   subroutine require_particle_cells_allocated()
      if (.not. allocated(cell_particles)) then
         error stop "LCM particle cells have not been allocated"
      endif
   end subroutine require_particle_cells_allocated

   subroutine require_grid_set()
      if (.not. grid_state%is_set) then
         error stop "LCM grid must be initialized before generating particle IDs"
      endif
   end subroutine require_grid_set

   subroutine validate_grid_cell_index(i, j, k)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k

      if (i < 1 .or. i > grid_state%grid%nx_local .or. &
          j < 1 .or. j > grid_state%grid%ny_local .or. &
          k < 1 .or. k > grid_state%grid%nz) then
         error stop "LCM particle ID cell index is outside the local grid"
      endif
   end subroutine validate_grid_cell_index

   subroutine validate_cell_index(i, j, k)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k

      if (i < 1 .or. i > size(cell_particles, 1) .or. &
          j < 1 .or. j > size(cell_particles, 2) .or. &
          k < 1 .or. k > size(cell_particles, 3)) then
         error stop "LCM particle cell index is outside the local grid"
      endif
   end subroutine validate_cell_index

   integer function next_capacity(current_capacity)
      integer, intent(in) :: current_capacity

      if (current_capacity == 0) then
         next_capacity = 1
      else
         next_capacity = 2 * current_capacity
      endif
   end function next_capacity

   subroutine resize_cell_storage(cell, new_capacity)
      type(lcm_particle_cell_t), intent(inout) :: cell
      integer, intent(in) :: new_capacity
      type(lcm_particle_t), allocatable :: old_particles(:)

      if (new_capacity <= cell%capacity) return

      if (allocated(cell%particles)) then
         ! Fortran allocatable arrays cannot be resized in place. Temporarily
         ! keep the active particles while replacing the storage array.
         allocate(old_particles(cell%count))
         old_particles = cell%particles(1:cell%count)
         deallocate(cell%particles)
      endif

      allocate(cell%particles(new_capacity))
      if (allocated(old_particles)) then
         cell%particles(1:size(old_particles)) = old_particles
      endif
      cell%capacity = new_capacity
   end subroutine resize_cell_storage

end module lcm_particle_state
