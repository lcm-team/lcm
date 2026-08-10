!> LCM-owned computational particle storage.
!!
!! Each local grid cell owns a dynamic list of computational particles.
!! cell_particles(i,j,k)%count is the number of active particles in that cell.
!! cell_particles(i,j,k)%capacity is the allocated storage length. When a cell
!! is full, lcm_add_particle resizes that cell's particle array and preserves
!! the existing particles.
module lcm_particle_state
   use iso_fortran_env, only: int64, real64
   use lcm_config_state, only: config_state
   use lcm_field_state, only: field_state
   use lcm_grid_state, only: grid_state
   use lcm_parameters, only: epsv, fac_cond, molecular_weight_water, pi, &
                             rho_liq, rv, surface_tension_water
   use lcm_particle_location, only: lcm_locate_particle_cell
   use lcm_exchange_directions, only: lcm_exchange_direction_from_offsets
   use lcm_particle_properties, only: lcm_particle_radius
   use lcm_types, only: lcm_particle_cell_t, lcm_particle_t
   use lcm_user, only: lcm_user_init_particles
   use lcm_velocity_interpolation, only: lcm_interpolate_particle_velocity
   use lcm_random, only: randomu

   implicit none

   private

   public :: cell_particles
   public :: lcm_allocate_particle_cells
   public :: lcm_add_particle
   public :: lcm_move_particle
   public :: lcm_remove_particle
   public :: lcm_reset_particle_cells
   public :: lcm_initialize_particle_id_counter
   public :: lcm_next_particle_id
   public :: lcm_initialize_computational_particles
   public :: lcm_prepare_particle_velocities
   public :: lcm_refine_particle_velocities
   public :: lcm_update_particle_positions
   public :: lcm_reflect_particles_at_top
   public :: lcm_apply_particles_at_bottom
   public :: lcm_extract_outgoing_particles
   public :: lcm_update_particle_cell_ownership
   public :: lcm_clear_outgoing_particles
   public :: set_particle_position
   public :: particle_random_seed_id
   public :: cell_tabs_lcm
   public :: cell_absolute_temperature
   public :: qsatw_liquid
   public :: inverse_error_function_approx
   public :: ice_freezing_temperature
   public :: lcm_outgoing_particle_t
   public :: outgoing_particles
   public :: n_outgoing_particles

   type :: lcm_outgoing_particle_t
      type(lcm_particle_t) :: particle
      integer :: direction = 0
   end type lcm_outgoing_particle_t

   type(lcm_particle_cell_t), allocatable, protected, save :: cell_particles(:,:,:)
   type(lcm_outgoing_particle_t), allocatable, protected, save :: outgoing_particles(:)
   type(lcm_particle_t), allocatable, save :: aerosol_initialization_template(:)
   integer, protected, save :: n_outgoing_particles = 0
   integer(int64), protected, save :: n_step_id = 0_int64
   integer(int64), protected, save :: n_id = 0_int64
   logical, save :: aerosol_initialization_template_is_set = .false.
   logical, protected, save :: particle_id_counter_is_initialized = .false.

   integer, parameter :: liquid_particle_type = 1
   integer, parameter :: ice_nucleating_particle_type = 2
   integer, parameter :: n_inv_cdf = 50000000

contains

   subroutine lcm_allocate_particle_cells(initial_capacity_per_cell, &
                                          initial_capacity_vertical_levels)
      integer, intent(in), optional :: initial_capacity_per_cell
      integer, intent(in), optional :: initial_capacity_vertical_levels
      integer :: capacity
      integer :: capacity_vertical_levels
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

      capacity_vertical_levels = grid_state%grid%nz
      if (present(initial_capacity_vertical_levels)) then
         if (initial_capacity_vertical_levels < 0 .or. &
             initial_capacity_vertical_levels > grid_state%grid%nz) then
            error stop "LCM initial-capacity levels must be within the vertical grid"
         endif
         capacity_vertical_levels = initial_capacity_vertical_levels
      endif

      ! Allocate one particle container per local grid cell. Particle counts
      ! are stored inside each cell and grow dynamically.
      allocate(cell_particles(grid_state%grid%nx_local, &
                              grid_state%grid%ny_local, &
                              grid_state%grid%nz))
      call lcm_clear_outgoing_particles()

      if (capacity > 0) then
         do k = 1, capacity_vertical_levels
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

      call append_particle_to_cell(cell_particles(i,j,k), particle)
   end subroutine lcm_add_particle

   subroutine lcm_remove_particle(i, j, k, particle_index)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: particle_index
      integer :: count

      call require_particle_cells_allocated()
      call validate_cell_index(i, j, k)

      count = cell_particles(i,j,k)%count
      call validate_particle_index(particle_index, count)

      ! Keep the remaining order stable for particle tracking and exchange.
      if (particle_index < count) then
         cell_particles(i,j,k)%particles(particle_index:count - 1) = &
            cell_particles(i,j,k)%particles(particle_index + 1:count)
      endif

      cell_particles(i,j,k)%particles(count) = lcm_particle_t()
      cell_particles(i,j,k)%count = count - 1
   end subroutine lcm_remove_particle

   subroutine lcm_move_particle(source_i, source_j, source_k, particle_index, &
                                destination_i, destination_j, destination_k)
      integer, intent(in) :: source_i
      integer, intent(in) :: source_j
      integer, intent(in) :: source_k
      integer, intent(in) :: particle_index
      integer, intent(in) :: destination_i
      integer, intent(in) :: destination_j
      integer, intent(in) :: destination_k
      type(lcm_particle_t) :: particle

      call require_particle_cells_allocated()
      call validate_cell_index(source_i, source_j, source_k)
      call validate_cell_index(destination_i, destination_j, destination_k)
      call validate_particle_index(particle_index, &
                                   cell_particles(source_i,source_j,source_k)%count)

      if (source_i == destination_i .and. source_j == destination_j .and. &
          source_k == destination_k) then
         return
      endif

      particle = cell_particles(source_i,source_j,source_k)%particles(particle_index)
      call lcm_add_particle(destination_i, destination_j, destination_k, particle)
      call lcm_remove_particle(source_i, source_j, source_k, particle_index)
   end subroutine lcm_move_particle

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

   subroutine lcm_clear_outgoing_particles()
      n_outgoing_particles = 0
   end subroutine lcm_clear_outgoing_particles

   subroutine lcm_initialize_computational_particles()
      integer :: i, j, k, n
      integer :: initialization_levels
      type(lcm_particle_t) :: particle

      call require_particle_cells_allocated()
      call lcm_reset_particle_cells()
      call lcm_initialize_particle_id_counter()
      call reset_aerosol_initialization_template()

      initialization_levels = min(config_state%nzl, grid_state%grid%nz)

      do k = 1, initialization_levels
         do j = 1, grid_state%grid%ny_local
            do i = 1, grid_state%grid%nx_local
               do n = 1, config_state%np_xy
                  particle = lcm_particle_t()
                  particle%id = lcm_next_particle_id()
                  call init_particle_properties(i, j, k, n, particle)
                  call lcm_add_particle(i, j, k, particle)
               enddo

               call lcm_init_particles(i, j, k)
            enddo
         enddo
      enddo
   end subroutine lcm_initialize_computational_particles

   subroutine lcm_initialize_particle_id_counter()
      integer(int64) :: n_subdomains
      integer(int64) :: subdomain_id

      call require_grid_set()

      n_subdomains = int(grid_state%grid%n_subdomains_x, int64) * &
                     int(grid_state%grid%n_subdomains_y, int64)
      subdomain_id = &
         int(grid_state%grid%subdomain_y_index, int64) * &
         int(grid_state%grid%n_subdomains_x, int64) + &
         int(grid_state%grid%subdomain_x_index, int64)

      n_step_id = huge(1_int64) / n_subdomains
      n_id = n_step_id * subdomain_id
      particle_id_counter_is_initialized = .true.
   end subroutine lcm_initialize_particle_id_counter

   integer(int64) function lcm_next_particle_id() result(particle_id)
      if (.not. particle_id_counter_is_initialized) then
         error stop "LCM particle ID counter has not been initialized"
      endif

      particle_id = n_id
      n_id = n_id + 1_int64
   end function lcm_next_particle_id

   subroutine lcm_prepare_particle_velocities()
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      logical :: is_local

      call require_particle_cells_allocated()

      do k = 1, size(cell_particles, 3)
         do j = 1, size(cell_particles, 2)
            do i = 1, size(cell_particles, 1)
               do n = 1, cell_particles(i,j,k)%count
                  call lcm_interpolate_particle_velocity( &
                     cell_particles(i,j,k)%particles(n)%x, &
                     cell_particles(i,j,k)%particles(n)%y, &
                     cell_particles(i,j,k)%particles(n)%z, &
                     cell_particles(i,j,k)%particles(n)%u, &
                     cell_particles(i,j,k)%particles(n)%v, &
                     cell_particles(i,j,k)%particles(n)%w, is_local)
                  if (.not. is_local) then
                     error stop "LCM particle velocity preparation found a nonlocal particle"
                  endif
               enddo
            enddo
         enddo
      enddo
   end subroutine lcm_prepare_particle_velocities

   subroutine lcm_refine_particle_velocities(time_step)
      real(real64), intent(in) :: time_step

      integer :: i
      integer :: j
      integer :: k
      integer :: n
      integer :: iteration
      logical :: is_local
      real(real64) :: predicted_u
      real(real64) :: predicted_v
      real(real64) :: predicted_w
      real(real64) :: unow
      real(real64) :: vnow
      real(real64) :: wnow
      real(real64) :: xp
      real(real64) :: yp
      real(real64) :: zp

      call require_particle_cells_allocated()

      if (config_state%n_uvw_iteration < 1) return

      do k = 1, size(cell_particles, 3)
         do j = 1, size(cell_particles, 2)
            do i = 1, size(cell_particles, 1)
               do n = 1, cell_particles(i,j,k)%count
                  unow = cell_particles(i,j,k)%particles(n)%u
                  vnow = cell_particles(i,j,k)%particles(n)%v
                  wnow = cell_particles(i,j,k)%particles(n)%w

                  do iteration = 1, config_state%n_uvw_iteration
                     xp = cell_particles(i,j,k)%particles(n)%x + &
                          cell_particles(i,j,k)%particles(n)%u * time_step
                     yp = cell_particles(i,j,k)%particles(n)%y + &
                          cell_particles(i,j,k)%particles(n)%v * time_step * &
                          real(grid_state%grid%y_dimension, real64)
                     zp = cell_particles(i,j,k)%particles(n)%z + &
                          cell_particles(i,j,k)%particles(n)%w * time_step

                     call lcm_interpolate_particle_velocity( &
                        xp, yp, zp, predicted_u, predicted_v, predicted_w, &
                        is_local, clamp_to_local=.true.)
                     cell_particles(i,j,k)%particles(n)%u = &
                        0.5_real64 * (unow + predicted_u)
                     cell_particles(i,j,k)%particles(n)%v = &
                        0.5_real64 * (vnow + predicted_v)
                     cell_particles(i,j,k)%particles(n)%w = &
                        0.5_real64 * (wnow + predicted_w)
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine lcm_refine_particle_velocities

   subroutine lcm_update_particle_positions(time_step, current_step)
      real(real64), intent(in) :: time_step
      integer, intent(in) :: current_step

      integer :: i
      integer :: j
      integer :: k
      integer :: n

      call require_particle_cells_allocated()

      do k = 1, size(cell_particles, 3)
         do j = 1, size(cell_particles, 2)
            do i = 1, size(cell_particles, 1)
               do n = 1, cell_particles(i,j,k)%count
                  cell_particles(i,j,k)%particles(n)%x = &
                     cell_particles(i,j,k)%particles(n)%x + &
                     cell_particles(i,j,k)%particles(n)%u * time_step
                  cell_particles(i,j,k)%particles(n)%y = &
                     cell_particles(i,j,k)%particles(n)%y + &
                     cell_particles(i,j,k)%particles(n)%v * time_step * &
                     real(grid_state%grid%y_dimension, real64)
                  cell_particles(i,j,k)%particles(n)%z = &
                     cell_particles(i,j,k)%particles(n)%z + &
                     cell_particles(i,j,k)%particles(n)%w * time_step

                  if (config_state%switch_sgs_velocities) then
                     cell_particles(i,j,k)%particles(n)%x = &
                        cell_particles(i,j,k)%particles(n)%x + &
                        cell_particles(i,j,k)%particles(n)%u_sgs * time_step
                     cell_particles(i,j,k)%particles(n)%y = &
                        cell_particles(i,j,k)%particles(n)%y + &
                        cell_particles(i,j,k)%particles(n)%v_sgs * time_step * &
                        real(grid_state%grid%y_dimension, real64)
                     cell_particles(i,j,k)%particles(n)%z = &
                        cell_particles(i,j,k)%particles(n)%z + &
                        cell_particles(i,j,k)%particles(n)%w_sgs * time_step
                  endif

                  if (config_state%switch_sedimentation .and. &
                      current_step >= config_state%start_sedimentation) then
                     cell_particles(i,j,k)%particles(n)%z = &
                        cell_particles(i,j,k)%particles(n)%z - &
                        cell_particles(i,j,k)%particles(n)%w_sedi * time_step
                  endif
               enddo
            enddo
         enddo
      enddo
   end subroutine lcm_update_particle_positions

   subroutine lcm_reflect_particles_at_top(ztop)
      real(real64), intent(in) :: ztop

      integer :: i
      integer :: j
      integer :: k
      integer :: n

      call require_particle_cells_allocated()

      do k = 1, size(cell_particles, 3)
         do j = 1, size(cell_particles, 2)
            do i = 1, size(cell_particles, 1)
               do n = 1, cell_particles(i,j,k)%count
                  if (cell_particles(i,j,k)%particles(n)%z > ztop) then
                     ! SAM-LCM micro_boundaries: ztop + (ztop - z), -abs(w_sgs).
                     cell_particles(i,j,k)%particles(n)%z = &
                        ztop + (ztop - cell_particles(i,j,k)%particles(n)%z)
                     cell_particles(i,j,k)%particles(n)%w_sgs = &
                        -abs(cell_particles(i,j,k)%particles(n)%w_sgs)
                  endif
               enddo
            enddo
         enddo
      enddo
   end subroutine lcm_reflect_particles_at_top

   subroutine lcm_apply_particles_at_bottom(zbottom, removal_radius, aerosol_density, &
                                            deleted_particle_count, deleted_water_mass, &
                                            deleted_aerosol_mass)
      real(real64), intent(in) :: zbottom
      real(real64), intent(in) :: removal_radius
      real(real64), intent(in) :: aerosol_density(:)
      integer(int64), intent(out) :: deleted_particle_count
      real(real64), intent(out) :: deleted_water_mass
      real(real64), intent(out) :: deleted_aerosol_mass

      integer :: i
      integer :: j
      integer :: k
      integer :: n
      real(real64) :: radius

      call require_particle_cells_allocated()

      deleted_particle_count = 0_int64
      deleted_water_mass = 0.0_real64
      deleted_aerosol_mass = 0.0_real64

      do k = 1, size(cell_particles, 3)
         do j = 1, size(cell_particles, 2)
            do i = 1, size(cell_particles, 1)
               n = 1
               do while (n <= cell_particles(i,j,k)%count)
                  if (cell_particles(i,j,k)%particles(n)%z < zbottom) then
                     radius = lcm_particle_radius( &
                        cell_particles(i,j,k)%particles(n), aerosol_density)
                     if (radius > removal_radius) then
                        deleted_particle_count = deleted_particle_count + 1_int64
                        deleted_water_mass = deleted_water_mass + &
                           cell_particles(i,j,k)%particles(n)%mass_water
                        deleted_aerosol_mass = deleted_aerosol_mass + &
                           cell_particles(i,j,k)%particles(n)%mass_aero
                        call lcm_remove_particle(i, j, k, n)
                     else
                        ! SAM-LCM: zbottom + (zbottom - z), abs(w_sgs).
                        cell_particles(i,j,k)%particles(n)%z = zbottom + &
                           (zbottom - cell_particles(i,j,k)%particles(n)%z)
                        cell_particles(i,j,k)%particles(n)%w_sgs = &
                           abs(cell_particles(i,j,k)%particles(n)%w_sgs)
                        n = n + 1
                     endif
                  else
                     n = n + 1
                  endif
               enddo
            enddo
         enddo
      enddo
   end subroutine lcm_apply_particles_at_bottom

   subroutine lcm_extract_outgoing_particles()
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      integer :: new_i
      integer :: new_j
      integer :: new_k
      integer :: x_direction
      integer :: y_direction
      logical :: is_local

      call require_particle_cells_allocated()
      call lcm_clear_outgoing_particles()

      do k = 1, size(cell_particles, 3)
         do j = 1, size(cell_particles, 2)
            do i = 1, size(cell_particles, 1)
               n = 1
               do while (n <= cell_particles(i,j,k)%count)
                  call lcm_locate_particle_cell( &
                     cell_particles(i,j,k)%particles(n)%x, &
                     cell_particles(i,j,k)%particles(n)%y, &
                     cell_particles(i,j,k)%particles(n)%z, &
                     new_i, new_j, new_k, is_local)

                  if (particle_leaves_to_neighbor_rank(new_i, new_j, new_k, &
                                                       x_direction, y_direction)) then
                     call append_outgoing_particle( &
                        cell_particles(i,j,k)%particles(n), x_direction, y_direction)
                     call lcm_remove_particle(i, j, k, n)
                  else
                     n = n + 1
                  endif
               enddo
            enddo
         enddo
      enddo
   end subroutine lcm_extract_outgoing_particles

   subroutine lcm_update_particle_cell_ownership()
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      integer :: new_i
      integer :: new_j
      integer :: new_k
      logical :: is_local

      call require_particle_cells_allocated()

      do k = 1, size(cell_particles, 3)
         do j = 1, size(cell_particles, 2)
            do i = 1, size(cell_particles, 1)
               n = 1
               do while (n <= cell_particles(i,j,k)%count)
                  call lcm_locate_particle_cell( &
                     cell_particles(i,j,k)%particles(n)%x, &
                     cell_particles(i,j,k)%particles(n)%y, &
                     cell_particles(i,j,k)%particles(n)%z, &
                     new_i, new_j, new_k, is_local)

                  if (is_local .and. (new_i /= i .or. new_j /= j .or. &
                                       new_k /= k)) then
                     call lcm_move_particle(i, j, k, n, new_i, new_j, new_k)
                  else
                     n = n + 1
                  endif
               enddo
            enddo
         enddo
      enddo
   end subroutine lcm_update_particle_cell_ownership

   subroutine init_particle_properties(i, j, k, cell_particle_index, particle)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: cell_particle_index
      type(lcm_particle_t), intent(inout) :: particle
      integer :: seed_id

      particle%cell_particle_index = cell_particle_index
      particle%last_i = grid_state%grid%global_i_start + i - 1
      particle%last_j = grid_state%grid%global_j_start + j - 1
      particle%last_k = k
      particle%spray_marker = -999
      call set_particle_tracking_marker(i, j, k, cell_particle_index, particle)

      particle%particle_type = liquid_particle_type
      seed_id = particle_random_seed_id()
      call set_particle_position(i, j, k, seed_id, particle)
      particle%w_sedi = 0.0_real64
      particle%temperature = cell_absolute_temperature(i, j, k, &
                                                       cell_tabs_lcm(i, j, k))
      particle%qv = field_state%qv(i,j,k)
   end subroutine init_particle_properties

   subroutine set_particle_position(i, j, k, seed_id, particle)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: seed_id
      type(lcm_particle_t), intent(inout) :: particle

      particle%x = (real(i - 1 + grid_state%nx_sb, real64) + &
                    real(randomu(seed_id), real64)) * grid_state%grid%dx
      if (grid_state%grid%y_dimension == 1) then
         particle%y = (real(j - 1 + grid_state%ny_sb, real64) + &
                       real(randomu(seed_id), real64)) * grid_state%grid%dy
      else
         particle%y = 0.5_real64 * grid_state%grid%dy
      endif
      particle%z = grid_state%grid%z_face(k) + &
                   real(randomu(seed_id), real64) * grid_state%grid%dz(k)
   end subroutine set_particle_position

   subroutine set_particle_tracking_marker(i, j, k, cell_particle_index, particle)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: cell_particle_index
      type(lcm_particle_t), intent(inout) :: particle

      integer :: marker_divisor
      integer :: particle_sequence
      integer :: tracking_period
      integer :: tracked_particles_per_subdomain

      particle%track_marker = -999
      if (.not. config_state%track_some_particles) return

      tracked_particles_per_subdomain = &
         max(1, config_state%n_track_particles / &
                grid_state%grid%n_subdomains_x / &
                grid_state%grid%n_subdomains_y)
      tracking_period = initial_particle_count_for_tracking() / &
                        tracked_particles_per_subdomain
      marker_divisor = initial_particle_count_for_tracking() / &
                       config_state%n_track_particles

      if (tracking_period <= 0 .or. marker_divisor <= 0) then
         error stop "LCM particle tracking requests more markers than initial particles"
      endif

      particle_sequence = initial_particle_sequence(i, j, k, cell_particle_index)
      if (mod(particle_sequence, tracking_period) == 0) then
         particle%track_marker = particle_subdomain_id() * 1000000 + &
                                 particle_sequence / marker_divisor
      endif
   end subroutine set_particle_tracking_marker

   integer function initial_particle_count_for_tracking() result(particle_count)
      particle_count = grid_state%grid%nx_local * grid_state%grid%ny_local * &
                       min(config_state%nzl, grid_state%grid%nz) * &
                       config_state%np_xy
   end function initial_particle_count_for_tracking

   integer function initial_particle_sequence(i, j, k, cell_particle_index) &
      result(sequence)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: cell_particle_index

      sequence = (((k - 1) * grid_state%grid%ny_local + (j - 1)) * &
                  grid_state%grid%nx_local + (i - 1)) * &
                 config_state%np_xy + cell_particle_index
   end function initial_particle_sequence

   subroutine lcm_init_particles(i, j, k)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k

      integer :: n, nn, m, m_max, i_mode, np_mode, np_done
      integer :: seed_id
      real(real64) :: dlogr, r_l, r_r, r_aero, n_dummy
      real(real64) :: f_ak, f_bk, f_ck, ak, bk, ck
      real(real64) :: p_random, p_usm, x_rand
      real(real64) :: in_typeprobs(4)
      real(real64) :: r_aero_min(4), r_aero_max(4)

      if (config_state%np_xy == 0) return

      if (uses_shared_aerosol_initialization_template() .and. &
          aerosol_initialization_template_is_set) then
         call copy_aerosol_initialization_template(i, j, k)
         call initialize_equilibrium_liquid_state(i, j, k)
         call lcm_user_init_particles(i, j, k, cell_particles(i,j,k))
         return
      endif

      seed_id = particle_random_seed_id()
      if (config_state%switch_ice) then
         call ice_type_probabilities(in_typeprobs)
      else
         in_typeprobs = 0.0_real64
      endif

      np_done = 0
      ! Initialize aerosols for every mode separately, preserving the
      ! SAM-LCM micro_init_particles loop structure.
      do i_mode = 1, 4
         np_mode = floor(real(config_state%np_xy, real64) * &
                         real(config_state%shares_aero(i_mode), real64) / &
                         sum(real(config_state%shares_aero, real64)))
         if (i_mode == 1) then
            n = 0
            do nn = 1, 4
               n = n + floor(real(config_state%np_xy, real64) * &
                             real(config_state%shares_aero(nn), real64) / &
                             sum(real(config_state%shares_aero, real64)))
            enddo
            np_mode = np_mode + config_state%np_xy - n
         endif

         if (np_mode == 0) cycle
         if (config_state%n_aero(i_mode) == 0.0) then
            error stop "LCM aerosol mode with zero concentration received particles"
         endif

         r_aero_min(i_mode) = 1.0e-9_real64
         r_aero_max(i_mode) = 1.0e-5_real64

         ak = r_aero_min(i_mode)
         f_ak = log_normal_cdf(real(config_state%n_aero(i_mode), real64), &
                               real(config_state%sigma_aero(i_mode), real64), &
                               real(config_state%rm_aero(i_mode), real64), &
                               ak) / real(config_state%n_aero(i_mode), real64) - &
                5.0e-5_real64

         bk = r_aero_max(i_mode)
         f_bk = log_normal_cdf(real(config_state%n_aero(i_mode), real64), &
                               real(config_state%sigma_aero(i_mode), real64), &
                               real(config_state%rm_aero(i_mode), real64), &
                               bk) / real(config_state%n_aero(i_mode), real64) - &
                5.0e-5_real64

         do
            ck = (ak * f_bk - bk * f_ak) / (f_bk - f_ak)
            f_ck = log_normal_cdf(real(config_state%n_aero(i_mode), real64), &
                                  real(config_state%sigma_aero(i_mode), real64), &
                                  real(config_state%rm_aero(i_mode), real64), &
                                  ck) / real(config_state%n_aero(i_mode), real64) - &
                   5.0e-5_real64

            if (f_ak * f_ck > 0.0_real64) then
               ak = ck
               f_ak = f_ck
            else
               bk = ck
               f_bk = f_ck
            endif

            if (abs(f_ck) < 1.0e-12_real64) exit
         enddo
         r_aero_min(i_mode) = ck

         ak = r_aero_min(i_mode)
         f_ak = log_normal_cdf(real(config_state%n_aero(i_mode), real64), &
                               real(config_state%sigma_aero(i_mode), real64), &
                               real(config_state%rm_aero(i_mode), real64), &
                               ak) / real(config_state%n_aero(i_mode), real64) - &
                (1.0_real64 - 5.0e-5_real64)

         bk = r_aero_max(i_mode)
         f_bk = log_normal_cdf(real(config_state%n_aero(i_mode), real64), &
                               real(config_state%sigma_aero(i_mode), real64), &
                               real(config_state%rm_aero(i_mode), real64), &
                               bk) / real(config_state%n_aero(i_mode), real64) - &
                (1.0_real64 - 5.0e-5_real64)

         do
            ck = (ak * f_bk - bk * f_ak) / (f_bk - f_ak)
            f_ck = log_normal_cdf(real(config_state%n_aero(i_mode), real64), &
                                  real(config_state%sigma_aero(i_mode), real64), &
                                  real(config_state%rm_aero(i_mode), real64), &
                                  ck) / real(config_state%n_aero(i_mode), real64) - &
                   (1.0_real64 - 5.0e-5_real64)

            if (f_ak * f_ck > 0.0_real64) then
               ak = ck
               f_ak = f_ck
            else
               bk = ck
               f_bk = f_ck
            endif

            if (abs(f_ck) < 1.0e-12_real64) exit
         enddo
         r_aero_max(i_mode) = ck

         if (trim(config_state%init_aero) == "bin") then
            dlogr = (log(r_aero_max(i_mode)) - log(r_aero_min(i_mode))) / &
                    real(np_mode, real64)

            do n = np_done + 1, np_done + np_mode
               nn = n - np_done

               r_l = exp(log(r_aero_min(i_mode)) + real(nn - 1, real64) * dlogr)
               r_r = exp(log(r_aero_min(i_mode)) + real(nn, real64) * dlogr)
               r_aero = sqrt(r_l * r_r)

               n_dummy = log_normal_pdf(real(config_state%n_aero(i_mode), real64), &
                                        real(config_state%sigma_aero(i_mode), real64), &
                                        real(config_state%rm_aero(i_mode), real64), &
                                        r_aero)
               if (trim(config_state%unit_aero) == "1/m3") then
                  cell_particles(i,j,k)%particles(n)%weight = &
                     int(n_dummy * r_aero * dlogr * real(grid_state%dV(k), real64), int64)
               else if (trim(config_state%unit_aero) == "1/kg") then
                  cell_particles(i,j,k)%particles(n)%weight = &
                     int(n_dummy * r_aero * dlogr * real(grid_state%dV(k), real64) * &
                         field_state%density(k), int64)
               else
                  error stop "LCM unknown aerosol initialization unit"
               endif
               if (config_state%switch_ice) then
                  cell_particles(i,j,k)%particles(n)%weight = &
                     int(real(cell_particles(i,j,k)%particles(n)%weight, real64) / &
                         (1.0_real64 - real(config_state%frac_ice, real64)), int64)
               endif

               call set_particle_aerosol_state(i, j, k, n, i_mode, r_aero)
               call apply_ice_initialization(i, j, k, n, i_mode, r_aero, &
                                             seed_id, in_typeprobs)
            enddo

            do n = np_done + 1, np_done + np_mode
               if (cell_particles(i,j,k)%particles(n)%weight == 0_int64) then
                  m_max = n
                  do m = np_done + 1, np_done + np_mode
                     if (cell_particles(i,j,k)%particles(m)%weight > &
                         cell_particles(i,j,k)%particles(m_max)%weight .and. &
                         cell_particles(i,j,k)%particles(m)%weight >= 2_int64) then
                        m_max = m
                     endif
                  enddo

                  if (cell_particles(i,j,k)%particles(m_max)%weight == 0_int64) then
                     error stop "LCM initialization of aerosols not possible"
                  endif

                  cell_particles(i,j,k)%particles(n)%weight = &
                     cell_particles(i,j,k)%particles(m_max)%weight / 2_int64
                  cell_particles(i,j,k)%particles(m_max)%weight = &
                     cell_particles(i,j,k)%particles(m_max)%weight - &
                     cell_particles(i,j,k)%particles(n)%weight

                  cell_particles(i,j,k)%particles(n)%mass_aero = &
                     cell_particles(i,j,k)%particles(m_max)%mass_aero / &
                     real(cell_particles(i,j,k)%particles(m_max)%weight + &
                          cell_particles(i,j,k)%particles(n)%weight, real64) * &
                     real(cell_particles(i,j,k)%particles(n)%weight, real64)
                  cell_particles(i,j,k)%particles(m_max)%mass_aero = &
                     cell_particles(i,j,k)%particles(m_max)%mass_aero - &
                     cell_particles(i,j,k)%particles(n)%mass_aero

                  cell_particles(i,j,k)%particles(n)%aerosol_portions = &
                     cell_particles(i,j,k)%particles(m_max)%aerosol_portions
                  cell_particles(i,j,k)%particles(n)%particle_type = &
                     cell_particles(i,j,k)%particles(m_max)%particle_type

                  if (cell_particles(i,j,k)%particles(n)%weight == 0_int64 .or. &
                      cell_particles(i,j,k)%particles(m_max)%weight == 0_int64) then
                     error stop "LCM initialization of aerosols not possible"
                  endif
               endif
            enddo
         else if (trim(config_state%init_aero) == "USM") then
            p_usm = 1.0_real64 / real(np_mode, real64)
            r_l = r_aero_min(i_mode)

            do n = np_done + 1, np_done + np_mode - 1
               nn = n - np_done
               r_r = aerosol_percentile_radius(i_mode, p_usm * real(nn, real64), &
                                                r_l, r_aero_max(i_mode))
               r_aero = sqrt(r_l * r_r)
               r_l = r_r

               call set_equal_mode_weight(i, j, k, n, i_mode, np_mode)
               call set_particle_aerosol_state(i, j, k, n, i_mode, r_aero)
               call apply_ice_initialization(i, j, k, n, i_mode, r_aero, &
                                             seed_id, in_typeprobs)
            enddo

            r_aero = sqrt(r_l * r_aero_max(i_mode))
            n = np_done + np_mode
            cell_particles(i,j,k)%particles(n)%weight = &
               int(real(grid_state%dV(k), real64) * &
                   sum(real(config_state%n_aero, real64)) / &
                   (real(config_state%np_xy, real64) * &
                    (1.0_real64 - real(config_state%frac_ice, real64))), int64)
            call set_particle_aerosol_state(i, j, k, n, i_mode, r_aero, &
                                            set_aerosol_portions=.false.)
            call apply_ice_initialization(i, j, k, n, i_mode, r_aero, &
                                          seed_id, in_typeprobs)
         else if (trim(config_state%init_aero) == "random") then
            p_random = 1.0_real64 / real(n_inv_cdf, real64)

            do n = np_done + 1, np_done + np_mode
               find_r: do nn = 1, 1000
                  x_rand = real(randomu(seed_id), real64)
                  r_aero = inverse_cdf_sample_radius(i_mode, x_rand, &
                                                     r_aero_min(i_mode), &
                                                     r_aero_max(i_mode), &
                                                     p_random)
                  n_dummy = equal_mode_weight(i_mode, k, np_mode)

                  if (n_dummy < log_normal_pdf(real(config_state%n_aero(i_mode), real64), &
                                               real(config_state%sigma_aero(i_mode), real64), &
                                               real(config_state%rm_aero(i_mode), real64), &
                                               r_aero)) then
                     exit find_r
                  endif
               enddo find_r

               cell_particles(i,j,k)%particles(n)%weight = int(n_dummy, int64)
               if (config_state%switch_ice) then
                  cell_particles(i,j,k)%particles(n)%weight = &
                     int(real(cell_particles(i,j,k)%particles(n)%weight, real64) / &
                         (1.0_real64 - real(config_state%frac_ice, real64)), int64)
               endif
               call set_particle_aerosol_state(i, j, k, n, i_mode, r_aero)
               call apply_ice_initialization(i, j, k, n, i_mode, r_aero, &
                                             seed_id, in_typeprobs)
            enddo
         else
            error stop "LCM unknown aerosol initialization mode"
         endif

         np_done = np_done + np_mode
      enddo

      if (uses_shared_aerosol_initialization_template()) then
         call save_aerosol_initialization_template(i, j, k)
      endif

      call initialize_equilibrium_liquid_state(i, j, k)

      call lcm_user_init_particles(i, j, k, cell_particles(i,j,k))
   end subroutine lcm_init_particles

   logical function uses_shared_aerosol_initialization_template() &
      result(uses_template)
      uses_template = trim(config_state%init_aero) == "bin" .or. &
                      trim(config_state%init_aero) == "USM"
   end function uses_shared_aerosol_initialization_template

   subroutine reset_aerosol_initialization_template()
      if (allocated(aerosol_initialization_template)) then
         deallocate(aerosol_initialization_template)
      endif
      aerosol_initialization_template_is_set = .false.
   end subroutine reset_aerosol_initialization_template

   subroutine save_aerosol_initialization_template(i, j, k)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k

      integer :: n

      if (allocated(aerosol_initialization_template)) then
         deallocate(aerosol_initialization_template)
      endif
      allocate(aerosol_initialization_template(config_state%np_xy))

      do n = 1, config_state%np_xy
         aerosol_initialization_template(n)%weight = &
            cell_particles(i,j,k)%particles(n)%weight
         aerosol_initialization_template(n)%mass_aero = &
            cell_particles(i,j,k)%particles(n)%mass_aero
         aerosol_initialization_template(n)%aerosol_portions = &
            cell_particles(i,j,k)%particles(n)%aerosol_portions
         aerosol_initialization_template(n)%particle_type = &
            cell_particles(i,j,k)%particles(n)%particle_type
         aerosol_initialization_template(n)%in_surface_area = &
            cell_particles(i,j,k)%particles(n)%in_surface_area
         aerosol_initialization_template(n)%freezing_temperature = &
            cell_particles(i,j,k)%particles(n)%freezing_temperature
      enddo

      aerosol_initialization_template_is_set = .true.
   end subroutine save_aerosol_initialization_template

   subroutine copy_aerosol_initialization_template(i, j, k)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k

      integer :: n

      if (.not. aerosol_initialization_template_is_set) then
         error stop "LCM aerosol initialization template has not been set"
      endif

      do n = 1, config_state%np_xy
         cell_particles(i,j,k)%particles(n)%weight = &
            aerosol_initialization_template(n)%weight
         cell_particles(i,j,k)%particles(n)%mass_aero = &
            aerosol_initialization_template(n)%mass_aero
         cell_particles(i,j,k)%particles(n)%aerosol_portions = &
            aerosol_initialization_template(n)%aerosol_portions
         cell_particles(i,j,k)%particles(n)%particle_type = &
            aerosol_initialization_template(n)%particle_type
         cell_particles(i,j,k)%particles(n)%in_surface_area = &
            aerosol_initialization_template(n)%in_surface_area
         cell_particles(i,j,k)%particles(n)%freezing_temperature = &
            aerosol_initialization_template(n)%freezing_temperature
      enddo
   end subroutine copy_aerosol_initialization_template

   integer function particle_random_seed_id() result(seed_id)
      seed_id = particle_subdomain_id() + 1
   end function particle_random_seed_id

   integer function particle_subdomain_id() result(subdomain_id)
      subdomain_id = grid_state%grid%subdomain_y_index * &
                     grid_state%grid%n_subdomains_x + &
                     grid_state%grid%subdomain_x_index
   end function particle_subdomain_id

   subroutine set_equal_mode_weight(i, j, k, n, i_mode, np_mode)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n
      integer, intent(in) :: i_mode
      integer, intent(in) :: np_mode

      cell_particles(i,j,k)%particles(n)%weight = &
         int(equal_mode_weight(i_mode, k, np_mode), int64)
   end subroutine set_equal_mode_weight

   real(real64) function equal_mode_weight(i_mode, k, np_mode) result(weight)
      integer, intent(in) :: i_mode
      integer, intent(in) :: k
      integer, intent(in) :: np_mode

      if (trim(config_state%unit_aero) == "1/m3") then
         weight = real(grid_state%dV(k), real64) * &
                  real(config_state%n_aero(i_mode), real64) / &
                  (real(np_mode, real64) * &
                   (1.0_real64 - real(config_state%frac_ice, real64)))
      else if (trim(config_state%unit_aero) == "1/kg") then
         weight = real(grid_state%dV(k), real64) * field_state%density(k) * &
                  real(config_state%n_aero(i_mode), real64) / &
                  (real(np_mode, real64) * &
                   (1.0_real64 - real(config_state%frac_ice, real64)))
      else
         error stop "LCM unknown aerosol initialization unit"
      endif
   end function equal_mode_weight

   subroutine set_particle_aerosol_state(i, j, k, n, i_mode, r_aero, &
                                         set_aerosol_portions)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n
      integer, intent(in) :: i_mode
      real(real64), intent(in) :: r_aero
      logical, intent(in), optional :: set_aerosol_portions
      logical :: should_set_aerosol_portions

      cell_particles(i,j,k)%particles(n)%mass_aero = &
         real(cell_particles(i,j,k)%particles(n)%weight, real64) * &
         4.0_real64 / 3.0_real64 * pi * r_aero**3 * &
         real(config_state%rho_aero(i_mode), real64)
      should_set_aerosol_portions = .true.
      if (present(set_aerosol_portions)) then
         should_set_aerosol_portions = set_aerosol_portions
      endif
      if (should_set_aerosol_portions) then
         cell_particles(i,j,k)%particles(n)%aerosol_portions = 0.0_real64
         cell_particles(i,j,k)%particles(n)%aerosol_portions(i_mode) = 1.0_real64
      endif
      cell_particles(i,j,k)%particles(n)%particle_type = liquid_particle_type
      cell_particles(i,j,k)%particles(n)%in_surface_area = 0.0_real64
      cell_particles(i,j,k)%particles(n)%freezing_temperature = 233.15_real64
   end subroutine set_particle_aerosol_state

   subroutine initialize_equilibrium_liquid_state(i, j, k)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k

      integer :: n
      integer :: nn
      real(real64) :: cell_qv
      real(real64) :: cell_temperature
      real(real64) :: dql
      real(real64) :: qsat
      real(real64) :: r_aero
      real(real64) :: ss_gb

      cell_qv = field_state%qv(i,j,k)
      cell_temperature = cell_tabs_lcm(i, j, k)

      ! Port of SAM-LCM micro_init_particles' five-iteration equilibrium
      ! initialization. The adjusted cell qv/T are local here because attached
      ! host fields are read-only in standalone LCM initialization.
      do nn = 1, 5
         qsat = qsatw_liquid(cell_temperature, field_state%pressure(k))
         ss_gb = (cell_qv - qsat) / qsat
         ss_gb = min(ss_gb, -0.01_real64)

         dql = 0.0_real64
         do n = 1, cell_particles(i,j,k)%count
            r_aero = particle_aerosol_radius(i, j, k, n)
            dql = dql - cell_particles(i,j,k)%particles(n)%mass_water
            cell_particles(i,j,k)%particles(n)%mass_water = &
               real(cell_particles(i,j,k)%particles(n)%weight, real64) * &
               4.0_real64 / 3.0_real64 * pi * rho_liq * &
               max(equilibrium_radius(ss_gb, cell_temperature, i, j, k, n), &
                   r_aero)**3
            dql = dql + cell_particles(i,j,k)%particles(n)%mass_water
         enddo

         dql = dql / (real(grid_state%dV(k), real64) * field_state%density(k))
         cell_qv = cell_qv - dql
         cell_temperature = cell_temperature + fac_cond * dql
      enddo

      do n = 1, cell_particles(i,j,k)%count
         cell_particles(i,j,k)%particles(n)%temperature = &
            cell_absolute_temperature(i, j, k, cell_temperature)
         cell_particles(i,j,k)%particles(n)%qv = cell_qv
         cell_particles(i,j,k)%particles(n)%eta = &
            cell_qv - qsatw_liquid(cell_temperature, field_state%pressure(k))
         cell_particles(i,j,k)%particles(n)%n_monomer = 1_int64
      enddo
   end subroutine initialize_equilibrium_liquid_state

   real(real64) function cell_tabs_lcm(i, j, k) result(tabs_lcm)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k

      ! SAM-LCM diagnoses tabs_LCM as liquid-water static energy minus gamaz.
      tabs_lcm = field_state%static_energy(i,j,k) - cell_gamaz(i, j, k)
   end function cell_tabs_lcm

   real(real64) function cell_absolute_temperature(i, j, k, tabs_lcm) result(temperature)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      real(real64), intent(in) :: tabs_lcm

      temperature = tabs_lcm + cell_gamaz(i, j, k)
   end function cell_absolute_temperature

   real(real64) function cell_gamaz(i, j, k) result(gamaz)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k

      gamaz = field_state%static_energy(i,j,k) - field_state%temperature(i,j,k)
   end function cell_gamaz

   real(real64) function particle_aerosol_radius(i, j, k, n) result(r_aero)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n

      r_aero = (particle_aerosol_volume(i, j, k, n) * &
                3.0_real64 / 4.0_real64 / pi)**(1.0_real64 / 3.0_real64)
   end function particle_aerosol_radius

   real(real64) function particle_aerosol_volume(i, j, k, n) result(volume)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n

      integer :: i_mode
      real(real64) :: weight

      volume = 0.0_real64
      weight = real(cell_particles(i,j,k)%particles(n)%weight, real64)
      if (weight <= 0.0_real64) return

      do i_mode = 1, 4
         if (config_state%rho_aero(i_mode) > 0.0 .and. &
             cell_particles(i,j,k)%particles(n)%aerosol_portions(i_mode) /= &
             0.0_real64) then
            volume = volume + &
               cell_particles(i,j,k)%particles(n)%aerosol_portions(i_mode) * &
               cell_particles(i,j,k)%particles(n)%mass_aero / &
               real(config_state%rho_aero(i_mode), real64) / weight
         endif
      enddo
   end function particle_aerosol_volume

   real(real64) function equilibrium_radius(ss_gb, temperature, i, j, k, n) &
      result(radius)
      real(real64), intent(in) :: ss_gb
      real(real64), intent(in) :: temperature
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n

      real(real64) :: afactor
      real(real64) :: bfactor
      real(real64) :: radius_old

      afactor = 2.0_real64 * surface_tension_water / &
                (rho_liq * rv * temperature)
      bfactor = particle_bfactor(i, j, k, n)

      radius_old = 1.0_real64
      radius = 1.0e-6_real64
      do while (abs((radius - radius_old) / radius_old) > 1.0e-20_real64)
         radius_old = radius
         radius = (bfactor / (afactor / radius - ss_gb))**(1.0_real64 / 3.0_real64)
      enddo
   end function equilibrium_radius

   real(real64) function particle_bfactor(i, j, k, n) result(bfactor)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n

      integer :: i_mode
      real(real64) :: weight

      bfactor = 0.0_real64
      weight = real(cell_particles(i,j,k)%particles(n)%weight, real64)
      if (weight <= 0.0_real64) return

      if (only_one_aerosol_mode()) then
         bfactor = bfactor_coefficient(1)
      else
         do i_mode = 1, 4
            if (config_state%molecular_weight_aero(i_mode) > 0.0 .and. &
                cell_particles(i,j,k)%particles(n)%aerosol_portions(i_mode) /= &
                0.0_real64) then
               bfactor = bfactor + bfactor_coefficient(i_mode) * &
                  cell_particles(i,j,k)%particles(n)%aerosol_portions(i_mode)
            endif
         enddo
      endif

      bfactor = bfactor * &
                cell_particles(i,j,k)%particles(n)%mass_aero / weight
   end function particle_bfactor

   logical function only_one_aerosol_mode() result(has_one_mode)
      has_one_mode = config_state%shares_aero(1) > 0.0 .and. &
                     all(config_state%shares_aero(2:4) == 0.0)
   end function only_one_aerosol_mode

   real(real64) function bfactor_coefficient(i_mode) result(coefficient)
      integer, intent(in) :: i_mode

      coefficient = &
         real(config_state%vanthoff_aero(i_mode), real64) / &
         real(config_state%molecular_weight_aero(i_mode), real64) * &
         molecular_weight_water / (rho_liq * 4.0_real64 * pi / 3.0_real64)
   end function bfactor_coefficient

   real(real64) function qsatw_liquid(temperature, pressure) result(qsatw)
      real(real64), intent(in) :: temperature
      real(real64), intent(in) :: pressure

      real(real64) :: e

      e = 610.94_real64 * &
          exp(17.625_real64 * (temperature - 273.15_real64) / &
              max(temperature - 30.11_real64, 1.0_real64))
      qsatw = epsv * e / max(pressure - e, 1.0_real64)
   end function qsatw_liquid

   subroutine ice_type_probabilities(in_typeprobs)
      real(real64), intent(out) :: in_typeprobs(4)
      integer :: m
      real(real64) :: n_ice_sum

      n_ice_sum = sum(real(config_state%n_ice, real64))
      do m = 1, 4
         in_typeprobs(m) = sum(real(config_state%n_ice(1:m), real64)) / n_ice_sum
      enddo
   end subroutine ice_type_probabilities

   subroutine apply_ice_initialization(i, j, k, n, i_mode, r_aero, &
                                       seed_id, in_typeprobs)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n
      integer, intent(in) :: i_mode
      integer, intent(in) :: seed_id
      real(real64), intent(in) :: r_aero
      real(real64), intent(in) :: in_typeprobs(4)

      integer :: intype
      real(real64) :: n_dummy
      real(real64) :: r_in
      real(real64) :: sfc_in
      real(real64) :: x_rand

      if (.not. config_state%switch_ice) return

      if (real(randomu(seed_id), real64) < real(config_state%frac_ice, real64)) then
         cell_particles(i,j,k)%particles(n)%particle_type = &
            ice_nucleating_particle_type

         x_rand = real(randomu(seed_id), real64)
         intype = minloc(in_typeprobs - x_rand, &
                         mask=in_typeprobs >= x_rand, dim=1)

         n_dummy = real(config_state%n_ice(intype), real64) * &
                   real(grid_state%dV(k), real64) / &
                   (real(config_state%np_xy, real64) * &
                    real(config_state%frac_ice, real64))
         if (trim(config_state%unit_aero) == "1/kg") then
            n_dummy = n_dummy * field_state%density(k)
         endif

         cell_particles(i,j,k)%particles(n)%mass_aero = &
            n_dummy * 4.0_real64 / 3.0_real64 * pi * &
            real(config_state%rho_aero(i_mode), real64) * r_aero**3
         cell_particles(i,j,k)%particles(n)%weight = int(n_dummy, int64)

         r_in = -1.0_real64
         do while (r_in <= 0.001e-6_real64 .or. r_in >= 10.0e-4_real64)
            x_rand = real(randomu(seed_id), real64)
            r_in = real(config_state%rm_in(intype), real64) * &
                   exp(sqrt(2.0_real64) * &
                       log(real(config_state%sigma_in(intype), real64)) * &
                       min(max(inverse_error_function_approx( &
                                  2.0_real64 * x_rand - 1.0_real64), &
                               -2.0_real64), 2.0_real64))
         enddo

         sfc_in = 4.0_real64 * pi * r_in**2
         cell_particles(i,j,k)%particles(n)%in_surface_area = 0.0_real64
         cell_particles(i,j,k)%particles(n)%in_surface_area(intype) = &
            n_dummy * sfc_in
         cell_particles(i,j,k)%particles(n)%freezing_temperature = &
            ice_freezing_temperature(intype, sfc_in, seed_id)
      else
         cell_particles(i,j,k)%particles(n)%particle_type = liquid_particle_type
         cell_particles(i,j,k)%particles(n)%in_surface_area = 0.0_real64
         cell_particles(i,j,k)%particles(n)%freezing_temperature = 233.15_real64
      endif
   end subroutine apply_ice_initialization

   real(real64) function ice_freezing_temperature(m, sfc, seed_id) result(t_frz)
      integer, intent(in) :: m
      integer, intent(in) :: seed_id
      real(real64), intent(in) :: sfc

      integer :: loc14
      integer :: tcounter
      real(real64) :: frzprob_dust(289)
      real(real64) :: frzprob_pollen(289)
      real(real64) :: frzprob_snowmax(289)
      real(real64) :: frzprob_soot(289)
      real(real64) :: ns_dust(289)
      real(real64) :: ns_pollen(289)
      real(real64) :: ns_snowmax(289)
      real(real64) :: ns_soot(289)
      real(real64) :: randnum
      real(real64) :: tempc(289)

      ! Formula ported from SAM-LCM micro_init_particles:t_frz.
      tempc = [(-38.0_real64 + (36.0_real64 / 288.0_real64) * &
                real(tcounter, real64), tcounter = 0, 288)]
      loc14 = minloc(abs(tempc + 14.0_real64), dim=1)

      ns_dust = exp(-0.517_real64 * &
                    max(-30.0_real64, min(tempc, -14.0_real64)) + &
                    8.781_real64)
      ns_soot = 7.463_real64 * &
                exp(-0.0101_real64 * &
                    max(-34.0_real64, min(tempc, -18.0_real64))**2 - &
                    0.8525_real64 * &
                    max(-34.0_real64, min(tempc, -18.0_real64)) + &
                    0.7667_real64)
      ns_snowmax = 4.53e11_real64 * &
                   (1.0_real64 - exp(-2.0e-10_real64 * &
                    exp(-max(-20.0_real64, min(tempc, -2.0_real64)) * &
                    2.34_real64)))
      ns_pollen(1:loc14) = &
         exp(14.701_real64 + 6.0842_real64 * &
             atan(-15.8840_real64 * &
                  max(min(tempc(1:loc14), -14.0_real64), -35.0_real64)) + &
             1.4573_real64 * &
             log(-0.45627_real64 * &
                 max(min(tempc(1:loc14), -14.0_real64), -35.0_real64)))
      ns_pollen(loc14+1:289) = &
         exp(-17.864_real64 - 6.9238_real64 * &
             max(min(tempc(loc14+1:289), -2.0_real64), -14.0_real64) - &
             0.59183_real64 * &
             max(min(tempc(loc14+1:289), -2.0_real64), -14.0_real64) - &
             1.7131e-2_real64 * &
             max(min(tempc(loc14+1:289), -2.0_real64), -14.0_real64))

      frzprob_dust = 1.0_real64 - exp(-sfc * ns_dust)
      frzprob_soot = 1.0_real64 - exp(-sfc * ns_soot)
      frzprob_snowmax = 1.0_real64 - exp(-sfc * ns_snowmax)
      frzprob_pollen = 1.0_real64 - exp(-sfc * ns_pollen)

      randnum = real(randomu(seed_id), real64)
      select case (m)
      case (1)
         t_frz = tempc(minloc(abs(frzprob_soot - randnum), dim=1))
      case (2)
         t_frz = tempc(minloc(abs(frzprob_dust - randnum), dim=1))
      case (3)
         t_frz = tempc(minloc(abs(frzprob_snowmax - randnum), dim=1))
      case (4)
         t_frz = tempc(minloc(abs(frzprob_pollen - randnum), dim=1))
      case default
         t_frz = -40.0_real64
      end select
      t_frz = max(t_frz, -40.0_real64)
      t_frz = t_frz + 273.15_real64
   end function ice_freezing_temperature

   real(real64) function inverse_error_function_approx(x) result(erfi)
      real(real64), intent(in) :: x
      real(real64), parameter :: a = 0.147_real64

      ! Winitzki approximation used by SAM-LCM micro_init_particles:ERFI.
      erfi = 2.0_real64 / (pi * a) + log(1.0_real64 - x**2) / 2.0_real64
      erfi = sign(1.0_real64, x) * &
             sqrt(sqrt(erfi**2 - log(1.0_real64 - x**2) / a) - erfi)
   end function inverse_error_function_approx

   real(real64) function inverse_cdf_sample_radius(i_mode, x_rand, &
                                                   r_aero_min, r_aero_max, &
                                                   p_random) result(r_aero)
      integer, intent(in) :: i_mode
      real(real64), intent(in) :: x_rand
      real(real64), intent(in) :: r_aero_min
      real(real64), intent(in) :: r_aero_max
      real(real64), intent(in) :: p_random

      integer :: cdf_index
      real(real64) :: left_percentile
      real(real64) :: right_percentile
      real(real64) :: r_l
      real(real64) :: r_r

      ! SAM-LCM stores inv_cdf(0:n_inv_cdf). Standalone computes the selected
      ! entry on demand to avoid allocating the full lookup table during tests.
      cdf_index = max(min(nint(x_rand * real(n_inv_cdf, real64)), n_inv_cdf), 0)

      if (cdf_index == 0) then
         r_aero = r_aero_min
      else if (cdf_index == n_inv_cdf) then
         r_aero = r_aero_max
      else
         right_percentile = p_random * real(cdf_index, real64)
         if (cdf_index == 1) then
            r_l = r_aero_min
         else
            left_percentile = p_random * real(cdf_index - 1, real64)
            r_l = aerosol_percentile_radius(i_mode, left_percentile, &
                                             r_aero_min, r_aero_max)
         endif
         r_r = aerosol_percentile_radius(i_mode, right_percentile, &
                                          r_l, r_aero_max)
         r_aero = sqrt(r_l * r_r)
      endif
   end function inverse_cdf_sample_radius

   real(real64) function aerosol_percentile_radius(i_mode, percentile, &
                                                   r_left, r_right) result(radius)
      integer, intent(in) :: i_mode
      real(real64), intent(in) :: percentile
      real(real64), intent(in) :: r_left
      real(real64), intent(in) :: r_right

      real(real64) :: ak
      real(real64) :: bk
      real(real64) :: ck
      real(real64) :: f_ak
      real(real64) :: f_bk
      real(real64) :: f_ck

      ak = r_left
      f_ak = log_normal_cdf(real(config_state%n_aero(i_mode), real64), &
                            real(config_state%sigma_aero(i_mode), real64), &
                            real(config_state%rm_aero(i_mode), real64), &
                            ak) / real(config_state%n_aero(i_mode), real64) - &
             percentile

      bk = r_right
      f_bk = log_normal_cdf(real(config_state%n_aero(i_mode), real64), &
                            real(config_state%sigma_aero(i_mode), real64), &
                            real(config_state%rm_aero(i_mode), real64), &
                            bk) / real(config_state%n_aero(i_mode), real64) - &
             percentile

      do
         ck = (ak * f_bk - bk * f_ak) / (f_bk - f_ak)
         f_ck = log_normal_cdf(real(config_state%n_aero(i_mode), real64), &
                               real(config_state%sigma_aero(i_mode), real64), &
                               real(config_state%rm_aero(i_mode), real64), &
                               ck) / real(config_state%n_aero(i_mode), real64) - &
                percentile

         if (f_ak * f_ck > 0.0_real64) then
            ak = ck
            f_ak = f_ck
         else
            bk = ck
            f_bk = f_ck
         endif

         if (abs(f_ck) < 1.0e-12_real64) exit
      enddo
      radius = ck
   end function aerosol_percentile_radius

   real(real64) function log_normal_pdf(n, s, m, x) result(pdf)
      real(real64), intent(in) :: n
      real(real64), intent(in) :: s
      real(real64), intent(in) :: m
      real(real64), intent(in) :: x

      pdf = 0.0_real64
      if (n /= 0.0_real64 .and. s /= 0.0_real64 .and. &
          m /= 0.0_real64 .and. x > 0.0_real64) then
         pdf = n / (sqrt(2.0_real64 * pi) * log(s) * x) * &
               exp(-1.0_real64 * log(x / m)**2 / &
                   (2.0_real64 * log(s)**2))
      endif
   end function log_normal_pdf

   real(real64) function log_normal_cdf(n, s, m, x) result(cdf)
      real(real64), intent(in) :: n
      real(real64), intent(in) :: s
      real(real64), intent(in) :: m
      real(real64), intent(in) :: x

      cdf = 0.0_real64
      if (n /= 0.0_real64 .and. s /= 0.0_real64 .and. &
          m /= 0.0_real64 .and. x > 0.0_real64) then
         cdf = n * 0.5_real64 * &
               (1.0_real64 + erf(log(x / m) / &
                (sqrt(2.0_real64) * log(s))))
      endif
   end function log_normal_cdf

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

   subroutine validate_particle_index(particle_index, count)
      integer, intent(in) :: particle_index
      integer, intent(in) :: count

      if (particle_index < 1 .or. particle_index > count) then
         error stop "LCM particle index is outside the cell container"
      endif
   end subroutine validate_particle_index

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

   subroutine append_particle_to_cell(cell, particle)
      type(lcm_particle_cell_t), intent(inout) :: cell
      type(lcm_particle_t), intent(in) :: particle

      if (cell%count == cell%capacity) then
         call resize_cell_storage(cell, next_capacity(cell%capacity))
      endif

      cell%count = cell%count + 1
      cell%particles(cell%count) = particle
   end subroutine append_particle_to_cell

   logical function particle_leaves_to_neighbor_rank(new_i, new_j, new_k, &
                                                     x_direction, y_direction)
      integer, intent(in) :: new_i
      integer, intent(in) :: new_j
      integer, intent(in) :: new_k
      integer, intent(out) :: x_direction
      integer, intent(out) :: y_direction

      x_direction = 0
      y_direction = 0
      particle_leaves_to_neighbor_rank = .false.

      if (new_k < 1 .or. new_k > grid_state%grid%nz) return

      if (new_i < 1) x_direction = -1
      if (new_i > grid_state%grid%nx_local) x_direction = 1
      if (new_j < 1) y_direction = -1
      if (new_j > grid_state%grid%ny_local) y_direction = 1

      if (x_direction == 0 .and. y_direction == 0) return
      if (x_direction < 0 .and. grid_state%grid%west_rank < 0) return
      if (x_direction > 0 .and. grid_state%grid%east_rank < 0) return
      if (y_direction < 0 .and. grid_state%grid%south_rank < 0) return
      if (y_direction > 0 .and. grid_state%grid%north_rank < 0) return

      particle_leaves_to_neighbor_rank = .true.
   end function particle_leaves_to_neighbor_rank

   subroutine append_outgoing_particle(particle, x_direction, y_direction)
      type(lcm_particle_t), intent(in) :: particle
      integer, intent(in) :: x_direction
      integer, intent(in) :: y_direction

      call ensure_outgoing_particle_capacity(n_outgoing_particles + 1)

      n_outgoing_particles = n_outgoing_particles + 1
      outgoing_particles(n_outgoing_particles)%particle = particle
      outgoing_particles(n_outgoing_particles)%direction = &
         lcm_exchange_direction_from_offsets(x_direction, y_direction)
   end subroutine append_outgoing_particle

   subroutine ensure_outgoing_particle_capacity(required_capacity)
      integer, intent(in) :: required_capacity
      type(lcm_outgoing_particle_t), allocatable :: previous_particles(:)
      integer :: new_capacity

      if (.not. allocated(outgoing_particles)) then
         allocate(outgoing_particles(max(1, required_capacity)))
         return
      endif
      if (required_capacity <= size(outgoing_particles)) return

      new_capacity = max(required_capacity, 2 * size(outgoing_particles))
      call move_alloc(outgoing_particles, previous_particles)
      allocate(outgoing_particles(new_capacity))
      outgoing_particles(1:n_outgoing_particles) = previous_particles(1:n_outgoing_particles)
   end subroutine ensure_outgoing_particle_capacity

end module lcm_particle_state
