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
   use lcm_types, only: lcm_particle_cell_t, lcm_particle_t
   use lcm_user, only: lcm_user_init_particles
   use micro_random, only: randomu

   implicit none

   private

   public :: cell_particles
   public :: lcm_allocate_particle_cells
   public :: lcm_add_particle
   public :: lcm_reset_particle_cells
   public :: lcm_initialize_particle_id_counter
   public :: lcm_next_particle_id
   public :: lcm_initialize_computational_particles
   public :: lcm_apply_particle_sources

   type(lcm_particle_cell_t), allocatable, protected, save :: cell_particles(:,:,:)
   integer(int64), protected, save :: n_step_id = 0_int64
   integer(int64), protected, save :: n_id = 0_int64
   logical, protected, save :: particle_id_counter_is_initialized = .false.

   real(real64), parameter :: pi = 3.141592653589793238462643383279502884_real64
   real(real64), parameter :: epsv = 0.622_real64
   real(real64), parameter :: fac_cond = 2.5e6_real64 / 1004.0_real64
   real(real64), parameter :: molecular_weight_water = 18.0_real64
   real(real64), parameter :: rho_liq = 1000.0_real64
   real(real64), parameter :: rv = 461.5_real64
   real(real64), parameter :: surface_tension_water = 0.0728_real64
   integer, parameter :: liquid_particle_type = 1
   integer, parameter :: ice_nucleating_particle_type = 2
   integer, parameter :: source_kind_aerosol = 1
   integer, parameter :: source_kind_in = 2
   integer, parameter :: n_inv_cdf = 50000000

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

   subroutine lcm_initialize_computational_particles()
      integer :: i, j, k, n
      integer :: initialization_levels
      type(lcm_particle_t) :: particle

      call require_particle_cells_allocated()
      call lcm_reset_particle_cells()
      call lcm_initialize_particle_id_counter()

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

   subroutine lcm_apply_particle_sources(time_step, model_step)
      real(real64), intent(in) :: time_step
      integer, intent(in) :: model_step

      integer :: i_mode
      real(real64) :: p_source_aero(4)
      real(real64) :: p_source_in(4)

      if (.not. config_state%switch_particle_source) return
      if (config_state%switch_moving_source) then
         error stop "LCM moving particle sources need standalone source-position state"
      endif

      p_source_aero = real(config_state%p_source_aero, real64)
      p_source_in = real(config_state%p_source_in, real64)

      if (source_is_inactive(model_step, config_state%n_start_source_aero, &
                             config_state%n_end_source_aero)) then
         p_source_aero = 0.0_real64
      endif
      if (source_is_inactive(model_step, config_state%n_start_source_in, &
                             config_state%n_end_source_in) .and. &
          config_state%switch_phase_sources .and. config_state%switch_ice) then
         p_source_in = 0.0_real64
      endif

      do i_mode = 1, 4
         if (max(p_source_aero(i_mode), p_source_in(i_mode)) <= 0.0_real64) exit

         if (p_source_aero(i_mode) > 0.0_real64) then
            call add_source_particles_for_region(i_mode, p_source_aero(i_mode), &
                                                 time_step, source_kind_aerosol)
         endif

         if (config_state%switch_phase_sources .and. config_state%switch_ice .and. &
             p_source_in(i_mode) > 0.0_real64) then
            call add_source_particles_for_region(i_mode, p_source_in(i_mode), &
                                                 time_step, source_kind_in)
         endif
      enddo
   end subroutine lcm_apply_particle_sources

   logical function source_is_inactive(model_step, n_start, n_end) result(is_inactive)
      integer, intent(in) :: model_step
      integer, intent(in) :: n_start
      integer, intent(in) :: n_end

      is_inactive = model_step < n_start .or. &
                    (model_step > n_end .and. n_end > 0) .or. &
                    n_start < 0
   end function source_is_inactive

   subroutine add_source_particles_for_region(i_mode, p_source, time_step, source_kind)
      integer, intent(in) :: i_mode
      real(real64), intent(in) :: p_source
      real(real64), intent(in) :: time_step
      integer, intent(in) :: source_kind

      integer :: i, j, k

      do k = 1, grid_state%grid%nz
         do j = 1, grid_state%grid%ny_local
            do i = 1, grid_state%grid%nx_local
               if (cell_is_in_source_region(i, j, k, i_mode, source_kind)) then
                  call add_source_particles_for_cell(i, j, k, i_mode, p_source, &
                                                     time_step, source_kind)
               endif
            enddo
         enddo
      enddo
   end subroutine add_source_particles_for_region

   logical function cell_is_in_source_region(i, j, k, i_mode, source_kind) &
      result(is_in_region)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: i_mode
      integer, intent(in) :: source_kind

      integer :: global_i
      integer :: global_j
      integer :: nx_w
      integer :: nx_e
      integer :: ny_s
      integer :: ny_n
      integer :: nz_b
      integer :: nz_t

      global_i = grid_state%grid%global_i_start + i - 1
      global_j = grid_state%grid%global_j_start + j - 1

      if (source_kind == source_kind_aerosol) then
         nx_w = config_state%nx_w_source_aero(i_mode)
         nx_e = config_state%nx_e_source_aero(i_mode)
         ny_s = config_state%ny_s_source_aero(i_mode)
         ny_n = config_state%ny_n_source_aero(i_mode)
         nz_b = config_state%nz_b_source_aero(i_mode)
         nz_t = config_state%nz_t_source_aero(i_mode)
      else
         nx_w = config_state%nx_w_source_in(i_mode)
         nx_e = config_state%nx_e_source_in(i_mode)
         ny_s = config_state%ny_s_source_in(i_mode)
         ny_n = config_state%ny_n_source_in(i_mode)
         nz_b = config_state%nz_b_source_in(i_mode)
         nz_t = config_state%nz_t_source_in(i_mode)
      endif

      is_in_region = k >= nz_b .and. k <= nz_t .and. &
                     global_j >= ny_s .and. global_j <= ny_n .and. &
                     global_i >= nx_w .and. global_i <= nx_e
   end function cell_is_in_source_region

   subroutine add_source_particles_for_cell(i, j, k, i_mode, p_source, &
                                            time_step, source_kind)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: i_mode
      real(real64), intent(in) :: p_source
      real(real64), intent(in) :: time_step
      integer, intent(in) :: source_kind

      integer :: n
      integer :: n_source_particles
      integer :: seed_id
      integer(int64) :: weight
      real(real64) :: x_rand

      if (p_source <= 0.0_real64) return

      seed_id = particle_random_seed_id()
      weight = source_particle_weight(i_mode, k, p_source, time_step, source_kind)
      if (weight == 0_int64) then
         error stop "LCM source produces particles with 0 weighting factor"
      endif

      x_rand = real(randomu(seed_id), real64)
      if (x_rand < p_source - floor(p_source)) then
         n_source_particles = ceiling(p_source)
      else
         n_source_particles = floor(p_source)
      endif

      do n = 1, n_source_particles
         call add_one_source_particle(i, j, k, i_mode, weight, source_kind, seed_id)
      enddo
   end subroutine add_source_particles_for_cell

   integer(int64) function source_particle_weight(i_mode, k, p_source, &
                                                  time_step, source_kind) &
      result(weight)
      integer, intent(in) :: i_mode
      integer, intent(in) :: k
      real(real64), intent(in) :: p_source
      real(real64), intent(in) :: time_step
      integer, intent(in) :: source_kind

      real(real64) :: source_strength

      if (source_kind == source_kind_aerosol) then
         source_strength = real(config_state%n_source_aero(i_mode), real64)
      else
         source_strength = real(config_state%n_source_in(i_mode), real64)
      endif

      if (trim(config_state%unit_source) == "1/m3/s") then
         weight = nint(source_strength * real(grid_state%dV(k), real64) * &
                       time_step / p_source, int64)
      else if (trim(config_state%unit_source) == "1/kg/s") then
         weight = nint(source_strength * real(grid_state%dV(k), real64) * &
                       field_state%density(k) * time_step / p_source, int64)
      else
         error stop "LCM unknown particle source unit"
      endif
   end function source_particle_weight

   subroutine add_one_source_particle(i, j, k, i_mode, weight, source_kind, seed_id)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: i_mode
      integer(int64), intent(in) :: weight
      integer, intent(in) :: source_kind
      integer, intent(in) :: seed_id

      type(lcm_particle_t) :: particle
      real(real64) :: r_aero
      real(real64) :: r_in
      real(real64) :: sfc_in
      real(real64) :: ss_gb

      particle = lcm_particle_t()
      particle%spray_marker = 100
      particle%cell_particle_index = cell_particles(i,j,k)%count + 1
      particle%last_i = grid_state%grid%global_i_start + i - 1
      particle%last_j = grid_state%grid%global_j_start + j - 1
      particle%last_k = k
      particle%id = lcm_next_particle_id()
      particle%cloudy_age = 0.0_real64
      call set_particle_position(i, j, k, seed_id, particle)
      particle%u = 0.0_real64
      particle%v = 0.0_real64
      particle%w = 0.0_real64
      particle%u_sgs = 0.0_real64
      particle%v_sgs = 0.0_real64
      particle%w_sgs = 0.0_real64
      particle%w_sedi = 0.0_real64

      ss_gb = source_cell_supersaturation(i, j, k)
      r_aero = source_aerosol_radius(i_mode, seed_id)
      call set_source_liquid_aerosol_state(i, j, k, i_mode, weight, r_aero, &
                                           ss_gb, particle)

      if (source_kind == source_kind_in) then
         r_in = source_in_radius(i_mode, seed_id, 100.0e-6_real64)
         sfc_in = 4.0_real64 * pi * r_in**2
         call set_source_in_state(i_mode, weight, sfc_in, seed_id, particle)
      else
         call maybe_set_embedded_source_in(i_mode, weight, seed_id, particle)
      endif

      call lcm_add_particle(i, j, k, particle)
   end subroutine add_one_source_particle

   real(real64) function source_cell_supersaturation(i, j, k) result(ss_gb)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k

      real(real64) :: qsat

      qsat = qsatw_liquid(field_state%temperature(i,j,k), field_state%pressure(k))
      ss_gb = (field_state%qv(i,j,k) - qsat) / qsat
      ss_gb = min(ss_gb, -0.01_real64)
   end function source_cell_supersaturation

   real(real64) function source_aerosol_radius(i_mode, seed_id) result(r_aero)
      integer, intent(in) :: i_mode
      integer, intent(in) :: seed_id

      r_aero = source_lognormal_radius(real(config_state%r_source_aero(i_mode), real64), &
                                       real(config_state%sigma_source_aero(i_mode), real64), &
                                       0.001e-6_real64, 10.0e-6_real64, seed_id)
   end function source_aerosol_radius

   real(real64) function source_in_radius(i_mode, seed_id, upper_bound) result(r_in)
      integer, intent(in) :: i_mode
      integer, intent(in) :: seed_id
      real(real64), intent(in) :: upper_bound

      r_in = source_lognormal_radius(real(config_state%r_source_in(i_mode), real64), &
                                     real(config_state%sigma_source_in(i_mode), real64), &
                                     0.001e-6_real64, upper_bound, seed_id)
   end function source_in_radius

   real(real64) function source_lognormal_radius(r_mean, sigma, lower_bound, &
                                                 upper_bound, seed_id) result(radius)
      real(real64), intent(in) :: r_mean
      real(real64), intent(in) :: sigma
      real(real64), intent(in) :: lower_bound
      real(real64), intent(in) :: upper_bound
      integer, intent(in) :: seed_id

      real(real64) :: ierf_dummy
      real(real64) :: x_rand

      radius = -1.0_real64
      do while (radius <= lower_bound .or. radius >= upper_bound)
         x_rand = real(randomu(seed_id), real64)
         ierf_dummy = min(max(inverse_error_function_approx( &
                          2.0_real64 * x_rand - 1.0_real64), &
                          -2.0_real64), 2.0_real64)
         radius = r_mean * exp(sqrt(2.0_real64) * log(sigma) * ierf_dummy)
      enddo
   end function source_lognormal_radius

   subroutine set_source_liquid_aerosol_state(i, j, k, i_mode, weight, r_aero, &
                                              ss_gb, particle)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: i_mode
      integer(int64), intent(in) :: weight
      real(real64), intent(in) :: r_aero
      real(real64), intent(in) :: ss_gb
      type(lcm_particle_t), intent(inout) :: particle

      particle%weight = weight
      particle%mass_aero = real(weight, real64) * 4.0_real64 / 3.0_real64 * &
                           pi * real(config_state%rho_aero(1), real64) * r_aero**3
      if (trim(config_state%type_source) == "equilibrium") then
         particle%mass_water = real(weight, real64) * 4.0_real64 / 3.0_real64 * &
                               pi * rho_liq * &
                               max(source_equilibrium_radius(ss_gb, &
                                     field_state%temperature(i,j,k), r_aero, &
                                     particle%mass_aero, weight), r_aero)**3
      else if (trim(config_state%type_source) == "seawater") then
         particle%mass_water = real(weight, real64) * 4.0_real64 / 3.0_real64 * &
                               pi * rho_liq * &
                               (r_aero**3 / (real(config_state%salinity_source, real64) * &
                                (rho_liq * (1.0_real64 + &
                                 real(config_state%salinity_source, real64))) / &
                                real(config_state%rho_aero(1), real64)))
      else
         error stop "LCM unknown particle source type"
      endif

      particle%aerosol_portions = 0.0_real64
      particle%aerosol_portions(1) = 1.0_real64
      particle%temperature = field_state%temperature(i,j,k)
      particle%qv = field_state%qv(i,j,k)
      particle%eta = field_state%qv(i,j,k) - &
                     qsatw_liquid(field_state%temperature(i,j,k), &
                                  field_state%pressure(k))
      particle%particle_type = liquid_particle_type
      particle%in_surface_area = 0.0_real64
      particle%freezing_temperature = 233.15_real64
      particle%n_monomer = 1_int64
   end subroutine set_source_liquid_aerosol_state

   real(real64) function source_equilibrium_radius(ss_gb, temperature, r_aero, &
                                                   mass_aero, weight) result(radius)
      real(real64), intent(in) :: ss_gb
      real(real64), intent(in) :: temperature
      real(real64), intent(in) :: r_aero
      real(real64), intent(in) :: mass_aero
      integer(int64), intent(in) :: weight

      real(real64) :: afactor
      real(real64) :: bfactor
      real(real64) :: radius_old

      afactor = 2.0_real64 * surface_tension_water / &
                (rho_liq * rv * temperature)
      bfactor = bfactor_coefficient(1) * mass_aero / real(weight, real64)

      radius_old = 1.0_real64
      radius = max(r_aero, 1.0e-6_real64)
      do while (abs((radius - radius_old) / radius_old) > 1.0e-20_real64)
         radius_old = radius
         radius = (bfactor / (afactor / radius - ss_gb))**(1.0_real64 / 3.0_real64)
      enddo
   end function source_equilibrium_radius

   subroutine maybe_set_embedded_source_in(i_mode, source_weight, seed_id, particle)
      integer, intent(in) :: i_mode
      integer(int64), intent(in) :: source_weight
      integer, intent(in) :: seed_id
      type(lcm_particle_t), intent(inout) :: particle

      real(real64) :: r_in
      real(real64) :: sfc_in
      real(real64) :: x_rand

      x_rand = real(randomu(seed_id), real64)
      if (x_rand < real(config_state%frac_in_source_particles(i_mode), real64) .and. &
          config_state%switch_ice .and. (.not. config_state%switch_phase_sources)) then
         r_in = source_in_radius(i_mode, seed_id, 10.0e-4_real64)
         sfc_in = 4.0_real64 * pi * r_in**2

         particle%weight = nint(real(source_weight, real64) * &
                                real(config_state%frac_in_source_weight(i_mode), real64), &
                                int64)
         particle%mass_aero = particle%mass_aero * &
                              real(config_state%frac_in_source_weight(i_mode), real64)
         particle%mass_water = particle%mass_water * &
                               real(config_state%frac_in_source_weight(i_mode), real64)
         call set_source_in_state(i_mode, particle%weight, sfc_in, seed_id, particle)
      endif
   end subroutine maybe_set_embedded_source_in

   subroutine set_source_in_state(i_mode, weight, sfc_in, seed_id, particle)
      integer, intent(in) :: i_mode
      integer(int64), intent(in) :: weight
      real(real64), intent(in) :: sfc_in
      integer, intent(in) :: seed_id
      type(lcm_particle_t), intent(inout) :: particle

      integer :: in_type

      in_type = config_state%type_source_in(i_mode)
      particle%particle_type = ice_nucleating_particle_type
      particle%in_surface_area = 0.0_real64
      particle%in_surface_area(in_type) = real(weight, real64) * sfc_in
      particle%freezing_temperature = ice_freezing_temperature(in_type, sfc_in, seed_id)
   end subroutine set_source_in_state

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
      particle%track_marker = -999

      particle%particle_type = liquid_particle_type
      seed_id = particle_random_seed_id()
      call set_particle_position(i, j, k, seed_id, particle)
      particle%w_sedi = 0.0_real64
      particle%temperature = field_state%temperature(i,j,k)
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

      call initialize_equilibrium_liquid_state(i, j, k)

      call lcm_user_init_particles(i, j, k, cell_particles(i,j,k))
   end subroutine lcm_init_particles

   integer function particle_random_seed_id() result(seed_id)
      integer :: subdomain_id

      subdomain_id = grid_state%grid%subdomain_y_index * &
                     grid_state%grid%n_subdomains_x + &
                     grid_state%grid%subdomain_x_index
      seed_id = subdomain_id + 1
   end function particle_random_seed_id

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
      cell_temperature = field_state%temperature(i,j,k)

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
         cell_particles(i,j,k)%particles(n)%temperature = cell_temperature
         cell_particles(i,j,k)%particles(n)%qv = cell_qv
         cell_particles(i,j,k)%particles(n)%eta = &
            cell_qv - qsatw_liquid(cell_temperature, field_state%pressure(k))
         cell_particles(i,j,k)%particles(n)%n_monomer = 1_int64
      enddo
   end subroutine initialize_equilibrium_liquid_state

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
