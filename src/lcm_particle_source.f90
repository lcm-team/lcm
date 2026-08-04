!> Source-particle injection for standalone LCM.
!!
!! This module owns the standalone counterpart of SAM-LCM micro_particle_source:
!! fixed and moving source regions, source-particle creation, and the IN
!! source profile diagnostic. Particle storage remains cell-owned through
!! lcm_particle_state.
module lcm_particle_source
   use iso_fortran_env, only: int64, real64
   use lcm_config_state, only: config_state
   use lcm_field_state, only: field_state
   use lcm_grid_state, only: grid_state
   use lcm_particle_state, only: cell_particles, cell_absolute_temperature, &
                                cell_tabs_lcm, ice_freezing_temperature, &
                                inverse_error_function_approx, lcm_add_particle, &
                                lcm_next_particle_id, particle_random_seed_id, &
                                qsatw_liquid, set_particle_position
   use lcm_parameters, only: molecular_weight_water, pi, rho_liq, rv, &
                             surface_tension_water
   use lcm_types, only: lcm_particle_t
   use lcm_random, only: randomu

   implicit none

   private

   public :: lcm_apply_particle_sources
   public :: nin_prof

   real(real64), allocatable, protected, save :: nin_prof(:)
   logical, save :: source_position_state_is_initialized = .false.
   real(real64), save :: x_w_source_aero(4) = 0.0_real64
   real(real64), save :: x_e_source_aero(4) = 0.0_real64
   real(real64), save :: y_s_source_aero(4) = 0.0_real64
   real(real64), save :: y_n_source_aero(4) = 0.0_real64
   real(real64), save :: z_b_source_aero(4) = 0.0_real64
   real(real64), save :: z_t_source_aero(4) = 0.0_real64
   real(real64), save :: x_w_source_in(4) = 0.0_real64
   real(real64), save :: x_e_source_in(4) = 0.0_real64
   real(real64), save :: y_s_source_in(4) = 0.0_real64
   real(real64), save :: y_n_source_in(4) = 0.0_real64
   real(real64), save :: z_b_source_in(4) = 0.0_real64
   real(real64), save :: z_t_source_in(4) = 0.0_real64

   integer, parameter :: liquid_particle_type = 1
   integer, parameter :: ice_nucleating_particle_type = 2
   integer, parameter :: source_kind_aerosol = 1
   integer, parameter :: source_kind_in = 2

contains

   subroutine lcm_apply_particle_sources(time_step, model_step)
      real(real64), intent(in) :: time_step
      integer, intent(in) :: model_step

      integer :: i_mode
      real(real64) :: p_source_aero(4)
      real(real64) :: p_source_in(4)

      if (.not. config_state%switch_particle_source) return
      call ensure_source_state()
      nin_prof = 0.0_real64

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

         if (config_state%switch_moving_source) then
            call move_source_positions(i_mode, p_source_aero(i_mode), &
                                       p_source_in(i_mode), time_step)
         endif

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

   subroutine ensure_source_state()
      if (.not. allocated(nin_prof)) then
         allocate(nin_prof(grid_state%grid%nz))
      endif

      if (.not. source_position_state_is_initialized) then
         call initialize_source_position_state()
         source_position_state_is_initialized = .true.
      endif
   end subroutine ensure_source_state

   subroutine initialize_source_position_state()
      integer :: i_mode

      do i_mode = 1, 4
         call initialize_source_position_for_kind(i_mode, source_kind_aerosol)
         call initialize_source_position_for_kind(i_mode, source_kind_in)
      enddo
   end subroutine initialize_source_position_state

   subroutine initialize_source_position_for_kind(i_mode, source_kind)
      integer, intent(in) :: i_mode
      integer, intent(in) :: source_kind

      integer :: nx_e
      integer :: nx_w
      integer :: ny_n
      integer :: ny_s
      integer :: nz_b
      integer :: nz_t

      call fixed_source_bounds(i_mode, source_kind, nx_w, nx_e, ny_s, ny_n, &
                               nz_b, nz_t)

      if (source_kind == source_kind_aerosol) then
         call initialize_source_position_from_bounds(nx_w, nx_e, ny_s, ny_n, &
                                                     nz_b, nz_t,             &
                                                     x_w_source_aero(i_mode),&
                                                     x_e_source_aero(i_mode),&
                                                     y_s_source_aero(i_mode),&
                                                     y_n_source_aero(i_mode),&
                                                     z_b_source_aero(i_mode),&
                                                     z_t_source_aero(i_mode))
      else
         call initialize_source_position_from_bounds(nx_w, nx_e, ny_s, ny_n, &
                                                     nz_b, nz_t,             &
                                                     x_w_source_in(i_mode),  &
                                                     x_e_source_in(i_mode),  &
                                                     y_s_source_in(i_mode),  &
                                                     y_n_source_in(i_mode),  &
                                                     z_b_source_in(i_mode),  &
                                                     z_t_source_in(i_mode))
      endif
   end subroutine initialize_source_position_for_kind

   subroutine initialize_source_position_from_bounds(nx_w, nx_e, ny_s, ny_n, &
                                                     nz_b, nz_t,             &
                                                     x_w, x_e, y_s, y_n,     &
                                                     z_b, z_t)
      integer, intent(in) :: nx_w
      integer, intent(in) :: nx_e
      integer, intent(in) :: ny_s
      integer, intent(in) :: ny_n
      integer, intent(in) :: nz_b
      integer, intent(in) :: nz_t
      real(real64), intent(out) :: x_w
      real(real64), intent(out) :: x_e
      real(real64), intent(out) :: y_s
      real(real64), intent(out) :: y_n
      real(real64), intent(out) :: z_b
      real(real64), intent(out) :: z_t

      x_w = real(nx_w - 1, real64) * grid_state%grid%dx
      x_e = source_upper_bound_position(nx_e, grid_state%grid%dx)
      y_s = real(ny_s - 1, real64) * grid_state%grid%dy
      y_n = source_upper_bound_position(ny_n, grid_state%grid%dy)

      if (nz_b >= 1 .and. nz_b <= grid_state%grid%nz + 1) then
         z_b = grid_state%grid%z_face(nz_b)
      else
         z_b = 0.0_real64
      endif

      if (nz_t >= 0 .and. nz_t < grid_state%grid%nz + 1) then
         z_t = source_vertical_upper_bound(nz_t)
      else
         z_t = 0.0_real64
      endif
   end subroutine initialize_source_position_from_bounds

   real(real64) function source_upper_bound_position(index, spacing) &
      result(position)
      integer, intent(in) :: index
      real(real64), intent(in) :: spacing

      if (index <= 0) then
         position = 0.0_real64
      else
         position = real(index, real64) * spacing - spacing * &
                    epsilon(1.0_real64)
      endif
   end function source_upper_bound_position

   real(real64) function source_vertical_upper_bound(nz_t) result(position)
      integer, intent(in) :: nz_t

      if (nz_t <= 0) then
         position = 0.0_real64
      else
         position = grid_state%grid%z_face(nz_t + 1) - &
                    max(1.0_real64, abs(grid_state%grid%z_face(nz_t + 1))) * &
                    epsilon(1.0_real64)
      endif
   end function source_vertical_upper_bound

   subroutine move_source_positions(i_mode, p_source_aero, p_source_in, time_step)
      integer, intent(in) :: i_mode
      real(real64), intent(in) :: p_source_aero
      real(real64), intent(in) :: p_source_in
      real(real64), intent(in) :: time_step

      if (p_source_aero > 0.0_real64) then
         call move_source_position_set(time_step,                                  &
                                       real(config_state%vx_source_aero(i_mode), real64), &
                                       real(config_state%vy_source_aero(i_mode), real64), &
                                       real(config_state%vz_source_aero(i_mode), real64), &
                                       x_w_source_aero(i_mode), x_e_source_aero(i_mode),  &
                                       y_s_source_aero(i_mode), y_n_source_aero(i_mode),  &
                                       z_b_source_aero(i_mode), z_t_source_aero(i_mode))
      endif

      if (p_source_in > 0.0_real64) then
         call move_source_position_set(time_step,                                  &
                                       real(config_state%vx_source_in(i_mode), real64), &
                                       real(config_state%vy_source_in(i_mode), real64), &
                                       real(config_state%vz_source_in(i_mode), real64), &
                                       x_w_source_in(i_mode), x_e_source_in(i_mode),    &
                                       y_s_source_in(i_mode), y_n_source_in(i_mode),    &
                                       z_b_source_in(i_mode), z_t_source_in(i_mode))
      endif
   end subroutine move_source_positions

   subroutine move_source_position_set(time_step, vx, vy, vz, x_w, x_e, y_s, &
                                       y_n, z_b, z_t)
      real(real64), intent(in) :: time_step
      real(real64), intent(in) :: vx
      real(real64), intent(in) :: vy
      real(real64), intent(in) :: vz
      real(real64), intent(inout) :: x_w
      real(real64), intent(inout) :: x_e
      real(real64), intent(inout) :: y_s
      real(real64), intent(inout) :: y_n
      real(real64), intent(inout) :: z_b
      real(real64), intent(inout) :: z_t

      x_e = mod(x_e + vx * time_step, source_domain_length_x())
      x_w = mod(x_w + vx * time_step, source_domain_length_x())
      y_n = y_n + vy * time_step
      y_s = y_s + vy * time_step
      if (grid_state%grid%y_dimension == 1) then
         y_n = mod(y_n, source_domain_length_y())
         y_s = mod(y_s, source_domain_length_y())
      endif
      z_t = z_t + vz * time_step
      z_b = z_b + vz * time_step
   end subroutine move_source_position_set

   real(real64) function source_domain_length_x() result(length)
      length = real(grid_state%grid%nx_global, real64) * grid_state%grid%dx
   end function source_domain_length_x

   real(real64) function source_domain_length_y() result(length)
      length = real(grid_state%grid%ny_global, real64) * grid_state%grid%dy
   end function source_domain_length_y

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

      call source_bounds(i_mode, source_kind, nx_w, nx_e, ny_s, ny_n, &
                         nz_b, nz_t)

      is_in_region = k >= nz_b .and. k <= nz_t .and. &
                     global_j >= ny_s .and. global_j <= ny_n .and. &
                     global_i >= nx_w .and. global_i <= nx_e
   end function cell_is_in_source_region

   subroutine source_bounds(i_mode, source_kind, nx_w, nx_e, ny_s, ny_n, &
                            nz_b, nz_t)
      integer, intent(in) :: i_mode
      integer, intent(in) :: source_kind
      integer, intent(out) :: nx_w
      integer, intent(out) :: nx_e
      integer, intent(out) :: ny_s
      integer, intent(out) :: ny_n
      integer, intent(out) :: nz_b
      integer, intent(out) :: nz_t

      if (config_state%switch_moving_source) then
         call moving_source_bounds(i_mode, source_kind, nx_w, nx_e, ny_s, &
                                   ny_n, nz_b, nz_t)
      else
         call fixed_source_bounds(i_mode, source_kind, nx_w, nx_e, ny_s, &
                                  ny_n, nz_b, nz_t)
      endif
   end subroutine source_bounds

   subroutine fixed_source_bounds(i_mode, source_kind, nx_w, nx_e, ny_s, ny_n, &
                                  nz_b, nz_t)
      integer, intent(in) :: i_mode
      integer, intent(in) :: source_kind
      integer, intent(out) :: nx_w
      integer, intent(out) :: nx_e
      integer, intent(out) :: ny_s
      integer, intent(out) :: ny_n
      integer, intent(out) :: nz_b
      integer, intent(out) :: nz_t

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
   end subroutine fixed_source_bounds

   subroutine moving_source_bounds(i_mode, source_kind, nx_w, nx_e, ny_s, ny_n, &
                                   nz_b, nz_t)
      integer, intent(in) :: i_mode
      integer, intent(in) :: source_kind
      integer, intent(out) :: nx_w
      integer, intent(out) :: nx_e
      integer, intent(out) :: ny_s
      integer, intent(out) :: ny_n
      integer, intent(out) :: nz_b
      integer, intent(out) :: nz_t

      if (source_kind == source_kind_aerosol) then
         call source_indices_from_positions(x_w_source_aero(i_mode), &
                                            x_e_source_aero(i_mode), &
                                            y_s_source_aero(i_mode), &
                                            y_n_source_aero(i_mode), &
                                            z_b_source_aero(i_mode), &
                                            z_t_source_aero(i_mode), &
                                            nx_w, nx_e, ny_s, ny_n, nz_b, nz_t)
      else
         call source_indices_from_positions(x_w_source_in(i_mode), &
                                            x_e_source_in(i_mode), &
                                            y_s_source_in(i_mode), &
                                            y_n_source_in(i_mode), &
                                            z_b_source_in(i_mode), &
                                            z_t_source_in(i_mode), &
                                            nx_w, nx_e, ny_s, ny_n, nz_b, nz_t)
      endif
   end subroutine moving_source_bounds

   subroutine source_indices_from_positions(x_w, x_e, y_s, y_n, z_b, z_t, &
                                             nx_w, nx_e, ny_s, ny_n, nz_b, nz_t)
      real(real64), intent(in) :: x_w
      real(real64), intent(in) :: x_e
      real(real64), intent(in) :: y_s
      real(real64), intent(in) :: y_n
      real(real64), intent(in) :: z_b
      real(real64), intent(in) :: z_t
      integer, intent(out) :: nx_w
      integer, intent(out) :: nx_e
      integer, intent(out) :: ny_s
      integer, intent(out) :: ny_n
      integer, intent(out) :: nz_b
      integer, intent(out) :: nz_t

      nx_w = int(x_w / grid_state%grid%dx) + 1
      nx_e = int(x_e / grid_state%grid%dx) + 1
      ny_s = int(y_s / grid_state%grid%dy) + 1
      ny_n = int(y_n / grid_state%grid%dy) + 1
      nz_b = search_sorted_z_face(z_b)
      nz_t = search_sorted_z_face(z_t)
   end subroutine source_indices_from_positions

   integer function search_sorted_z_face(z) result(index)
      real(real64), intent(in) :: z

      integer :: finish
      integer :: mid
      integer :: range
      integer :: start

      start = 1
      finish = grid_state%grid%nz + 1
      range = finish - start
      mid = (start + finish) / 2

      do while (range > 1)
         if (z > grid_state%grid%z_face(mid)) then
            start = mid
         else
            finish = mid
         endif
         range = finish - start
         mid = (start + finish) / 2
      enddo

      index = start
   end function search_sorted_z_face

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
         call add_nin_profile(k, weight)
      else
         call maybe_set_embedded_source_in(i_mode, k, weight, seed_id, particle)
      endif

      call lcm_add_particle(i, j, k, particle)
   end subroutine add_one_source_particle

   real(real64) function source_cell_supersaturation(i, j, k) result(ss_gb)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k

      real(real64) :: qsat

      qsat = qsatw_liquid(cell_tabs_lcm(i, j, k), field_state%pressure(k))
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
                                     cell_tabs_lcm(i, j, k), r_aero, &
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
      particle%temperature = cell_absolute_temperature(i, j, k, &
                                                       cell_tabs_lcm(i, j, k))
      particle%qv = field_state%qv(i,j,k)
      particle%eta = field_state%qv(i,j,k) - &
                     qsatw_liquid(cell_tabs_lcm(i, j, k), &
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
      bfactor = source_bfactor_coefficient(1) * mass_aero / real(weight, real64)

      radius_old = 1.0_real64
      radius = max(r_aero, 1.0e-6_real64)
      do while (abs((radius - radius_old) / radius_old) > 1.0e-20_real64)
         radius_old = radius
         radius = (bfactor / (afactor / radius - ss_gb))**(1.0_real64 / 3.0_real64)
      enddo
   end function source_equilibrium_radius

   subroutine maybe_set_embedded_source_in(i_mode, k, source_weight, seed_id, &
                                           particle)
      integer, intent(in) :: i_mode
      integer, intent(in) :: k
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
         call add_nin_profile(k, particle%weight)
      endif
   end subroutine maybe_set_embedded_source_in

   subroutine add_nin_profile(k, weight)
      integer, intent(in) :: k
      integer(int64), intent(in) :: weight

      if (.not. allocated(nin_prof)) then
         error stop "LCM source injection profile is not allocated"
      endif

      nin_prof(k) = nin_prof(k) + real(weight, real64)
   end subroutine add_nin_profile

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


   real(real64) function source_bfactor_coefficient(i_mode) result(coefficient)
      integer, intent(in) :: i_mode

      coefficient = &
         real(config_state%vanthoff_aero(i_mode), real64) / &
         real(config_state%molecular_weight_aero(i_mode), real64) * &
         molecular_weight_water / (rho_liq * 4.0_real64 * pi / 3.0_real64)
   end function source_bfactor_coefficient

end module lcm_particle_source
