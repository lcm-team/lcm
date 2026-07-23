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

   implicit none

   private

   public :: cell_particles
   public :: lcm_allocate_particle_cells
   public :: lcm_add_particle
   public :: lcm_reset_particle_cells
   public :: lcm_initialize_particle_id_counter
   public :: lcm_next_particle_id
   public :: lcm_initialize_computational_particles

   type(lcm_particle_cell_t), allocatable, protected, save :: cell_particles(:,:,:)
   integer(int64), protected, save :: n_step_id = 0_int64
   integer(int64), protected, save :: n_id = 0_int64
   logical, protected, save :: particle_id_counter_is_initialized = .false.

   real(real64), parameter :: pi = 3.141592653589793238462643383279502884_real64
   integer, parameter :: liquid_particle_type = 1

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

      particle%cell_particle_index = cell_particle_index
      particle%last_i = grid_state%grid%global_i_start + i - 1
      particle%last_j = grid_state%grid%global_j_start + j - 1
      particle%last_k = k
      particle%spray_marker = -999
      particle%track_marker = -999

      particle%particle_type = liquid_particle_type
      particle%w_sedi = 0.0_real64
      particle%temperature = field_state%temperature(i,j,k)
      particle%qv = field_state%qv(i,j,k)
   end subroutine init_particle_properties

   subroutine lcm_init_particles(i, j, k)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k

      integer :: n, nn, m, m_max, i_mode, np_mode, np_done
      real(real64) :: dlogr, r_l, r_r, r_aero, n_dummy
      real(real64) :: f_ak, f_bk, f_ck, ak, bk, ck
      real(real64) :: r_aero_min(4), r_aero_max(4)

      if (config_state%np_xy == 0) return

      if (config_state%switch_ice) then
         error stop "LCM ice initialization must be ported from SAM-LCM micro_init_particles"
      endif

      if (trim(config_state%init_aero) /= "bin") then
         error stop "LCM random and USM aerosol initialization must be ported from SAM-LCM micro_init_particles"
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

            cell_particles(i,j,k)%particles(n)%mass_aero = &
               real(cell_particles(i,j,k)%particles(n)%weight, real64) * &
               4.0_real64 / 3.0_real64 * pi * r_aero**3 * &
               real(config_state%rho_aero(i_mode), real64)

            cell_particles(i,j,k)%particles(n)%aerosol_portions = 0.0_real64
            cell_particles(i,j,k)%particles(n)%aerosol_portions(i_mode) = 1.0_real64
            cell_particles(i,j,k)%particles(n)%particle_type = liquid_particle_type
            cell_particles(i,j,k)%particles(n)%in_surface_area = 0.0_real64
            cell_particles(i,j,k)%particles(n)%freezing_temperature = 233.15_real64
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

         np_done = np_done + np_mode
      enddo

      do n = 1, cell_particles(i,j,k)%count
         cell_particles(i,j,k)%particles(n)%temperature = field_state%temperature(i,j,k)
         cell_particles(i,j,k)%particles(n)%qv = field_state%qv(i,j,k)
      enddo
   end subroutine lcm_init_particles

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
