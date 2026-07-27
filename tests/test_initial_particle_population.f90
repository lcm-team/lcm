program test_initial_particle_population
   use iso_fortran_env, only: int64, real64
   use lcm_host_interface, only: lcm_attach_fields, lcm_config_t, lcm_grid_t, &
                                 lcm_init_config, lcm_init_grid,             &
                                 lcm_initialize, lcm_set_config, lcm_set_grid
   use lcm_initialization_state, only: lcm_initial_particle_count_global, &
                                       lcm_initial_particle_count_local
   use lcm_particle_state, only: cell_particles

   implicit none

   integer, parameter :: nx_local = 2
   integer, parameter :: ny_local = 2
   integer, parameter :: nz = 3
   integer, parameter :: nzl = 2
   integer, parameter :: np_xy = 8
   real(real64), parameter :: pi = 3.141592653589793238462643383279502884_real64
   integer(int64), parameter :: expected_count = int(nx_local, int64) * &
                                                 int(ny_local, int64) * &
                                                 int(nzl, int64) * &
                                                 int(np_xy, int64)

   type(lcm_config_t) :: config
   type(lcm_grid_t) :: grid
   real(real64), target :: u(nx_local, ny_local, nz)
   real(real64), target :: v(nx_local, ny_local, nz)
   real(real64), target :: w(nx_local, ny_local, nz)
   real(real64), target :: static_energy(nx_local, ny_local, nz)
   real(real64), target :: temperature(nx_local, ny_local, nz)
   real(real64), target :: qv(nx_local, ny_local, nz)
   real(real64), target :: pressure(nz)
   real(real64), target :: density(nz)
   real(real64) :: dz(nz)
   real(real64) :: z_center(nz)
   real(real64) :: z_face(nz + 1)
   character(len=32) :: case_name

   call get_command_argument(1, case_name)
   if (len_trim(case_name) == 0) case_name = "one_mode"

   call lcm_init_config(config)
   call configure_case(trim(case_name), config)
   call lcm_set_config(config)

   dz = [100.0_real64, 100.0_real64, 100.0_real64]
   z_center = [50.0_real64, 150.0_real64, 250.0_real64]
   z_face = [0.0_real64, 100.0_real64, 200.0_real64, 300.0_real64]

   call lcm_init_grid(grid, nx_local, ny_local, nz, nx_local, ny_local, &
                      1, nx_local, 1, ny_local,                         &
                      50.0_real64, 50.0_real64, z_center, z_face, dz,   &
                      0, 0, 0, 0, 0, 1, 1, 1, 1,                        &
                      -1, -1, -1, -1, .true., .true., 1)
   call lcm_set_grid(grid)

   u = 0.0_real64
   v = 0.0_real64
   w = 0.0_real64
   static_energy = 300.0_real64
   temperature = 280.0_real64
   qv = 0.010_real64
   pressure = [100000.0_real64, 99000.0_real64, 98000.0_real64]
   density = [1.2_real64, 1.1_real64, 1.0_real64]

   call lcm_attach_fields(u=u, v=v, w=w, static_energy=static_energy, &
                          temperature=temperature, qv=qv,             &
                          pressure=pressure, density=density)
   call lcm_initialize()

   call assert_true(lcm_initial_particle_count_local == expected_count, &
                    "local initial particle count")
   call assert_true(lcm_initial_particle_count_global == expected_count, &
                    "global initial particle count")
   call assert_particle_cells()

contains

   subroutine configure_case(case_name, config)
      character(len=*), intent(in) :: case_name
      type(lcm_config_t), intent(inout) :: config

      config%np_xy = np_xy
      config%nzl = nzl

      select case (case_name)
      case ("default_aerosol")
         continue
      case ("one_mode")
         config%init_aero = "bin"
         config%unit_aero = "1/kg"
         config%switch_ice = .false.
         config%n_aero = [100.0e6, 0.0, 0.0, 0.0]
         config%rm_aero = [0.04e-6, 0.0, 0.0, 0.0]
         config%sigma_aero = [1.4, 0.0, 0.0, 0.0]
         config%shares_aero = [1.0, 0.0, 0.0, 0.0]
         config%rho_aero = [2170.0, 0.0, 0.0, 0.0]
      case ("multi_mode")
         config%init_aero = "bin"
         config%unit_aero = "1/kg"
         config%switch_ice = .false.
         config%n_aero = [100.0e6, 50.0e6, 25.0e6, 0.0]
         config%rm_aero = [0.04e-6, 0.08e-6, 0.12e-6, 0.0]
         config%sigma_aero = [1.4, 1.6, 1.8, 0.0]
         config%shares_aero = [1.0, 2.0, 1.0, 0.0]
         config%rho_aero = [2170.0, 1800.0, 1500.0, 0.0]
      case ("zero_weight_split")
         config%init_aero = "bin"
         config%unit_aero = "1/kg"
         config%switch_ice = .false.
         config%n_aero = [0.1, 0.0, 0.0, 0.0]
         config%rm_aero = [0.04e-6, 0.0, 0.0, 0.0]
         config%sigma_aero = [1.4, 0.0, 0.0, 0.0]
         config%shares_aero = [1.0, 0.0, 0.0, 0.0]
         config%rho_aero = [2170.0, 0.0, 0.0, 0.0]
      case ("usm_aerosol")
         config%init_aero = "USM"
         config%unit_aero = "1/kg"
         config%switch_ice = .false.
         config%n_aero = [100.0e6, 0.0, 0.0, 0.0]
         config%rm_aero = [0.04e-6, 0.0, 0.0, 0.0]
         config%sigma_aero = [1.4, 0.0, 0.0, 0.0]
         config%shares_aero = [1.0, 0.0, 0.0, 0.0]
         config%rho_aero = [2170.0, 0.0, 0.0, 0.0]
      case ("random_aerosol")
         config%init_aero = "random"
         config%unit_aero = "1/kg"
         config%switch_ice = .false.
         config%n_aero = [100.0e6, 0.0, 0.0, 0.0]
         config%rm_aero = [0.04e-6, 0.0, 0.0, 0.0]
         config%sigma_aero = [1.4, 0.0, 0.0, 0.0]
         config%shares_aero = [1.0, 0.0, 0.0, 0.0]
         config%rho_aero = [2170.0, 0.0, 0.0, 0.0]
      case ("bin_ice")
         config%init_aero = "bin"
         config%unit_aero = "1/kg"
         config%switch_ice = .true.
         config%frac_ice = 0.999999
         config%n_aero = [100.0e6, 0.0, 0.0, 0.0]
         config%rm_aero = [0.04e-6, 0.0, 0.0, 0.0]
         config%sigma_aero = [1.4, 0.0, 0.0, 0.0]
         config%shares_aero = [1.0, 0.0, 0.0, 0.0]
         config%rho_aero = [2170.0, 0.0, 0.0, 0.0]
         config%n_ice = [0.0, 0.0, 0.001e6, 0.0]
         config%rm_in = [1.0e-6, 0.2e-6, 1.0e-6, 0.0]
         config%sigma_in = [1.4, 1.4, 1.4, 0.0]
      case default
         error stop "unknown initial particle population test case"
      end select
   end subroutine configure_case

   subroutine assert_particle_cells()
      integer :: i, j, k, n
      integer :: total_count

      total_count = 0
      do k = 1, nz
         do j = 1, ny_local
            do i = 1, nx_local
               if (k <= nzl) then
                  call assert_true(cell_particles(i,j,k)%count == np_xy, &
                                   "initialized cell particle count")
                  call assert_true(cell_particles(i,j,k)%capacity >= &
                                   cell_particles(i,j,k)%count, &
                                   "cell-local active particle range")
                  do n = 1, cell_particles(i,j,k)%count
                     total_count = total_count + 1
                     call assert_particle(i, j, k, n)
                  enddo
               else
                  call assert_true(cell_particles(i,j,k)%count == 0, &
                                   "uninitialized cell particle count")
               endif
            enddo
         enddo
      enddo

      call assert_true(total_count == expected_count, "summed particle count")
      call assert_unique_ids()
   end subroutine assert_particle_cells

   subroutine assert_particle(i, j, k, n)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n
      integer(int64) :: expected_weight(np_xy)
      real(real64) :: expected_aerosol_mass(np_xy)
      real(real64) :: expected_portions(4, np_xy)

      expected_portions = 0.0_real64

      select case (trim(config%init_aero))
      case ("bin")
         if (config%switch_ice) then
            expected_portions(1, :) = 1.0_real64
         else
            call reference_bin_aerosol(k, expected_weight, expected_aerosol_mass, &
                                       expected_portions)
            call assert_true(cell_particles(i,j,k)%particles(n)%weight == &
                             expected_weight(n), "SAM-LCM bin aerosol weight")
            call assert_near_relative(cell_particles(i,j,k)%particles(n)%mass_aero, &
                                      expected_aerosol_mass(n), &
                                      "SAM-LCM bin aerosol mass")
         endif
      case ("USM")
         call reference_usm_aerosol(k, expected_weight, expected_aerosol_mass, &
                                    expected_portions)
         call assert_true(cell_particles(i,j,k)%particles(n)%weight == &
                          expected_weight(n), "SAM-LCM USM aerosol weight")
         call assert_near_relative(cell_particles(i,j,k)%particles(n)%mass_aero, &
                                   expected_aerosol_mass(n), &
                                   "SAM-LCM USM aerosol mass")
      case ("random")
         call reference_random_aerosol(k, expected_weight, expected_portions)
         call assert_true(cell_particles(i,j,k)%particles(n)%weight == &
                          expected_weight(n), "SAM-LCM random aerosol weight")
         call assert_true(cell_particles(i,j,k)%particles(n)%mass_aero > 0.0_real64, &
                          "SAM-LCM random aerosol mass")
      case default
         error stop "unknown aerosol initialization mode"
      end select

      call assert_true(cell_particles(i,j,k)%particles(n)%weight > 0_int64, &
                       "particle weight")
      call assert_true(cell_particles(i,j,k)%particles(n)%mass_aero > 0.0_real64, &
                       "particle aerosol mass")
      if (config%switch_ice) then
         call assert_true(cell_particles(i,j,k)%particles(n)%particle_type == 2, &
                          "ice-nucleating particle type")
         call assert_true(cell_particles(i,j,k)%particles(n)%in_surface_area(3) > &
                          0.0_real64, "ice-nucleating particle surface area")
         call assert_near(maxval(abs(cell_particles(i,j,k)%particles(n)%in_surface_area([1,2,4]))), &
                          0.0_real64, "inactive IN surface areas")
         call assert_true(cell_particles(i,j,k)%particles(n)%freezing_temperature >= &
                          233.15_real64, "ice freezing temperature lower bound")
         call assert_true(cell_particles(i,j,k)%particles(n)%freezing_temperature <= &
                          273.15_real64, "ice freezing temperature upper bound")
      else
         call assert_true(cell_particles(i,j,k)%particles(n)%particle_type == 1, &
                          "particle type")
         call assert_near(cell_particles(i,j,k)%particles(n)%freezing_temperature, &
                          233.15_real64, "default homogeneous freezing temperature")
         call assert_near(maxval(abs(cell_particles(i,j,k)%particles(n)%in_surface_area)), &
                          0.0_real64, "default IN surface area")
      endif
      call assert_near(maxval(abs(cell_particles(i,j,k)%particles(n)%aerosol_portions - &
                                  expected_portions(:, n))), &
                       0.0_real64, "SAM-LCM aerosol portions")
      call assert_near(cell_particles(i,j,k)%particles(n)%temperature, &
                       temperature(i,j,k), "particle temperature")
      call assert_near(cell_particles(i,j,k)%particles(n)%qv, qv(i,j,k), &
                       "particle qv")
   end subroutine assert_particle

   subroutine reference_bin_aerosol(k, expected_weight, expected_aerosol_mass, &
                                    expected_portions)
      integer, intent(in) :: k
      integer(int64), intent(out) :: expected_weight(np_xy)
      real(real64), intent(out) :: expected_aerosol_mass(np_xy)
      real(real64), intent(out) :: expected_portions(4, np_xy)

      integer :: i_mode
      integer :: m
      integer :: m_max
      integer :: n
      integer :: nn
      integer :: np_done
      integer :: np_mode
      real(real64) :: ak
      real(real64) :: bk
      real(real64) :: ck
      real(real64) :: dlogr
      real(real64) :: f_ak
      real(real64) :: f_bk
      real(real64) :: f_ck
      real(real64) :: n_dummy
      real(real64) :: r_aero
      real(real64) :: r_aero_max(4)
      real(real64) :: r_aero_min(4)
      real(real64) :: r_l
      real(real64) :: r_r
      real(real64) :: volume

      expected_weight = 0_int64
      expected_aerosol_mass = 0.0_real64
      expected_portions = 0.0_real64
      volume = 50.0_real64 * 50.0_real64 * dz(k)

      np_done = 0
      do i_mode = 1, 4
         np_mode = floor(real(np_xy, real64) * &
                         real(config%shares_aero(i_mode), real64) / &
                         sum(real(config%shares_aero, real64)))
         if (i_mode == 1) then
            n = 0
            do nn = 1, 4
               n = n + floor(real(np_xy, real64) * &
                             real(config%shares_aero(nn), real64) / &
                             sum(real(config%shares_aero, real64)))
            enddo
            np_mode = np_mode + np_xy - n
         endif

         if (np_mode == 0) cycle
         call assert_true(config%n_aero(i_mode) /= 0.0, &
                          "reference aerosol mode concentration")

         r_aero_min(i_mode) = 1.0e-9_real64
         r_aero_max(i_mode) = 1.0e-5_real64

         ak = r_aero_min(i_mode)
         f_ak = reference_log_normal_cdf(real(config%n_aero(i_mode), real64), &
                                         real(config%sigma_aero(i_mode), real64), &
                                         real(config%rm_aero(i_mode), real64), ak) / &
                real(config%n_aero(i_mode), real64) - 5.0e-5_real64

         bk = r_aero_max(i_mode)
         f_bk = reference_log_normal_cdf(real(config%n_aero(i_mode), real64), &
                                         real(config%sigma_aero(i_mode), real64), &
                                         real(config%rm_aero(i_mode), real64), bk) / &
                real(config%n_aero(i_mode), real64) - 5.0e-5_real64

         do
            ck = (ak * f_bk - bk * f_ak) / (f_bk - f_ak)
            f_ck = reference_log_normal_cdf(real(config%n_aero(i_mode), real64), &
                                            real(config%sigma_aero(i_mode), real64), &
                                            real(config%rm_aero(i_mode), real64), ck) / &
                   real(config%n_aero(i_mode), real64) - 5.0e-5_real64

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
         f_ak = reference_log_normal_cdf(real(config%n_aero(i_mode), real64), &
                                         real(config%sigma_aero(i_mode), real64), &
                                         real(config%rm_aero(i_mode), real64), ak) / &
                real(config%n_aero(i_mode), real64) - &
                (1.0_real64 - 5.0e-5_real64)

         bk = r_aero_max(i_mode)
         f_bk = reference_log_normal_cdf(real(config%n_aero(i_mode), real64), &
                                         real(config%sigma_aero(i_mode), real64), &
                                         real(config%rm_aero(i_mode), real64), bk) / &
                real(config%n_aero(i_mode), real64) - &
                (1.0_real64 - 5.0e-5_real64)

         do
            ck = (ak * f_bk - bk * f_ak) / (f_bk - f_ak)
            f_ck = reference_log_normal_cdf(real(config%n_aero(i_mode), real64), &
                                            real(config%sigma_aero(i_mode), real64), &
                                            real(config%rm_aero(i_mode), real64), ck) / &
                   real(config%n_aero(i_mode), real64) - &
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

            n_dummy = reference_log_normal_pdf(real(config%n_aero(i_mode), real64), &
                                               real(config%sigma_aero(i_mode), real64), &
                                               real(config%rm_aero(i_mode), real64), &
                                               r_aero)
            if (trim(config%unit_aero) == "1/m3") then
               expected_weight(n) = int(n_dummy * r_aero * dlogr * volume, int64)
            else if (trim(config%unit_aero) == "1/kg") then
               expected_weight(n) = int(n_dummy * r_aero * dlogr * volume * &
                                        density(k), int64)
            else
               error stop "unknown reference aerosol unit"
            endif

            expected_aerosol_mass(n) = real(expected_weight(n), real64) * &
                                       4.0_real64 / 3.0_real64 * pi * &
                                       r_aero**3 * &
                                       real(config%rho_aero(i_mode), real64)
            expected_portions(:, n) = 0.0_real64
            expected_portions(i_mode, n) = 1.0_real64
         enddo

         do n = np_done + 1, np_done + np_mode
            if (expected_weight(n) == 0_int64) then
               m_max = n
               do m = np_done + 1, np_done + np_mode
                  if (expected_weight(m) > expected_weight(m_max) .and. &
                      expected_weight(m) >= 2_int64) then
                     m_max = m
                  endif
               enddo

               if (expected_weight(m_max) == 0_int64) then
                  error stop "reference initialization of aerosols not possible"
               endif

               expected_weight(n) = expected_weight(m_max) / 2_int64
               expected_weight(m_max) = expected_weight(m_max) - expected_weight(n)
               expected_aerosol_mass(n) = expected_aerosol_mass(m_max) / &
                                          real(expected_weight(m_max) + &
                                               expected_weight(n), real64) * &
                                          real(expected_weight(n), real64)
               expected_aerosol_mass(m_max) = expected_aerosol_mass(m_max) - &
                                              expected_aerosol_mass(n)
               expected_portions(:, n) = expected_portions(:, m_max)

               if (expected_weight(n) == 0_int64 .or. &
                   expected_weight(m_max) == 0_int64) then
                  error stop "reference initialization of aerosols not possible"
               endif
            endif
         enddo

         np_done = np_done + np_mode
      enddo
   end subroutine reference_bin_aerosol

   subroutine reference_usm_aerosol(k, expected_weight, expected_aerosol_mass, &
                                    expected_portions)
      integer, intent(in) :: k
      integer(int64), intent(out) :: expected_weight(np_xy)
      real(real64), intent(out) :: expected_aerosol_mass(np_xy)
      real(real64), intent(out) :: expected_portions(4, np_xy)

      integer :: i_mode
      integer :: n
      integer :: nn
      integer :: np_done
      integer :: np_mode
      real(real64) :: p_usm
      real(real64) :: r_aero
      real(real64) :: r_aero_max
      real(real64) :: r_aero_min
      real(real64) :: r_l
      real(real64) :: r_r
      real(real64) :: volume

      expected_weight = 0_int64
      expected_aerosol_mass = 0.0_real64
      expected_portions = 0.0_real64
      volume = 50.0_real64 * 50.0_real64 * dz(k)

      np_done = 0
      do i_mode = 1, 4
         call aerosol_mode_particle_count(i_mode, np_mode)
         if (np_mode == 0) cycle
         call assert_true(config%n_aero(i_mode) /= 0.0, &
                          "reference USM aerosol mode concentration")

         call reference_aerosol_radius_bounds(i_mode, r_aero_min, r_aero_max)
         p_usm = 1.0_real64 / real(np_mode, real64)
         r_l = r_aero_min

         do n = np_done + 1, np_done + np_mode - 1
            nn = n - np_done
            r_r = reference_percentile_radius(i_mode, p_usm * real(nn, real64), &
                                              r_l, r_aero_max)
            r_aero = sqrt(r_l * r_r)
            r_l = r_r

            if (trim(config%unit_aero) == "1/m3") then
               expected_weight(n) = int(volume * real(config%n_aero(i_mode), real64) / &
                                        real(np_mode, real64), int64)
            else if (trim(config%unit_aero) == "1/kg") then
               expected_weight(n) = int(volume * density(k) * &
                                        real(config%n_aero(i_mode), real64) / &
                                        real(np_mode, real64), int64)
            else
               error stop "unknown reference aerosol unit"
            endif

            expected_aerosol_mass(n) = real(expected_weight(n), real64) * &
                                       4.0_real64 / 3.0_real64 * pi * &
                                       r_aero**3 * &
                                       real(config%rho_aero(i_mode), real64)
            expected_portions(i_mode, n) = 1.0_real64
         enddo

         n = np_done + np_mode
         r_aero = sqrt(r_l * r_aero_max)
         expected_weight(n) = int(volume * sum(real(config%n_aero, real64)) / &
                                  real(np_xy, real64), int64)
         expected_aerosol_mass(n) = real(expected_weight(n), real64) * &
                                    4.0_real64 / 3.0_real64 * pi * &
                                    r_aero**3 * &
                                    real(config%rho_aero(i_mode), real64)
         np_done = np_done + np_mode
      enddo
   end subroutine reference_usm_aerosol

   subroutine reference_random_aerosol(k, expected_weight, expected_portions)
      integer, intent(in) :: k
      integer(int64), intent(out) :: expected_weight(np_xy)
      real(real64), intent(out) :: expected_portions(4, np_xy)

      integer :: i_mode
      integer :: n
      integer :: np_done
      integer :: np_mode
      real(real64) :: volume

      expected_weight = 0_int64
      expected_portions = 0.0_real64
      volume = 50.0_real64 * 50.0_real64 * dz(k)

      np_done = 0
      do i_mode = 1, 4
         call aerosol_mode_particle_count(i_mode, np_mode)
         if (np_mode == 0) cycle
         call assert_true(config%n_aero(i_mode) /= 0.0, &
                          "reference random aerosol mode concentration")

         do n = np_done + 1, np_done + np_mode
            if (trim(config%unit_aero) == "1/m3") then
               expected_weight(n) = int(volume * real(config%n_aero(i_mode), real64) / &
                                        real(np_mode, real64), int64)
            else if (trim(config%unit_aero) == "1/kg") then
               expected_weight(n) = int(volume * density(k) * &
                                        real(config%n_aero(i_mode), real64) / &
                                        real(np_mode, real64), int64)
            else
               error stop "unknown reference aerosol unit"
            endif
            expected_portions(i_mode, n) = 1.0_real64
         enddo
         np_done = np_done + np_mode
      enddo
   end subroutine reference_random_aerosol

   subroutine aerosol_mode_particle_count(i_mode, np_mode)
      integer, intent(in) :: i_mode
      integer, intent(out) :: np_mode

      integer :: n
      integer :: nn

      np_mode = floor(real(np_xy, real64) * &
                      real(config%shares_aero(i_mode), real64) / &
                      sum(real(config%shares_aero, real64)))
      if (i_mode == 1) then
         n = 0
         do nn = 1, 4
            n = n + floor(real(np_xy, real64) * &
                          real(config%shares_aero(nn), real64) / &
                          sum(real(config%shares_aero, real64)))
         enddo
         np_mode = np_mode + np_xy - n
      endif
   end subroutine aerosol_mode_particle_count

   subroutine reference_aerosol_radius_bounds(i_mode, r_aero_min, r_aero_max)
      integer, intent(in) :: i_mode
      real(real64), intent(out) :: r_aero_min
      real(real64), intent(out) :: r_aero_max

      r_aero_min = reference_percentile_radius(i_mode, 5.0e-5_real64, &
                                               1.0e-9_real64, 1.0e-5_real64)
      r_aero_max = reference_percentile_radius(i_mode, 1.0_real64 - 5.0e-5_real64, &
                                               r_aero_min, 1.0e-5_real64)
   end subroutine reference_aerosol_radius_bounds

   real(real64) function reference_percentile_radius(i_mode, percentile, &
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
      f_ak = reference_log_normal_cdf(real(config%n_aero(i_mode), real64), &
                                      real(config%sigma_aero(i_mode), real64), &
                                      real(config%rm_aero(i_mode), real64), ak) / &
             real(config%n_aero(i_mode), real64) - percentile

      bk = r_right
      f_bk = reference_log_normal_cdf(real(config%n_aero(i_mode), real64), &
                                      real(config%sigma_aero(i_mode), real64), &
                                      real(config%rm_aero(i_mode), real64), bk) / &
             real(config%n_aero(i_mode), real64) - percentile

      do
         ck = (ak * f_bk - bk * f_ak) / (f_bk - f_ak)
         f_ck = reference_log_normal_cdf(real(config%n_aero(i_mode), real64), &
                                         real(config%sigma_aero(i_mode), real64), &
                                         real(config%rm_aero(i_mode), real64), ck) / &
                real(config%n_aero(i_mode), real64) - percentile

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
   end function reference_percentile_radius

   real(real64) function reference_log_normal_pdf(n, s, m, x) result(pdf)
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
   end function reference_log_normal_pdf

   real(real64) function reference_log_normal_cdf(n, s, m, x) result(cdf)
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
   end function reference_log_normal_cdf

   subroutine assert_unique_ids()
      integer :: i, j, k, n
      integer :: idx
      integer :: other
      integer(int64) :: ids(expected_count)

      idx = 0
      do k = 1, nzl
         do j = 1, ny_local
            do i = 1, nx_local
               do n = 1, cell_particles(i,j,k)%count
                  idx = idx + 1
                  ids(idx) = cell_particles(i,j,k)%particles(n)%id
               enddo
            enddo
         enddo
      enddo

      do idx = 1, int(expected_count)
         do other = idx + 1, int(expected_count)
            call assert_true(ids(idx) /= ids(other), "unique particle IDs")
         enddo
      enddo
   end subroutine assert_unique_ids

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         error stop "failed assertion: " // label
      endif
   end subroutine assert_true

   subroutine assert_near(actual, expected, label)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: label

      if (abs(actual - expected) > 1.0e-12_real64) then
         error stop "failed assertion: " // label
      endif
   end subroutine assert_near

   subroutine assert_near_relative(actual, expected, label)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: label
      real(real64) :: scale

      scale = max(1.0e-30_real64, abs(expected))
      if (abs(actual - expected) / scale > 1.0e-10_real64) then
         error stop "failed assertion: " // label
      endif
   end subroutine assert_near_relative

end program test_initial_particle_population
