program test_particle_source
   use iso_fortran_env, only: int64, real64
   use lcm_host_interface, only: lcm_advance, lcm_attach_fields, lcm_config_t, &
                                 lcm_grid_t, lcm_init_config, lcm_init_grid,  &
                                 lcm_initialize, lcm_set_config, lcm_set_grid
   use lcm_particle_state, only: cell_particles

   implicit none

   integer, parameter :: nx_local = 2
   integer, parameter :: ny_local = 2
   integer, parameter :: nz = 3
   integer, parameter :: nzl = 2
   integer, parameter :: np_xy = 8

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
   if (len_trim(case_name) == 0) case_name = "aerosol"

   call configure_model(trim(case_name))
   call lcm_initialize()
   call lcm_advance(1.0_real64)

   select case (trim(case_name))
   case ("aerosol")
      call assert_aerosol_source()
   case ("in")
      call assert_in_source()
   case default
      error stop "unknown particle source test case"
   end select

contains

   subroutine configure_model(case_name)
      character(len=*), intent(in) :: case_name

      call lcm_init_config(config)
      call configure_base_config(config)

      select case (case_name)
      case ("aerosol")
         call configure_aerosol_source(config)
      case ("in")
         call configure_in_source(config)
      case default
         error stop "unknown particle source test case"
      end select

      call lcm_set_config(config)
      call configure_grid_and_fields()
   end subroutine configure_model

   subroutine configure_base_config(config)
      type(lcm_config_t), intent(inout) :: config

      config%np_xy = np_xy
      config%nzl = nzl
      config%init_aero = "bin"
      config%unit_aero = "1/kg"
      config%n_aero = [100.0e6, 0.0, 0.0, 0.0]
      config%rm_aero = [0.04e-6, 0.0, 0.0, 0.0]
      config%sigma_aero = [1.4, 0.0, 0.0, 0.0]
      config%shares_aero = [1.0, 0.0, 0.0, 0.0]
      config%rho_aero = [2170.0, 0.0, 0.0, 0.0]
      config%switch_particle_source = .true.
      config%unit_source = "1/m3/s"
      config%r_source_aero = [0.0625e-6, 0.0, 0.0, 0.0]
      config%sigma_source_aero = [2.0, 0.0, 0.0, 0.0]
   end subroutine configure_base_config

   subroutine configure_aerosol_source(config)
      type(lcm_config_t), intent(inout) :: config

      config%switch_ice = .false.
      config%switch_phase_sources = .false.
      config%n_start_source_aero = config%nstep_lcm_init
      config%n_end_source_aero = -1
      config%type_source = "equilibrium"
      config%n_source_aero = [100.0, 0.0, 0.0, 0.0]
      config%p_source_aero = [1.0, 0.0, 0.0, 0.0]
      config%nx_w_source_aero = [1, 1, 1, 1]
      config%nx_e_source_aero = [1, 0, 0, 0]
      config%ny_s_source_aero = [1, 1, 1, 1]
      config%ny_n_source_aero = [1, 0, 0, 0]
      config%nz_b_source_aero = [1, 1, 1, 1]
      config%nz_t_source_aero = [1, 0, 0, 0]
   end subroutine configure_aerosol_source

   subroutine configure_in_source(config)
      type(lcm_config_t), intent(inout) :: config

      config%switch_ice = .true.
      config%frac_ice = 0.0
      config%switch_phase_sources = .true.
      config%n_start_source_in = config%nstep_lcm_init
      config%n_end_source_in = -1
      config%n_source_in = [100.0, 0.0, 0.0, 0.0]
      config%p_source_in = [1.0, 0.0, 0.0, 0.0]
      config%r_source_in = [1.0e-6, 0.0, 0.0, 0.0]
      config%sigma_source_in = [2.0, 0.0, 0.0, 0.0]
      config%type_source_in = [3, 3, 3, 3]
      config%nx_w_source_in = [1, 1, 1, 1]
      config%nx_e_source_in = [1, 0, 0, 0]
      config%ny_s_source_in = [1, 1, 1, 1]
      config%ny_n_source_in = [1, 0, 0, 0]
      config%nz_b_source_in = [1, 1, 1, 1]
      config%nz_t_source_in = [1, 0, 0, 0]
   end subroutine configure_in_source

   subroutine configure_grid_and_fields()
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
   end subroutine configure_grid_and_fields

   subroutine assert_aerosol_source()
      integer :: expected_source_weight
      integer :: source_index

      expected_source_weight = nint(100.0_real64 * 50.0_real64 * &
                                    50.0_real64 * dz(1))
      source_index = np_xy + 1

      call assert_source_cell_counts()
      call assert_true(cell_particles(1,1,1)%particles(source_index)%spray_marker == 100, &
                       "source particle spray marker")
      call assert_true(cell_particles(1,1,1)%particles(source_index)%weight == &
                       int(expected_source_weight, int64), &
                       "SAM-LCM source particle weight")
      call assert_true(cell_particles(1,1,1)%particles(source_index)%particle_type == 1, &
                       "source particle liquid type")
      call assert_true(cell_particles(1,1,1)%particles(source_index)%mass_aero > 0.0_real64, &
                       "source particle aerosol mass")
      call assert_true(cell_particles(1,1,1)%particles(source_index)%mass_water > 0.0_real64, &
                       "source particle liquid mass")
      call assert_near(cell_particles(1,1,1)%particles(source_index)%temperature, &
                       static_energy(1,1,1), "source particle absolute temperature")
      call assert_true(cell_particles(1,1,1)%particles(source_index)%aerosol_portions(1) == &
                       1.0_real64, "source particle aerosol portion")
      call assert_unique_ids()
   end subroutine assert_aerosol_source

   subroutine assert_in_source()
      integer :: expected_source_weight
      integer :: source_index

      expected_source_weight = nint(100.0_real64 * 50.0_real64 * &
                                    50.0_real64 * dz(1))
      source_index = np_xy + 1

      call assert_source_cell_counts()
      call assert_true(cell_particles(1,1,1)%particles(source_index)%spray_marker == 100, &
                       "IN source particle spray marker")
      call assert_true(cell_particles(1,1,1)%particles(source_index)%weight == &
                       int(expected_source_weight, int64), &
                       "SAM-LCM IN source particle weight")
      call assert_true(cell_particles(1,1,1)%particles(source_index)%particle_type == 2, &
                       "IN source particle type")
      call assert_true(cell_particles(1,1,1)%particles(source_index)%in_surface_area(3) > &
                       0.0_real64, "IN source particle surface area")
      call assert_true(cell_particles(1,1,1)%particles(source_index)%freezing_temperature > &
                       0.0_real64, "IN source particle freezing temperature")
      call assert_near(cell_particles(1,1,1)%particles(source_index)%temperature, &
                       static_energy(1,1,1), "IN source particle absolute temperature")
      call assert_unique_ids()
   end subroutine assert_in_source

   subroutine assert_source_cell_counts()
      call assert_true(cell_particles(1,1,1)%count == np_xy + 1, &
                       "source cell particle count")
      call assert_true(cell_particles(2,1,1)%count == np_xy, &
                       "non-source x-neighbor particle count")
      call assert_true(cell_particles(1,2,1)%count == np_xy, &
                       "non-source y-neighbor particle count")
      call assert_true(cell_particles(1,1,2)%count == np_xy, &
                       "non-source z-neighbor particle count")
   end subroutine assert_source_cell_counts

   subroutine assert_unique_ids()
      integer :: i, j, k, n
      integer :: idx
      integer :: other
      integer :: total_count
      integer(int64), allocatable :: ids(:)

      total_count = 0
      do k = 1, nzl
         do j = 1, ny_local
            do i = 1, nx_local
               total_count = total_count + cell_particles(i,j,k)%count
            enddo
         enddo
      enddo

      allocate(ids(total_count))
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

      do idx = 1, total_count
         do other = idx + 1, total_count
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
end program test_particle_source
