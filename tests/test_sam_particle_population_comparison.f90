program test_sam_particle_population_comparison
   use iso_fortran_env, only: int64, real64
   use lcm_host_interface, only: lcm_attach_fields, lcm_config_t, lcm_grid_t, &
                                 lcm_init_config, lcm_init_grid,             &
                                 lcm_initialize, lcm_set_config, lcm_set_grid
   use lcm_particle_state, only: cell_particles

   implicit none

   integer, parameter :: np_xy = 8
   integer, parameter :: nx_local = 1
   integer, parameter :: ny_local = 1
   integer, parameter :: nz = 1
   integer, parameter :: nzl = 1

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
   integer(int64) :: sam_weight(np_xy)
   real(real64) :: sam_mass_aero(np_xy)
   real(real64) :: sam_aerosol_portions(4, np_xy)
   real(real64) :: sam_in_surface_area(4, np_xy)
   real(real64) :: sam_freezing_temperature(np_xy)
   integer :: sam_particle_type(np_xy)
   character(len=1024) :: helper_path
   character(len=1024) :: output_path
   character(len=32) :: aerosol_mode
   character(len=32) :: ice_mode
   character(len=2048) :: command
   integer :: command_status
   integer :: exit_status

   call get_command_argument(1, helper_path)
   if (len_trim(helper_path) == 0) then
      error stop "SAM-LCM particle initialization probe path is required"
   endif
   call get_command_argument(2, aerosol_mode)
   if (len_trim(aerosol_mode) == 0) aerosol_mode = "bin"
   call get_command_argument(3, ice_mode)

   output_path = "sam_" // trim(aerosol_mode) // &
                 trim(ice_mode) // "_particle_population_reference.dat"
   command = trim(helper_path) // " " // trim(output_path) // " " // &
             trim(aerosol_mode) // " " // trim(ice_mode)
   call execute_command_line(trim(command), exitstat=exit_status, &
                             cmdstat=command_status)
   if (command_status /= 0 .or. exit_status /= 0) then
      error stop "SAM-LCM particle initialization probe failed"
   endif

   call read_sam_reference(trim(output_path), sam_weight, sam_mass_aero, &
                           sam_aerosol_portions, sam_particle_type,     &
                           sam_in_surface_area, sam_freezing_temperature)
   call initialize_standalone_lcm_case(trim(aerosol_mode), trim(ice_mode))

   call assert_true(cell_particles(1,1,1)%count == np_xy, &
                    "standalone particle count")
   call assert_matches_sam_reference(sam_weight, sam_mass_aero, &
                                     sam_aerosol_portions, sam_particle_type, &
                                     sam_in_surface_area,                    &
                                     sam_freezing_temperature)

contains

   subroutine initialize_standalone_lcm_case(aerosol_mode, ice_mode)
      character(len=*), intent(in) :: aerosol_mode
      character(len=*), intent(in) :: ice_mode

      call lcm_init_config(config)
      config%np_xy = np_xy
      config%nzl = nzl
      config%init_aero = aerosol_mode
      config%unit_aero = "1/kg"
      config%switch_ice = trim(ice_mode) == "ice"
      if (config%switch_ice) then
         config%frac_ice = 0.999999
      else
         config%frac_ice = 0.0
      endif
      config%n_aero = [100.0e6, 0.0, 0.0, 0.0]
      config%rm_aero = [0.04e-6, 0.0, 0.0, 0.0]
      config%sigma_aero = [1.4, 0.0, 0.0, 0.0]
      config%shares_aero = [1.0, 0.0, 0.0, 0.0]
      config%rho_aero = [2170.0, 0.0, 0.0, 0.0]
      config%n_ice = [0.0, 0.0, 0.001e6, 0.0]
      config%rm_in = [1.0e-6, 0.2e-6, 1.0e-6, 0.0]
      config%sigma_in = [1.4, 1.4, 1.4, 0.0]
      call lcm_set_config(config)

      dz = [100.0_real64]
      z_center = [50.0_real64]
      z_face = [0.0_real64, 100.0_real64]

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
      pressure = [100000.0_real64]
      density = [1.2_real64]

      call lcm_attach_fields(u=u, v=v, w=w, static_energy=static_energy, &
                             temperature=temperature, qv=qv,             &
                             pressure=pressure, density=density)
      call lcm_initialize()
   end subroutine initialize_standalone_lcm_case

   subroutine read_sam_reference(path, weight, mass_aero, aerosol_portions, &
                                 particle_type, in_surface_area,            &
                                 freezing_temperature)
      character(len=*), intent(in) :: path
      integer(int64), intent(out) :: weight(np_xy)
      real(real64), intent(out) :: mass_aero(np_xy)
      real(real64), intent(out) :: aerosol_portions(4, np_xy)
      integer, intent(out) :: particle_type(np_xy)
      real(real64), intent(out) :: in_surface_area(4, np_xy)
      real(real64), intent(out) :: freezing_temperature(np_xy)

      integer :: i
      integer :: record_index
      integer :: input_unit

      open(newunit=input_unit, file=path, status="old", action="read")
      do i = 1, np_xy
         read(input_unit, *) record_index, weight(i), mass_aero(i), &
                            aerosol_portions(1,i), aerosol_portions(2,i), &
                            aerosol_portions(3,i), aerosol_portions(4,i), &
                            particle_type(i), in_surface_area(1,i),       &
                            in_surface_area(2,i), in_surface_area(3,i),   &
                            in_surface_area(4,i), freezing_temperature(i)
         call assert_true(record_index == i, "SAM reference record index")
      enddo
      close(input_unit)
   end subroutine read_sam_reference

   subroutine assert_matches_sam_reference(weight, mass_aero, aerosol_portions, &
                                           particle_type, in_surface_area,   &
                                           freezing_temperature)
      integer(int64), intent(in) :: weight(np_xy)
      real(real64), intent(in) :: mass_aero(np_xy)
      real(real64), intent(in) :: aerosol_portions(4, np_xy)
      integer, intent(in) :: particle_type(np_xy)
      real(real64), intent(in) :: in_surface_area(4, np_xy)
      real(real64), intent(in) :: freezing_temperature(np_xy)

      integer :: n

      do n = 1, np_xy
         call assert_true(cell_particles(1,1,1)%particles(n)%weight == &
                          weight(n), "SAM aerosol weight")
         call assert_near_relative(cell_particles(1,1,1)%particles(n)%mass_aero, &
                                   mass_aero(n), 1.0e-10_real64, &
                                   "SAM aerosol mass")
         call assert_near_relative( &
            cell_particles(1,1,1)%particles(n)%aerosol_portions(1), &
            aerosol_portions(1,n), 1.0e-10_real64, &
            "SAM aerosol soluble aerosol portion")
         call assert_near_relative( &
            cell_particles(1,1,1)%particles(n)%aerosol_portions(2), &
            aerosol_portions(2,n), 1.0e-10_real64, &
            "SAM aerosol insoluble aerosol portion")
         call assert_near_relative( &
            cell_particles(1,1,1)%particles(n)%aerosol_portions(3), &
            aerosol_portions(3,n), 1.0e-10_real64, &
            "SAM aerosol mixed aerosol portion")
         call assert_near_relative( &
            cell_particles(1,1,1)%particles(n)%aerosol_portions(4), &
            aerosol_portions(4,n), 1.0e-10_real64, &
            "SAM aerosol immersion-freezing portion")
         call assert_true(cell_particles(1,1,1)%particles(n)%particle_type == &
                          particle_type(n), "SAM aerosol particle type")
         if (maxval(abs(in_surface_area(:,n))) > 0.0_real64) then
            call assert_true(cell_particles(1,1,1)%particles(n)%in_surface_area(3) > &
                             0.0_real64, "SAM active IN surface area")
            call assert_near_relative( &
               maxval(abs(cell_particles(1,1,1)%particles(n)%in_surface_area([1,2,4]))), &
               0.0_real64, 1.0e-10_real64, "SAM inactive IN surface areas")
            call assert_true(cell_particles(1,1,1)%particles(n)%freezing_temperature >= &
                             233.15_real64, "SAM freezing temperature lower bound")
            call assert_true(cell_particles(1,1,1)%particles(n)%freezing_temperature <= &
                             273.15_real64, "SAM freezing temperature upper bound")
         else
            call assert_near_relative( &
               maxval(abs(cell_particles(1,1,1)%particles(n)%in_surface_area)), &
               0.0_real64, 1.0e-10_real64, "SAM IN surface area")
            call assert_near_relative( &
               cell_particles(1,1,1)%particles(n)%freezing_temperature, &
               freezing_temperature(n), 1.0e-10_real64, &
               "SAM freezing temperature")
         endif
      enddo
   end subroutine assert_matches_sam_reference

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_true

   subroutine assert_near_relative(actual, expected, tolerance, message)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      real(real64), intent(in) :: tolerance
      character(len=*), intent(in) :: message
      real(real64) :: scale

      scale = max(1.0e-30_real64, abs(expected))
      if (abs(actual - expected) > tolerance * scale) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_near_relative

end program test_sam_particle_population_comparison
