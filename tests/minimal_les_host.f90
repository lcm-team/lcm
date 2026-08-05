!> Minimal host-model fixture for exercising the public LCM host interface.
!!
!! The arrays in this module represent host-owned LES fields. They are module
!! variables with the TARGET attribute because lcm_attach_fields stores pointers
!! to them; tests must keep the host fields alive for the full LCM run.
module minimal_les_host
   use iso_fortran_env, only: real64
   use lcm_host_interface, only: lcm_attach_fields, lcm_config_t, lcm_grid_t, &
                                 lcm_init_config, lcm_init_grid,             &
                                 lcm_initialize, lcm_set_config, lcm_set_grid

   implicit none

   private

   public :: minimal_host_config
   public :: minimal_host_grid
   public :: minimal_host_nx_local
   public :: minimal_host_ny_local
   public :: minimal_host_nz
   public :: minimal_host_nzl
   public :: minimal_host_np_xy
   public :: minimal_host_dx
   public :: minimal_host_dy
   public :: minimal_host_u
   public :: minimal_host_v
   public :: minimal_host_w
   public :: minimal_host_static_energy
   public :: minimal_host_temperature
   public :: minimal_host_qv
   public :: minimal_host_pressure
   public :: minimal_host_density
   public :: minimal_host_dz
   public :: minimal_host_z_center
   public :: minimal_host_z_face
   public :: setup_minimal_les_host

   integer, parameter :: minimal_host_nx_local = 2
   integer, parameter :: minimal_host_ny_local = 2
   integer, parameter :: minimal_host_nz = 2
   integer, parameter :: minimal_host_nzl = 1
   integer, parameter :: minimal_host_np_xy = 4
   real(real64), parameter :: minimal_host_dx = 50.0_real64
   real(real64), parameter :: minimal_host_dy = 50.0_real64

   type(lcm_config_t), save :: minimal_host_config
   type(lcm_grid_t), save :: minimal_host_grid

   real(real64), target, save :: minimal_host_u( &
      minimal_host_nx_local, minimal_host_ny_local, minimal_host_nz)
   real(real64), target, save :: minimal_host_v( &
      minimal_host_nx_local, minimal_host_ny_local, minimal_host_nz)
   real(real64), target, save :: minimal_host_w( &
      minimal_host_nx_local, minimal_host_ny_local, minimal_host_nz)
   real(real64), target, save :: minimal_host_static_energy( &
      minimal_host_nx_local, minimal_host_ny_local, minimal_host_nz)
   real(real64), target, save :: minimal_host_temperature( &
      minimal_host_nx_local, minimal_host_ny_local, minimal_host_nz)
   real(real64), target, save :: minimal_host_qv( &
      minimal_host_nx_local, minimal_host_ny_local, minimal_host_nz)
   real(real64), target, save :: minimal_host_pressure(minimal_host_nz)
   real(real64), target, save :: minimal_host_density(minimal_host_nz)

   real(real64), save :: minimal_host_dz(minimal_host_nz)
   real(real64), save :: minimal_host_z_center(minimal_host_nz)
   real(real64), save :: minimal_host_z_face(minimal_host_nz + 1)

contains

   subroutine setup_minimal_les_host(case_name)
      character(len=*), intent(in), optional :: case_name

      call configure_minimal_host_config(case_name)
      call configure_minimal_host_grid()
      call configure_minimal_host_fields()

      call lcm_attach_fields(u=minimal_host_u, v=minimal_host_v, &
                             w=minimal_host_w,                   &
                             static_energy=minimal_host_static_energy, &
                             temperature=minimal_host_temperature, &
                             qv=minimal_host_qv,                  &
                             pressure=minimal_host_pressure,      &
                             density=minimal_host_density)

      call lcm_initialize()
   end subroutine setup_minimal_les_host

   subroutine configure_minimal_host_config(case_name)
      character(len=*), intent(in), optional :: case_name
      character(len=32) :: selected_case

      call lcm_init_config(minimal_host_config)

      minimal_host_config%np_xy = minimal_host_np_xy
      minimal_host_config%nzl = minimal_host_nzl
      minimal_host_config%init_aero = "bin"
      minimal_host_config%unit_aero = "1/kg"
      minimal_host_config%switch_ice = .false.
      minimal_host_config%n_aero = [100.0e6, 0.0, 0.0, 0.0]
      minimal_host_config%rm_aero = [0.04e-6, 0.0, 0.0, 0.0]
      minimal_host_config%sigma_aero = [1.4, 0.0, 0.0, 0.0]
      minimal_host_config%shares_aero = [1.0, 0.0, 0.0, 0.0]
      minimal_host_config%rho_aero = [2170.0, 0.0, 0.0, 0.0]

      selected_case = "bin"
      if (present(case_name)) selected_case = case_name
      select case (trim(selected_case))
      case ("bin")
         minimal_host_config%init_aero = "bin"
      case ("USM")
         minimal_host_config%init_aero = "USM"
      case ("random")
         minimal_host_config%init_aero = "random"
      case ("bin_ice")
         minimal_host_config%init_aero = "bin"
         minimal_host_config%switch_ice = .true.
         minimal_host_config%frac_ice = 0.999999
         minimal_host_config%n_ice = [0.0, 0.0, 0.001e6, 0.0]
         minimal_host_config%rm_in = [1.0e-6, 0.2e-6, 1.0e-6, 0.0]
         minimal_host_config%sigma_in = [1.4, 1.4, 1.4, 0.0]
      case default
         error stop "unknown minimal LES host initialization case"
      end select

      call lcm_set_config(minimal_host_config)
   end subroutine configure_minimal_host_config

   subroutine configure_minimal_host_grid()
      minimal_host_dz = [100.0_real64, 120.0_real64]
      minimal_host_z_face = [0.0_real64, 100.0_real64, 220.0_real64]
      minimal_host_z_center = [50.0_real64, 160.0_real64]

      call lcm_init_grid(minimal_host_grid,                         &
                         minimal_host_nx_local, minimal_host_ny_local, &
                         minimal_host_nz,                           &
                         minimal_host_nx_local, minimal_host_ny_local, &
                         1, minimal_host_nx_local,                  &
                         1, minimal_host_ny_local,                  &
                         minimal_host_dx, minimal_host_dy,          &
                         minimal_host_z_center, minimal_host_z_face,&
                         minimal_host_dz,                           &
                         0, 0, 0, 0, 0, 1, 1, 1, 1,                &
                         -1, -1, -1, -1, .true., .true., 1)
      call lcm_set_grid(minimal_host_grid)
   end subroutine configure_minimal_host_grid

   subroutine configure_minimal_host_fields()
      minimal_host_u = 0.10_real64
      minimal_host_v = -0.05_real64
      minimal_host_w = 0.02_real64

      minimal_host_temperature(:,:,1) = 280.0_real64
      minimal_host_temperature(:,:,2) = 278.0_real64
      minimal_host_static_energy(:,:,1) = 300.0_real64
      minimal_host_static_energy(:,:,2) = 299.0_real64
      minimal_host_qv(:,:,1) = 0.010_real64
      minimal_host_qv(:,:,2) = 0.008_real64

      minimal_host_pressure = [100000.0_real64, 98000.0_real64]
      minimal_host_density = [1.20_real64, 1.05_real64]
   end subroutine configure_minimal_host_fields

end module minimal_les_host
