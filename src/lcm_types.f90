module lcm_types
   use iso_fortran_env, only: int64, real64

   implicit none

   private

   public :: lcm_grid_t
   public :: lcm_config_t
   public :: lcm_field_state_t
   public :: lcm_particle_t
   public :: particle_type_liquid_no_in
   public :: particle_type_liquid_with_in
   public :: particle_type_ice_no_in
   public :: particle_type_ice_with_in
   public :: lcm_particle_mpi_message_t
   public :: lcm_particle_cell_t
   public :: lcm_init_grid
   public :: lcm_init_config
   public :: timer_tic_ifc
   public :: timer_toc_ifc
   public :: lcm_timer_tic
   public :: lcm_timer_toc

   ! SAM-LCM iTYPE encoding: sign is phase, magnitude records IN membership.
   integer, parameter :: particle_type_liquid_no_in = 1
   integer, parameter :: particle_type_liquid_with_in = 2
   integer, parameter :: particle_type_ice_no_in = -1
   integer, parameter :: particle_type_ice_with_in = -2

   type :: lcm_grid_t
      ! Physical structured-grid geometry only. MPI rank topology is derived
      ! from the host communicator and stored in lcm_parallel_t.
      integer :: nx_local = 0
      integer :: ny_local = 0
      integer :: nz = 0

      integer :: nx_global = 0
      integer :: ny_global = 0

      !> Interior indices in the host model's local arrays.  For example,
      !> a DALES adapter can pass i=2:i1 and j=2:j1.  These indices describe
      !> the host array layout; they are not LCM particle-cell indices.
      integer :: i_start = 0
      integer :: i_end = 0
      integer :: j_start = 0
      integer :: j_end = 0

      real(real64) :: dx = 0.0_real64
      real(real64) :: dy = 0.0_real64

      real(real64), allocatable :: z_center(:)
      real(real64), allocatable :: z_face(:)
      real(real64), allocatable :: dz(:)

      integer :: halo_x = 0
      integer :: halo_y = 0
      integer :: halo_z = 0

      !> One-based global indices of LCM local cell (1,1).  Together with
      !> nx_local and ny_local, these define the physical MPI subdomain used
      !> for particle location and outgoing-particle classification.
      integer :: global_i_start = 0
      integer :: global_j_start = 0

      integer :: y_dimension = 1

   end type lcm_grid_t

   type :: lcm_config_t
      integer :: np_xy = 50
      integer :: nzl = 50
      integer :: nstep_lcm_init = 3
      integer :: nstep_cut_supersat = 0
      integer :: n_uvw_iteration = 1
      integer :: start_sedimentation = 1
      integer :: start_collection = 1
      integer :: breakup_conserve_mode = 1
      integer :: n_start_source_aero = -1
      integer :: n_end_source_aero = -1
      integer :: n_start_source_in = -1
      integer :: n_end_source_in = -1
      integer :: start_ice = 0
      integer :: spec_nsave = 30
      integer :: spec_nsavestart = 30
      integer :: n_track_particles = 1003
      integer :: particle_ntrack = 100
      integer :: particle_ntrackstart = 300
      integer :: particle_ntrackstop = 999999999
      integer :: nz_b_source_aero(4) = [1, 1, 1, 1]
      integer :: nz_t_source_aero(4) = [0, 0, 0, 0]
      integer :: ny_s_source_aero(4) = [1, 1, 1, 1]
      integer :: ny_n_source_aero(4) = [0, 0, 0, 0]
      integer :: nx_w_source_aero(4) = [1, 1, 1, 1]
      integer :: nx_e_source_aero(4) = [0, 0, 0, 0]
      integer :: nz_b_source_in(4) = [1, 1, 1, 1]
      integer :: nz_t_source_in(4) = [0, 0, 0, 0]
      integer :: ny_s_source_in(4) = [1, 1, 1, 1]
      integer :: ny_n_source_in(4) = [0, 0, 0, 0]
      integer :: nx_w_source_in(4) = [1, 1, 1, 1]
      integer :: nx_e_source_in(4) = [0, 0, 0, 0]
      integer :: type_source_in(4) = [3, 3, 3, 3]
      logical :: switch_sgs_velocities = .true.
      logical :: switch_sedimentation = .true.
      logical :: switch_collection = .true.
      logical :: switch_linear_sampling_collection = .true.
      logical :: switch_breakup = .true.
      logical :: switch_analytical_micro = .false.
      logical :: switch_radiation_micro = .false.
      logical :: switch_lem_micro = .true.
      logical :: switch_koehler_micro = .true.
      logical :: switch_ffc_micro = .false.
      logical :: switch_steady_pichamber = .false.
      logical :: switch_cloud_tracking = .false.
      logical :: switch_phase_sources = .true.
      logical :: switch_particle_source = .false.
      logical :: switch_walls_pichamber = .false.
      logical :: switch_reflective_walls_pichamber = .false.
      logical :: switch_ice = .false.
      logical :: switch_habit_ice = .true.
      logical :: switch_icebreak = .true.
      logical :: switch_rs = .false.
      logical :: switch_ds = .true.
      logical :: switch_moving_source = .false.
      logical :: switch_diskice_sizethresh = .true.
      logical :: switch_steady_aerosol = .false.
      logical :: switch_spec_aero = .false.
      logical :: track_some_particles = .false.
      logical :: horiz_average_spec = .true.
      logical :: start_from_restart = .false.
      real(real64) :: bottom_removal_radius = 1.0e-6_real64
      real :: r_separate_micro = 25.0e-6
      real :: supersat_threshold = 0.0
      real :: alpha_spec = 0.5
      real :: r_start_spec = 0.001e-6
      real :: r_end_spec = 6000.0e-6
      real :: frac_ice = 0.0
      real :: salinity_source = 0.037
      real :: sizethresh_chen = 10.0e-6
      real :: sizethresh_diskice = 5.0e-6
      real :: inject_bottom = 300.0
      real :: inject_top = 500.0
      real :: n_aero(4) = [100.0e6, 0.0, 0.0, 0.0]
      real :: rm_aero(4) = [0.04e-6, 0.0, 0.0, 0.0]
      real :: sigma_aero(4) = [1.4, 0.0, 0.0, 0.0]
      real :: shares_aero(4) = [1.0, 0.0, 0.0, 0.0]
      real :: n_ice(4) = [0.0, 0.0, 0.001e6, 0.0]
      real :: rm_in(4) = [1.0e-6, 0.2e-6, 1.0e-6, 0.0]
      real :: sigma_in(4) = [1.4, 1.4, 1.4, 0.0]
      real :: n_source_aero(4) = [0.0, 0.0, 0.0, 0.0]
      real :: n_source_in(4) = [0.0, 0.0, 0.0, 0.0]
      real :: r_source_aero(4) = [0.0625e-6, 0.0625e-6, 0.0625e-6, 0.0625e-6]
      real :: r_source_in(4) = [1.0e-6, 0.0, 0.0, 0.0]
      real :: sigma_source_aero(4) = [2.0, 2.0, 2.0, 2.0]
      real :: sigma_source_in(4) = [2.0, 2.0, 2.0, 2.0]
      real :: p_source_aero(4) = [0.0, 0.0, 0.0, 0.0]
      real :: p_source_in(4) = [0.0, 0.0, 0.0, 0.0]
      real :: frac_in_source_particles(4) = [0.0, 0.0, 0.0, 0.0]
      real :: frac_in_source_weight(4) = [0.001, 0.001, 0.001, 0.001]
      real :: vx_source_aero(4) = [0.0, 0.0, 0.0, 0.0]
      real :: vy_source_aero(4) = [0.0, 0.0, 0.0, 0.0]
      real :: vz_source_aero(4) = [0.0, 0.0, 0.0, 0.0]
      real :: vx_source_in(4) = [0.0, 0.0, 0.0, 0.0]
      real :: vy_source_in(4) = [0.0, 0.0, 0.0, 0.0]
      real :: vz_source_in(4) = [0.0, 0.0, 0.0, 0.0]
      real :: rho_aero(4) = [2170.0, 0.0, 0.0, 0.0]
      real :: molecular_weight_aero(4) = [58.443, 0.0, 0.0, 0.0]
      real :: vanthoff_aero(4) = [2.0, 0.0, 0.0, 0.0]
      character(len=6) :: init_aero = 'bin'
      character(len=12) :: unit_aero = '1/kg'
      character(len=12) :: unit_source = '1/kg/s'
      character(len=12) :: type_source = 'equilibrium'
      character(len=10) :: habitparam = 'Chen'
      character(len=10) :: ds_param = 'Phillips'
      character(len=8) :: icebreak_param = 'Phillips'
   end type lcm_config_t

   type :: lcm_field_state_t
      real(real64), pointer :: u(:,:,:) => null()
      real(real64), pointer :: v(:,:,:) => null()
      real(real64), pointer :: w(:,:,:) => null()
      real(real64), pointer :: static_energy(:,:,:) => null()
      real(real64), pointer :: temperature(:,:,:) => null()
      real(real64), pointer :: qv(:,:,:) => null()
      real(real64), pointer :: cloud_liquid(:,:,:) => null()
      real(real64), pointer :: pressure(:) => null()
      real(real64), pointer :: density(:) => null()

      logical :: u_is_attached = .false.
      logical :: v_is_attached = .false.
      logical :: w_is_attached = .false.
      logical :: static_energy_is_attached = .false.
      logical :: temperature_is_attached = .false.
      logical :: qv_is_attached = .false.
      logical :: cloud_liquid_is_attached = .false.
      logical :: pressure_is_attached = .false.
      logical :: density_is_attached = .false.
   end type lcm_field_state_t

   type :: lcm_particle_t
      ! SAM-LCM ivar(iA), ivar(iID), and ivar(iTYPE).
      integer(int64) :: weight = 0_int64
      integer(int64) :: id = 0_int64
      integer :: particle_type = particle_type_liquid_no_in

      ! SAM-LCM motion var(iX:iWsgs). Cell ownership remains the primary
      ! location, but these fields preserve the current velocity-advection path.
      real(real64) :: x = 0.0_real64
      real(real64) :: y = 0.0_real64
      real(real64) :: z = 0.0_real64
      real(real64) :: u = 0.0_real64
      real(real64) :: v = 0.0_real64
      real(real64) :: w = 0.0_real64
      real(real64) :: u_sgs = 0.0_real64
      real(real64) :: v_sgs = 0.0_real64
      real(real64) :: w_sgs = 0.0_real64

      ! SAM-LCM microphysical var(iM), var(iN), var(iWsedi), var(iT),
      ! var(iq), and var(ieta).
      real(real64) :: mass_water = 0.0_real64
      real(real64) :: mass_aero = 0.0_real64
      real(real64) :: w_sedi = 0.0_real64
      real(real64) :: temperature = 0.0_real64
      real(real64) :: qv = 0.0_real64
      real(real64) :: eta = 0.0_real64

      ! SAM-LCM ice and aerosol-composition state:
      ! var(iAXS), var(iCXS), var(iDNS), var(iMrime), var(iTFZ),
      ! var(iSC:iSP), and var(iAero_portions).
      real(real64) :: axis_a = 0.0_real64
      real(real64) :: axis_c = 0.0_real64
      real(real64) :: density = 0.0_real64
      real(real64) :: rimed_mass = 0.0_real64
      real(real64) :: freezing_temperature = 0.0_real64
      real(real64) :: in_surface_area(4) = 0.0_real64
      real(real64) :: aerosol_portions(4) = 0.0_real64

      ! SAM-LCM timing var(itime).
      real(real64) :: cloudy_age = 0.0_real64

      ! SAM-LCM integer metadata ivar(iLEM), ivar(iilst:iklst),
      ! ivar(iNmono), ivar(iSPR), and ivar(iTRACK).
      integer :: cell_particle_index = 0
      integer :: last_i = 0
      integer :: last_j = 0
      integer :: last_k = 0
      integer(int64) :: n_monomer = 0_int64
      integer :: spray_marker = -999
      integer :: track_marker = -999
   end type lcm_particle_t

   !> Temporary particle-exchange payload.
   !>
   !> global_i, global_j, and k identify the receiving cell. They are message
   !> metadata only: lcm_particle_t remains free of permanent ownership fields.
   type :: lcm_particle_mpi_message_t
      integer :: global_i = 0
      integer :: global_j = 0
      integer :: k = 0
      type(lcm_particle_t) :: particle
   end type lcm_particle_mpi_message_t

   type :: lcm_particle_cell_t
      integer :: count = 0
      integer :: capacity = 0
      type(lcm_particle_t), allocatable :: particles(:)
   end type lcm_particle_cell_t

   ! matches the host model (dales) timer routines, with its parameters
   abstract interface
      subroutine timer_tic_ifc(timer_name, opt_nvtx_id)
         character(*), intent(in)      :: timer_name
         integer, intent(in), optional :: opt_nvtx_id
      end subroutine

      subroutine timer_toc_ifc(timer_name)
         character(*), intent(in) :: timer_name
      end subroutine
   end interface

   procedure(timer_tic_ifc), pointer :: lcm_timer_tic => timer_tic_noop
   procedure(timer_toc_ifc), pointer :: lcm_timer_toc => timer_toc_noop

contains

   ! these two routines do nothing and only serve as a placeholder routine to point at
   subroutine timer_tic_noop(timer_name, opt_nvtx_id)
      character(*), intent(in)      :: timer_name
      integer, intent(in), optional :: opt_nvtx_id
   end subroutine

   subroutine timer_toc_noop(timer_name)
      character(*), intent(in) :: timer_name
   end subroutine

   subroutine lcm_init_config(config)
      type(lcm_config_t), intent(out) :: config

      config = lcm_config_t()
   end subroutine lcm_init_config

   subroutine lcm_init_grid(grid, nx_local, ny_local, nz, nx_global, ny_global, &
                            i_start, i_end, j_start, j_end, dx, dy, z_center, &
                            z_face, dz, halo_x, halo_y, halo_z,                &
                            global_i_start, global_j_start, y_dimension)
      type(lcm_grid_t), intent(out) :: grid
      integer, intent(in) :: nx_local
      integer, intent(in) :: ny_local
      integer, intent(in) :: nz
      integer, intent(in) :: nx_global
      integer, intent(in) :: ny_global
      integer, intent(in) :: i_start
      integer, intent(in) :: i_end
      integer, intent(in) :: j_start
      integer, intent(in) :: j_end
      real(real64), intent(in) :: dx
      real(real64), intent(in) :: dy
      real(real64), intent(in) :: z_center(:)
      real(real64), intent(in) :: z_face(:)
      real(real64), intent(in) :: dz(:)
      integer, intent(in) :: halo_x
      integer, intent(in) :: halo_y
      integer, intent(in) :: halo_z
      integer, intent(in) :: global_i_start
      integer, intent(in) :: global_j_start
      integer, intent(in), optional :: y_dimension

      grid%nx_local = nx_local
      grid%ny_local = ny_local
      grid%nz = nz

      grid%nx_global = nx_global
      grid%ny_global = ny_global

      grid%i_start = i_start
      grid%i_end = i_end
      grid%j_start = j_start
      grid%j_end = j_end

      grid%dx = dx
      grid%dy = dy

      allocate(grid%z_center(size(z_center)))
      allocate(grid%z_face(size(z_face)))
      allocate(grid%dz(size(dz)))

      grid%z_center = z_center
      grid%z_face = z_face
      grid%dz = dz

      grid%halo_x = halo_x
      grid%halo_y = halo_y
      grid%halo_z = halo_z

      grid%global_i_start = global_i_start
      grid%global_j_start = global_j_start

      if (present(y_dimension)) then
         grid%y_dimension = y_dimension
      endif
   end subroutine lcm_init_grid

end module lcm_types
