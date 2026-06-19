 MODULE micro_vars
!
!-- Microphysical variables
    USE micro_user_vars, ONLY: nvU, nivU

    USE params,          ONLY: pi
    USE domain,          ONLY: YES3D
#ifdef SAMCHEM
    USE domain,          ONLY: naqspc
#endif /*SAMCHEM*/
    USE grid,            ONLY: nx, ny, nz, nzm, nsubdomains, rank

    IMPLICIT NONE
    PUBLIC
    SAVE

!
!-- Variables and arrays to store SAM variables. Store SAM values so that:
!--   - no modification will be made in SAM.
!--   - conflict with different real kind between SAM and LCM is avoided.
    INTEGER, PARAMETER                  :: nx_LCM = nx, ny_LCM = ny, &
                                           nz_LCM = nzm
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: u_LCM, v_LCM, w_LCM, tabs_LCM,      &
                                           tke_LCM, diss_LCM, tk_LCM, smix_LCM
    REAL, DIMENSION(nz_LCM)             :: pres_LCM, rho_air, dz_LCM
    REAL                                :: dtn_LCM, dx_LCM, dy_LCM
!
!-- Precision of the LCM real variables: either 4 or 8 byte
    INTEGER, PARAMETER :: lcmp = 8
!
!-- Computational particles
    INTEGER :: np_xy = 50         ! Number of particles placed in the same grid box
    INTEGER :: nzl = 50           ! Highest grid level for placing particles
    INTEGER :: nstep_LCM_init = 3 ! nstep_LCM_init is the time step at which the LCM is initialized
    INTEGER :: nstep_cut_supersat = 0 ! number of steps (at beginning) to cut supersaturation to supersat_threshold (only for liquid now since I assume ice is always started later on)
    REAL    :: supersat_threshold = 0.0 ! threshold to cut supersaturation until nstep_cut_supersat
!
!-- Particle motion
    INTEGER :: n_uvw_iteration = 1      ! Number of iteration to find velocity at t=t+1/2. If set 0, the particle position is updated with the forward scheme: x(t+1)=x(t)+u(t)*dt
    LOGICAL :: switch_sgs_velocities = .TRUE. ! Switch to indicate use of sgs velocities
!
!-- Pi champer controls
    LOGICAL :: switch_steady_pichamber = .FALSE.           ! Pi chamber boundary condition (re-injection of absorbed droplets as aerosols)
    LOGICAL :: switch_reflective_walls_pichamber = .FALSE. ! Pi chamber reflective (.TRUE.) or absorbing (.FALSE.) walls
    LOGICAL :: switch_walls_pichamber = .FALSE.            ! Pi chamber walls
    REAL    :: wall_x_w_pichamber = 0.001, wall_x_e_pichamber = 1.999
    REAL    :: wall_y_s_pichamber = 0.001, wall_y_n_pichamber = 1.999
!
!-- Cloud tracking
    LOGICAL                 :: switch_cloud_tracking = .FALSE.
    REAL                    :: tau_SCLR = 1800.0  ! decay timescale [60.0 s in Dawe and Austin (2012), 1800.0 s in Heus and Seifert (2013)]
    REAL, DIMENSION(nz_LCM) :: mean_SCLR_prof     ! OUTPUT profile of mean of plume-detecting scalar
    REAL, DIMENSION(nz_LCM) :: std_SCLR_prof      ! OUTPUT profile of standard deviation of plume-detecting scalar
!
!-- Sedimentation
    LOGICAL :: switch_sedimentation = .TRUE. ! Switch for sedimentation of cloud droplets
    INTEGER :: start_sedimentation = 1! # of timesteps to initiate sedimentation after model starts
    LOGICAL :: switch_steady_aerosol = .FALSE. !Whether to re-inject aerosol from particles hitting the surface
    REAL    :: inject_bottom         = 300.0 !m, minimum height at which to re-inject aerosol if switch_steady_aerosol is .TRUE.
    REAL    :: inject_top            = 500.0 !m, maximum height at which to re-inject aerosol if switch_steady_aerosol is .TRUE.
!
!-- Collisional growth
    LOGICAL          :: switch_collection = .TRUE. ! Switch for collection growth
    INTEGER          :: start_collection = 1  ! # of timesteps to initiate collection after model starts
    LOGICAL          :: switch_linear_sampling_collection = .TRUE. ! LSM method to speed up collision process (Shima et al. 2009)
!
!-- Breakup
    LOGICAL          :: switch_breakup = .TRUE. !Whether to do collisional breakup in the simulation
    INTEGER          :: breakup_conserve_mode = 1 ! Set to 1 to attempt to maintain existing superdroplet diameters and to reassign fragments
                                                  ! to the closest-size superdroplet within the gridbox. Set to 2 to follow the approach of
                                                  ! De Jong et al. (2023) and set the smaller weighting-factor superdroplet to be all fragments
!
!-- Diffusional growth
    LOGICAL :: switch_analytical_micro = .FALSE. ! Switch for using analytical condensational growth
    LOGICAL :: switch_radiation_micro  = .FALSE. ! Switch for radiation effects of condensation growth
    LOGICAL :: switch_LEM_micro   = .TRUE.       ! Switch for supersaturation sgs model
    LOGICAL :: switch_koehler_micro = .TRUE.     ! Switch for Köhler theory
    LOGICAL :: switch_FFC_micro = .FALSE.        ! Switch for film-forming compounds effect on diffusional growth (Feingold and Chuang 2002)

    REAL, DIMENSION(4) :: n_aero = (/ 100.0E6, 0.0, 0.0, 0.0 /)              ! aerosol/droplet number mixing ratio (#/m3)
    REAL, DIMENSION(4) :: rm_aero = (/ 0.04E-6, 0.0, 0.0, 0.0 /) ! Geometric mean radius of aerosol size distribution (m)
    REAL, DIMENSION(4) :: sigma_aero = (/ 1.4, 0.0, 0.0, 0.0 /)              ! Geometric stddev of aerosol size distribution (unitless)
    REAL, DIMENSION(4) :: shares_aero = (/ 1, 1, 1, 1 /)                  ! shares of superparticles used to create each mode (like cocktail-mixing, "how many parts does each mode get?")
    LOGICAL :: switch_inmass          = .FALSE.             ! Use the ice nucleation scheme from the particle model, which keeps track of IN mass

    REAL, DIMENSION(4) :: rho_aero = (/ 2170.0, 0.0, 0.0, 0.0 /) ! mass density aerosol (kg/m3)
    REAL, DIMENSION(4) :: molecular_weight_aero = (/ 58.443, 0.0, 0.0, 0.0 /) ! molecular weight aerosol (g/mol)
    REAL, DIMENSION(4) :: vanthoff_aero = (/ 2.0, 0.0, 0.0, 0.0  /) ! van‘t Hoff factor aerosol
    REAL :: sizethresh_chen = 10.0E-6       ! size threshold below which to force isometric ice growth for the Chen & Lamb (1994) parameterization
    REAL:: sizethresh_diskice = 5.0E-6     ! size threshold below which to force isometric ice growth for the DiSKICE (Zhang and Harrington, 2014) parameterization
    LOGICAl :: switch_diskice_sizethresh = .TRUE.  ! force the diskice scheme to use a minimum size threshold for habit development. Default is true

    CHARACTER (LEN=6)  :: init_aero = 'bin'    ! 'bin' or 'random' aerosol initialization
    CHARACTER (LEN=12) :: unit_aero = '1/kg'   ! '1/kg' or '1/m3' aerosol initialization units
    CHARACTER (LEN=10) :: habitparam = "Chen" !Options are currently 'Chen' 'DISKICE'. Will add 'Weiss' soon: B. Ascher 15.01.2025

    INTEGER, PARAMETER   :: n_inv_cdf = 50000000 ! number of nodes to discretize inverse aerosol cdf
    REAL,    ALLOCATABLE :: inv_cdf(:, :)

    REAL, PARAMETER   :: r_act = 1.0E-6  ! radius above which particles are considered activated (in meters)
!
!-- Parameters and constants not available from SAM:
    REAL, PARAMETER :: rho_liq = 1000.0              ! Density of water (kg/m3)
    REAL, PARAMETER :: sigma = 0.0728                ! Suface tension of water (J/m2) (Straka: 0.0761 - 0.000155 * ( tabs(i,j,k) - 273.15 ) )
    REAL, PARAMETER :: molecular_weight_water = 18.0 ! Molecular weight of water (g/mol)
    REAL, PARAMETER :: N_Avogadro = 6.02214076E23    ! Avogadro constant (1/mol)
    REAL, PARAMETER :: g_asym_param_LCM = 0.85       ! asymmetry parameter (Lacis and Hansen, 1974, JAS)

    REAL, PARAMETER :: nv = 3.35*10.0**(28) !Number concentration of water molecules, #/m^3
    REAL, PARAMETER :: boltzk = 1.38*10.0**(-23) !Boltzmann's constant, J/K
    REAL, PARAMETER :: planckh = 6.626*10.0**(-34) !Planck's constant, J*s
    REAL, PARAMETER :: av = -2.323*10.0**(-18) !Fit parameter 1, J - from Kuhn (2011)
    REAL, PARAMETER :: bv = -1.075*10.0**(-20) !Fit parameter 2, J/K - from Kuhn (2011)
    REAL, PARAMETER :: dw = -2.13*10.0**(-20) !Perturbation in free energy of freezing, J - from Kuhn (2011)
    REAL, PARAMETER :: srfthick = 2.5*10.0**(-9) !Thickness of "droplet surface layer", m - from Kuhn (2011)
!
!-- Ice parameters
    REAL, PARAMETER :: A_ice            = 9.99E-10   ! Coefficient nucleation scheme (1/s), Hartmann et al. 2013 (ACP)
    REAL, PARAMETER :: B_ice            = -2.34      ! Coefficient nucleation scheme (1/°C), Hartmann et al. 2013 (ACP)
    REAL, PARAMETER :: k_ice            = 1.0        ! Number of ice nucleating protein complexes
    REAL, PARAMETER :: rho_ice          = 916.7      ! Mass desity ice (kg/m3)
!
!-- B. Ascher 13.05.2025 - Add different modes of IN. They have different mean sizes and different temperatures of nucleation
!-- Following Shima et al. (2020), we'll use these different types and masses of IN to assign a "freezing temperature" to each droplet
!-- When the ambient temperature (or droplet temperature if using the LEM) falls below the freezing tempearature,
!-- the droplet will freeze into an ice crystal
!-- Mode 1: soot; Mode 2: Mineral Dust; Mode 3: Snowmax bacteria; Mode 4: Pollen
    REAL, DIMENSION(4) :: n_ice         = (/0.0, 0.0, 0.001E6, 0./)    ! concentration ice nuclei (1/m3) for each IN mode
    REAL, DIMENSION(4) :: rm_in         = (/1E-6, 0.2E-6, 1E-6, 0./)      ! Default radius of ice nuclei (m) of each mode (m)
    REAL, DIMENSION(4) :: sigma_in      = (/1.4, 1.4, 1.4, 0./)           ! Geometric std. deviation of ice particle size (-)
    REAL            :: frac_ice         = 0.0        ! fraction of superdroplets which contain IN (-)
    LOGICAL         :: switch_ice       = .FALSE.    ! Ice switch
    LOGICAL         :: switch_habit_ice = .TRUE.     ! Ice habit switch
    ! LOGICAL         :: switch_deterministic_ice = .TRUE. ! Whether the transition from quasi-spherical to habit-bearing ice crystals
                                                     ! is handled deterministically or probablistically - currently inactive
    INTEGER         :: start_ice        = 0          ! start ice physics at timestep start_ice

!-- Secondary ice production (SIP)
    LOGICAL         :: switch_icebreak = .TRUE. !Allow SIP from fragments produced by ice-ice collisions
    CHARACTER(LEN=8):: icebreak_param  = "Phillips" !Which collisional breakup parameterization to use for SIP. Current options are 'Phillips' or 'Sullivan', for the 
                                                    !approaches of Phillips et al., (2017) and Sullivan et al., (2018), respectively
    LOGICAL         :: switch_rs       = .FALSE. !Allow SIP from rime-splintering (Hallet-Mossop process)
    LOGICAL         :: switch_ds       = .TRUE. !Allow SIP from droplet shattering during freezing
    CHARACTER(LEN=10):: ds_param   = "Phillips" !Which shattering parameterization to use. Current options are "Phillips", "Lachapelle", or "Lauber", for the parameterizations
                                                !of Phillips et al. (2018), Lachapelle et al. (2024), or Lauber et al. (2021), respectively
!
!-- Minimum cloud water to identify clouds
    REAL, PARAMETER :: qnmin_LCM = 0.01E-3
!
!-- Separation radius between cloud and rain drops (in meters)
    REAL            :: r_separate_micro = 25.0e-6
!
!-- Wavelength for optical thickness, albedo
    REAL, PARAMETER :: lambda_micro = 500.0E-9
!
!-- Add a particle source

    LOGICAL :: switch_moving_source = .FALSE.
    LOGICAL :: switch_no_background = .TRUE.
    LOGICAL :: switch_particle_source = .FALSE.
    LOGICAL :: switch_phase_sources   = .TRUE. !Whether to use separate CCN and IN source volumes with their own emission rates
    !INTEGER :: n_ptcl_sprayed  = 1
    INTEGER :: n_start_source_aero = -1, n_end_source_aero = -1 
!
!-- Timesteps at which to begin and end spraying CCN (and IN if switch_phase_sources is .FALSE. and frac_in_source_particles is > 0)
!-- Defaults are -1, meaning that spraying will never occur unless set in prm
    INTEGER :: n_start_source_in = -1, n_end_source_in = -1
!
!-- Timesteps at which to begin and end spraying IN if switch_phase_sources if .TRUE.
!-- Defaults are -1, meaning that IN spraying will never occur unless set in prm
!
!-- Location of source in grid point space
    INTEGER, DIMENSION(4) :: nz_b_source_aero = (/1,1,1,1/), nz_t_source_aero = (/nzm,nzm,nzm,nzm/),             &
                             ny_s_source_aero = (/1,1,1,1/), ny_n_source_aero = (/ny_LCM,ny_LCM,ny_LCM,ny_LCM/), &
                             nx_w_source_aero = (/1,1,1,1/), nx_e_source_aero = (/nx_LCM,nx_LCM,nx_LCM,nx_LCM/), &
                             nz_b_source_in = (/1,1,1,1/), nz_t_source_in = (/nzm,nzm,nzm,nzm/),                 &
                             ny_s_source_in = (/1,1,1,1/), ny_n_source_in = (/ny_LCM,ny_LCM,ny_LCM,ny_LCM/),     &
                             nx_w_source_in = (/1,1,1,1/), nx_e_source_in = (/nx_LCM,nx_LCM,nx_LCM,nx_LCM/)
!
!-- Location of source in physical space
    REAL, DIMENSION(4)    :: z_b_source_aero = (/0.0,0.0,0.0,0.0/), z_t_source_aero = (/1.0,1.0,1.0,1.0/), &
                             y_s_source_aero = (/0.0,0.0,0.0,0.0/), y_n_source_aero = (/1.0,1.0,1.0,1.0/), &
                             x_w_source_aero = (/0.0,0.0,0.0,0.0/), x_e_source_aero = (/1.0,1.0,1.0,1.0/), &
                             z_b_source_in = (/0.0,0.0,0.0,0.0/), z_t_source_in = (/1.0,1.0,1.0,1.0/), &
                             y_s_source_in = (/0.0,0.0,0.0,0.0/), y_n_source_in = (/1.0,1.0,1.0,1.0/), &
                             x_w_source_In = (/0.0,0.0,0.0,0.0/), x_e_source_in = (/1.0,1.0,1.0,1.0/)

    REAL, DIMENSION(4)    :: n_source_aero          = (/0.0,0.0,0.0,0.0/)      ! CCN source strength in m^-3 s^-1
    REAL, DIMENSION(4)    :: n_source_in      = (/0.0,0.0,0.0,0.0/)           ! IN source strength (m^-3 s^-1)
    REAL, DIMENSION(4)    :: r_source_aero    = (/0.0625E-6,0.0625E-6,0.0625E-6,0.0625E-6/) ! geometric mean of aerosol source (m)
    REAL, DIMENSION(4)    :: r_source_in      = (/1.0E-6,0.0,0.0,0.0/) !Geometric mean of IN (m)
    REAL, DIMENSION(4)    :: sigma_source_aero = (/2.0,2.0,2.0,2.0/)      ! geometric standard deviation of CCN source
    REAL, DIMENSION(4)    :: sigma_source_in = (/2.0,2.0,2.0,2.0/)      ! geometric standard deviation of IN source
    INTEGER, DIMENSION(4) :: type_source_in  = (/3, 3, 3, 3/)           ! The type of IN emitted by each source
                                                                        ! 1: soot; 2: dust; 3: Snowmax; 4: pollen
    REAL, DIMENSION(4)    :: p_source_aero     = (/0.0,0.0,0.0,0.0/)      ! probability with which new CCN particles are created
    REAL, DIMENSION(4)    :: p_source_in       = (/0.0,0.0,0.0,0.0/)      ! probability with which new IN particles are created
    REAL, DIMENSION(4)    :: vx_source_aero    = (/0.0,0.0,0.0,0.0/), &   ! CCN source velocity in x direction (m/s)
                             vy_source_aero    = (/0.0,0.0,0.0,0.0/), &   ! CCN source velocity in y direction (m/s)
                             vz_source_aero    = (/0.0,0.0,0.0,0.0/), &   ! CCN source velocity in z direction (m/s)
                             vx_source_in    = (/0.0,0.0,0.0,0.0/), &   ! IN source velocity in x direction (m/s)
                             vy_source_in    = (/0.0,0.0,0.0,0.0/), &   ! IN source velocity in y direction (m/s)
                             vz_source_in    = (/0.0,0.0,0.0,0.0/)      ! IN source velocity in z direction (m/s)


    REAL, DIMENSION(4)    :: frac_in_source_particles   = (/0.0,0.0,0.0,0.0/)      ! fraction of source superdroplets containing IN (if switch_phase_sources = .FALSE.)
    REAL, DIMENSION(4)    :: frac_in_source_weight      = (/0.001,0.001,0.001,0.001/) !Factor by which to multiply the number concentration of source particles considered ice
    CHARACTER (LEN=12)    :: unit_source       = '1/kg/s'   ! '1/kg/s' or '1/m3/s' source units
    CHARACTER (LEN=12)    :: type_source       = 'equilibrium' ! type of source: 'equilibrium' or 'seawater'
    REAL                  :: salinity_source   = 0.037 ! salinity of source if type_source = 'seawater'
!
!-- ID counter for each subdomain
    INTEGER(8) :: n_step_ID = HUGE(1_8) / nsubdomains
    INTEGER(8) :: n_ID
!
!-- OUTPUT quantities (profiles)
    REAL, DIMENSION(nz_LCM) :: nsd_prof,   nac_prof,  ncc_prof,   nrc_prof,  cnd_prof,   &
                               evp_prof,   aut_prof,  acc_prof,   spm_prof,  sps_prof,   &
                               sfx_prof,   act_prof,  ninc_prof,  ar_prof,   am_prof,    &
                               dea_prof,   T_prof,    qv_prof,    tau_prof,  ref_prof,   &
                               rsd_prof,   qa_prof,   qc_prof,    qr_prof,   nin_prof,   &
                               afx_prof,   sgfx_prof, refs_prof,  r_prof,    nic_prof,   &
                               qi_prof,    imm_prof,  hom_prof,   hzn_prof,  nuc_prof,   &
                               rhice_prof, rmf_prof, nml_prof,  sub_prof,                &
                               dep_prof,   frz_prof,  mel_prof,   iREF_prof, iR_prof,    &
                               iPHI_prof,  iRSD_prof, iaavg_prof, icavg_prof,            &
                               fl_prof,    fi_prof,   fa_prof,   fin_prof,   dummy_prof, &
                               iREFS_prof, agg_prof,  rim_prof, coal_prof, brk_prof,     &
                               qs_prof,    nsc_prof,   qg_prof,   ngc_prof,              &
                               rho_app_prof, Nmono_prof, nsub_prof, t_max_prof,          &
                               t_mean_prof, t_count_prof, t_mean_prof_liq,               &
                               t_count_prof_liq, t_max_prof_liq, t_mean_prof_ice,        &
                               t_count_prof_ice, t_max_prof_ice,                         &
                               tfrz_prof, ssa_prof, dsa_prof, bsa_prof, psa_prof,        &
                               rs_prof, ds_prof, icebrk_prof, tfrzc_prof, tfrzr_prof,    &
                               tfrzi_prof, tfrzs_prof, tfrzg_prof, fc_prof, fr_prof,    &
                               fp_prof, fs_prof, fg_prof
!
!-- OUTPUT of droplet spectra
    INTEGER :: spec_nsave = 30, spec_nsavestart = 30
    INTEGER :: n_bins_spec
    REAL    :: alpha_spec = 0.5
    REAL    :: r_start_spec = 0.001E-6, r_end_spec = 6000.0E-6
    REAL    :: x_end_spec, x_start_spec
    ! rug spectra (aerosol and pH)
    LOGICAL :: switch_spec_aero = .FALSE.
#ifdef SAMCHEM
    LOGICAL :: switch_spec_pH = .FALSE.
    REAL    :: pH_start_spec = -3.0, pH_end_spec = 10.0
    REAL, DIMENSION(:), ALLOCATABLE :: pHl_spec, pHm_spec, pHr_spec
    REAL(KIND=4), DIMENSION(:,:,:,:), ALLOCATABLE :: spec_pH
#endif /*SAMCHEM*/

    ! switch to track particles
    LOGICAL :: track_some_particles = .FALSE.
    INTEGER :: n_track_particles    = 1003,   & ! number of particles to track
             & particle_ntrack      = 100,    & ! turnus of saving particle tracking data in steps
             & particle_ntrackstart = 300,    & ! number of steps after which to start particle tracking
             & particle_ntrackstop  = 999999999 ! number of steps after which to start particle tracking

    REAL        , DIMENSION(:),       ALLOCATABLE :: rl_spec, rm_spec, rr_spec
    REAL(KIND=4), DIMENSION(:,:,:,:), ALLOCATABLE :: spec_liq, spec_ice, spec_aero
    REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: spec_global

    CHARACTER (LEN=128) :: format_spec
    CHARACTER (LEN=128) :: format_particle_tracking

    LOGICAL :: horiz_average_spec = .TRUE.

!
!-- Freezing function for IN
    REAL, EXTERNAL :: t_frz

!
!-- flags for output time steps
!-- (only calculate output data if needed at current time step, initialized true for safety)
    LOGICAL :: t_out_stat       = .TRUE., & ! true if any output to OUT_STAT  at current time step
             & t_out_field      = .TRUE., & ! true if any output to OUT_xD    at current time step
             & t_out_spec       = .TRUE., & ! true if any output of spectra   at current time step
             & t_out_particle   = .TRUE., & ! true if any output to particles at current time step
             & t_out_any        = .TRUE.    ! true if any output              at current time step

!===============================================================================
!                 NO EDIT BELOW (unless you know what you are doing)
!===============================================================================
!
!-- Horizontal additional grid points for u, v, w, and scalar arrays.
!-- Parcels can be out of the current subdomain for the interpolation at half
!-- time step or guessed new position for resolved scale velocity iteration.
!--
!-- SAM's horizontal grid index orientation (s=scalar):
!--
!--         ---------------v(i=1,j=2)----------------
!--        |                                         |
!--   u(i=1,j=1)     w(i=1,j=1), s(i=1,j=1)     u(i=2,j=1)
!--        |                                         |
!--         ---------------v(i=1,j=1)----------------
!--
!-- All additional grid points have to be more than or equal to 1 for MPI exchange routine
!--
!-- THESE ADDITIONAL POINTS ARE FOR 5TH LAGRANGIAN POLYNOMIAL INTERPOLATION.
!-- IT WORKS FOR ALL AVAILABLE INTERPOLATION SCHEMES.
       INTEGER, PARAMETER :: ux1 = 3
       INTEGER, PARAMETER :: ux2 = 4
       INTEGER, PARAMETER :: uy1 = 3 * YES3D + 1
       INTEGER, PARAMETER :: uy2 = 3 * YES3D + 1
       INTEGER, PARAMETER :: vx1 = 3 + YES3D
       INTEGER, PARAMETER :: vx2 = 4
       INTEGER, PARAMETER :: vy1 = 2 * YES3D + 1
       INTEGER, PARAMETER :: vy2 = 3 * YES3D + 1
       INTEGER, PARAMETER :: wx1 = 4
       INTEGER, PARAMETER :: wx2 = 4
       INTEGER, PARAMETER :: wy1 = 3 * YES3D + 1
       INTEGER, PARAMETER :: wy2 = 3 * YES3D + 1
       INTEGER, PARAMETER :: sx1 = 4
       INTEGER, PARAMETER :: sx2 = 4
       INTEGER, PARAMETER :: sy1 = 3 * YES3D + 1
       INTEGER, PARAMETER :: sy2 = 3 * YES3D + 1
!
!-- Horizontal indices for u, v, w, and scalar arrays. Index for y direction
!-- is 1 for 2D run.
       INTEGER, PARAMETER :: x1u =  1 - ux1
       INTEGER, PARAMETER :: x2u = nx_LCM + ux2
       INTEGER, PARAMETER :: y1u =  1 - uy1 * YES3D
       INTEGER, PARAMETER :: y2u = ny_LCM + uy2 * YES3D
       INTEGER, PARAMETER :: x1v =  1 - vx1
       INTEGER, PARAMETER :: x2v = nx_LCM + vx2
       INTEGER, PARAMETER :: y1v =  1 - vy1 * YES3D
       INTEGER, PARAMETER :: y2v = ny_LCM + vy2 * YES3D
       INTEGER, PARAMETER :: x1w =  1 - wx1
       INTEGER, PARAMETER :: x2w = nx_LCM + wx2
       INTEGER, PARAMETER :: y1w =  1 - wy1 * YES3D
       INTEGER, PARAMETER :: y2w = ny_LCM + wy2 * YES3D
       INTEGER, PARAMETER :: x1s =  1 - sx1
       INTEGER, PARAMETER :: x2s = nx_LCM + sx2
       INTEGER, PARAMETER :: y1s =  1 - sy1 * YES3D
       INTEGER, PARAMETER :: y2s = ny_LCM + sy2 * YES3D
!
!-- Subdomain grid position information
    INTEGER :: nx_sb, ny_sb   ! total number of grid point to the current subdomain from (0,0)
!
!-- Grid center
    REAL :: gpxc(-1:nx_LCM+1)
    REAL :: gpyc(-1:ny_LCM+1)
    REAL :: gpzc(nz_LCM)
!
!-- Grid side (wall)
    REAL :: gpxs(0:nx_LCM+1)   ! 1 & nx+1 are at subdomain boundary
    REAL :: gpys(0:ny_LCM+1)
    REAL :: gpzs(nz_LCM+1)     ! bottom of grid box
!
!-- Horizontal width of whole domain
    REAL :: xwidth
    REAL :: ywidth
!
!-- Inverse of grid space
    REAL :: idx
    REAL :: idy
!
!-- Globals start id's of i and j
    INTEGER :: i_global_LCM, j_global_LCM
!
!-- Volume of a grid box
    REAL :: dV(nz_LCM)
!
!-- Limits of the particle domain
    REAL, PARAMETER :: boundary_distance = 0.001 ! standard distance to boundaries below which particles will be deleted, reflected, ...
    REAL            :: ztop, zbottom
!
!-- MPI parameters
    INTEGER               :: newtype
    INTEGER, DIMENSION(8) :: n_out
    INTEGER               :: np_out, np_in, np_deleted
!
!-- Lagrangian particle information
    INTEGER :: n_init            ! # of particles for subdomain at initial time step
    INTEGER :: n_total           ! # of particles for whole domain
    INTEGER :: n_particles       ! # of particles for subdomain at current time step
    INTEGER :: n_particles_alloc ! # of particles allocated for subdomain at current time step
!
!-- Number of Lagrangian particle features
    INTEGER, PARAMETER :: nv0  = 9  ! # of default variables: position, resolved velocities
#ifdef SAMCHEM
    INTEGER, PARAMETER :: nv1  = 20 + naqspc  ! # of microphysical variables
#else
    INTEGER, PARAMETER :: nv1  = 20 ! # of microphysical variables
#endif /*SAMCHEM*/
    INTEGER, PARAMETER :: niv0 = 10  ! # of integer variables: tag, sgs tag
!
!-- Indices of Lagrangian particle features
!-- Motion:
    INTEGER, PARAMETER :: iX    = 1 ! index of SD's x location
    INTEGER, PARAMETER :: iY    = 2 ! index of SD's y location
    INTEGER, PARAMETER :: iZ    = 3 ! index of SD's z location
    INTEGER, PARAMETER :: iU    = 4 ! index of SD's u velocity
    INTEGER, PARAMETER :: iV    = 5 ! index of SD's v velocity
    INTEGER, PARAMETER :: iW    = 6 ! index of SD's w velocity
    INTEGER, PARAMETER :: iUsgs = 7 ! index of SD's SGS u velocity
    INTEGER, PARAMETER :: iVsgs = 8 ! index of SD's SGS v velocity
    INTEGER, PARAMETER :: iWsgs = 9 ! index of SD's SGS w velocity
!
!-- Microphysics (general):
    INTEGER, PARAMETER :: iM     = ( 9 ) +  1 ! index of SD's total mass
    INTEGER, PARAMETER :: iN     = ( 9 ) +  2 ! index of SD's total aerosol mass
    INTEGER, PARAMETER :: iWsedi = ( 9 ) +  3 ! index of SD's terminal velocity
!
!-- Microphysics (ice):
    INTEGER, PARAMETER :: iAXS  = ( 9 ) +  4 ! index of SD's a axis
    INTEGER, PARAMETER :: iCXS  = ( 9 ) +  5 ! index of SD's c axis
    INTEGER, PARAMETER :: iDNS  = ( 9 ) +  6 ! index of SD's apparent density
    INTEGER, PARAMETER :: iMrime  = ( 9 ) +  7 ! index of SD's rimed mass
    INTEGER, PARAMETER :: iSC   = ( 9 ) +  8 ! Index of SD's black carbon (soot) IN surface area
    INTEGER, PARAMETER :: iSD   = ( 9 ) +  9 ! Index of SD's dust IN surface area
    INTEGER, PARAMETER :: iSB   = ( 9 ) + 10 ! Index of SD's Snowmax (bacterial) IN surface area
    INTEGER, PARAMETER :: iSP   = ( 9 ) + 11 ! Index of SD's pollen IN surface area
    INTEGER, PARAMETER :: iTFZ   = ( 9 ) +  12 ! Index of SD's freezing temperature
!
!-- Microphysics (SGS model):
    INTEGER, PARAMETER :: iT    = ( 9 ) +  13 ! index of SD's temperature
    INTEGER, PARAMETER :: iq    = ( 9 ) +  14 ! index of SD's humidity
    INTEGER, PARAMETER :: ieta  = ( 9 ) +  15 ! index of SD's last experienced absolute supersaturation

!
!-- Microphysics (timing):
    INTEGER, PARAMETER :: itime = ( 9 ) + 16 ! index of the superdroplet's cloudy lifetime
!
!-- Aerosol:
    INTEGER, DIMENSION(4), PARAMETER :: iAero_portions = ( 9 ) + (/ 17, 18, 19, 20 /) ! mass portions the aerosol has from each mode (0<values<1, sum(values)=1)
#ifdef SAMCHEM
!
!-- Chemistry (aqueous chemical species)
    INTEGER :: iiii
    INTEGER, DIMENSION(naqspc), PARAMETER :: iaqspc  =  ( 9 ) + nv1 - naqspc + (/ (iiii, iiii = 1, naqspc) /) ! indices of aqueous chemical species
#endif /*SAMCHEM*/

!
!-- Integers:
    INTEGER, PARAMETER :: iA    = 1 ! index of SD's weighting factor

    INTEGER, PARAMETER :: iLEM  = 2 ! index of SGS tag
    INTEGER, PARAMETER :: iilst = 3 ! index of last i location
    INTEGER, PARAMETER :: ijlst = 4 ! index of last j location
    INTEGER, PARAMETER :: iklst = 5 ! index of last k location

    INTEGER, PARAMETER :: iID   = 6 ! unique particle ID

    INTEGER, PARAMETER :: iTYPE = 7   ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
    INTEGER, PARAMETER :: iNmono  = 8 ! number of monomers
    INTEGER, PARAMETER :: iSPR  = 9   ! marking the sprayed aerosols

!
!-- Integer denoting the index one before the IN surface area variables - used to update/access IN information
!-- IMPORTANT! CHANGE THIS IF PARTICLE QUANTITIES ARE ADDED/REMOVED BETWEEN iMRIME and iSC!!!!
    INTEGER, PARAMETER :: in_startind = iSC - 1
!
    INTEGER, PARAMETER :: iTRACK  = 10 ! marking the tracked aerosols
!
!
!-- Derived type pointer array for particles. Number of particles changes with
!-- time, so pointer-arrays and pointers are used in stead of actual array.
!-- This is probably efficient because the linked-lists for outgoing particles
!-- can be made while counting outgoing particles right before MPI exchages.
    TYPE particle_pointer
       INTEGER(8)      :: ivar(niv0+nivU)   ! unique integer for each particle
       REAL(KIND=lcmp) :: var(nv0+nv1+nvU)
       TYPE(particle_pointer), POINTER :: next
    END TYPE particle_pointer

    TYPE(particle_pointer), POINTER, DIMENSION(:) :: particle
    TYPE(particle_pointer), POINTER, DIMENSION(:) :: particle_in
    TYPE(particle_pointer), POINTER               :: particle_head, particle_tail

    TYPE(particle_pointer), POINTER :: particle_out_nn_head, particle_out_nn_tail, &
                                       particle_out_ne_head, particle_out_ne_tail, &
                                       particle_out_ee_head, particle_out_ee_tail, &
                                       particle_out_se_head, particle_out_se_tail, &
                                       particle_out_ss_head, particle_out_ss_tail, &
                                       particle_out_sw_head, particle_out_sw_tail, &
                                       particle_out_ww_head, particle_out_ww_tail, &
                                       particle_out_nw_head, particle_out_nw_tail
!
!-- Derived type particle array for restart I/O and MPI
    TYPE particle_array
       INTEGER(8)      :: ivar(niv0+nivU)
       REAL(KIND=lcmp) :: var(nv0+nv1+nvU)
    END TYPE particle_array
!
!-- Particles are sorted by their grid box for efficient calculation of
!-- microphysical processes. The following fields indicate the start/end index
!-- of particles inside grid box i, j, k
    INTEGER :: n_start(nx_LCM,ny_LCM,nz_LCM)
    INTEGER :: n_end(nx_LCM,ny_LCM,nz_LCM)
!
!-- flag to indicate initialization of LCM has been carried out
    LOGICAL :: flag_LCM_init  = .FALSE.
!
!--  Output information. File name and file extension:
    CHARACTER (LEN=200) :: filename_LCM, filename_LCM_ext
!
!-- Current number of total history length:
    INTEGER             :: history_length !RELEVANT?

!
! -- handy constant to reduce redundant multiplications
     REAL, DIMENSION(4) :: bfactor_coefficient=0.0
     LOGICAL :: only_one_aerosol_mode=.FALSE.
     INTEGER :: n_aerosol_modes=4

    CONTAINS

    ! this function calculates a bfactor for Köhler theory needed for the solute effect
    ! the solute effect should, provided the result of this function, be calculated like bfactor/r_water**3 (no r_aero!)
    FUNCTION get_bfactor(n) RESULT(bfactor)
      ! input
      INTEGER, INTENT(IN) :: n ! index of particle to get bfactor of
      ! output
      REAL :: bfactor
      ! internal
      INTEGER :: iMode

      IF (only_one_aerosol_mode) THEN
          bfactor = bfactor_coefficient(1)
      ELSE
        bfactor = 0.0
        DO iMode = 1, n_aerosol_modes
          !bfactor = SUM(vanthoff_aero*particle(n)%var(iAero_portions)/molecular_weight_aero)*(particle(n)%var(iN)/particle(n)%ivar(iA)) * molecular_weight_water / (rho_liq*4*pi/3)
          bfactor = bfactor + bfactor_coefficient(iMode)*particle(n)%var(iAero_portions(iMode))
        END DO
      END IF

      bfactor = bfactor * particle(n)%var(iN)/particle(n)%ivar(iA)
    END FUNCTION get_bfactor

    ! this function calculates the equivalent radius of an aerosol contained in a particle
    ! (equivalent because it is assumed to be dissolved in the water the particle contains)
    FUNCTION get_r_aero(n) RESULT(r_aero)
      USE params, ONLY: pi
      ! input
      INTEGER, INTENT(IN) :: n ! index of particle to calculate the aerosol radius of
      ! output
      REAL :: r_aero

      r_aero = (get_V_aero(n) * 3.0 / 4.0 / pi)**(0.3333333333333)
    END FUNCTION get_r_aero

    ! this function calculates the radius of a particle
    FUNCTION get_r_particle(n) RESULT(r_particle)
      USE params, ONLY: pi
      ! input
      INTEGER, INTENT(IN) :: n ! index of particle to calculate the radius of
      ! output
      REAL :: r_particle
      ! internal
      REAL :: V_particle

      ! aqueous volume
      V_particle = particle(n)%var(iM)/particle(n)%ivar(iA)/rho_liq
      ! if aerosol mass exceeds 10% of water mass consider aerosol volume as well
      ! (makes an error in order of 1% (1.05**0.333) in radius for density of aerosol twice the density of water)
      ! otherwise ignore aerosol portion (doesn't make a difference probably..)
      IF (particle(n)%var(iN)/particle(n)%var(iM)>0.1) THEN
        V_particle = V_particle + get_V_aero(n)
      END IF
      r_particle = (V_particle * 3.0 / 4.0 / pi)**0.3333333333333
    END FUNCTION get_r_particle

    ! this function calculates the equivalent volume of an aerosol contained in a particle
    ! (equivalent because it is dissolved in the water the particle contains)
    FUNCTION get_V_aero(n) RESULT(V_aero)
      ! input
      INTEGER, INTENT(IN) :: n ! index of particle to calculate the aerosol volume of
      ! output
      REAL :: V_aero
      ! internal
      INTEGER :: i

      V_aero = 0.0
      DO i = 1, n_aerosol_modes
        V_aero = V_aero + particle(n)%var(iAero_portions(i))*particle(n)%var(iN) / rho_aero(i) / particle(n)%ivar(iA)
      END DO
    END FUNCTION get_V_aero


 REAL FUNCTION r_equi( S, T, r_aero, n )
!
!-- Calculates equilibrium radius using an iterative solver
    USE params,     ONLY: rv
    REAL, INTENT(IN)  :: S, T, r_aero
    INTEGER, OPTIONAL :: n
    REAL             :: afactor, bfactor, C, r_equi_0
!
!-- Curvature effect (afactor)
    afactor = 2.0 * sigma / ( rho_liq * rv * T )
!
!-- Solute effect (bfactor)
    IF (PRESENT(n)) THEN
      bfactor = get_bfactor(n)
    ELSE
      bfactor = vanthoff_aero(1) * rho_aero(1) * molecular_weight_water &
              / ( rho_liq * molecular_weight_aero(1) ) * r_aero**3
    END IF
!
!-- Iterative solver ( 0 = S - A / r + B / r^3 => r = ( B / ( A / r - S ) )^(1/3) )
    r_equi_0 = 1.0
    r_equi   = 1.0E-6
    DO WHILE ( ABS( ( r_equi - r_equi_0 ) / r_equi_0 ) .GT. 1.0E-20 )
       r_equi_0 = r_equi
       r_equi   = ( bfactor / ( afactor / r_equi - S ) )**(1.0/3.0)
    END DO

 END FUNCTION r_equi

 FUNCTION search_sorted_array(array, val, N) RESULT(ind)
   ! perform a binary search on a sorted array, way faster than minloc (bc using sortedness)
   ! returns index i where array(i)<val<array(i+1) (because we need it to get k index of grid boxes)
   REAL, DIMENSION(:), INTENT(IN) :: array
   REAL              , INTENT(IN) :: val
   INTEGER           , INTENT(IN) :: N

   INTEGER            :: ind

   INTEGER :: start, mid, finish, range

   start = 1
   finish = N
   range=finish-start
   mid = (start + finish)/2 ! this is integer division (without remainder)

  DO WHILE (range>1)

     IF (val > array(mid)) then
       start = mid
     ELSE
       finish = mid
     END IF
     range = finish - start
     mid = (start + finish)/2 ! this is integer division (without remainder)

  END DO

  ind = start
 END FUNCTION search_sorted_array


 END MODULE micro_vars
