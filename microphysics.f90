!TODO: make it possible to start LCM at first time step: pull preparations infront of
!      the actual microphysics

 MODULE microphysics

    USE grid,   ONLY: nx, ny, nz, nzm, dimx1_s, dimx2_s, dimy1_s, dimy2_s,     &
                     compute_reffc, compute_reffi, z
    USE params, ONLY: doprecip

    USE micro_vars

#ifdef SAMCHEM
    USE chemistry, ONLY: chem_supply_caretaking, chem_set_initial_concentrations, chem_monster_proc, flag_chem_init, nstep_chem_init
    USE domain,    ONLY: naqspc
#endif /*SAMCHEM*/

    IMPLICIT none
!
!-- Dummy arguments for cloud optics option added by P. Blossey for SAM6.10.0
!-- Effective radius - reffc should be filled during micro_proc and refereed in Get_reffc
    REAL, PUBLIC, DIMENSION(nx,ny,nzm)          :: reffc ! microns
    REAL, PUBLIC, DIMENSION(nx,ny,nzm) :: reffi
!-- B. Ascher 12.05.2025 - Previously, reffi was set to be an allocatable array. However, this caused
!-- crashes when using RAD_CAM. As far as I'm aware, since this is a dummy argument that doesn't impact the rest
!-- of the code, it shouldn't matter if we just set it to an array here, aside from consuming very slightly more memory
!
!-- Flags and variables related to M2005 cloud optics routines, included here so that RRTM will
!-- compile against this microphysics as well.
    LOGICAL, PUBLIC :: dosnow_radiatively_active = .FALSE.
    LOGICAL, PUBLIC :: dorrtm_cloud_optics_from_effrad_LegacyOption = .TRUE.


    REAL, PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: CloudLiquidMassMixingRatio,     &
                                                   CloudLiquidGammaExponent,       &
                                                   CloudLiquidLambda,              &
                                                   CloudIceMassMixingRatio,        &
                                                   SnowMassMixingRatio, reffs
!
!-- Total number of prognostic variables
    INTEGER, PARAMETER :: nmicro_fields = 6
!
!-- Number of additional diagnostic variables to micro_field, which will be used for statistics and
!-- 3D output. Currently 4 additional variables: qc, qr, nc, nr
#ifdef SAMCHEM
   INTEGER, PARAMETER :: nmicro_fields_diag = 6+6+2+2+4+1+3+1+2+2+9+3+4+5+6+6+5+3+6 !start counting from ina onwards
#else
   INTEGER, PARAMETER :: nmicro_fields_diag = 6+6+2+2+4+1+3+1+2+2+9+3+4+5+6+6+5+3 !start counting from ina onwards
#endif /*SAMCHEM*/
!
!-- Index for variables used within this module
    INTEGER, PARAMETER :: imv      = 1  ! water vapor mixing ratio              (12)
    INTEGER, PARAMETER :: iTlcm    = 2  ! LCM static energy                     (13)
    INTEGER, PARAMETER :: iTllcm   = 3  ! LCM liquid water static energy        (14)

    INTEGER, PARAMETER :: ilastT   = 4 ! advected temperature  (for forcings)   (15)
    INTEGER, PARAMETER :: ilastp   = 5 ! advected pressure     (for forcings)   (16)

    INTEGER, PARAMETER :: iSCLR    = 6  ! scalar for detecting plumes/thermals  (17)

    INTEGER, PARAMETER :: ina      = 7  ! aerosol number mixing ratio           (18)
    INTEGER, PARAMETER :: iqa      = 8  ! aerosol water mass mixing ratio       (19)
    INTEGER, PARAMETER :: inc      = 9  ! cloud droplet number mixing ratio     (20)
    INTEGER, PARAMETER :: iqc      = 10  ! cloud water mass mixing ratio        (21)
    INTEGER, PARAMETER :: inr      = 11 ! rain drop number mixing ratio         (22)
    INTEGER, PARAMETER :: iqr      = 12 ! rain water mass mixing ratio          (23)

    INTEGER, PARAMETER :: ini      = 13 ! cloud ice number mixing ratio         (24)
    INTEGER, PARAMETER :: iqi      = 14 ! cloud ice mass mixing ratio           (25)
    INTEGER, PARAMETER :: ins      = 15 ! snow number mixing ratio              (26)
    INTEGER, PARAMETER :: iqs      = 16 ! snow water mass mixing ratio          (27)
    INTEGER, PARAMETER :: ing      = 17 ! graupel/hail number mixing ratio      (28)
    INTEGER, PARAMETER :: iqg      = 18 ! graupel/hail number mixing ratio      (29)

    INTEGER, PARAMETER :: ipl      = 19 ! liquid precipitation flux             (30)
    INTEGER, PARAMETER :: ipi      = 20 ! frozen precipitation flux             (31)

    INTEGER, PARAMETER :: iRsd     = 21 ! effective radius standard deviation   (32)
    INTEGER, PARAMETER :: iiREF    = 22 ! ice effective radius                  (33)

    INTEGER, PARAMETER :: isLES    = 23 ! LES supersaturation                   (34)
    INTEGER, PARAMETER :: issp     = 24 ! standard deviation supersaturation    (35)
    INTEGER, PARAMETER :: imsp     = 25 ! mean supersaturation                  (36)
    INTEGER, PARAMETER :: iiss     = 26 ! supersaturation w.r.t ice             (37)

    INTEGER, PARAMETER :: idiss    = 27 ! dissipation rate                      (38)

    INTEGER, PARAMETER :: iPlume     = 28 ! plume grid box                      (39)
    INTEGER, PARAMETER :: iCondensed = 29 ! condesed grid box                   (40)
    INTEGER, PARAMETER :: iCore      = 30 ! core grid box                       (41)

    INTEGER, PARAMETER :: ininc      = 31 ! ice nuclei concentration            (42)

    INTEGER, PARAMETER :: iasp      = 32 ! number of sprayed particles               (43)
    INTEGER, PARAMETER :: icsp      = 33 ! number of droplets from sprayed particle  (44)
    
    INTEGER, PARAMETER :: iact      = 34 ! activation rate                           (45)
    INTEGER, PARAMETER :: idea      = 35 ! deactivation rate                         (46)
    
    INTEGER, PARAMETER :: itfrz     = 36 ! freezing temperature                      (47)
    INTEGER, PARAMETER :: issa      = 37 ! Soot IN surface area                      (48)
    INTEGER, PARAMETER :: idsa      = 38 ! Dust IN surface area                      (49)
    INTEGER, PARAMETER :: ibsa      = 39 ! Snowmax (bacterial) IN surface area       (50)
    INTEGER, PARAMETER :: ipsa      = 40 ! Pollen IN surface area (WIP)              (51)
    INTEGER, PARAMETER :: iimm      = 41 !Immersion freezing rate                    (52)
    INTEGER, PARAMETER :: ihzn      = 42 !Haze/Deposition nucleation rate            (53)
    INTEGER, PARAMETER :: ihom      = 43 !Homogeneous nucleation rate                (54)
    INTEGER, PARAMETER :: inuc      = 44 !Total ice nucleation rate                  (55)

    INTEGER, PARAMETER :: irim      = 45 !Riming rate                                (56)
    INTEGER, PARAMETER :: ifrz      = 46 !Total freezing rate (includes riming)      (57)
    INTEGER, PARAMETEr :: imel      = 47 !Total melting rate                         (58)

    INTEGER, PARAMETER :: icon      = 48 !Condensation rate                          (59)
    INTEGER, PARAMETER :: ievp      = 49 !Evaporation rate                           (60)
    INTEGER, PARAMETER :: idep      = 50 !Deposition rate                            (61)
    INTEGER, PARAMETER :: isub      = 51 !Sublimation rate                           (62)
    
    INTEGER, PARAMETER :: iphi      = 52 !Ice effective aspect ratio                 (63)
    INTEGER, PARAMETER :: iaavg     = 53 !Ice average a-axis length                  (64)
    INTEGER, PARAMETER :: icavg     = 54 !Ice average c-axis length                  (65)
    INTEGER, PARAMETER :: irhoapp   = 55 !Ice apparent density                       (66)
    INTEGER, PARAMETER :: irmf      = 56 !Rime mass fraction                         (67)

    INTEGER, PARAMETER :: itmean    = 57 !Mean particle lifetime                     (68)
    INTEGER, PARAMETER :: itmax     = 58 !Max particle lifetime                      (69)
    INTEGER, PARAMETER :: iltmean   = 59 !Mean liquid particle lifetime              (70)
    INTEGER, PARAMETER :: iltmax    = 60 !Max liquid particle lifetime               (71)
    INTEGER, PARAMETER :: iitmean   = 61 !Mean ice particle lifetime                 (72)
    INTEGER, PARAMETER :: iitmax    = 62 !Max ice particle lifetime                  (73)

    INTEGER, PARAMETER :: icoal     = 63 !Rate of coalescences                       (74)
    INTEGER, PARAMETER :: ibreak    = 64 !Rate of droplet collisional breakups       (75)
    INTEGER, PARAMETER :: iagg      = 65 !Rate of aggregation                        (76)
    INTEGER, PARAMETER :: irs       = 66 !Rate of rime splintering                   (77)
    INTEGER, PARAMETER :: ids       = 67 !Rate of droplet shattering                 (78)
    INTEGER, PARAMETER :: iicebreak = 68 !Rate of ice collisional breakup            (79)

    INTEGER, PARAMETER :: itfrzc    = 69 !Mean freezing temperature of cloud droplets(80)
    INTEGER, PARAMETER :: itfrzr    = 70 !Mean freezing temperature of rain          (81)
    INTEGER, PARAMETER :: itfrzi    = 71 !Mean freezing temperature of pristine ice  (82)
    INTEGER, PARAMETER :: itfrzs    = 72 !Mean freezing temperature of snow          (83)
    INTEGER, PARAMETER :: itfrzg    = 73 !Mean freezing temperature of graupel/hail  (84)

    INTEGER, PARAMETER :: iss       = 74 ! current supersaturation                   (85)
    INTEGER, PARAMETER :: isr       = 75 ! aerosol radius                            (86)
    INTEGER, PARAMETER :: ism       = 76 ! aerosol mass                              (87)
#ifdef SAMCHEM
    INTEGER, PARAMETER :: ipH       = 77 ! current pH in grid boxes                  (88)
    INTEGER, PARAMETER :: idelR     = 78 ! current pH in grid boxes                  (89)
    INTEGER, PARAMETER :: iSbyO3t    = 79 ! total O3->S(VI) in grid boxes            (90)
    INTEGER, PARAMETER :: iSbyO3m    = 80 ! mean O3->S(VI) in grid boxes             (91)
    INTEGER, PARAMETER :: iSbyH2O2t  = 81 ! total H2O2->S(VI) in grid boxes          (92)
    INTEGER, PARAMETER :: iSbyH2O2m  = 82 ! mean H2O2->S(VI) in grid boxes           (93)
#endif /*SAMCHEM*/

    REAL,    PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: micro_field
!
!-- Flag for water mixing ratio, regardless of water species
!-- Use a simple array of flags with 1 marking the water mass variable
    INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       :: flag_wmass
!
!-- Flag for number concentration
    INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       :: flag_number
!
!-- Flag for variables associated with precipipitation
    INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       :: flag_precip
!
!-- Index for variable that contains water vapor
    INTEGER, PARAMETER :: index_water_vapor = 1
!
!-- Index for cloud ice (sedimentation)
    INTEGER, PARAMETER :: index_cloud_ice = -1
!
!-- Here number concetnrations are added. They are output if docloud.
    INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: flag_micro3Dout
!
!-- The following arrays are needed to hold the turbulent surface and domain-top fluxes for the
!-- microphysics prognostic variables:
    REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: fluxbmk ! surface fluxes
    REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: fluxtmk ! top boundary fluxes
#ifdef PICHAMBER
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: fluxlmk, fluxrmk !left/right fluxes
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: fluxqmk, fluxhmk !front/back fluxes
#endif /*PICHAMBER*/
!
!-- These arrays are needed for output statistics from advection and diffusion routines
    REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: mkwle   ! resolved vertical flux
    REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: mkwsb   ! SGS vertical flux
    REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: mkadv   ! tendency due to vertical advection
    REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: mkdiff  ! tendency due to vertical diffusion
    REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: mklsadv ! tendency due to subsidence
!
!-- Arrays with names/units for microphysical outputs in statistics.
    CHARACTER(LEN=10),  PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: mkname
    CHARACTER(LEN=80), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: mklongname
    CHARACTER(LEN=10), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: mkunits
    REAL,              PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: mkoutputscale

    CONTAINS

 SUBROUTINE micro_setparm()
!
!-- Read microphysics options from prm file
    USE micro_vars, ONLY: np_xy, nzl, nstep_LCM_init, n_uvw_iteration,         &
                          switch_sgs_velocities, r_separate_micro,             &
                          n_init, rm_aero, sigma_aero,                         &
                          switch_sedimentation, switch_collection,             &
                          start_collection, start_sedimentation, init_aero,    &
                          switch_radiation_micro, switch_analytical_micro,     &
                          switch_habit_ice, switch_diskice_sizethresh,          &
                          sizethresh_diskice, sizethresh_chen

    USE micro_random, ONLY: randomu

    USE grid,       ONLY: case, masterproc, rank
    use params

    IMPLICIT NONE

   INTEGER :: ierr, ios, ios_missing_namelist, place_holder, iMode, i
   CHARACTER(len=1000) :: line
   CHARACTER(LEN=10), DIMENSION(3), PARAMETER :: habitparamlist = (/"Chen      ", "DISKICE   ", "Weiss     "/)
!
!-- List namelist variables
   NAMELIST /MICRO_LAGRANGE/ &
   np_xy, nzl, nstep_LCM_init, n_uvw_iteration, switch_sgs_velocities,     &
   r_separate_micro, n_aero, rm_aero, shares_aero, sigma_aero,             &
   n_ice, rm_in, sigma_in,                                                 &
   switch_sedimentation, spec_nsave, spec_nsavestart, switch_LEM_micro,    &
   switch_collection, start_collection, start_sedimentation, init_aero,    &
   switch_radiation_micro, switch_analytical_micro,                        &
   switch_steady_pichamber, switch_cloud_tracking, switch_phase_sources,   &
   switch_particle_source, nz_b_source_aero,                               &
   nz_t_source_aero, ny_s_source_aero, ny_n_source_aero,                   &
   nx_w_source_aero, nx_e_source_aero, nz_b_source_in, nz_t_source_in,     &
   ny_s_source_in, ny_n_source_in, nx_w_source_in, nx_e_source_in,         &
   n_start_source_aero, n_end_source_aero, n_source_aero,                  &
   n_start_source_in, n_end_source_in, n_source_in, r_source_aero,         &
   r_source_in, alpha_spec, r_start_spec, switch_walls_pichamber,          &
   switch_reflective_walls_pichamber, p_source_aero, p_source_in,          &
   switch_koehler_micro, switch_linear_sampling_collection,                &
   sigma_source_aero, sigma_source_in, type_source_in,                     &
   switch_FFC_micro, horiz_average_spec,                                   &
   frac_ice, frac_in_source_particles, frac_in_source_weight,              &
   switch_ice, switch_habit_ice, habitparam, start_ice, switch_icebreak,   &
   switch_rs, switch_ds, ds_param, icebreak_param,                         &
   unit_source, unit_aero, switch_moving_source,                           &
   vx_source_aero, vy_source_aero, vz_source_aero,                         &
   vx_source_in, vy_source_in, vz_source_in, salinity_source, type_source, &
   r_end_spec, switch_diskice_sizethresh, sizethresh_chen, sizethresh_diskice, &
   switch_breakup, breakup_conserve_mode, switch_steady_aerosol,           &
   inject_bottom, inject_top, switch_spec_aero,                            &
#ifdef SAMCHEM
   switch_spec_pH,                                                         &
#endif /*SAMCHEM*/
   track_some_particles, n_track_particles, particle_ntrack,               &
   particle_ntrackstart, rho_aero, vanthoff_aero, molecular_weight_aero,   &
   particle_ntrackstop, nstep_cut_supersat, supersat_threshold

   ! NAMELIST /PARAMETERS/ &
   !   caseid, nrestart, caseid_restart, case_restart, restart_sep, nrestart_skip, &
   !    nstop, nelapse, nelapsemin, day0, ncycle_max, nprint, nstat, nstatfreq, &
   !    doSAMconditionals, dosatupconditionals, doisccp, domodis, domisr,          &
   !    dosimfileout, output_sep, nsave2D, nsave2Dstart, nsave2Dend, save2Dbin,    &
   !    save2Dsep, save2Davg, dogzip2D, nsave3D, nsave3Dstart, nsave3Dend,         &
   !    save3Dbin, save3Dsep, dogzip3D, nstatmom, nstatmomstart, nstatmomend,      &
   !    savemombin, savemomsep, nmovie, nmoviestart, nmovieend, dx, dy, dz,        &
   !    dz_constant, dt, LES, CRM, LAND, OCEAN, soil_wetness, z0, dowally, dowallx,&
   !    rundatadir, dosgs, doscalar, dodamping, doupperbound, docloud, doprecip,   &
   !    longitude0, latitude0, docoriolis, docoriolisz, dofplane, fcor,            &
   !    dolongwave, doshortwave, nrad, doperpetual, dosolarconstant,               &
   !    zenith_angle, doseasons, doradforcing, doradhomo, doradlon, roradlat,      &
   !    nxco2, dotracers, dosmoke, doradsimple, dosurface, dolargescale,           &
   !    dosfcforcing, ocean_type, doseawater, tabs_s, delta_sst, dodynamic,        &
   !    timesimpleocean, depth_slab_ocean, Szero, deltas, donudging_uv,            &
   !    donudging_tauz_uv, donudging_tq, donudging_tauz_tq, donudging_t,           &
   !    donudging_q, nudging_uv_z1, nudging_t_z1, nudging_q_z1, tauls, tautqls,    &
   !    doensemble, nensemble, SFC_FLX_FXD, SFC_TAU_FLX, SFC_TAU_FXD, tau0,        &
   !    fluxt0, fluxq0, ug, vg, perturb_type, bubble_x0, bubble_y0, bubble_z0,     &
   !    bubble_radius_hor, bubble_radius_ver, bubble_dtemp, bubble_dq,             &
   !    nudging_tauz_z1, nudging_tauz_z2, tauz_tau, nudging_tauz_z1_uv,            &
   !    nudging_tauz_z2_uv, tauz_tau_uv
!-- Create dummy namelist, so that we can figure out error code for a missing namelist.
!-- This lets us differentiate between missing namelists and those with an error within
!-- the namelist.
    NAMELIST /BNCUIODSBJCB/ place_holder
!
!-- Read MICRO_LAGRANGE namelist from prm file
#ifdef SAMFABIAN
!B. Ascher, 05.05.2025: permit use of user-chosen prm file names
   open(55,file='./'//trim(case)//'/'//trim(prmfile), status='old',form='formatted')
#else
   open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')
#endif /*SAMFABIAN*/
!
!-- Get error code for missing namelist (by giving the name for a namelist that
!-- doesn't exist in the prm file).
    READ( UNIT=55, NML=BNCUIODSBJCB, IOSTAT=ios_missing_namelist )
    REWIND( UNIT=55 )

    READ( UNIT=55, NML=MICRO_LAGRANGE, IOSTAT=ios )
    print*, "Missing namelist code:", ios_missing_namelist, "Actual namelist error flag: ", ios
!
!-- Check namelist for errors
    IF ( ios /= 0 ) THEN
       IF ( ios /= ios_missing_namelist ) THEN
          PRINT*, '****** LCM: bad specification in MICRO_LAGRANGE namelist ******'
          BACKSPACE(55)
          READ(55,fmt='(A)') line
          WRITE(*,'(A)') &
          "Invalid line in LCM namelist _prm_ : " //TRIM(line)
          CALL task_abort()
       ELSE IF ( masterproc ) THEN
          PRINT*, '*****************************************************'
          PRINT*, '******* no MICRO_LAGRANGE namelist in prm file ******'
          PRINT*, '*****************************************************'
       ENDIF
    ENDIF
      
    CLOSE(UNIT=55)

    ! avoid trying to double-allocate when restart=2
    IF (ALLOCATED(flag_wmass)) THEN
      DEALLOCATE( flag_wmass, flag_number, flag_precip, flag_micro3Dout, fluxbmk, fluxtmk, &
#ifdef PICHAMBER
       fluxlmk, fluxrmk, fluxqmk, fluxhmk, &
#endif /*PICHAMBER*/
                  mkwle, mkwsb, mkadv, mkdiff, mklsadv, mkname, mklongname, mkunits, mkoutputscale )
    END IF

!
!-- Allocate all necessary arrays
    ALLOCATE( &
!
!--    Flags
       flag_wmass( nmicro_fields ),                         &
       flag_number( nmicro_fields+nmicro_fields_diag ),     &
       flag_precip( nmicro_fields ),                        &
       flag_micro3Dout( nmicro_fields+nmicro_fields_diag ), &
!
!--    Fluxes at the botoom and top boundary
       fluxbmk( nx, ny, nmicro_fields ), &
       fluxtmk( nx, ny, nmicro_fields ), &
#ifdef PICHAMBER
       fluxlmk( ny, nzm, nmicro_fields ), fluxrmk( ny, nzm, nmicro_fields ), &
       fluxqmk( nx, nzm, nmicro_fields ), fluxhmk( nx, nzm, nmicro_fields ), &
#endif /*PICHAMBER*/
!
!--    Tendencies for statistics output
       mkwle( nz, nmicro_fields+nmicro_fields_diag ),   &
       mkwsb( nz, nmicro_fields+nmicro_fields_diag ),   &
       mkadv( nz, nmicro_fields+nmicro_fields_diag ),   &
       mkdiff( nz, nmicro_fields+nmicro_fields_diag ),  &
       mklsadv( nz, nmicro_fields+nmicro_fields_diag ), &
!
!--    Output information
       mkname( nmicro_fields+nmicro_fields_diag ),        &
       mklongname( nmicro_fields+nmicro_fields_diag ),    &
       mkunits( nmicro_fields+nmicro_fields_diag ),       &
       mkoutputscale( nmicro_fields+nmicro_fields_diag ), &
       STAT=ierr )

    IF ( ierr /= 0 ) THEN
       WRITE(*,*) 'LCM: Failed to allocate microphysical arrays'
       CALL task_abort()
    ENDIF

    ! avoid trying to double-allocate when restart=2
    ! also, don't deallocatem because for nrestart=2, the last micro_field is read 
    ! from the restart file before here
    IF (.NOT. ALLOCATED(micro_field)) THEN
      ALLOCATE( &
         micro_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nz_LCM,nmicro_fields+nmicro_fields_diag), &
         STAT=ierr )
    END IF

    IF ( ierr /= 0 ) THEN
       WRITE(*,*) 'LCM: Failed to allocate microphysical arrays'
       CALL task_abort()
    ENDIF
    IF (.NOT.(ANY(habitparamlist == habitparam))) THEN
      WRITE(*,*) "Ice habit parameterization must be one of 'Chen', 'Weiss', or 'Harrington'"
      CALL task_abort()
    ENDIF
!
!-- Turn on calculation of effective radii in the LCM
    compute_reffc = .TRUE.
    compute_reffi = .TRUE.
!
!-- Which PROGNOSTIC micro_field quantity represents water mass, a number
!-- concentration, or precipitation?
    flag_wmass(:)   = 0
    flag_wmass(imv) = 1 ! water vapor
    flag_number(:)  = 0 ! could be extended to diagnostic quantities ... but I prefer mixing ratios
    flag_precip(:)  = 0
!
!-- Set output quantities (automatic output of QV, QN, and QP by SAM)
    flag_micro3Dout(:) = 0

    flag_micro3Dout( iss )   = 1
    flag_micro3Dout( ina )   = 1
    flag_micro3Dout( iqa )   = 1
    flag_micro3Dout( inc )   = 1
    flag_micro3Dout( iqc )   = 1
    flag_micro3Dout( inr )   = 1
    flag_micro3Dout( iqr )   = 1
    flag_micro3Dout( ipl )   = 1
    IF ( switch_ice )  THEN
       flag_micro3Dout( iiss )  = 1
       flag_micro3Dout( ini )   = 1
       flag_micro3Dout( iqi )   = 1
       flag_micro3Dout( ins )   = 1
       flag_micro3Dout( iqs )   = 1
       flag_micro3Dout( ing )   = 1
       flag_micro3Dout( iqg )   = 1
       flag_micro3Dout( ininc ) = 1
       flag_micro3Dout( ipi )   = 1
    ENDIF

    flag_micro3Dout( iRsd )  = 1

    flag_micro3Dout( isLES ) = 1
    IF ( switch_LEM_micro )  THEN
       flag_micro3Dout( issp )  = 1
       flag_micro3Dout( imsp )  = 1
    ENDIF

    flag_micro3Dout( idiss ) = 1

    flag_micro3Dout( itmax ) = 1
    flag_micro3Dout( itmean ) = 1
    flag_micro3Dout( iltmax) = 1
    flag_micro3Dout( iltmean) = 1
    IF (switch_ice) THEN
      flag_micro3Dout( iitmax) = 1
      flag_micro3Dout( iitmean) = 1
    ENDIF
    IF (switch_ice) flag_micro3Dout( iiref ) = 1

    IF ( switch_cloud_tracking )  THEN
       flag_micro3Dout( iSCLR )      = 1
       flag_micro3Dout( iPlume )     = 1
       flag_micro3Dout( iCondensed ) = 1
       flag_micro3Dout( iCore )      = 1
    ENDIF

    IF (switch_particle_source) THEN
      flag_micro3Dout( iasp )      = 1
      flag_micro3Dout( icsp )      = 1
    ENDIF

    flag_micro3Dout( iact ) = 1
    flag_micro3Dout( idea ) = 1

    IF (switch_ice) THEN
      flag_micro3Dout( itfrz) = 1
      flag_micro3Dout( itfrzc) = 1
      flag_micro3Dout( itfrzr) = 1
      flag_micro3Dout( itfrzi) = 1
      flag_micro3Dout( itfrzs) = 1
      flag_micro3Dout( itfrzg) = 1
      flag_micro3Dout( issa ) = 1
      flag_micro3Dout( idsa ) = 1
      flag_micro3Dout( ibsa ) = 1
      flag_micro3Dout( ipsa ) = 1

      flag_micro3Dout( iimm ) = 1
         !flag_micro3Dout( ihzn ) = 1 
         !Disable haze/deposition output for now, until we implement that back into our IN scheme - B. Ascher 26.05.2025
      flag_micro3Dout( inuc ) = 1
      flag_micro3Dout( ihom ) = 1

      flag_micro3Dout( irim ) = 1
      flag_micro3Dout( ifrz ) = 1
      flag_micro3Dout( imel ) = 1

    ENDIF
      flag_micro3Dout( icon ) = 1
      flag_micro3Dout( ievp ) = 1
    IF (switch_ice) THEN
      flag_micro3Dout( idep ) = 1
      flag_micro3Dout( isub ) = 1

      IF (switch_habit_ice) THEN
         flag_micro3Dout( iphi ) = 1
         flag_micro3Dout( iaavg ) = 1
         flag_micro3Dout( icavg ) = 1
      ENDIF
      flag_micro3Dout( irmf ) = 1
      flag_micro3Dout( irhoapp) = 1
    ENDIF

    IF (switch_collection) THEN
      flag_micro3Dout( icoal ) = 1
      flag_micro3Dout( ibreak ) = 1
      IF (switch_ice) THEN
         flag_micro3Dout( iagg ) = 1
         flag_micro3Dout( iicebreak ) = 1
         flag_micro3Dout( ids ) = 1
         flag_micro3Dout( irs ) = 1
      ENDIF
    ENDIF

    flag_micro3Dout( isr )  = 1
    flag_micro3Dout( ism )  = 1

#ifdef SAMCHEM
    if (dochemistry) flag_micro3Dout( ipH )       = 1
    if (dochemistry) flag_micro3Dout( idelR )     = 1
    if (dochemistry) flag_micro3Dout( iSbyO3t )   = 1
    if (dochemistry) flag_micro3Dout( iSbyO3m )   = 1
    if (dochemistry) flag_micro3Dout( iSbyH2O2t ) = 1
    if (dochemistry) flag_micro3Dout( iSbyH2O2m ) = 1
#endif /*SAMCHEM*/

    ! check for errors in aerosol specification 
    ! determine active modes (non-zero modes)
    DO iMode = 1, 4
      IF (n_aero(iMode)==0.0) THEN
        ! no shares for inactive mode
        shares_aero(iMode) = 0.0
      ELSE
        ! check for values that need to be given (rm_aero, sigma_aero)
        IF (rm_aero(iMode)==0.0) THEN
            WRITE(*,*) 'Error in MICRO_LAGRANGE namelist. Bad aerosol specification: please give mean radii for all specified n_aero.'
            CALL task_abort
        END IF
        IF (sigma_aero(iMode)==0.0) THEN
            WRITE(*,*) 'Error in MICRO_LAGRANGE namelist. Bad aerosol specification: please give standard deviation for all specified n_aero.'
            CALL task_abort
        END IF
        ! check for values that don't need to be given (rho_aero, molecular_weight_aero, vanthoff_aero)
        ! if not given, use first value in array for all modes to either apply user-given or default values
        IF (rho_aero(iMode)==0.0) THEN
          rho_aero(iMode) = rho_aero(1)
        END IF
        IF (molecular_weight_aero(iMode)==0.0) THEN
          molecular_weight_aero(iMode) = molecular_weight_aero(1)
        END IF
        IF (vanthoff_aero(iMode)==0.0) THEN
          vanthoff_aero(iMode) = vanthoff_aero(1)
        END IF
      END IF
    END DO

    ! move all non-zero aerosol modes to the start of the array
    ! (just to make the constant bfactor work in case only one mode which is not the first is specified)
    DO i = 1, 3
      DO iMode = 1, 3
        IF (shares_aero(iMode)==0.0 .AND. shares_aero(iMode+1)/=0.0) THEN
          shares_aero(iMode)           = shares_aero(iMode+1)
          n_aero(iMode)                = n_aero(iMode+1)
          rm_aero(iMode)               = rm_aero(iMode+1)
          sigma_aero(iMode)            = sigma_aero(iMode+1)
          rho_aero(iMode)              = rho_aero(iMode+1)
          molecular_weight_aero(iMode) = molecular_weight_aero(iMode+1)
          vanthoff_aero(iMode)         = vanthoff_aero(iMode+1)

          shares_aero(iMode+1)           = 0
          n_aero(iMode+1)                = 0.0
          rm_aero(iMode+1)               = 0.0
          sigma_aero(iMode+1)            = 0.0
          rho_aero(iMode+1)              = 0.0
          molecular_weight_aero(iMode+1) = 0.0
          vanthoff_aero(iMode+1)         = 0.0
        END IF
      END DO
    END DO

    ! catch error or n_aero=0.0 everywhere
    IF (shares_aero(1)==0.0) THEN
        WRITE(*,*) 'No aerosol mode specified correctly. Please check aerosol parameters.'
        call task_abort
    END IF

    ! save how many aerosol modes there are and explicitly if there is only one to avoid loop
    DO iMode = 2, 4
      IF (shares_aero(iMode) == 0.0) THEN
        IF (iMode==2) only_one_aerosol_mode = .TRUE.
        n_aerosol_modes = iMode-1
        EXIT
      END IF
    END DO

    ! calculate once the constant, input-dependent bfactor coefficients
    DO iMode = 1, n_aerosol_modes
      bfactor_coefficient(iMode) = (vanthoff_aero(iMode)/molecular_weight_aero(iMode)) * molecular_weight_water / (rho_liq*4*pi/3)
    END DO

 END SUBROUTINE micro_setparm

 SUBROUTINE micro_init()
!
!-- Initialize microphysics: this one is guaranteed to be called by SAM at the
!-- beginning of each run, initial or restart:
    USE micro_vars, ONLY: r_separate_micro, rho_liq, spec_liq, spec_ice, spec_aero
    USE vars,       ONLY: qsatw, tabs0, q0, nrestart, pres
    USE params,     ONLY: docloud, pi
    USE grid,       ONLY: masterproc, rank, nrestart, caseid, case, &
                          nsubdomains_x, nsubdomains_y, dx, dy, dz, dt, zi
#ifdef SAMCHEM
    USE micro_vars, ONLY: spec_pH
    USE params,     ONLY: dochemistry
    USE chemistry,  ONLY: aqspc_molmasses, aqspc_names
#endif /*SAMCHEM*/


    IMPLICIT NONE

    INTEGER                 :: k, stat, n
    CHARACTER (LEN=7), SAVE :: rank_char
    CHARACTER (LEN=8), SAVE :: char_n_bins_spec
    CHARACTER(1000) :: aq_names

    REAL :: xl, xr, xm
!
!-- Output of droplet spectra using variable bins
    x_end_spec   = 4.0 / 3.0 * pi * rho_liq * ( r_end_spec )**3
    x_start_spec = 4.0 / 3.0 * pi * rho_liq * ( r_start_spec )**3
    n_bins_spec  = CEILING( LOG( x_end_spec / x_start_spec ) / LOG( 2.0**alpha_spec ) )

    ALLOCATE( rl_spec(n_bins_spec), rr_spec(n_bins_spec), rm_spec(n_bins_spec))
#ifdef SAMCHEM
    ALLOCATE( pHl_spec(n_bins_spec), pHr_spec(n_bins_spec), pHm_spec(n_bins_spec)  )
#endif /*SAMCHEM*/

    ALLOCATE( spec_global(nzl,n_bins_spec) )
    IF ( horiz_average_spec )  THEN
      IF ( switch_ice )  THEN
         ALLOCATE( spec_ice(1,1,nzl,n_bins_spec) )
      ENDIF
      ALLOCATE( spec_liq(1,1,nzl,n_bins_spec) )
      IF (switch_spec_aero) ALLOCATE( spec_aero(1,1,nzl,n_bins_spec) )
#ifdef SAMCHEM
      IF (switch_spec_pH)   ALLOCATE( spec_pH(1,1,nzl,n_bins_spec) )
#endif /*SAMCHEM*/
    ELSE
      IF ( switch_ice )  THEN
         ALLOCATE( spec_ice(nx_LCM*nsubdomains_x,ny_LCM*nsubdomains_y,nzl,n_bins_spec) )
      ENDIF
      ALLOCATE( spec_liq(nx_LCM*nsubdomains_x,ny_LCM*nsubdomains_y,nzl,n_bins_spec) )
      IF (switch_spec_aero) ALLOCATE( spec_aero(nx_LCM*nsubdomains_x,ny_LCM*nsubdomains_y,nzl,n_bins_spec) )
#ifdef SAMCHEM
      IF (switch_spec_pH)   ALLOCATE( spec_pH(nx_LCM*nsubdomains_x,ny_LCM*nsubdomains_y,nzl,n_bins_spec) )
#endif /*SAMCHEM*/
    ENDIF

    DO n = 1, n_bins_spec
       xl = x_start_spec * (2.0**alpha_spec)**(n-1)
       xr = x_start_spec * (2.0**alpha_spec)**n
       xm = SQRT( xl * xr )
       rm_spec(n) = ( xm / ( 4.0 / 3.0 * pi * rho_liq ) )**(1.0/3.0)
       rl_spec(n) = ( xl / ( 4.0 / 3.0 * pi * rho_liq ) )**(1.0/3.0)
       rr_spec(n) = ( xr / ( 4.0 / 3.0 * pi * rho_liq ) )**(1.0/3.0)
    ENDDO
#ifdef SAMCHEM
    DO n = 1, n_bins_spec
       pHm_spec(n) = pH_start_spec + (n-0.5) * (pH_end_spec-pH_start_spec)/n_bins_spec
       pHl_spec(n) = pH_start_spec + (n-1.0) * (pH_end_spec-pH_start_spec)/n_bins_spec
       pHr_spec(n) = pH_start_spec + n       * (pH_end_spec-pH_start_spec)/n_bins_spec
    ENDDO
#endif /*SAMCHEM*/

!
!-- Convert rank number into string
    WRITE(rank_char,'(I4)') rank
!
!-- Convert number of aqueous species into string
    format_particle_tracking = '(I9,1X,I16,1X,3(I5,1X),I13,1X,2(E9.3,1X)'
#ifdef SAMCHEM
    IF ( dochemistry ) THEN
      WRITE(char_n_bins_spec,'(I8)') naqspc
      char_n_bins_spec = ADJUSTL(char_n_bins_spec) ! remove leading blanks
      format_particle_tracking = TRIM(format_particle_tracking)//",1X,"//TRIM(char_n_bins_spec)//'(E9.3,1X))'
    ELSE
#endif /*SAMCHEM*/
      format_particle_tracking = TRIM(format_particle_tracking)//")"
#ifdef SAMCHEM
    END IF
#endif /*SAMCHEM*/

!
!-- Convert number of bins into string
    WRITE(char_n_bins_spec,'(I8)') n_bins_spec
    char_n_bins_spec = ADJUSTL(char_n_bins_spec) ! remove leading blanks
    ! format_spec = '(E12.5,1X,E12.5,1X,'//TRIM(char_n_bins_spec)//'(E12.5,1X))'
    ! format_spec = '(E9.3,1X,E9.3,1X,'//TRIM(char_n_bins_spec)//'(E9.3,1X))'
    format_spec = '(I5,1X,I5,1X,I5,1X,I9,1X,'//TRIM(char_n_bins_spec)//'(E9.3,1X))'

    IF ( masterproc )  THEN
       IF ( switch_ice )  THEN
          IF ( nrestart .NE. 1 )  THEN
            OPEN(UNIT=98,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_ice_header',STATUS='REPLACE',FORM='FORMATTED')

             WRITE(98,*) "time;", "horiz. x;", "horiz. y;", "height;", "geometric mean bin radius;", "left boundary of bin;", "right boundary of bin;", "number concentration per bin;"
             WRITE(98,*) "TIME;", "X;", "Y;", "Z;", "RM;", "RL;", "RR;", "N(RM);"
             WRITE(98,*) "s ", "m ", "m ", "m ", "m ", "m ", "m ", "1/m3/m "
             WRITE(98,*) dx, dy, dz, dt

             WRITE(98,*) n_bins_spec
             IF ( horiz_average_spec )  THEN
                WRITE(98,*) 1
                WRITE(98,*) 1
                WRITE(98,*) nzl
             ELSE
                WRITE(98,*) nx_LCM*nsubdomains_x
                WRITE(98,*) ny_LCM*nsubdomains_y
                WRITE(98,*) nzl
             ENDIF

             WRITE(98,format_spec) -1, -1, -1, -1, ( rm_spec(n), n = 1, n_bins_spec )
             WRITE(98,format_spec) -1, -1, -1, -1, ( rl_spec(n), n = 1, n_bins_spec )
             WRITE(98,format_spec) -1, -1, -1, -1, ( rr_spec(n), n = 1, n_bins_spec )
          ELSE
             OPEN(UNIT=98,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_ice_header',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND',IOSTAT=stat)
             IF ( stat .NE. 0 )  THEN
                OPEN(UNIT=98,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_ice_header',STATUS='NEW',FORM='FORMATTED')

                WRITE(98,*) "time;", "horiz. x;", "horiz. y;", "height;", "geometric mean bin radius;", "left boundary of bin;", "right boundary of bin;", "number concentration per bin;"
                WRITE(98,*) "TIME;", "X;", "Y;", "Z;", "RM;", "RL;", "RR;", "N(RM);"
                WRITE(98,*) "s ", "m ", "m ", "m ", "m ", "m ", "m ", "1/m3/m "

                WRITE(98,*) n_bins_spec
                IF ( horiz_average_spec )  THEN
                   WRITE(98,*) 1
                   WRITE(98,*) 1
                   WRITE(98,*) nzl
                ELSE
                   WRITE(98,*) nx_LCM*nsubdomains_x
                   WRITE(98,*) ny_LCM*nsubdomains_y
                   WRITE(98,*) nzl
                ENDIF

                WRITE(98,format_spec) -1, -1, -1, -1, ( rm_spec(n), n = 1, n_bins_spec )
                WRITE(98,format_spec) -1, -1, -1, -1, ( rl_spec(n), n = 1, n_bins_spec )
                WRITE(98,format_spec) -1, -1, -1, -1, ( rr_spec(n), n = 1, n_bins_spec )
             ENDIF
          ENDIF
          CLOSE(98)
       ENDIF
#ifdef SAMCHEM
       IF ( switch_spec_pH ) THEN
           IF ( nrestart .NE. 1 )  THEN
              OPEN(UNIT=96,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_pH_header',STATUS='REPLACE',FORM='FORMATTED')

              WRITE(96,*) "time;", "horiz. x;", "horiz. y;", "height;", "mean bin pH;", "left boundary of bin;", "right boundary of bin;", "number concentration per bin;"
              WRITE(96,*) "TIME;", "X;", "Y;", "Z;", "pHM;", "pHL;", "pHR;", "N(pHM);"
              WRITE(96,*) "s ", "m ", "m ", "m ", "m ", "m ", "m ", "1/m3/m "
              WRITE(96,*) dx, dy, dz, dt

              WRITE(96,*) n_bins_spec
              IF ( horiz_average_spec )  THEN
                 WRITE(96,*) 1
                 WRITE(96,*) 1
                 WRITE(96,*) nzl
              ELSE
                 WRITE(96,*) nx_LCM*nsubdomains_x
                 WRITE(96,*) ny_LCM*nsubdomains_y
                 WRITE(96,*) nzl
              ENDIF

              WRITE(96,format_spec) -1, -1, -1, -1, ( pHm_spec(n), n = 1, n_bins_spec )
              WRITE(96,format_spec) -1, -1, -1, -1, ( pHl_spec(n), n = 1, n_bins_spec )
              WRITE(96,format_spec) -1, -1, -1, -1, ( pHr_spec(n), n = 1, n_bins_spec )
           ELSE
              OPEN(UNIT=96,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_pH_header',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND',IOSTAT=stat)
              IF ( stat .NE. 0 )  THEN
                 OPEN(UNIT=96,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_pH_header',STATUS='NEW',FORM='FORMATTED')

                 WRITE(96,*) "time;", "horiz. x;", "horiz. y;", "height;", "mean bin pH;", "left boundary of bin;", "right boundary of bin;", "litres per bin;"
                 WRITE(96,*) "TIME;", "X;", "Y;", "Z;", "pHM;", "pHL;", "pHR;", "Litres(pHM);"
                 WRITE(96,*) "s ", "m ", "m ", "m ", "m ", "m ", "m ", "l/m3/pH "
                 WRITE(96,*) dx, dy, dz, dt

                 WRITE(96,*) n_bins_spec
                 IF ( horiz_average_spec )  THEN
                    WRITE(96,*) 1
                    WRITE(96,*) 1
                    WRITE(96,*) nzl
                 ELSE
                    WRITE(96,*) nx_LCM*nsubdomains_x
                    WRITE(96,*) ny_LCM*nsubdomains_y
                    WRITE(96,*) nzl
                 ENDIF

                 WRITE(96,format_spec) -1, -1, -1, -1, ( pHm_spec(n), n = 1, n_bins_spec )
                 WRITE(96,format_spec) -1, -1, -1, -1, ( pHl_spec(n), n = 1, n_bins_spec )
                 WRITE(96,format_spec) -1, -1, -1, -1, ( pHr_spec(n), n = 1, n_bins_spec )
              ENDIF
           ENDIF
           CLOSE(96)
       END IF
#endif /*SAMCHEM*/
       IF (switch_spec_aero) THEN
           IF ( nrestart .NE. 1 )  THEN
              OPEN(UNIT=97,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_aero_header',STATUS='REPLACE',FORM='FORMATTED')

              WRITE(97,*) "time;", "horiz. x;", "horiz. y;", "height;", "geometric mean bin radius;", "left boundary of bin;", "right boundary of bin;", "number concentration per bin;"
              WRITE(97,*) "TIME;", "X;", "Y;", "Z;", "RM;", "RL;", "RR;", "N(RM);"
              WRITE(97,*) "s ", "m ", "m ", "m ", "m ", "m ", "m ", "1/m3/m "
              WRITE(97,*) dx, dy, dz, dt

              WRITE(97,*) n_bins_spec
              IF ( horiz_average_spec )  THEN
                 WRITE(97,*) 1
                 WRITE(97,*) 1
                 WRITE(97,*) nzl
              ELSE
                 WRITE(97,*) nx_LCM*nsubdomains_x
                 WRITE(97,*) ny_LCM*nsubdomains_y
                 WRITE(97,*) nzl
              ENDIF

              WRITE(97,format_spec) -1, -1, -1, -1, ( rm_spec(n), n = 1, n_bins_spec )
              WRITE(97,format_spec) -1, -1, -1, -1, ( rl_spec(n), n = 1, n_bins_spec )
              WRITE(97,format_spec) -1, -1, -1, -1, ( rr_spec(n), n = 1, n_bins_spec )
           ELSE
              OPEN(UNIT=97,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_aero_header',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND',IOSTAT=stat)
              IF ( stat .NE. 0 )  THEN
                 OPEN(UNIT=97,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_aero_header',STATUS='NEW',FORM='FORMATTED')

                 WRITE(97,*) "time;", "horiz. x;", "horiz. y;", "height;", "geometric mean bin radius;", "left boundary of bin;", "right boundary of bin;", "number concentration per bin;"
                 WRITE(97,*) "TIME;", "X;", "Y;", "Z;", "RM;", "RL;", "RR;", "N(RM);"
                 WRITE(97,*) "s ", "m ", "m ", "m ", "m ", "m ", "m ", "1/m3/m "
                 WRITE(97,*) dx, dy, dz, dt

                 WRITE(97,*) n_bins_spec
                 IF ( horiz_average_spec )  THEN
                    WRITE(97,*) 1
                    WRITE(97,*) 1
                    WRITE(97,*) nzl
                 ELSE
                    WRITE(97,*) nx_LCM*nsubdomains_x
                    WRITE(97,*) ny_LCM*nsubdomains_y
                    WRITE(97,*) nzl
                 ENDIF

                 WRITE(97,format_spec) -1, -1, -1, -1, ( rm_spec(n), n = 1, n_bins_spec )
                 WRITE(97,format_spec) -1, -1, -1, -1, ( rl_spec(n), n = 1, n_bins_spec )
                 WRITE(97,format_spec) -1, -1, -1, -1, ( rr_spec(n), n = 1, n_bins_spec )
              ENDIF
           ENDIF
           CLOSE(97)
       ENDIF

       IF ( nrestart .NE. 1 )  THEN
          OPEN(UNIT=99,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_liq_header',STATUS='REPLACE',FORM='FORMATTED')

          WRITE(99,*) "time;", "horiz. x;", "horiz. y;", "height;", "geometric mean bin radius;", "left boundary of bin;", "right boundary of bin;", "number concentration per bin;"
          WRITE(99,*) "TIME;", "X;", "Y;", "Z;", "RM;", "RL;", "RR;", "N(RM);"
          WRITE(99,*) "s ", "m ", "m ", "m ", "m ", "m ", "m ", "1/m3/m "
          WRITE(99,*) dx, dy, dz, dt

          WRITE(99,*) n_bins_spec
          IF ( horiz_average_spec )  THEN
             WRITE(99,*) 1
             WRITE(99,*) 1
             WRITE(99,*) nzl
          ELSE
             WRITE(99,*) nx_LCM*nsubdomains_x
             WRITE(99,*) ny_LCM*nsubdomains_y
             WRITE(99,*) nzl
          ENDIF

          WRITE(99,format_spec) -1, -1, -1, -1, ( rm_spec(n), n = 1, n_bins_spec )
          WRITE(99,format_spec) -1, -1, -1, -1, ( rl_spec(n), n = 1, n_bins_spec )
          WRITE(99,format_spec) -1, -1, -1, -1, ( rr_spec(n), n = 1, n_bins_spec )
       ELSE
          OPEN(UNIT=99,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_liq_header',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND',IOSTAT=stat)
          IF ( stat .NE. 0 )  THEN
             OPEN(UNIT=99,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_liq_header',STATUS='NEW',FORM='FORMATTED')

             WRITE(99,*) "time;", "horiz. x;", "horiz. y;", "height;", "geometric mean bin radius;", "left boundary of bin;", "right boundary of bin;", "number concentration per bin;"
             WRITE(99,*) "TIME;", "X;", "Y;", "Z;", "RM;", "RL;", "RR;", "N(RM);"
             WRITE(99,*) "s ", "m ", "m ", "m ", "m ", "m ", "m ", "1/m3/m "

             WRITE(99,*) n_bins_spec
             IF ( horiz_average_spec )  THEN
                WRITE(99,*) 1
                WRITE(99,*) 1
                WRITE(99,*) nzl
             ELSE
                WRITE(99,*) nx_LCM*nsubdomains_x
                WRITE(99,*) ny_LCM*nsubdomains_y
                WRITE(99,*) nzl
             ENDIF

             WRITE(99,format_spec) -1, -1, -1, -1, ( rm_spec(n), n = 1, n_bins_spec )
             WRITE(99,format_spec) -1, -1, -1, -1, ( rl_spec(n), n = 1, n_bins_spec )
             WRITE(99,format_spec) -1, -1, -1, -1, ( rr_spec(n), n = 1, n_bins_spec )
          ENDIF
       ENDIF
       CLOSE(99)
    ENDIF
    IF (track_some_particles) THEN
       IF ( nrestart .NE. 1 )  THEN
         OPEN(UNIT=126,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.tracked_particles_'//TRIM(ADJUSTL(rank_char)),STATUS='REPLACE',FORM='FORMATTED',POSITION='APPEND')
         WRITE(126,*) '-------------------------------------------'
         WRITE(126,*) 'Data-file for tracked Lagrangian particles.'
         WRITE(126,*) '-------------------------------------------'
         WRITE(126,*) ''
         WRITE(126,'(A,I4)') 'FROM RANK (processor ID): ', rank
         WRITE(126,*) ''
         WRITE(126,*) '------------------------'
         WRITE(126,*) 'Meta info BEGIN'
         WRITE(126,*) 'rho_aero (kg/m3)'
         WRITE(126,*) rho_aero
#ifdef SAMCHEM
         IF (dochemistry) THEN
           WRITE(126,*) 'aqspc_molmasses (g/mol)'
           WRITE(126,*) aqspc_molmasses
         END IF
#endif /*SAMCHEM*/
         WRITE(126,*) 'Meta info END'
         WRITE(126,*) '------------------------'
         WRITE(126,*) ''
         WRITE(126,*) 'Data is stored as follows:'
         WRITE(126,*) 'nstep, ID, i, j, k, multiplicity (iA) (#), water mass (iM) (kg), aerosol mass (iN) (kg)'
#ifdef SAMCHEM
         IF (dochemistry) THEN
           aq_names = ""
           DO k = 1, naqspc
             aq_names = TRIM(ADJUSTL(aq_names))//TRIM(aqspc_names(k))//","
           END DO
           WRITE(126,'(A,A)') 'and in molec: ', TRIM(aq_names)
         ELSE
         END IF
#endif /*SAMCHEM*/
         WRITE(126,*) ''
         WRITE(126,*) ''
         WRITE(126,*) 'Start tracking ....'
         WRITE(126,*) ''
         WRITE(126,*) ''
         WRITE(126,*) ''
       ELSE
         OPEN(UNIT=126,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.tracked_particles_'//TRIM(ADJUSTL(rank_char)),STATUS='OLD',FORM='FORMATTED',POSITION='APPEND')
       END IF
      !CLOSE(126) ! never close like in horiz_average_spec
     END IF



    IF ( horiz_average_spec )  THEN
       IF ( switch_ice )  THEN
          IF ( nrestart .NE. 1 )  THEN
             OPEN(UNIT=98,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_ice',STATUS='REPLACE',FORM='FORMATTED',POSITION='APPEND',IOSTAT=stat)
          ELSE
             OPEN(UNIT=98,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_ice',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND',IOSTAT=stat)
          ENDIF
          CLOSE(98)
       ENDIF
#ifdef SAMCHEM
       IF ( switch_spec_pH ) THEN
          IF ( nrestart .NE. 1 )  THEN
             OPEN(UNIT=96,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_pH',STATUS='REPLACE',FORM='FORMATTED',POSITION='APPEND',IOSTAT=stat)
          ELSE
             OPEN(UNIT=96,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_pH',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND',IOSTAT=stat)
          ENDIF
          CLOSE(96)
       END IF
#endif /*SAMCHEM*/
       IF ( switch_spec_aero ) THEN
          IF ( nrestart .NE. 1 )  THEN
             OPEN(UNIT=97,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_aero',STATUS='REPLACE',FORM='FORMATTED',POSITION='APPEND',IOSTAT=stat)
          ELSE
             OPEN(UNIT=97,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_aero',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND',IOSTAT=stat)
          ENDIF
          CLOSE(97)
       END IF
       IF ( nrestart .NE. 1 )  THEN
          OPEN(UNIT=99,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_liq',STATUS='REPLACE',FORM='FORMATTED',POSITION='APPEND',IOSTAT=stat)
       ELSE
          OPEN(UNIT=99,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_liq',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND',IOSTAT=stat)
       ENDIF
       CLOSE(99)
    ENDIF
!
!-- Initialize user output quantities
    CALL micro_user_init
!
!-- First, setup the underlying Lagrangian particle model
    CALL micro_set
!
!-- Second, setup micro_fields and surface/top fluxes (fluxbmk, fluxtmk)
    IF ( nrestart .eq. 0 )  THEN
       micro_field = 0.0
       DO k = 1, nz_LCM
          micro_field(:,:,k,imv) = q0(k)
          micro_field(:,:,k,iss) = q0(k) - qsatw( tabs0(k), pres(k) )
          micro_field(:,:,k,ilastT) = tabs0(k)
          micro_field(:,:,k,ilastp) = pres(k)

       ENDDO
       fluxbmk = 0.0
       fluxtmk = 0.0
    ELSE
       fluxbmk = 0.0
       fluxtmk = 0.0
    ENDIF
#ifdef PICHAMBER
    fluxlmk = 0.0
    fluxrmk = 0.0
    fluxqmk = 0.0
    fluxhmk = 0.0
#endif /*PICHAMBER*/

    fluxbmk(:,:,iSCLR) = 1.0 ! surface flux of plume-detecting scalar is always one
!
!-- Third, setup output quantities
    mkwle  = 0.0
    mkwsb  = 0.0
    mkadv  = 0.0
    mkdiff = 0.0

    mkname(imv) = 'QQ'
    mklongname(imv) = 'WATER VAPOR'
    mkunits(imv) = 'g/kg'
    mkoutputscale(imv) = 1.0E3

    mkname(iTlcm) = 'S_LCM'
    mklongname(iTlcm) = 'static energy'
    mkunits(iTlcm) = 'K'
    mkoutputscale(iTlcm) = 1.0

    mkname(iTllcm) = 'SL_LCM'
    mklongname(iTllcm) = 'liquid water static energy'
    mkunits(iTllcm) = 'K'
    mkoutputscale(iTllcm) = 1.0



    mkname(iss) = 'ETA'
    mklongname(iss) = 'QV-QSAT, SUPERSATURATION'
    mkunits(iss) = 'g/kg'
    mkoutputscale(iss) = 1.0E3

    mkname(iSCLR) = 'SCLR'
    mklongname(iSCLR) = 'plume-detecting scalar (unitless)'
    mkunits(iSCLR) = '-'
    mkoutputscale(iSCLR) = 1.0

    mkname(ina) = 'NA'
    mklongname(ina) = 'AEROSOL CONCENTRATION'
    mkunits(ina) = '#/mg'
    mkoutputscale(ina) = 1.0E-6

    mkname(iqa) = 'QA'
    mklongname(iqa) = 'AEROSOL WATER MIXING RATIO'
    mkunits(iqa) = 'g/kg'
    mkoutputscale(iqa) = 1.0E3

    mkname(inc) = 'NC'
    mklongname(inc) = 'CLOUD DROPLET CONCENTRATION'
    mkunits(inc) = '#/mg'
    mkoutputscale(inc) = 1.0E-6

    mkname(iqc) = 'QC'
    mklongname(iqc) = 'CLOUD WATER'
    mkunits(iqc) = 'g/kg'
    mkoutputscale(iqc) = 1.0E3

    mkname(inr) = 'NR'
    mklongname(inr) = 'RAIN DROPLET NUMBER MIXING RATIO'
    mkunits(inr) = '#/mg'
    mkoutputscale(inr) = 1.0E-6

    mkname(iqr) = 'QR'
    mklongname(iqr) = 'RAIN WATER'
    mkunits(iqr) = 'g/kg'
    mkoutputscale(iqr) = 1.0E3



    mkname(ini) = 'NI'
    mklongname(ini) = 'ICE NUMBER MIXING RATIO'
    mkunits(ini) = '#/mg'
    mkoutputscale(ini) = 1.0E-6

    mkname(iqi) = 'QI'
    mklongname(iqi) = 'ICE WATER'
    mkunits(iqi) = 'g/kg'
    mkoutputscale(iqi) = 1.0E3

    mkname(ins) = 'NS'
    mklongname(ins) = 'SNOW NUMBER MIXING RATIO'
    mkunits(ins) = '#/mg'
    mkoutputscale(ins) = 1.0E-6

    mkname(iqs) = 'QS'
    mklongname(iqs) = 'SNOW WATER'
    mkunits(iqs) = 'g/kg'
    mkoutputscale(iqs) = 1.0E3

    mkname(ing) = 'NG'
    mklongname(ing) = 'GRAUPEL NUMBER MIXING RATIO'
    mkunits(ing) = '#/mg'
    mkoutputscale(ing) = 1.0E-6

    mkname(iqg) = 'QG'
    mklongname(iqg) = 'GRAUPEL WATER'
    mkunits(iqg) = 'g/kg'
    mkoutputscale(iqg) = 1.0E3

    mkname(ipl) = "P_LIQ"
    mklongname(ipl) = "LIQUID PRECIP. FLUX"
    mkunits(ipl) = 'kg/m^2/s'
    mkoutputscale(ipl) = 1.0

    IF (switch_ice) THEN
      mkname(ipi) = "P_ICE"
      mklongname(ipi) = "LIQUID PRECIP. FLUX"
      mkunits(ipi) = 'kg/m^2/s'
      mkoutputscale(ipi) = 1.0
    ENDIF

    mkname(ininc) = 'NIN'
    mklongname(ininc) = 'ICE NUCLEI MIXING RATIO (part of NA, NC, NR)'
    mkunits(ininc) = '#/mg'
    mkoutputscale(ininc) = 1.0E-6



    mkname(iRsd) = 'RSD'
    mklongname(iRsd) = 'eff. radius SD'
    mkunits(iRsd) = 'um'
    mkoutputscale(iRsd) = 1.0E6

    mkname(iiREF) = "iREF"
    mklongname(iiREF) = "ice effective radius"
    mkunits(iiREF) = "um"
    mkoutputscale(iiREF) = 1.0E6

    mkname(isLES) = 'SLS'
    mklongname(isLES) = 'S LES'
    mkunits(isLES) = '%'
    mkoutputscale(isLES) = 100.0

    mkname(issp) = 'S_SP'
    mklongname(issp) = 'standard deviation of supersaturation perturbation'
    mkunits(issp) = '%'
    mkoutputscale(issp) = 100.0

    mkname(imsp) = 'M_SP'
    mklongname(imsp) = 'mean of supersaturation perturbation'
    mkunits(imsp) = '%'
    mkoutputscale(imsp) = 100.0

    mkname(iiss) = "ETA_ICE"
    mklongname(iiss) = "saturation vapor mixing ratio deficit over ice"
    mkunits(iiss) = 'g/kg'
    mkoutputscale(iiss) = 1.0E3



    mkname(idiss) = 'DISS'
    mklongname(idiss) = 'kinetic energy dissipation rate'
    mkunits(idiss) = 'cm2 s-3'
    mkoutputscale(idiss) = 10000.0



    mkname(itmax) = 'TMX'
    mklongname(itmax) = 'maximum residence time'
    mkunits(itmax) = 's'
    mkoutputscale(itmax) = 1.0

    mkname(itmean) = 'TMN'
    mklongname(itmean) = 'mean residence time'
    mkunits(itmean) = 's'
    mkoutputscale(itmean) = 1.0

    mkname(iltmax) = "TMX_LIQ"
    mklongname(iltmax) = "liquid particle maximum residence time"
    mkunits(iltmax) = 's'
    mkoutputscale(iltmax) = 1.0

    mkname(iltmean) = "TMN_LIQ"
    mklongname(iltmean) = "liquid particle mean residence time"
    mkunits(iltmean) = 's'
    mkoutputscale(iltmean) = 1.0

    mkname(iitmax) = "TMX_ICE"
    mklongname(iitmax) = "ice particle maximum residence time"
    mkunits(iitmax) = 's'
    mkoutputscale(iitmax) = 1.0

    mkname(iitmean) = "TMN_ICE"
    mklongname(iitmean) = "ice particle mean residence time"
    mkunits(iitmean) = 's'
    mkoutputscale(iitmean) = 1.0

    mkname(iPlume) = 'plume'
    mklongname(iPlume) = 'plume grid box for cloud tracking'
    mkunits(iPlume) = '-'
    mkoutputscale(iPlume) = 1.0

    mkname(iCondensed) = 'condensed'
    mklongname(iCondensed) = 'condensed grid box for cloud tracking'
    mkunits(iCondensed) = '-'
    mkoutputscale(iCondensed) = 1.0

    mkname(iCore) = 'core'
    mklongname(iCore) = 'core grid box for cloud tracking'
    mkunits(iCore) = '-'
    mkoutputscale(iCore) = 1.0

    mkname(iasp) = 'ASP'
    mklongname(iasp) = 'Sprayed particle number conentration'
    mkunits(iasp) = '#/mg'
    mkoutputscale(iasp) = 1.0E-6

    mkname(icsp) = 'CSP'
    mklongname(icsp) = 'Cloud droplets from sprayed particle number concentration'
    mkunits(icsp) = '#/mg'
    mkoutputscale(icsp) = 1.0E-6

    mkname(iACT) = 'NACT'
    mklongname(iACT) = 'activation rate'
    mkunits(iACT) = '#/mg/s'
    mkoutputscale(iACT) = 1.0E-6

    mkname(iDEA) = 'NDEA'
    mklongname(iDEA) = 'deactivation rate'
    mkunits(iDEA) = '#/mg/s'
    mkoutputscale(iDEA) = 1.0E-6

    mkname(itfrz) = "TFRZ"
    mklongname(itfrz) = "Average droplet freezing temperature"
    mkunits(itfrz) = "K"
    mkoutputscale(itfrz) = 1.0

    mkname(issa) = "SSA"
    mklongname(issa) = "Average soot IN surface area per particle"
    mkunits(issa) = "um2"
    mkoutputscale(issa) = 1.0E12

    mkname(idsa) = "DSA"
    mklongname(idsa) = "Average dust IN surface area per particle"
    mkunits(idsa) = "um2"
    mkoutputscale(idsa) = 1.0E12

    mkname(ibsa) = "BSA"
    mklongname(ibsa) = "Average Snowmax IN surface area per particle"
    mkunits(ibsa) = "um2"
    mkoutputscale(ibsa) = 1.0E12

    mkname(ipsa) = "PSA"
    mklongname(ipsa) = "Average pollen IN surface area per particle"
    mkunits(ipsa) = "um2"
    mkoutputscale(ipsa) = 1.0E12

    mkname(iimm) = "IMM"
    mklongname(iimm) = "immersion freezing rate"
    mkunits(iimm) = "#/g/s"  !-- Note different units; this is because I'm not sure how it will get handled in floats and I don't want
    mkoutputscale(iimm) = 1.0E-3 !-- these fields to get truncated in the output

    mkname(ihzn) = "HZN"
    mklongname(ihzn) = "deposition nucleation rate"
    mkunits(ihzn) = "#/g/s"
    mkoutputscale(ihzn) = 1.0E-3

    mkname(ihom) = "HOM"
    mklongname(ihom) = "homogeneous nucleation rate"
    mkunits(ihom) = "#/g/s"
    mkoutputscale(ihom) = 1.0E-3

    mkname(inuc) = "NUC"
    mklongname(inuc) = "total ice nucleation rate"
    mkunits(inuc) = "#/g/s"
    mkoutputscale(inuc) = 1.0E-3


    mkname(irim) = "RIMR"
    mklongname(irim) = "total riming rate"
    mkunits(irim) = "mg/kg/s"
    mkoutputscale(irim) = 1.0E6

    mkname(ifrz) = "FRZ"
    mklongname(ifrz) = "total freezing rate"
    mkunits(ifrz) = "mg/kg/s"
    mkoutputscale(ifrz) = 1.0E6

    mkname(imel) = "MEL"
    mklongname(imel) = "total melting rate"
    mkunits(imel) = "mg/kg/s"
    mkoutputscale(imel) = 1.0E6

    mkname(icon) = "CON"
    mklongname(icon) = "net condensation rate"
    mkunits(icon) = "mg/kg/s"
    mkoutputscale(icon) = 1.0E6

    mkname(ievp) = "EVAP"
    mklongname(ievp) = "net evaporation rate"
    mkunits(ievp) = "mg/kg/s"
    mkoutputscale(ievp) = 1.0E6

    mkname(idep) = "DEP"
    mklongname(idep) = "net deposition rate"
    mkunits(idep) = "mg/kg/s"
    mkoutputscale(idep) = 1.0E6

    mkname(isub) = "SUB"
    mklongname(isub) = "net sublimation rate"
    mkunits(isub) = "mg/kg/s"
    mkoutputscale(isub) = 1.0E6


    mkname(iphi) = "PHI"
    mklongname(iphi) = "ice aspect ratio"
    mkunits(iphi) = "unitless"
    mkoutputscale(iphi) = 1.0

    mkname(iaavg) = "AAXIS"
    mklongname(iaavg) = "ice mass-weighted a-axis length"
    mkunits(iaavg) = "um"
    mkoutputscale(iaavg) = 1.0E6

    mkname(icavg) = "CAXIS"
    mklongname(icavg) = "ice mass-weighted c-axis length"
    mkunits(icavg) = "um"
    mkoutputscale(icavg) = 1.0E6

    mkname(irmf) = "RIMEF"
    mklongname(irmf) = "rime mass fraction"
    mkunits(irmf) = "unitless"
    mkoutputscale(irmf) = 1.0

    mkname(irhoapp) = "RHO_APP"
    mklongname(irhoapp) = "ice apparent density"
    mkunits(irhoapp) = "kg m-3"
    mkoutputscale(irhoapp) = 1.0

    mkname(icoal) = "COAL_RATE"
    mklongname(icoal) = "coalescence rate"
    mkunits(icoal) = "mg/kg/s"
    mkoutputscale(icoal) = 1.0E6

    mkname(ibreak) = "BRKUP_RATE"
    mklongname(ibreak) = "breakup rate"
    mkunits(ibreak) = "mg/kg/s"
    mkoutputscale(ibreak) = 1.0E6

    mkname(iagg) = "AGG_RATE"
    mklongname(iagg) = "aggregation rate"
    mkunits(iagg) = "mg/kg/s"
    mkoutputscale(iagg) = 1.0E6

    mkname(iicebreak) = "ICEBRK_RAT"
    mklongname(iicebreak) = "ice breakup rate"
    mkunits(iicebreak) = "1/kg/s"
    mkoutputscale(iicebreak) = 1.0

    mkname(ids) = "DS_RATE"
    mklongname(ids) = "droplet shattering rate"
    mkunits(ids) = "1/kg/s"
    mkoutputscale(ids) = 1.0

    mkname(irs) = "RS_RATE"
    mklongname(irs) = "rime splintering rate"
    mkunits(irs) = "1/kg/s"
    mkoutputscale(irs) = 1.0

    mkname(itfrzc) = "TFZ_C"
    mklongname(itfrzc) = "mass-weighted cloud droplet freezing temperature"
    mkunits(itfrzc) = "K"
    mkoutputscale(itfrzc) = 1.0

    mkname(itfrzr) = "TFZ_R"
    mklongname(itfrzr) = "mass-weighted rain drop freezing temperature"
    mkunits(itfrzr) = "K"
    mkoutputscale(itfrzr) = 1.0

    mkname(itfrzi) = "TFZ_I"
    mklongname(itfrzi) = "mass-weighted cloud ice freezing temperature"
    mkunits(itfrzi) = "K"
    mkoutputscale(itfrzi) = 1.0

    mkname(itfrzs) = "TFZ_S"
    mklongname(itfrzs) = "mass-weighted snow freezing temperature"
    mkunits(itfrzs) = "K"
    mkoutputscale(itfrzs) = 1.0

    mkname(itfrzg) = "TFZ_G"
    mklongname(itfrzg) = "mass-weighted graupel freezing temperature"
    mkunits(itfrzg) = "K"
    mkoutputscale(itfrzg) = 1.0

    mkname(isr) = "AERO_R"
    mklongname(isr) = "mean aerosol radius"
    mkunits(isr) = "nm"
    mkoutputscale(isr) = 1.0E9

    mkname(ism) = "AEROMIX"
    mklongname(ism) = "aerosol mass mixing ratio"
    mkunits(ism) = "mg/kg"
    mkoutputscale(ism) = 1.0E6


#ifdef SAMCHEM
    mkname(ipH) = 'PH'
    mklongname(ipH) = 'volume-mean pH-value of droplets with chemistry on (-log10([Hp]))'
    mkunits(ipH) = '-'
    mkoutputscale(ipH) = 1.0

    mkname(idelR) = "delR"
    mklongname(idelR) = 'absolute growth of aerosol mass in the grid box'
    mkunits(idelR) = 'kg'
    mkoutputscale(idelR) = 1.0

    mkname(iSbyO3t) = "SbyO3t"
    mklongname(iSbyO3t) = 'sum of O3->S(VI) reaction rates in the grid box'
    mkunits(iSbyO3t) = 'mol/m3/s'
    mkoutputscale(iSbyO3t) = 1.0

    mkname(iSbyO3m) = "SbyO3m"
    mklongname(iSbyO3m) = 'mean O3->S(VI) reaction rates in the grid box'
    mkunits(iSbyO3m) = 'mol/l/s'
    mkoutputscale(iSbyO3m) = 1.0

    mkname(iSbyH2O2t) = "SbyH2O2t"
    mklongname(iSbyH2O2t) = 'sum of H2O2->S(VI) reaction rates in the grid box'
    mkunits(iSbyH2O2t) = 'mol/m3/s'
    mkoutputscale(iSbyH2O2t) = 1.0

    mkname(iSbyH2O2m) = "SbyH2O2m"
    mklongname(iSbyH2O2m) = 'mean of H2O2->S(VI) reaction rates in the grid box'
    mkunits(iSbyH2O2m) = 'mol/l/s'
    mkoutputscale(iSbyH2O2m) = 1.0
#endif /*SAMCHEM*/

!
!-- Calculate sources in physical space from sources in grid space
!

    IF (.NOT. switch_phase_sources .OR. (MAXVAL(p_source_aero) .GT. 0.0_8)) THEN !Only bother calculating physical locations of ccn sources if we have them
      z_b_source_aero = REAL(zi(nz_b_source_aero))
      z_t_source_aero = REAL(zi(nz_t_source_aero + 1))
      y_s_source_aero = ( ny_s_source_aero - 1.0 ) * dy_LCM
      y_n_source_aero = ( ny_n_source_aero - 1.0e-20 ) * dy_LCM
      x_w_source_aero = ( nx_w_source_aero - 1.0 ) * dx_LCM
      x_e_source_aero = ( nx_e_source_aero - 1.0e-20 ) * dx_LCM
    ENDIF
    IF (SWITCH_PHASE_SOURCES .AND. (MAXVAL(p_source_in) .GT. 0.0_8)) THEN !Only bother calculating physical locations of in sources if we have them
      z_b_source_in = real(zi(nz_b_source_in))
      z_t_source_in = real(zi(nz_t_source_in + 1))
      y_s_source_in = ( ny_s_source_in - 1.0 ) * dy_LCM
      y_n_source_in = ( ny_n_source_in - 1.0e-20 ) * dy_LCM
      x_w_source_in = ( nx_w_source_in - 1.0 ) * dx_LCM
      x_e_source_in = ( nx_e_source_in - 1.0e-20 ) * dx_LCM
    ENDIF
!
!-- Finally, initial diagnose of microphysical quantities
    IF ( docloud )  CALL micro_diagnose()

 END SUBROUTINE micro_init

 SUBROUTINE micro_flux()
!
!-- surface and top boundary fluxes of water vapor
#ifdef PICHAMBER

    USE vars, ONLY: fluxbq, fluxtq, fluxlq, fluxrq, fluxqq, fluxhq

    fluxbmk(:,:,:) = 0. ! initialize all fluxes at surface to zero
    fluxtmk(:,:,:) = 0. ! initialize all fluxes at top of domain to zero
    fluxlmk(:,:,:) = 0. ! initialize all fluxes at left wall to zero
    fluxrmk(:,:,:) = 0. ! initialize all fluxes at right wall of domain to zero
    fluxqmk(:,:,:) = 0. ! initialize all fluxes at front wall to zero
    fluxhmk(:,:,:) = 0. ! initialize all fluxes at back wall of domain to zero
    fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
    fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)
    fluxlmk(:,:,index_water_vapor) = fluxlq(:,:)
    fluxrmk(:,:,index_water_vapor) = fluxrq(:,:)
    fluxqmk(:,:,index_water_vapor) = fluxqq(:,:)
    fluxhmk(:,:,index_water_vapor) = fluxhq(:,:)

#else

    USE vars, ONLY: fluxbq, fluxtq

    fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
    fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)

#endif /*PICHAMBER*/

 END SUBROUTINE micro_flux

 SUBROUTINE micro_proc()
!
!-- Compute local microphysics processes (beyond advection and SGS diffusion):
!--
!-- This is the place where the condensation/sublimation, accretion, coagulation, freezing,
!-- melting, etc., are done, that is  all the microphysics processes except for the spatial
!-- transport and mixing.
!
!-- IMPORTANT: For consistancy, you need to use thermodynamic constants like specific heat,
!-- specific heat of condensation, gas constant, etc, the same as in file params.f90.
!-- Also, you should assume that the conservative thermodynamic variable during these
!-- proceses is the liquid/ice water static energy: t = tabs + gz - Lc (qc+qr) - Ls (qi+qs+qg)
!-- It should not change during all of your point microphysical processes!
    USE micro_vars,         ONLY: n_particles, flag_LCM_init, n_start, n_end, &
                                  nstep_LCM_init, tabs_LCM, pres_LCM,         &
                                  switch_sedimentation, gpzs, dtn_LCM, dV, rho_air,&
                                  tau_SCLR, vanthoff_aero, molecular_weight_aero,  &
                                  bfactor_coefficient, molecular_weight_water, rho_liq, &
                                  & spec_nsave, spec_nsavestart, particle_ntrack, particle_ntrackstart

    USE micro_restart,      ONLY: micro_restart_write
    USE micro_random,       ONLY: randomu
#ifdef SAMTAK
USE params,             ONLY: docloud, fac_cond, pi, dorestart_last_only
#else
USE params,             ONLY: docloud, fac_cond, pi
#endif /*SAMTAK*/
    USE vars,               ONLY: total_water_prec, qsatw, t, tabs, gamaz, precsfc, prec_xy
    USE grid,               ONLY: nstep, nstat, nrestart_skip, nstop, nelapse, icycle, ncycle, masterproc, rank, &
                                & nsave3D_mod, nsave3Dstart_mod, nsave3Dend_mod, nsave2D_mod, nsave2Dstart_mod,  &
                                & nsave2Dend_mod, nsave2D, nsave2Dstart, nsave2Dend, nsave3D, nsave3Dstart,      &
                                & nsave3Dend, nstatis

#ifdef SAMCHEM
    USE params,             ONLY: dochemistry
    USE chemistry,          ONLY: flag_chem_init
#endif /*SAMCHEM*/

    IMPLICIT NONE

    INTEGER :: i, j, k, n, m
#ifdef SAMCHEM
    REAL, ALLOCATABLE :: chem_new_radii_all(:)
#endif /*SAMCHEM*/


    CALL t_startf ('micro_LCM')

    ! check if output values need to be calculated
    t_out_field = ( MOD(nstep,nsave3D_mod) == 0 .AND. nstep >= nsave3Dstart_mod .AND. nstep <= nsave3Dend_mod ) .OR. &
                & (mod(nstep,nsave3D).eq.0.and.nstep.ge.nsave3Dstart.and.nstep.le.nsave3Dend )                  .OR. &
                & (mod(nstep,nsave2D).eq.0.and.nstep.ge.nsave2Dstart .and.nstep.le.nsave2Dend)                  .OR. &
                & ( MOD(nstep,nsave2D_mod) == 0 .AND. nstep >= nsave2Dstart_mod .AND. nstep <= nsave2Dend_mod )

    t_out_stat = (MOD(nstep,nstatis) .EQ. 0) ! condition to call statistics and micro_statistics afterwards, copied from stepout.f90
    t_out_spec = ( MOD(nstep,spec_nsave) .EQ. 0  .AND.  nstep .GE. spec_nsavestart  .AND.  icycle .EQ. ncycle )
    t_out_particle =  ( MOD(nstep,particle_ntrack) .EQ. 0 .AND. &
                      & nstep .GE. particle_ntrackstart   .AND. &
                      & nstep .LE. particle_ntrackstop    .AND. &
                       icycle .EQ. ncycle                       )
    t_out_any = (t_out_field .OR. t_out_stat .OR. t_out_spec .OR. t_out_particle)

!
!-- Initialization of LCM
    IF ( .NOT. flag_LCM_init  .AND.  nstep .EQ. nstep_LCM_init )  THEN
!
!--    Initialize particle IDs, value of n_ID depends on rank to assure individual IDs across the domain
       n_ID = n_step_ID * rank
!
!--    Initialize temperature fields
       micro_field(:,:,:,iTlcm)  = t(:,:,:)
       micro_field(:,:,:,iTllcm) = t(:,:,:)
!
!--    Exchange LES fields
       CALL micro_coupling
!
!--    Place and initialize particles
       CALL micro_place_particle
!
!--    Calculate release/depletion of water vapor/absolute supersaturation
       DO k = 1, nz_LCM
           DO j = 1, ny_LCM
              DO i = 1, nx_LCM
                 micro_field(i,j,k,iss) = micro_field(i,j,k,imv) - qsatw( tabs_LCM(i,j,k), pres_LCM(k) )
              ENDDO
           ENDDO
           micro_field(1:nx_LCM, 1:ny_LCM, k, ilastp) = pres_LCM(k) !advection of pressure will be calculated explicitly
       ENDDO
       micro_field(1:nx_LCM, 1:ny_LCM, 1:nz_LCM, ilastT) = tabs_LCM(1:nx_LCM, 1:ny_LCM ,1:nz_LCM) !advection of temperature will be calculated explicitly
!
!--    flag_LCM_init has to be true, now
       flag_LCM_init = .TRUE.

       IF (masterproc) THEN
          PRINT*, '************************************************************************************'
          PRINT*, '*** MICRO_LAGRANGE: initialization completed, Lagrangian cloud physics start now ***'
          PRINT*, '************************************************************************************'

       ENDIF

    ENDIF

#ifdef SAMCHEM
    IF ( .NOT. flag_chem_init ) THEN
      IF ( nstep .EQ. nstep_chem_init .AND. dochemistry ) THEN
        ! call chem_init, split chem_init into chem_setparm and read mechanism here and allocate particles and give concentrations..
        ! dont forget to add concentrations in particle_source routine as well

        IF (masterproc) WRITE(*,*) 'Chemistry :: Starting!'

        CALL chem_set_initial_concentrations

        flag_chem_init = .TRUE.
      END IF
    END IF
#endif /*SAMCHEM*/

    IF ( flag_LCM_init ) THEN
!
      !--    Set rate fields to zero
      micro_field(:,:,:,iACT)   = 0.0
      micro_field(:,:,:,iDEA)   = 0.0

!
!--    Set diagnostic microphysical fieds to zero
       cnd_prof  = 0.0
       evp_prof  = 0.0
       aut_prof  = 0.0
       acc_prof  = 0.0
       sfx_prof  = 0.0
       afx_prof  = 0.0
       sgfx_prof = 0.0
       act_prof  = 0.0
       imm_prof  = 0.0
       hzn_prof  = 0.0
       hom_prof  = 0.0
       nuc_prof  = 0.0
       nml_prof  = 0.0
       dea_prof  = 0.0
       ref_prof  = 0.0
       rsd_prof  = 0.0
       nin_prof  = 0.0
       frz_prof  = 0.0
       mel_prof  = 0.0
       dep_prof  = 0.0
       sub_prof  = 0.0
       nsub_prof = 0.0
       agg_prof  = 0.0
       icebrk_prof  = 0.0
       ds_prof  = 0.0
       rs_prof  = 0.0
       rim_prof  = 0.0
       coal_prof = 0.0
       brk_prof = 0.0

       std_SCLR_prof  = 0.0
       mean_SCLR_prof = 0.0

       precsfc = 0.0
       prec_xy = 0.0
!
!--    Exchange LES fields
       IF ( nstep .GT. nstep_LCM_init ) CALL micro_coupling ! already done for nstep .EQ. nstep_LCM_init
!
!--    Compute motion of particles
       DO k = 1, nz_LCM
           DO j = 1, ny_LCM
              DO i = 1, nx_LCM

                 IF ( n_end(i,j,k) - n_start(i,j,k) .LT. 0 )  CYCLE !if there are no particles inside a gridbox, do nothing
                 IF ( switch_sedimentation  .AND.  nstep .GE. start_sedimentation )  THEN
                    CALL micro_terminal_velocity(i,j,k)               ! Compute terminal velocity
                 ENDIF
                 CALL micro_uvw(i,j,k)                                ! Compute resolved velocity at t=t+1/2
                 IF ( switch_sgs_velocities )  CALL micro_sgs_uvw(i,j,k)    ! Compute SGS velocity
              ENDDO
          ENDDO
       ENDDO
!
!--    Update particle x, y, z position
       CALL micro_xyz
!
!--    Boundary conditions (bottom/top)
       CALL micro_boundaries
!
!--    Parcel exchange across processors
       CALL micro_exchange_particles
!
!--    Sort all particles by their new grid box (this makes them easily accessible
!--    via n_start and n_end in a 3D loop structure)
       CALL micro_sort_by_grid_box
!
!--    Add additional particles to domain
       IF ( switch_particle_source  .AND.  nstep .GE. MIN(n_start_source_aero, n_start_source_in) &
          .AND.  (nstep .LE. MAX(n_end_source_aero, n_end_source_in ) .OR. (MIN(n_end_source_aero, n_end_source_in) .LT. 0)) .AND. (n_start_source_aero .GT. 0 .OR. n_start_source_in .GT. 0) &
          )  CALL micro_particle_source
#ifdef SAMCHEM
!
!--    take care of the supplies (arrays) needed for chemistry numerics
       ALLOCATE(chem_new_radii_all(n_particles))
       IF (flag_chem_init) THEN
         ! HERE could be your monster-chem ...
         chem_new_radii_all = 0.0
         !nTimestepsZ = 0
!--    call chemistry for all grid boxes at once (currently not used but important option)
         !CALL chem_monster_proc(chem_new_radii_all, micro_field(:,:,:,imv))
       END IF
#endif /*SAMCHEM*/

!--    Compute microphysics

       DO k = 1, nz_LCM
           DO j = 1, ny_LCM
              DO i = 1, nx_LCM

                 IF ( n_end(i,j,k) - n_start(i,j,k) .LT. 0 )  CYCLE !if there are no particles inside a gridbox, do nothing
!
!--              Exponential decay of plume-detecting scalar
                 micro_field(i,j,k,iSCLR) = micro_field(i,j,k,iSCLR) - micro_field(i,j,k,iSCLR) / tau_SCLR * dtn_LCM
                 IF ( switch_LEM_micro )  CALL micro_sgs_mixing(i,j,k) ! Compute SGS mixing
!
!--              Compute condensation/evaporation with possibility to get new radii from chemistry
#ifdef SAMCHEM
                 CALL micro_condensation(i,j,k, chem_new_radii_all(n_start(i,j,k):n_end(i,j,k)))
#else
                 CALL micro_condensation(i,j,k)
#endif /*SAMCHEM*/
!
                 IF ( switch_collection  .AND.  nstep .GE. start_collection )  THEN
                    CALL micro_collision(i,j,k) ! Compute collision/coalescence
                 ENDIF
!
!--              Prepare calculation of particle motion for next time step
                 CALL micro_uvw_prep(i,j,k)

              ENDDO
          ENDDO
       ENDDO
       
!
!--    Diagnose microphysical quantities
       CALL micro_diagnose()
!
!--    If required, save restart data
       IF ( icycle == ncycle ) THEN
          IF ( MOD( nstep, nstat * ( 1 + nrestart_skip ) ) .EQ. 0  &
               .OR.  nstep .EQ. nstop  .OR.  nelapse .EQ. 0 ) THEN

#ifdef SAMTAK
             IF ( .NOT. dorestart_last_only  .OR.  ( dorestart_last_only .AND. &
                  ( nstep .EQ. nstop  .OR.  nelapse .EQ. 0 ) ) )  CALL micro_restart_write
#else
             CALL micro_restart_write
#endif /*SAMTAK*/

          ENDIF
       ENDIF

    ENDIF

    CALL t_stopf ('micro_LCM')

 END SUBROUTINE micro_proc

 SUBROUTINE micro_diagnose()
!
!-- Diagnostics for several parameters.
!-- tabs, qv, qcl, qpl are essential for SAM's dynamical core
    USE micro_vars, ONLY: tabs_LCM, rho_air, dtn_LCM, gpzs, dz_LCM, nac_prof,  &
                          ncc_prof, nrc_prof, nsd_prof, rho_liq, iZ, ref_prof, &
                          rsd_prof, qa_prof, diss_LCM, qnmin_LCM,              &
                          iREF_prof, iR_prof, iPHI_prof, iAXS, iCXS, iRSD_prof,&
                          fl_prof, fi_prof, fin_prof, fc_prof, fr_prof,        &
                          fp_prof, fs_prof, fg_prof, dummy_prof, iREFS_prof,   &
                          itime, t_mean_prof, t_max_prof, t_count_prof,        &
                          t_mean_prof_liq, t_max_prof_liq, t_count_prof_liq,   &
                          t_mean_prof_ice, t_max_prof_ice, t_count_prof_ice,   &
                          get_r_particle, t_out_stat, t_out_field, t_out_spec, &
                          t_out_particle, t_out_any

    USE micro_user_vars

    USE vars
! #ifdef SAMTAK
!     USE vars,   ONLY: qv, qcl, qci, qpl, t, tabs, opdqcl, opdqpl, qclopd_xy, qplopd_xy, &
!                       ncl_xy, npl_xy, qsatw, w, gamaz
! #else
!     USE vars,   ONLY: qv, qcl, qci, qpl, t, tabs, qsatw, w, gamaz
! #endif /*SAMTAK*/

    USE grid,   ONLY: dtfactor, masterproc, time, nstep, nstop, ncycle, icycle, &
                      rank, nsubdomains_x, nsubdomains_y, z
    USE params, ONLY: pi, rv, fac_cond, epsv, ggr, cp, fac_sub

    USE domain, ONLY: nx_gl, ny_gl

#ifdef SAMCHEM
    USE chemistry, ONLY: calc_ionic_strength, iHp, iOHm, ionic_strength_threshold
    USE params,    ONLY: dochemistry
#endif /*SAMCHEM*/

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    INTEGER :: i, j, k, l, n, ierr, ir, n_bin
    REAL    :: optdepth, radius, r_eff, r0_ijk, r1_ijk, r2_ijk, r3_ijk, r4_ijk,          &
               qsatw_pre, mass, mean_SCLR_local(nz_LCM), std_SCLR_local(nz_LCM),         &
               min_SCLR, T_local(nz_LCM), qv_local(nz_LCM), ql_local(nz_LCM),            &
               T_global(nz_LCM), qv_global(nz_LCM), ql_global(nz_LCM), B_global(nz_LCM), &
               B_local, S_sgs, mean2_dummy, T_dummy, qv_dummy, particle_tabs, Q_ext,     &
               dV_spec, beta_dummy, ir0_ijk, ir1_ijk, ir2_ijk,                           &
               ir3_ijk, ir4_ijk, PHI, Nmono, rimemass, icemass_total, phi_mom0, phi_mom1,&
               rho_mom1, ias_sum, ics_sum, incount, ssasum, dsasum, bsasum, psasum,      &
               mass_aero, mean_dummy, av0_ijk, av1_ijk, am_ijk,                          &
               cldmass_total, rainmass_total, prismass_total, snowmass_total,            &
               graupelmass_total

    REAL, EXTERNAL :: capacitance

    LOGICAL :: cloudy_flag

#ifdef SAMCHEM
    REAL :: pH_particle, ionic_strength
    LOGICAL :: calc_pH_flag
#endif /*SAMCHEM*/

    CHARACTER(LEN=9) :: nstep_char

    REAL, PARAMETER :: mass_1 = 4.0 / 3.0 * pi * rho_liq * ( 3.125E-6 / 2.0 )**3.0 ! upper boundary of first TAU bin

    CALL t_startf ('micro_diagnose')

    ! output-only variables
    ! to avoid untangling dependence-trees of output variables, all output is diagnosed if any output is required in the current time step
    IF ( flag_LCM_init &
            !.AND. t_out_any &
        )  THEN
!
!--    Create user output quantities
       CALL micro_user_diagnose
!
!--    Recalculate diagnostic fields for aerosol number, cloud water ...
       micro_field(:,:,:,ina)    = 0.0
       micro_field(:,:,:,iqa)    = 0.0
       micro_field(:,:,:,inc)    = 0.0
       micro_field(:,:,:,iqc)    = 0.0
       micro_field(:,:,:,inr)    = 0.0
       micro_field(:,:,:,iqr)    = 0.0

       micro_field(:,:,:,ini)    = 0.0
       micro_field(:,:,:,iqi)    = 0.0
       micro_field(:,:,:,ins)    = 0.0
       micro_field(:,:,:,iqs)    = 0.0
       micro_field(:,:,:,ing)    = 0.0
       micro_field(:,:,:,iqg)    = 0.0

       micro_field(:,:,:,ininc)  = 0.0

       micro_field(:,:,:,iRSD)   = 0.0
       micro_field(:,:,:,iIREF)   = 0.0

       micro_field(:,:,:,iCXS)   = 0.0
       micro_field(:,:,:,iAXS)   = 0.0

       micro_Field(:,:,:,itfrz)  = 0.0
       micro_field(:,:,:,itfrzc) = 0.0
       micro_field(:,:,:,itfrzr) = 0.0
       micro_field(:,:,:,itfrzi) = 0.0
       micro_field(:,:,:,itfrzs) = 0.0
       micro_field(:,:,:,itfrzg) = 0.0

       micro_field(:,:,:,imsp)   = 0.0
       micro_field(:,:,:,issp)   = 0.0
       micro_field(:,:,:,iiss)   = 0.0

       micro_field(:,:,:,itmean) = 0.0
       micro_field(:,:,:,itmax)  = 0.0
       micro_field(:,:,:,iltmean) = 0.0
       micro_field(:,:,:,iltmax)  = 0.0
       micro_field(:,:,:,iitmean) = 0.0
       micro_field(:,:,:,iitmax)  = 0.0
!
!--    sprayed particle quantities
       micro_field(:,:,:,iasp)    = 0.0
       micro_field(:,:,:,icsp)    = 0.0

       micro_field(:,:,:,iphi)    = 0.0
       micro_field(:,:,:,irmf)    = 0.0
       micro_field(:,:,:,irhoapp) = 0.0

       IF ( switch_ice )  THEN
          spec_ice = 0.0
       ENDIF
#ifdef SAMCHEM
       IF (switch_spec_pH) THEN
         spec_pH   = 0.0
       ENDIF
#endif /*SAMCHEM*/
       IF (switch_spec_aero) THEN
         spec_aero = 0.0
       ENDIF
       spec_liq = 0.0

       nsd_prof = 0.0

       nac_prof  = 0.0
       ncc_prof  = 0.0
       nrc_prof  = 0.0
       nic_prof  = 0.0
       ninc_prof = 0.0

       spm_prof = 0.0
       sps_prof = 0.0
       T_prof   = 0.0
       qv_prof  = 0.0

       ref_prof  = 0.0
       refs_prof = 0.0
       r_prof    = 0.0
       rsd_prof  = 0.0
       ar_prof   = 0.0
       am_prof   = 0.0

       t_mean_prof = 0.0
       t_max_prof = 0.0
       t_count_prof = 0.0
       t_mean_prof_liq = 0.0
       t_max_prof_liq = 0.0
       t_count_prof_liq = 0.0
       t_mean_prof_ice = 0.0
       t_max_prof_ice = 0.0
       t_count_prof_ice = 0.0

       qa_prof  = 0.0
       qc_prof  = 0.0
       qr_prof  = 0.0

       qi_prof  = 0.0
       qs_prof  = 0.0
       qg_prof  = 0.0

       rhice_prof = 0.0

       iREF_prof  = 0.0
       iREFS_prof = 0.0
       iR_prof    = 0.0
       iRSD_prof  = 0.0
       iaavg_prof  = 0.0
       icavg_prof  = 0.0
       rho_app_prof = 0.0

       IF (switch_ice) THEN
         rmf_prof = 0.0
         tfrz_prof = 0.0
         tfrzc_prof = 0.0
         tfrzr_prof = 0.0
         tfrzi_prof = 0.0
         tfrzs_prof = 0.0
         tfrzg_prof = 0.0
         ssa_prof = 0.0
         dsa_prof = 0.0
         bsa_prof = 0.0
         psa_prof = 0.0
       ENDIF

       iPHI_prof  = 0.0
       Nmono_prof = 0.0

       fl_prof = 0.0
       fa_prof = 0.0
       IF (switch_ice) THEN
         fi_prof = 0.0
         fin_prof = 0.0
         fc_prof = 0.0
         fr_prof = 0.0
         fp_prof = 0.0
         fs_prof = 0.0
         fg_prof = 0.0
       ENDIF

#ifdef SAMTAK
       opdqcl = 0.0
       opdqpl = 0.0
#endif /*SAMTAK*/

       tau_prof    = 0.0

       DO k = 1, nz_LCM

          DO j = 1, ny_LCM
             DO i = 1, nx_LCM
!
!--             Get relative humidity over ice
                rhice_prof(k) = rhice_prof(k) + micro_field(i,j,k,imv) / qsati( tabs_LCM(i,j,k), pres_LCM(k) )
!
!--             Get supersaturation over ice
                micro_field(i,j,k,iiss) = micro_field(i,j,k,imv) - qsati( tabs_LCM(i,j,k), pres_LCM(k) )
!
!--             Get kinetic energy dissipation rate
                micro_field(i,j,k,idiss) = diss_LCM(i,j,k)
!
!--             Compute diagnostic fields
                IF ( n_end(i,j,k) - n_start(i,j,k) .LT. 0 )  CYCLE !if there are no particles inside that gridbox, do nothing
!
!--             A priori estimate for cloudy grid box
                IF ( SUM(particle(n_start(i,j,k):n_end(i,j,k))%var(iM), mask=particle(n_start(i,j,k):n_end(i,j,k))%ivar(iTYPE) .GT. 0) / ( dV(k) * rho_air(k) ) .GT. qnmin_LCM )  THEN
                   cloudy_flag = .TRUE.
                ELSE
                   cloudy_flag = .FALSE.
                ENDIF

                r0_ijk = 0.0
                r1_ijk = 0.0
                r2_ijk = 0.0
                r3_ijk = 0.0
                r4_ijk = 0.0

                av0_ijk = 0.0
                av1_ijk = 0.0
                am_ijk  = 0.0

                ir0_ijk = 0.0
                ir1_ijk = 0.0
                ir2_ijk = 0.0
                ir3_ijk = 0.0
                ir4_ijk = 0.0

                PHI  = 0.0
                Nmono = 0.0

                rimemass = 0.0
                icemass_total = 0.0
                cldmass_total = 0.0
                rainmass_total = 0.0
                prismass_total = 0.0
                snowmass_total = 0.0
                graupelmass_total = 0.0
                phi_mom0 = 0.0
                phi_mom1 = 0.0
               !  rho_mom0 = 0.0
                rho_mom1 = 0.0
                ias_sum = 0.0
                ics_sum = 0.0
                IF (switch_ice) THEN
                  ssasum = 0.0
                  dsasum = 0.0
                  bsasum = 0.0
                  psasum = 0.0
                ENDIF

                incount = 0
!
!--             number of super-droplets profile
                nsd_prof(k) = nsd_prof(k) + ( n_end(i,j,k) - n_start(i,j,k) + 1.0 )
!
!--             sgs mixing model profiles
                IF ( switch_LEM_micro )  THEN

                   mean_dummy  = 0.0
                   mean2_dummy = 0.0
                   T_dummy     = 0.0
                   qv_dummy    = 0.0
                   DO n = n_start(i,j,k), n_end(i,j,k)
                      particle_tabs = particle(n)%var(iT) - particle(n)%var(iZ) * ggr / cp
                      S_sgs = particle(n)%var(iq) / qsatw( particle_tabs, pres_LCM(k) ) &
                              - micro_field(i,j,k,imv) / qsatw( tabs_LCM(i,j,k), pres_LCM(k) )
                      mean_dummy  = mean_dummy  + S_sgs
                      mean2_dummy = mean2_dummy + S_sgs**2
                      T_dummy     = T_dummy     + particle_tabs
                      qv_dummy    = qv_dummy    + particle(n)%var(iq)
                   END DO
                   mean_dummy  = mean_dummy  / ( n_end(i,j,k) - n_start(i,j,k) + 1 )
                   mean2_dummy = mean2_dummy / ( n_end(i,j,k) - n_start(i,j,k) + 1 )
                   T_dummy     = T_dummy     / ( n_end(i,j,k) - n_start(i,j,k) + 1 )
                   qv_dummy    = qv_dummy    / ( n_end(i,j,k) - n_start(i,j,k) + 1 )

                   spm_prof(k) = spm_prof(k) + mean_dummy
                   sps_prof(k) = sps_prof(k) + SQRT( MAX( 0.0, mean2_dummy - mean_dummy**2 ) )
                   T_prof(k)   = T_prof(k)   + T_dummy
                   qv_prof(k)  = qv_prof(k)  + qv_dummy

                   micro_field(i,j,k,imsp) = mean_dummy
                   micro_field(i,j,k,issp) = SQRT( MAX( 0.0, mean2_dummy - mean_dummy**2 ) )

                ENDIF

                DO  n = n_start(i,j,k), n_end(i,j,k)

                  IF (switch_ice) THEN
!
!-- Get average freezing temperature in a gridbox
                     micro_field(i,j,k,itfrz) = micro_field(i,j,k,itfrz) + particle(n)%var(iTFZ)*particle(n)%var(iM)
                  ENDIF

!
!--                Get particle residence times
                  IF ( cloudy_flag .and. (particle(n)%ivar(iTYPE) .GT. 0) )  THEN
                     particle(n)%var(itime) = particle(n)%var(itime) + dtn_LCM

                     IF ( particle(n)%var(itime) .GT. micro_field(i,j,k,itmax) )  THEN
                        micro_field(i,j,k,itmax) = particle(n)%var(itime)
                     ENDIF
                     micro_field(i,j,k,itmean) = micro_field(i,j,k,itmean) + particle(n)%var(itime) * particle(n)%ivar(iA)

                     IF ( particle(n)%var(itime) .GT. t_max_prof(k) )  THEN
                        t_max_prof(k) = particle(n)%var(itime)
                     ENDIF
                     t_mean_prof(k)  = t_mean_prof(k) + particle(n)%var(itime) * particle(n)%ivar(iA)
                     t_count_prof(k) = t_count_prof(k) + particle(n)%ivar(iA)

                     IF ( particle(n)%var(itime) .GT. micro_field(i,j,k,iltmax) )  THEN
                        micro_field(i,j,k,iltmax) = particle(n)%var(itime)
                     ENDIF
                     micro_field(i,j,k,iltmean) = micro_field(i,j,k,iltmean) + particle(n)%var(itime) * particle(n)%ivar(iA)

                     IF ( particle(n)%var(itime) .GT. t_max_prof_liq(k) )  THEN
                        t_max_prof_liq(k) = particle(n)%var(itime)
                     ENDIF
                     t_mean_prof_liq(k)  = t_mean_prof_liq(k) + particle(n)%var(itime) * particle(n)%ivar(iA)
                     t_count_prof_liq(k) = t_count_prof_liq(k) + particle(n)%ivar(iA)
                  ELSEIF ((micro_field(i,j,k,iiss) .GT. 0.0) .and. (particle(n)%ivar(iTYPE) .LT. 0)) THEN
                     particle(n)%var(itime) = particle(n)%var(itime) + dtn_LCM

                     IF ( particle(n)%var(itime) .GT. micro_field(i,j,k,itmax) )  THEN
                        micro_field(i,j,k,itmax) = particle(n)%var(itime)
                     ENDIF
                     micro_field(i,j,k,itmean) = micro_field(i,j,k,itmean) + particle(n)%var(itime) * particle(n)%ivar(iA)

                     IF ( particle(n)%var(itime) .GT. t_max_prof(k) )  THEN
                        t_max_prof(k) = particle(n)%var(itime)
                     ENDIF
                     t_mean_prof(k)  = t_mean_prof(k) + particle(n)%var(itime) * particle(n)%ivar(iA)
                     t_count_prof(k) = t_count_prof(k) + particle(n)%ivar(iA)
                     
                     IF ( particle(n)%var(itime) .GT. micro_field(i,j,k,iitmax) )  THEN
                        micro_field(i,j,k,iitmax) = particle(n)%var(itime)
                     ENDIF
                     micro_field(i,j,k,iitmean) = micro_field(i,j,k,iitmean) + particle(n)%var(itime) * particle(n)%ivar(iA)

                     IF ( particle(n)%var(itime) .GT. t_max_prof_ice(k) )  THEN
                        t_max_prof_ice(k) = particle(n)%var(itime)
                     ENDIF
                     t_mean_prof_ice(k)  = t_mean_prof_ice(k) + particle(n)%var(itime) * particle(n)%ivar(iA)
                     t_count_prof_ice(k) = t_count_prof_ice(k) + particle(n)%ivar(iA)
                  ELSE
                     particle(n)%var(itime) = 0.0_lcmp
                  ENDIF
!
!--                sprayed aerosol estimation
                   IF (particle(n)%ivar(iSPR) .GT. 0) THEN
                      micro_field(i,j,k,iasp)    = micro_field(i,j,k,iasp) + particle(n)%ivar(iA)
                   ENDIF

                   IF (track_some_particles) THEN
                       IF ( MOD(nstep,particle_ntrack) .EQ. 0 .AND. &
                          & nstep .GE. particle_ntrackstart   .AND. &
                          & nstep .LE. particle_ntrackstop    .AND. &
                          & icycle .EQ. ncycle                      ) THEN

                         IF (particle(n)%ivar(iTRACK) .GE. 0) THEN
                           !OPEN(UNIT=126,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.tracked_particles',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND')
#ifdef SAMCHEM
                           IF (dochemistry) THEN
                             WRITE(126, format_particle_tracking) nstep, particle(n)%ivar(iTRACK),       &
                                                                & particle(n)%var(iX),                   &
                                                                & particle(n)%var(iY),                   &
                                                                & particle(n)%var(iZ),                   &
                                                                & particle(n)%ivar(iA),                  &
                                                                & particle(n)%var(iM),                   &
                                                                & particle(n)%var(iN),                   &
                                                                & micro_field(i,j,k,isLES),              &
                                                                & particle(n)%var(iaqspc)
                             FLUSH(126)
                           ELSE
#endif /*SAMCHEM*/
                             WRITE(126, format_particle_tracking) nstep, particle(n)%ivar(iTRACK),       &
                                                                & particle(n)%var(iX),                   &
                                                                & particle(n)%var(iY),                   &
                                                                & particle(n)%var(iZ),                   &
                                                                & particle(n)%ivar(iA),                  &
                                                                & particle(n)%var(iM),                   &
                                                                & particle(n)%var(iN),                   &
                                                                & micro_field(i,j,k,isLES)
                             FLUSH(126)
#ifdef SAMCHEM
                           END IF
#endif /*SAMCHEM*/
                           !CLOSE(126)
                         END IF
                       END IF
                   END IF

#ifdef SAMCHEM
                   IF ( switch_spec_pH .AND. flag_chem_init ) THEN
                       ionic_strength = calc_ionic_strength(particle(n)%var(iaqspc)/N_avogadro, particle(n)%var(iM), iHp, iOHm)
                       calc_pH_flag = ionic_strength < ionic_strength_threshold
                       IF (calc_pH_flag) THEN
                         ! convert from                    molec   to    mol   to  mol/l
                         pH_particle = -LOG10( particle(n)%var(iaqspc(iHp))/N_avogadro /   particle(n)%var(iM))
                       ELSE
                         ! set pH to zero because chemistry is off in particle
                         pH_particle = 0.0
                       END IF

                       n_bin = CEILING( (pH_particle-pH_start_spec) / (pH_end_spec-pH_start_spec) * n_bins_spec )
                       IF ( n_bin .GE. 1  .AND.  n_bin .LE. n_bins_spec  .AND.  k .LE. nzl )  THEN
                          IF ( horiz_average_spec )  THEN
                             !spec_pH(1,1,k,n_bin) = spec_pH(1,1,k,n_bin) + particle(n)%ivar(iA)
                             ! i think adding the litres is better to know what budget of reactions we have, i.e.,
                             ! knowing how many droplets have a certain pH is not really useful without knowing their radius,
                             ! i think saving the litres is a better measure
                             spec_pH(1,1,k,n_bin) = spec_pH(1,1,k,n_bin) + particle(n)%var(iM)
                          ELSE
                             !spec_pH(i_global_LCM+i,j_global_LCM+j,k,n_bin) = spec_pH(i_global_LCM+i,j_global_LCM+j,k,n_bin) + particle(n)%ivar(iA)
                             ! why save iM? see comments above
                             spec_pH(i_global_LCM+i,j_global_LCM+j,k,n_bin) = spec_pH(i_global_LCM+i,j_global_LCM+j,k,n_bin) + particle(n)%var(iM)
                          ENDIF
                       ENDIF
                   END IF
#endif /*SAMCHEM*/

                   IF (switch_spec_aero) THEN
!--                    Calculate aerosol spectrum
                       ! the following is not applicable here as aerosols have different mass (and densities might change with time
                       ! due to collisions) and therefore mass bin <-> radius bin connection doesnt work as it's created with water density
                       !mass_aero = particle(n)%var(iN) / particle(n)%ivar(iA)
                       !n_bin       = CEILING( LOG( mass_aero / x_start_spec ) / LOG(2.0**alpha_spec) )
                       ! therefore, search bin directly
                       n_bin = search_sorted_array(rl_spec, get_r_aero(n), n_bins_spec)
                       IF ( n_bin .GE. 1  .AND.  n_bin .LE. n_bins_spec  .AND.  k .LE. nzl )  THEN
                          IF ( horiz_average_spec )  THEN
                             spec_aero(1,1,k,n_bin) = spec_aero(1,1,k,n_bin) + particle(n)%ivar(iA)
                          ELSE
                             spec_aero(i_global_LCM+i,j_global_LCM+j,k,n_bin) = spec_aero(i_global_LCM+i,j_global_LCM+j,k,n_bin) + particle(n)%ivar(iA)
                          ENDIF
                       ENDIF
                   ENDIF

                   mass   = particle(n)%var(iM) / particle(n)%ivar(iA)
!
!--                Get statistics of haze, cloud drops, rain drops
                   IF ( particle(n)%ivar(iTYPE) .GT. 0 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN

                      radius = get_r_particle(n)
                      mass   = (particle(n)%var(iM)+particle(n)%var(iN)) / particle(n)%ivar(iA)
!
!--                   optical thickness with Mie scattering
                      beta_dummy = Q_ext(radius) * pi * radius**2 * particle(n)%ivar(iA) / ( dx_LCM * dy_LCM )
                      tau_prof(k)                = tau_prof(k)                + beta_dummy
!
!--                   Calculate spectrum
                      n_bin = CEILING( LOG( mass / x_start_spec ) / LOG(2.0**alpha_spec) )
                      IF ( n_bin .GE. 1  .AND.  n_bin .LE. n_bins_spec  .AND.  k .LE. nzl )  THEN
                         IF ( horiz_average_spec )  THEN
                            spec_liq(1,1,k,n_bin) = spec_liq(1,1,k,n_bin) + particle(n)%ivar(iA)
                         ELSE
                            spec_liq(i_global_LCM+i,j_global_LCM+j,k,n_bin) = spec_liq(i_global_LCM+i,j_global_LCM+j,k,n_bin) + particle(n)%ivar(iA)
                         ENDIF
                      ENDIF

!
!--                   Get aerosol mass and radius
                      IF (particle(n)%var(iN) .gt. 0) THEN
                         av0_ijk = av0_ijk + particle(n)%ivar(iA)
                         av1_ijk = av1_ijk + SUM((3 * particle(n)%var(iAero_portions) * particle(n)%var(iN) / (4 * pi * particle(n)%ivar(iA) * rho_aero)))**(1./3.) * particle(n)%ivar(iA)
                         am_ijk = am_ijk + particle(n)%var(iN) * particle(n)%ivar(iA)

                      ENDIF
!
!--                   Get IN among aerosols, cloud droplets, rain drops
                      IF ( particle(n)%ivar(iTYPE) .EQ. 2 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
                         micro_field(i,j,k,ininc) = micro_field(i,j,k,ininc) + particle(n)%ivar(iA)
                         incount = incount + particle(n)%ivar(iA)
                         ssasum = ssasum + particle(n)%var(iSC)
                         dsasum = dsasum + particle(n)%var(iSD)
                         bsasum = bsasum + particle(n)%var(iSB)
                         psasum = psasum + particle(n)%var(iSP)
                      ENDIF
!
!--                   Get data on haze, cloud, and rain particles
                      IF ( radius .LT. r_act )  THEN
                         micro_field(i,j,k,iqa) = micro_field(i,j,k,iqa) + particle(n)%var(iM)  !also consider mass of aerosols here
                         micro_field(i,j,k,ina) = micro_field(i,j,k,ina) + particle(n)%ivar(iA)
                         qa_prof(k)             = qa_prof(k) + particle(n)%var(iM)

                      ELSE
!
!--                      moments of radius profile
                         r0_ijk = r0_ijk + radius**0 * particle(n)%ivar(iA)
                         r1_ijk = r1_ijk + radius**1 * particle(n)%ivar(iA)
                         r2_ijk = r2_ijk + radius**2 * particle(n)%ivar(iA)
                         r3_ijk = r3_ijk + radius**3 * particle(n)%ivar(iA)
                         r4_ijk = r4_ijk + radius**4 * particle(n)%ivar(iA)

                         IF ( radius .LT. r_separate_micro )  THEN
#ifdef SAMTAK
                            opdqcl(i,j) = opdqcl(i,j) + particle(n)%ivar(iA) * radius**2            !optical thickness cloud droplets
#endif /*SAMTAK*/
                            micro_field(i,j,k,iqc) = micro_field(i,j,k,iqc) + particle(n)%var(iM)   !mass cloud droplets
                            micro_field(i,j,k,inc) = micro_field(i,j,k,inc) + particle(n)%ivar(iA)  !number cloud droplets
                            qc_prof(k)             = qc_prof(k) + particle(n)%var(iM)
!
!-- Get freezing temperature of cloud droplets
                            IF (switch_ice) THEN
                              micro_field(i,j,k,itfrzc) = micro_field(i,j,k,itfrzc) + particle(n)%var(iTFZ)*particle(n)%var(iM)
                              cldmass_total = cldmass_total + particle(n)%var(iM)
                            ENDIF
!
!--                         Sprayed aerosol estimation
                            IF (particle(n)%ivar(iSPR) .GT. 0) THEN
                              micro_field(i,j,k,icsp)    = micro_field(i,j,k,icsp) + particle(n)%ivar(iA)
                            ENDIF

                         ELSE
#ifdef SAMTAK
                            opdqpl(i,j) = opdqpl(i,j) + particle(n)%ivar(iA) * radius**2            !optical thickness rain drops
#endif /*SAMTAK*/
                            micro_field(i,j,k,iqr) = micro_field(i,j,k,iqr) + particle(n)%var(iM)   !mass rain drops
                            micro_field(i,j,k,inr) = micro_field(i,j,k,inr) + particle(n)%ivar(iA)  !number rain drops
                            qr_prof(k)             = qr_prof(k) + particle(n)%var(iM)
!
!-- Get freezizng temperature of rain drops
                            IF (switch_ice) THEN
                              micro_field(i,j,k,itfrzr) = micro_field(i,j,k,itfrzr) + particle(n)%var(iTFZ) * particle(n)%var(iM)
                              rainmass_total = rainmass_total + particle(n)%var(iM)
                            ENDIF
                         ENDIF

                      ENDIF
!
!--                Get statistics of ice crystals
                   ELSEIF  ( particle(n)%ivar(iTYPE) .LT. 0 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
                      radius = MAX( particle(n)%var(iCXS), particle(n)%var(iAXS) )
                      icemass_total = icemass_total + particle(n)%var(iM)
                      rimemass = rimemass + particle(n)%var(iMrime)
                      phi_mom0 = phi_mom0 + particle(n)%ivar(iA) * capacitance(particle(n)%var(iAXS), particle(n)%var(iCXS))
                      phi_mom1 = phi_mom1 + (particle(n)%var(iCXS)/particle(n)%var(iAXS)) * &
                      particle(n)%ivar(iA) * capacitance(particle(n)%var(iAXS), particle(n)%var(iCXS))
                     !  rho_mom0 = rho_mom0 + particle(n)%var(iM)
                      rho_mom1 = rho_mom1 + particle(n)%var(iM) * particle(n)%var(iDNS)
                      ias_sum = ias_sum + particle(n)%var(iAXS) * particle(n)%ivar(iA)
                      ics_sum = ics_sum + particle(n)%var(iCXS) * particle(n)%ivar(iA)
                      
                      IF (particle(n)%ivar(iTYPE) .EQ. -2) THEN !Get IN in ice crystals
                        incount = incount + particle(n)%ivar(iA)
                        ssasum = ssasum + particle(n)%var(iSC)
                        dsasum = dsasum + particle(n)%var(iSD)
                        bsasum = bsasum + particle(n)%var(iSB)
                        psasum = psasum + particle(n)%var(iSP)
                      ENDIF

!
!--                   optical thickness with Mie scattering
                      tau_prof(k) = tau_prof(k) + Q_ext(radius) * pi * radius**2 * particle(n)%ivar(iA) / ( dx_LCM * dy_LCM )
!
!--                   Calculate spectrum
                      n_bin = CEILING( LOG( mass / x_start_spec ) / LOG(2.0**alpha_spec) )
                      IF ( n_bin .GE. 1  .AND.  n_bin .LE. n_bins_spec  .AND.  k .LE. nzl .and. switch_ice)  THEN
                         IF ( horiz_average_spec )  THEN
                            spec_ice(1,1,k,n_bin) = spec_ice(1,1,k,n_bin) + particle(n)%ivar(iA)
                         ELSE
                            spec_ice(i_global_LCM+i,j_global_LCM+j,k,n_bin) = spec_ice(i_global_LCM+i,j_global_LCM+j,k,n_bin) + particle(n)%ivar(iA)
                         ENDIF
                      ENDIF
!
!--                   moments of radius profile
                      ir0_ijk = ir0_ijk + radius**0 * particle(n)%ivar(iA)
                      ir1_ijk = ir1_ijk + radius**1 * particle(n)%ivar(iA)
                      ir2_ijk = ir2_ijk + radius**2 * particle(n)%ivar(iA)
                      ir3_ijk = ir3_ijk + radius**3 * particle(n)%ivar(iA)
                      ir4_ijk = ir4_ijk + radius**4 * particle(n)%ivar(iA)

                      IF ( particle(n)%var(iMrime) / particle(n)%var(iM) .GT. 0.3 )  THEN ! graupel, see Shima et al. (2020), Sec. 7.1
                         micro_field(i,j,k,iqg) = micro_field(i,j,k,iqg) + particle(n)%var(iM)   !mass ice
                         micro_field(i,j,k,ing) = micro_field(i,j,k,ing) + particle(n)%ivar(iA)  !number ice
                         qg_prof(k)             = qg_prof(k) + particle(n)%var(iM)
!
!--                   freezing temperature of graupel
                         micro_field(i,j,k,itfrzg) = micro_field(i,j,k,itfrzg) + particle(n)%var(iTFZ) * particle(n)%var(iM)
                         graupelmass_total = graupelmass_total + particle(n)%var(iM)
                      ELSEIF ( particle(n)%ivar(iNmono) .GT. 10 )  THEN ! snow, see Shima et al. (2020), Sec. 7.1
                         micro_field(i,j,k,iqs) = micro_field(i,j,k,iqs) + particle(n)%var(iM)   !mass ice
                         micro_field(i,j,k,ins) = micro_field(i,j,k,ins) + particle(n)%ivar(iA)  !number ice
                         qs_prof(k)             = qs_prof(k) + particle(n)%var(iM)
!
!--                   freezing temperature of snow
                         micro_field(i,j,k,itfrzs) = micro_field(i,j,k,itfrzs) + particle(n)%var(iTFZ)*particle(n)%var(iM)
                         snowmass_total = snowmass_total + particle(n)%var(iM)
                      ELSE ! ice, see Shima et al. (2020), Sec. 7.1
                         micro_field(i,j,k,iqi) = micro_field(i,j,k,iqi) + particle(n)%var(iM)   !mass ice
                         micro_field(i,j,k,ini) = micro_field(i,j,k,ini) + particle(n)%ivar(iA)  !number ice
                         qi_prof(k)             = qi_prof(k) + particle(n)%var(iM)
!
!--                   freezing temperature of 'pristine' ice
                         micro_field(i,j,k,itfrzi) = micro_field(i,j,k,itfrzi) + particle(n)%var(iTFZ)*particle(n)%var(iM)
                         prismass_total = prismass_total + particle(n)%var(iM)
                      ENDIF

                      PHI  = PHI  + particle(n)%var(iCXS) / ( particle(n)%var(iAXS) + 1.0E-20 ) * particle(n)%ivar(iA)
                      Nmono = Nmono + particle(n)%ivar(iNmono) * particle(n)%ivar(iA)

                   ENDIF

                END DO

                IF ((switch_ice) .and. (incount .gt. 0.0)) THEN
                  fin_prof(k) = fin_prof(k) + 1.0 ! Count of columns with IN
                  ssa_prof(k) = ssa_prof(k) + ssasum
                  dsa_prof(k) = dsa_prof(k) + dsasum
                  bsa_prof(k) = bsa_prof(k) + bsasum
                  psa_prof(k) = psa_prof(k) + psasum
                  micro_field(i,j,k,issa) = ssasum ! Average soot surface area within the gridbox
                  micro_field(i,j,k,idsa) = dsasum ! Average dust surface area within the gridbox
                  micro_field(i,j,k,ibsa) = bsasum ! Average Snowmax surface area within the gridbox
                  micro_field(i,j,k,ipsa) = psasum ! Average pollen surface area within the gridbox
                ENDIF
                IF (switch_ice) THEN
                  IF (cldmass_total .gt. 0.0) fc_prof(k) = fc_prof(k) + 1 !Count of columns with cloud droplets
                  IF (rainmass_total .gt. 0.0) fr_prof(k) = fr_prof(k) + 1 !Count of columns with cloud droplet
                  IF (prismass_total .gt. 0.0) fp_prof(k) = fp_prof(k) + 1 !Count of columns with 'pristine' ice
                  IF (snowmass_total .gt. 0.0) fs_prof(k) = fs_prof(k) + 1 !Count of columns with snow (aggregates)
                  IF (graupelmass_total .gt. 0.0) fg_prof(k) = fg_prof(k) + 1 !Count of columns with graupel
                ENDIF

                micro_field(i,j,k,itmean) = micro_field(i,j,k,itmean) / SUM( particle(n_start(i,j,k):n_end(i,j,k))%ivar(iA) )
                micro_field(i,j,k,iltmean) = micro_field(i,j,k,iltmean) / SUM(particle(n_start(i,j,k):n_end(i,j,k))%ivar(iA), mask = particle(n_start(i,j,k):n_end(i,j,k))%ivar(iTYPE) .gt. 0)
                micro_field(i,j,k,iitmean) = micro_field(i,j,k,iitmean) / SUM(particle(n_start(i,j,k):n_end(i,j,k))%ivar(iA), mask = particle(n_start(i,j,k):n_end(i,j,k))%ivar(iTYPE) .lt. 0)
!
!--             Freezing temperatures
                IF (switch_ice) THEN
                  micro_field(i,j,k,itfrz) = micro_field(i,j,k,itfrz) / SUM( particle(n_start(i,j,k):n_end(i,j,k))%var(iM))
                  micro_field(i,j,k,itfrzc) = micro_field(i,j,k,itfrzc) / cldmass_total
                  micro_field(i,j,k,itfrzr) = micro_field(i,j,k,itfrzr) / rainmass_total
                  micro_field(i,j,k,itfrzi) = micro_field(i,j,k,itfrzi) / prismass_total
                  micro_field(i,j,k,itfrzs) = micro_field(i,j,k,itfrzs) / snowmass_total
                  micro_field(i,j,k,itfrzg) = micro_field(i,j,k,itfrzg) / graupelmass_total
                ENDiF

!
!--             Effective radius standard deviation
                IF ( r2_ijk .GT. 0.0 )  THEN
                   micro_field(i,j,k,iRsd) = SQRT( MAX( r4_ijk - 2.0 * r3_ijk / r2_ijk * r3_ijk + ( r3_ijk / r2_ijk )**2 * r2_ijk, 0.0 ) / r2_ijk )
                ENDIF
!
!--             Mean and standard deviation of radius profiles
                IF ( r0_ijk .GT. 0.0 )  THEN
                   r_prof(k)    = r_prof(k)    + r1_ijk / r0_ijk  ! droplet radius mean
                   rsd_prof(k)  = rsd_prof(k)  + SQRT( MAX( r2_ijk / r0_ijk - ( r1_ijk / r0_ijk )**2, 0.0 ) ) ! droplet radius standard deviation
                   ref_prof(k)  = ref_prof(k)  + r3_ijk / r2_ijk ! droplet effective radius
                   refs_prof(k) = refs_prof(k) + SQRT( MAX( r4_ijk - 2.0 * r3_ijk / r2_ijk * r3_ijk + ( r3_ijk / r2_ijk )**2 * r2_ijk, 0.0 ) / r2_ijk ) ! droplet effective radius standard deviation
                   fl_prof(k)   = fl_prof(k)   + 1.0  ! droplet count
                ENDIF

!
!--             Mean of aerosol radius and Koehler equilibrium radius
                IF (av0_ijk .GT. 0.0 ) THEN
                  micro_field(i,j,k,isr) = av1_ijk/av0_ijk
                  micro_field(i,j,k,ism) = am_ijk
                  ar_prof(k) = ar_prof(k) + av1_ijk/av0_ijk
                  am_prof(k) = am_prof(k) + am_ijk
                  fa_prof(k) = fa_prof(k) + 1.0 !Aerosol count
                ENDIF

!
!--             Ice mean and standard deviation of radius profiles
                IF ( ir0_ijk .GT. 0.0 )  THEN
                   iPHI_prof(k)  = iPHI_prof(k)  + PHI / ir0_ijk ! ice aspect ratio
                   Nmono_prof(k) = Nmono_prof(k) + Nmono / ir0_ijk ! number of monomers
                   iR_prof(k)    = iR_prof(k)    + ir1_ijk / ir0_ijk ! ice radius mean
                   iaavg_prof(k) = iaavg_prof(k) + ias_sum / ir0_ijk ! number-weighted mean a-axis length
                   icavg_prof(k) = icavg_prof(k) + ics_sum / ir0_ijk ! number-weighted mean c-axis length
                   iRSD_prof(k)  = iRSD_prof(k)  + SQRT( MAX( ir2_ijk / ir0_ijk - ( ir1_ijk / ir0_ijk )**2, 0.0 ) ) ! ice radius standard deviation
                   iREF_prof(k)  = iREF_prof(k)  + ir3_ijk / ir2_ijk ! ice effective radius
                   iREFS_prof(k) = iREFS_prof(k) + SQRT( MAX( ir4_ijk - 2.0 * ir3_ijk / ir2_ijk * ir3_ijk + ( ir3_ijk / ir2_ijk )**2 * ir2_ijk, 0.0 ) / ir2_ijk ) ! ice effective radius standard deviation
                   rho_app_prof(k) = rho_mom1 / icemass_total ! Mass-weighted mean ice apparent density
                   fi_prof(k)    = fi_prof(k)    + 1.0 ! ice count
                   micro_field(i,j,k,iphi) = phi_mom1 / phi_mom0 ! Capacitance-weighted aspect ratio of ice in the gridbox
                   micro_field(i,j,k,iaavg) = ias_sum / ir0_ijk ! Average a-axis length of ice crystals in the gridbox
                   micro_field(i,j,k,icavg) = ics_sum / ir0_ijk ! Average c-axis length of ice crystals in the gridbox
                   micro_field(i,j,k,irmf) = rimemass/icemass_total ! Fraction of ice mass which is due to riming
                   micro_field(i,j,k,iiref) = ir3_ijk / ir2_ijk ! Ice effective radius
                   micro_field(i,j,k,irhoapp) = rho_mom1 / icemass_total ! Mass-weighted apparent density of ice
                ENDIF

             ENDDO
          ENDDO

          micro_field(:,:,k,ina)   = micro_field(:,:,k,ina)   / ( dV(k) * rho_air(k) )
          micro_field(:,:,k,inc)   = micro_field(:,:,k,inc)   / ( dV(k) * rho_air(k) )
          micro_field(:,:,k,inr)   = micro_field(:,:,k,inr)   / ( dV(k) * rho_air(k) )
          micro_field(:,:,k,ini)   = micro_field(:,:,k,ini)   / ( dV(k) * rho_air(k) )
          micro_field(:,:,k,ins)   = micro_field(:,:,k,ins)   / ( dV(k) * rho_air(k) )
          micro_field(:,:,k,ing)   = micro_field(:,:,k,ing)   / ( dV(k) * rho_air(k) )

          micro_field(:,:,k,ipl)   = micro_field(:,:,k,ipl)   / dV(k)
          micro_field(:,:,k,ipi)   = micro_field(:,:,k,ipi)   / dV(k)
!
!--       Sprayed aerosol quantites
          micro_field(:,:,k,iasp)  = micro_field(:,:,k,iasp)  / ( dV(k) * rho_air(k) )
          micro_field(:,:,k,icsp)  = micro_field(:,:,k,icsp)  / ( dV(k) * rho_air(k) )

          micro_field(:,:,k,iqa)   = micro_field(:,:,k,iqa)   / ( dV(k) * rho_air(k) )
          micro_field(:,:,k,iqc)   = micro_field(:,:,k,iqc)   / ( dV(k) * rho_air(k) )
          micro_field(:,:,k,iqr)   = micro_field(:,:,k,iqr)   / ( dV(k) * rho_air(k) )
          micro_field(:,:,k,iqi)   = micro_field(:,:,k,iqi)   / ( dV(k) * rho_air(k) )
          micro_field(:,:,k,iqs)   = micro_field(:,:,k,iqs)   / ( dV(k) * rho_air(k) )
          micro_field(:,:,k,iqg)   = micro_field(:,:,k,iqg)   / ( dV(k) * rho_air(k) )

          micro_field(:,:,k,ininc) = micro_field(:,:,k,ininc) / ( dV(k) * rho_air(k) )

          micro_field(:,:,k,iACT)  = micro_field(:,:,k,iACT) / ( dV(k) * rho_air(k) * dtn_LCM )
          micro_field(:,:,k,iDEA)  = micro_field(:,:,k,iDEA) / ( dV(k) * rho_air(k) * dtn_LCM )

          micro_field(:,:,k,inuc)  = micro_field(:,:,k,inuc) / ( dv(K) * rho_air(k) * dtn_LCM )
          micro_field(:,:,k,iimm)  = micro_field(:,:,k,iimm) / ( dv(K) * rho_air(k) * dtn_LCM )
          micro_field(:,:,k,ihzn)  = micro_field(:,:,k,ihzn) / ( dv(K) * rho_air(k) * dtn_LCM )
          micro_field(:,:,k,ihom)  = micro_field(:,:,k,ihom) / ( dv(K) * rho_air(k) * dtn_LCM )

          micro_field(:,:,k,irim)  = micro_field(:,:,k,irim) / ( dv(K) * rho_air(k) * dtn_LCM )
          micro_field(:,:,k,ifrz)  = micro_field(:,:,k,ifrz) / ( dv(K) * rho_air(k) * dtn_LCM )
          micro_field(:,:,k,imel)  = micro_field(:,:,k,imel) / ( dv(K) * rho_air(k) * dtn_LCM )

          micro_field(:,:,k,icon)  = micro_field(:,:,k,icon) / ( dv(K) * rho_air(k) * dtn_LCM )
          micro_field(:,:,k,ievp)  = micro_field(:,:,k,ievp) / ( dv(K) * rho_air(k) * dtn_LCM )
          micro_field(:,:,k,idep)  = micro_field(:,:,k,idep) / ( dv(K) * rho_air(k) * dtn_LCM )
          micro_field(:,:,k,isub)  = micro_field(:,:,k,isub) / ( dv(K) * rho_air(k) * dtn_LCM )

          micro_field(:,:,k,icoal) = micro_field(:,:,k,icoal) / ( dv(K) * rho_air(k) * dtn_LCM )
          micro_field(:,:,k,ibreak) = micro_field(:,:,k, ibreak) / ( dv(K) * rho_air(k) * dtn_LCM )
          micro_field(:,:,k,iagg) = micro_field(:,:,k,iagg) / ( dv(K) * rho_air(K) * dtn_LCM )
          micro_field(:,:,k,iicebreak) = micro_field(:,:,k,iicebreak) / ( dv(K) * rho_air(K) * dtn_LCM )
          micro_field(:,:,k,ids) = micro_field(:,:,k,ids) / ( dv(K) * rho_air(K) * dtn_LCM )
          micro_field(:,:,k,irs) = micro_field(:,:,k,irs) / ( dv(K) * rho_air(K) * dtn_LCM )
!
!--       Vertically integrated droplet number concentration
#ifdef SAMTAK
          ncl_xy(1:nx_LCM,1:ny_LCM) = ncl_xy(1:nx_LCM,1:ny_LCM) + micro_field(1:nx_LCM,1:ny_LCM,k,inc) * ( gpzs(k+1) - gpzs(k) ) * dtfactor * rho_air(k)
          npl_xy(1:nx_LCM,1:ny_LCM) = npl_xy(1:nx_LCM,1:ny_LCM) + micro_field(1:nx_LCM,1:ny_LCM,k,inr) * ( gpzs(k+1) - gpzs(k) ) * dtfactor * rho_air(k)
#endif /*SAMTAK*/
!
!--       Profiles of aerosol, cloud droplet, and rain drop concentration as well as freezing temperature
          nac_prof(k)  = SUM(micro_field(:,:,k,ina))
          ncc_prof(k)  = SUM(micro_field(:,:,k,inc))
          nrc_prof(k)  = SUM(micro_field(:,:,k,inr))
          IF (switch_ice) THEN
            rmf_prof(k)  = SUM(micro_field(:,:,k,irmf))
            nic_prof(k)  = SUM(micro_field(:,:,k,ini))
            nsc_prof(k)  = SUM(micro_field(:,:,k,ins))
            ngc_prof(k)  = SUM(micro_field(:,:,k,ing))
            ninc_prof(k) = SUM(micro_field(:,:,k,ininc))
            tfrz_prof(k) = SUM(micro_field(:,:,k,itfrz))
            tfrzc_prof(k) = SUM(micro_field(:,:,k,itfrzc))
            tfrzr_prof(k) = SUM(micro_field(:,:,k,itfrzr))
            tfrzi_prof(k) = SUM(micro_field(:,:,k,itfrzi))
            tfrzs_prof(k) = SUM(micro_field(:,:,k,itfrzs))
            tfrzg_prof(k) = SUM(micro_field(:,:,k,itfrzg))
          ENDIF

       ENDDO

       CALL MPI_ALLREDUCE(fl_prof,dummy_prof,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       fl_prof = dummy_prof + 1.0E-20
       CALL MPI_ALLREDUCE(fi_prof,dummy_prof,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       fi_prof = dummy_prof + 1.0E-20
       CALL MPI_ALLREDUCE(fa_prof,dummy_prof,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       fa_prof = fa_prof + 1.0E-20
       CALL MPI_ALLREDUCE(fin_prof,dummy_prof,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       fin_prof = fin_prof + 1.0E-20
       CALL MPI_ALLREDUCE(fc_prof,dummy_prof,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       fc_prof = fc_prof + 1.0E-20
       CALL MPI_ALLREDUCE(fr_prof,dummy_prof,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       fr_prof = fr_prof + 1.0E-20
       CALL MPI_ALLREDUCE(fp_prof,dummy_prof,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       fp_prof = fp_prof + 1.0E-20
       CALL MPI_ALLREDUCE(fs_prof,dummy_prof,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       fs_prof = fs_prof + 1.0E-20
       CALL MPI_ALLREDUCE(fg_prof,dummy_prof,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       fg_prof = fg_prof + 1.0E-20

       DO k = 1, nz_LCM

          r_prof(k)     = r_prof(k)     / fl_prof(k) * nsubdomains_x * nsubdomains_y
          rsd_prof(k)   = rsd_prof(k)   / fl_prof(k) * nsubdomains_x * nsubdomains_y
          ref_prof(k)   = ref_prof(k)   / fl_prof(k) * nsubdomains_x * nsubdomains_y
          refs_prof(k)  = refs_prof(k)  / fl_prof(k) * nsubdomains_x * nsubdomains_y

          iR_prof(k)    = iR_prof(k)    / fi_prof(k) * nsubdomains_x * nsubdomains_y
          iRSD_prof(k)  = iRSD_prof(k)  / fi_prof(k) * nsubdomains_x * nsubdomains_y
          iREF_prof(k)  = iREF_prof(k)  / fi_prof(k) * nsubdomains_x * nsubdomains_y
          iREFS_prof(k) = iREFS_prof(k) / fi_prof(k) * nsubdomains_x * nsubdomains_y

          iaavg_prof(k) = iaavg_prof(k) / fi_prof(k) * nsubdomains_x * nsubdomains_y
          icavg_prof(k) = icavg_prof(k) / fi_prof(k) * nsubdomains_x * nsubdomains_y
          rho_app_prof(k) = rho_app_prof(k) / fi_prof(k) * nsubdomains_x * nsubdomains_y

          iPHI_prof(k)  = iPHI_prof(k)  / fi_prof(k) * nsubdomains_x * nsubdomains_y
          rmf_prof(k)   = rmf_prof(k) / fi_prof(k) * nsubdomains_x * nsubdomains_y

          Nmono_prof(k) = Nmono_prof(k) / fi_prof(k) * nsubdomains_x * nsubdomains_y

          ar_prof(k) = ar_prof(k) / fa_prof(k) * nsubdomains_x * nsubdomains_y

          ssa_prof(k) = ssa_prof(k) / fin_prof(k) * nsubdomains_x * nsubdomains_y
          dsa_prof(k) = dsa_prof(k) / fin_prof(k) * nsubdomains_x * nsubdomains_y
          bsa_prof(k) = bsa_prof(k) / fin_prof(k) * nsubdomains_x * nsubdomains_y
          psa_prof(k) = psa_prof(k) / fin_prof(k) * nsubdomains_x * nsubdomains_y

          tfrzc_prof(k) = tfrzc_prof(k) / fc_prof(k) * nsubdomains_x * nsubdomains_y
          tfrzr_prof(k) = tfrzr_prof(k) / fr_prof(k) * nsubdomains_x * nsubdomains_y
          tfrzi_prof(k) = tfrzi_prof(k) / fp_prof(k) * nsubdomains_x * nsubdomains_y
          tfrzs_prof(k) = tfrzs_prof(k) / fs_prof(k) * nsubdomains_x * nsubdomains_y
          tfrzg_prof(k) = tfrzg_prof(k) / fg_prof(k) * nsubdomains_x * nsubdomains_y

       ENDDO

#ifdef SAMTAK
       opdqcl = opdqcl * 2.0 * pi / ( dx_LCM * dy_LCM )  ! 2.0 is the assumed extinction coefficient, see Twomey (1977, JAS)
       opdqpl = opdqpl * 2.0 * pi / ( dx_LCM * dy_LCM )

       qclopd_xy = qclopd_xy + opdqcl * dtfactor
       qplopd_xy = qplopd_xy + opdqpl * dtfactor
#endif /*SAMTAK*/

       DO k = 2, nz_LCM
          tau_prof(k) = tau_prof(k) + tau_prof(k-1) !optical thickness with Mie scattering
       ENDDO
!
!--    Store particle spectra
       IF ( MOD(nstep,spec_nsave) .EQ. 0  .AND.  nstep .GE. spec_nsavestart  .AND.  icycle .EQ. ncycle ) THEN
          IF ( rank .EQ. 0 )  THEN
             PRINT*, "****** LCM: output of spectra ******"
          ENDIF
          IF ( horiz_average_spec )  THEN
             IF ( switch_ice )  THEN
                CALL MPI_ALLREDUCE(spec_ice(1,1,:,:),spec_global,nzl*n_bins_spec,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,ierr)
                DO n = 1, n_bins_spec
                   spec_global(:,n) = spec_global(:,n) / ( rr_spec(n) - rl_spec(n) )
                ENDDO
                do k = 1, nzl
                  dv_spec     = dz_LCM(k) * xwidth * ywidth
                  spec_global(k,:) = spec_global(k,:) / dv_spec
                end do
                IF ( masterproc )  THEN
                  OPEN(UNIT=98,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_ice',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND')
                  DO k = 1, nzl
                    WRITE(98,format_spec) 1, 1, k, nstep, ( spec_global(k,n), n = 1, n_bins_spec )
                    FLUSH(98)
                  ENDDO
                  CLOSE(98)
                ENDIF
                spec_global = 0.0
             ENDIF

             CALL MPI_ALLREDUCE(spec_liq(1,1,:,:),spec_global,nzl*n_bins_spec,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,ierr)
             DO n = 1, n_bins_spec
                spec_global(:,n) = spec_global(:,n) / ( rr_spec(n) - rl_spec(n) )
             ENDDO
             IF ( masterproc )  THEN  
               OPEN(UNIT=99,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_liq',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND')
               DO k = 1, nzl
                dv_spec     = dz_LCM(k) * xwidth * ywidth
                spec_global(k,:) = spec_global(k,:) / dv_spec
                WRITE(99,format_spec) 1, 1, k, nstep, ( spec_global(k,n), n = 1, n_bins_spec )
                FLUSH(99)
               ENDDO
               CLOSE(99)
             ENDIF
             spec_global = 0.0

             IF (switch_spec_aero) THEN
                 ! do aerosol spectra
                 CALL MPI_ALLREDUCE(spec_aero(1,1,:,:),spec_global,nzl*n_bins_spec,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,ierr)
                 DO n = 1, n_bins_spec
                    spec_global(:,n) = spec_global(:,n) / ( rr_spec(n) - rl_spec(n) )
                 ENDDO
                 DO k=1,nzl
                   dV_spec          = dz_LCM(k) * xwidth * ywidth
                   spec_global(k,:) = spec_global(k,:) / dV_spec
                 ENDDO
                 IF ( masterproc )  THEN
                    OPEN(UNIT=97,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_aero',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND')
                    DO k = 1, nzl
                       WRITE(97,format_spec) 1, 1, k, nstep, ( spec_global(k,n), n = 1, n_bins_spec )
                       FLUSH(97)
                    ENDDO
                    CLOSE(97)
                 ENDIF
                 spec_global = 0.0
             END IF

#ifdef SAMCHEM
             IF ( switch_spec_pH ) THEN
                 CALL MPI_ALLREDUCE(spec_pH(1,1,:,:),spec_global,nzl*n_bins_spec,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,ierr)
                 DO n = 1, n_bins_spec
                    spec_global(:,n) = spec_global(:,n) / ( pHr_spec(n) - pHl_spec(n) )
                 ENDDO
                 DO k=1,nzl
                   dV_spec          = dz_LCM(k) * xwidth * ywidth
                   spec_global(k,:) = spec_global(k,:) / dV_spec
                 ENDDO
                 IF ( masterproc )  THEN
                    OPEN(UNIT=96,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_pH',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND')
                    DO k = 1, nzl
                       WRITE(96,format_spec) 1, 1, k, nstep, ( spec_global(k,n), n = 1, n_bins_spec )
                       FLUSH(96)
                    ENDDO
                    CLOSE(96)
                 ENDIF
                 spec_global = 0.0
             END IF
#endif /*SAMCHEM*/
          ELSE
            IF ( masterproc )  THEN
               WRITE(nstep_char,'(I9)') nstep
               nstep_char = ADJUSTL(nstep_char)
               OPEN(UNIT=99,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_liq.'//TRIM(nstep_char),STATUS='NEW',FORM='FORMATTED')
               IF ( switch_ice )  THEN
                  OPEN(UNIT=98,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_ice.'//TRIM(nstep_char),STATUS='NEW',FORM='FORMATTED')
               ENDIF
               IF ( switch_spec_aero ) OPEN(UNIT=97,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_aero.'//TRIM(nstep_char),STATUS='NEW',FORM='FORMATTED')
#ifdef SAMCHEM
               IF ( switch_spec_pH   ) OPEN(UNIT=96,FILE='./OUT_USER/'//TRIM(case)//'_'//TRIM(caseid)//'.spec_pH.'//TRIM(nstep_char),STATUS='NEW',FORM='FORMATTED')
#endif /*SAMCHEM*/
            ENDIF
             DO i = 1, nx_LCM*nsubdomains_x
                DO j = 1, ny_LCM*nsubdomains_y
                   IF ( switch_ice )  THEN
                      CALL MPI_ALLREDUCE(spec_ice(i,j,:,:),spec_global,nzl*n_bins_spec,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,ierr)
                      DO n = 1, n_bins_spec
                         spec_global(:,n) = spec_global(:,n) / ( rr_spec(n) - rl_spec(n) )
                      ENDDO
                     IF ( masterproc )  THEN
                        DO k = 1, nzl
                            dv_spec     = dz_LCM(k) * dx_LCM * dy_LCM
                            spec_global(k,:) = spec_global(k,:) / dv_spec
                            WRITE(98,format_spec) i, j, k, nstep, ( spec_global(k,n), n = 1, n_bins_spec )
                            FLUSH(98)
                        ENDDO
                     ENDIF
                      spec_global = 0.0
                   ENDIF
                   CALL MPI_ALLREDUCE(spec_liq(i,j,:,:),spec_global,nzl*n_bins_spec,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,ierr)
                   DO n = 1, n_bins_spec
                      spec_global(:,n) = spec_global(:,n) / ( rr_spec(n) - rl_spec(n) )
                   ENDDO
                   IF ( masterproc )  THEN
                     DO k=1, nzl
                         dv_spec     = dz_LCM(k) * dx_LCM * dy_LCM
                         spec_global(k,:) = spec_global(k,:) / dv_spec
                         WRITE(99,format_spec) i, j, k, nstep, ( spec_global(k,n), n = 1, n_bins_spec )
                         FLUSH(99)
                     ENDDO
                   ENDIF
                   spec_global = 0.0
                   IF (switch_spec_aero) THEN
                       ! do aerosol spectra
                       CALL MPI_ALLREDUCE(spec_aero(i,j,:,:),spec_global,nzl*n_bins_spec,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,ierr)
                       DO n = 1, n_bins_spec
                          spec_global(:,n) = spec_global(:,n) / ( rr_spec(n) - rl_spec(n) )
                       ENDDO
                       IF ( masterproc )  THEN
                          DO k = 1, nzl
                             dV_spec     = dz_LCM(k) * dx_LCM * dy_LCM
                             spec_global(k,:) = spec_global(k,:) / dV_spec
                             WRITE(97,format_spec) i, j, k, nstep, ( spec_global(k,n), n = 1, n_bins_spec )
                             FLUSH(97)
                          ENDDO
                       ENDIF
                       spec_global = 0.0
                   ENDIF
#ifdef SAMCHEM
                   IF ( switch_spec_pH ) THEN
                       CALL MPI_ALLREDUCE(spec_pH(i,j,:,:),spec_global,nzl*n_bins_spec,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,ierr)
                       DO n = 1, n_bins_spec
                          spec_global(:,n) = spec_global(:,n) / ( pHr_spec(n) - pHl_spec(n) )
                       ENDDO
                       IF ( masterproc )  THEN
                          DO k = 1, nzl
                             dV_spec     = dz_LCM(k) * dx_LCM * dy_LCM
                             spec_global(k,:) = spec_global(k,:) / dV_spec
                             WRITE(96,format_spec) i, j, k, nstep, ( spec_global(k,n), n = 1, n_bins_spec )
                             FLUSH(96)
                          ENDDO
                       ENDIF
                       spec_global = 0.0
                   ENDIF
#endif /*SAMCHEM*/
                ENDDO
             ENDDO
             IF ( masterproc )  THEN
               CLOSE(99)
               IF ( switch_ice ) THEN
                  CLOSE(98)
               ENDIF
#ifdef SAMCHEM
               IF (switch_spec_pH)   CLOSE(96)
#endif /*SAMCHEM*/
               IF (switch_spec_aero) CLOSE(97)
            ENDIF
          ENDIF
       ENDIF
!
!--    Get profiles of the mean and standard deviation of the plume-detecting scalar
       IF ( switch_cloud_tracking )  THEN

          DO k = 1, nz_LCM
             std_SCLR_local(k)  = SUM( micro_field(1:nx_LCM,1:ny_LCM,k,iSCLR)**2 )
             mean_SCLR_local(k) = SUM( micro_field(1:nx_LCM,1:ny_LCM,k,iSCLR) )
             T_local(k)         = SUM( tabs_LCM(1:nx_LCM,1:ny_LCM,k) )
             qv_local(k)        = SUM( micro_field(1:nx_LCM,1:ny_LCM,k,imv) )
             ql_local(k)        = SUM( micro_field(1:nx_LCM,1:ny_LCM,k,iqc) + micro_field(1:nx_LCM,1:ny_LCM,k,iqr) )
          ENDDO

          CALL MPI_ALLREDUCE(std_SCLR_local,std_SCLR_prof,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
          CALL MPI_ALLREDUCE(mean_SCLR_local,mean_SCLR_prof,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
          CALL MPI_ALLREDUCE(T_local,T_global,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
          CALL MPI_ALLREDUCE(qv_local,qv_global,nz_LCM,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

          mean_SCLR_prof = mean_SCLR_prof / FLOAT( nx_gl * ny_gl )
          std_SCLR_prof  = std_SCLR_prof / FLOAT( nx_gl * ny_gl )

          std_SCLR_prof = SQRT( MAX( std_SCLR_prof - mean_SCLR_prof**2, 0.0 ) )

          T_global  = T_global  / FLOAT( nx_gl * ny_gl )
          qv_global = qv_global / FLOAT( nx_gl * ny_gl )
          ql_global = ql_global / FLOAT( nx_gl * ny_gl )

          B_global = T_global * ( 1.0 + epsv * qv_global - ql_global )

          micro_field(:,:,:,iPlume)     = 0.0
          micro_field(:,:,:,iCondensed) = 0.0
          micro_field(:,:,:,iCore)      = 0.0
          DO k = 1, nz_LCM

             min_SCLR = SUM(std_SCLR_prof(1:k)) / FLOAT( k ) * 0.05 !+ 1.0E-20

             DO j = 1, ny_LCM
                DO i = 1, nx_LCM

                   IF ( micro_field(i,j,k,iqc) .GT. qnmin_LCM )  THEN
                      micro_field(i,j,k,iCondensed) = 1.0

                      B_local = tabs_LCM(i,j,k) * ( 1.0 + epsv * micro_field(i,j,k,imv) - ( micro_field(i,j,k,iqc) + micro_field(i,j,k,iqr) ) )
                      IF ( w(i,j,k) .GT. 0.0  .AND.  B_local - B_global(k) .GT. 0.0 )  THEN
                         micro_field(i,j,k,iCore) = 1.0
                      ENDIF

                   ENDIF

                   IF ( micro_field(i,j,k,iSCLR) - mean_SCLR_prof(k) .GT. MAX( min_SCLR, std_SCLR_prof(k) ) )  THEN
                       micro_field(i,j,k,iPlume) = 1.0
                   ENDIF

                ENDDO
             ENDDO
          ENDDO

       ENDIF


    ENDIF


    ! not output-only variables (always diagnose)
    IF ( flag_LCM_init )  THEN
!
!--    Return LCM arrays to SAM
       qv(1:nx,1:ny,1:nzm)     = micro_field(1:nx_LCM,1:ny_LCM,1:nz_LCM,imv)
       qcl(1:nx,1:ny,1:nzm)    = micro_field(1:nx_LCM,1:ny_LCM,1:nz_LCM,iqa)   &
                                 + micro_field(1:nx_LCM,1:ny_LCM,1:nz_LCM,iqc)
       IF ( doprecip ) THEN
          qpl(1:nx,1:ny,1:nzm) = micro_field(1:nx_LCM,1:ny_LCM,1:nz_LCM,iqr)
       ELSE
          qpl(1:nx,1:ny,1:nzm) = 0.0
       ENDIF

       IF ( switch_ice )  THEN
          qci(1:nx,1:ny,1:nzm) = micro_field(1:nx_LCM,1:ny_LCM,1:nz_LCM,iqi)
          IF (doprecip) THEN
            qpi(1:nx,1:ny,1:nzm) =  micro_field(1:nx_LCM,1:ny_LCM,1:nz_LCM,iqs) + micro_field(1:nx_LCM,1:ny_LCM,1:nz_LCM,iqg)
          ENDIF
       ELSE
          qci(1:nx,1:ny,1:nzm) = 0.0 ! no ice
       ENDIF
!
!--    Set temperatures for next timestep
       DO k = 1, nzm
          micro_field(1:nx,1:ny,k,iTlcm)  = tabs_LCM(1:nx,1:ny,k) + gamaz(k)
          micro_field(1:nx,1:ny,k,iTllcm) = micro_field(1:nx,1:ny,k,iTlcm) - ( qcl(1:nx,1:ny,k) + qpl(1:nx,1:ny,k) ) * fac_cond - (qci(1:nx,1:ny,k) + qpi(1:nx,1:ny,k)) * fac_sub 
          t(1:nx,1:ny,k)                  = micro_field(1:nx,1:ny,k,iTllcm)
          tabs(1:nx,1:ny,k)               = tabs_LCM(1:nx,1:ny,k)
       ENDDO
       CALL task_exchange(micro_field(:,:,:,iTlcm),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4)
       CALL task_exchange(micro_field(:,:,:,iTllcm),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4)
       CALL task_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4)

    ENDIF



    CALL t_stopf ('micro_diagnose')

 END SUBROUTINE micro_diagnose

 SUBROUTINE micro_precip_fall()
!
!-- REQUIRED BY SAM. NOT USED HERE.
 END SUBROUTINE micro_precip_fall

 SUBROUTINE micro_print()
!
!-- This one is called when stepout() is called.
!-- REQUIRED BY SAM. NOT USED HERE.
 END SUBROUTINE micro_print

 REAL FUNCTION total_water()
!
!-- Function that computes total water in the domain. REQUIRED BY SAM.
    USE micro_vars, ONLY: particle, n_particles, iM, dx_LCM, dy_LCM, dz_LCM, rho_air, gpzs

    USE vars,       ONLY: nstep, nprint

    IMPLICIT NONE

    INTEGER k, m

    total_water = 0.0
    IF ( MOD( nstep, nprint ) /= 0 ) RETURN
    DO k = 1, nz_LCM
       total_water = total_water + SUM( micro_field(1:nx_LCM,1:ny_LCM,k,imv) ) * ( gpzs(k+1) - gpzs(k) ) * rho_air(k)
    ENDDO
    total_water = total_water + SUM( particle(1:n_particles)%var(iM) ) / ( dx_LCM * dy_LCM )

 END FUNCTION total_water

 FUNCTION Get_reffc()
!
!-- Function to compute liquid water effective radius (in microns).
!-- REQUIRED BY SAM.
    USE micro_vars, ONLY: rho_liq, particle, n_start, n_end, iM, iA, iN,       &
                          rho_aero, vanthoff_aero, molecular_weight_water,     &
                          molecular_weight_aero, sigma, tabs_LCM, iTYPE,       &
                          get_r_particle
    USE params,     ONLY: pi, rv

    IMPLICIT NONE

    REAL, DIMENSION(nx,ny,nzm) :: Get_reffc
    REAL                       :: numerator, denominator, radius
    INTEGER                    :: i, j, k, n

    Get_reffc(:,:,:) = 0.0
    DO k = 1, nz_LCM
       DO j = 1, ny_LCM
          DO i = 1, nx_LCM
             IF ( n_end(i,j,k) - n_start(i,j,k) .LT. 0 )  CYCLE !if there are no particles inside that gridbox, do nothing

             numerator   = 0
             denominator = 0
             DO n = n_start(i,j,k), n_end (i,j,k)
                IF ( particle(n)%ivar(iTYPE) .GT. 0 )  THEN  ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
                   radius      = get_r_particle(n)
                   numerator   = numerator   + radius**3 * particle(n)%ivar(iA)
                   denominator = denominator + radius**2 * particle(n)%ivar(iA)
                ENDIF
             ENDDO

             IF ( numerator > 0.0 .AND. denominator > 0.0 ) THEN
                Get_reffc(i,j,k) = 1.0E6 * numerator / denominator
             ENDIF

          ENDDO
       ENDDO
    ENDDO

 END FUNCTION Get_reffc

 FUNCTION Get_reffi() ! ice
!
!-- Function that computes ice water effective radius.
!-- REQUIRED BY SAM.
    USE micro_vars, ONLY: rho_liq, particle, n_start, n_end, iM, iA, iN,       &
                          rho_aero, vanthoff_aero, molecular_weight_water,     &
                          iTYPE, iAXS, iCXS
    USE params,     ONLY: pi, rv

    IMPLICIT NONE

    REAL, DIMENSION(nx,ny,nzm) :: Get_reffi
    REAL                       :: numerator, denominator, radius
    INTEGER                    :: i, j, k, n

    Get_reffi(:,:,:) = 0.0
    DO k = 1, nz_LCM
       DO j = 1, ny_LCM
          DO i = 1, nx_LCM
             IF ( n_end(i,j,k) - n_start(i,j,k) .LT. 0 )  CYCLE !if there are no particles inside that gridbox, do nothing

             numerator   = 0
             denominator = 0
             DO n = n_start(i,j,k), n_end (i,j,k)
                IF ( particle(n)%ivar(iTYPE) .LT. 0 )  THEN  ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
                   radius      = MAX( particle(n)%var(iAXS), particle(n)%var(iCXS) )
                   numerator   = numerator   + radius**3 * particle(n)%ivar(iA)
                   denominator = denominator + radius**2 * particle(n)%ivar(iA)
                ENDIF
             ENDDO

             IF ( numerator > 0.0 .AND. denominator > 0.0 ) THEN
                Get_reffi(i,j,k) = 1.0E6 * numerator / denominator
             ENDIF

          ENDDO
       ENDDO
    ENDDO

 END FUNCTION Get_reffi

 FUNCTION Get_nca()
!
!-- This function is required when doradaerosimple = .TRUE.
!-- This function returns aerosol number concentration in #/cm3
!-- Make blank function if aerosol number concentration is unavailable.
!-- This is equivalent to doradaerosimple = .FALSE.
!-- REQUIRED BY SAM.
    USE vars, ONLY: rho

    IMPLICIT NONE

    REAL, DIMENSION(nx,ny,nzm) :: Get_nca
    INTEGER                    :: k

    DO k = 1, nzm
       Get_nca(1:nx,1:ny,k) = micro_field(1:nx_LCM,1:ny_LCM,k,ina) * rho_air(k)
    ENDDO

 END FUNCTION Get_nca

 FUNCTION micro_scheme_name()
!
!-- Return the scheme name. REQUIRED BY SAM.
    IMPLICIT NONE

    CHARACTER(LEN=32) :: micro_scheme_name

    micro_scheme_name = "LAGRANGIAN"

 END FUNCTION micro_scheme_name

 SUBROUTINE micro_hbuf_init(namelist, deflist, unitlist, status, average_type, count, microcount)
!
!-- Initialize the list of microphysics statistics that will be outputted to
!-- *.stat statistics file. REQUIRED BY SAM.
    USE micro_vars, ONLY: switch_LEM_micro

    USE grid,       ONLY: masterproc

    IMPLICIT NONE
!
!-- Input/Output
    CHARACTER(*) :: namelist(*), deflist(*), unitlist(*)
    INTEGER      :: status(*), average_type(*), count, microcount
!
!-- Local
    CHARACTER (LEN=8)  :: name
    CHARACTER (LEN=80) :: longname
    CHARACTER (LEN=10) :: units
    INTEGER            :: n, ncond

    microcount = 0

    name = 'lag_NSD'
    longname = 'number of super-droplets in grid-box'
    units = '-'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_NAC'
    longname = 'aerosol concentration'
    units = 'mg-1'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_NCC'
    longname = 'cloud droplet concentration'
    units = 'mg-1'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_NRC'
    longname = 'rain drop concentration'
    units = 'mg-1'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_CND'
    longname = 'condensation rate'
    units = 'g/kg/s'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_EVP'
    longname = 'evaporation rate'
    units = 'g/kg/s'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_AUT'
    longname = 'autoconversion rate'
    units = 'g/kg/s'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_ACC'
    longname = 'accretion rate'
    units = 'g/kg/s'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_pQLF'
    longname = 'precip. (sedimentation) downward liquid water flux'
    units = 'W/m2'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_rQLF'
    longname = 'resolved downward liquid water flux'
    units = 'W/m2'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_sQLF'
    longname = 'SGS downward liquid water flux'
    units = 'W/m2'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_ACT'
    longname = 'particle activation rate'
    units = '1/mg/s'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_DEA'
    longname = 'particle deactivation rate'
    units = '1/mg/s'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_TAU'
    longname = 'optical thickness (aerosol + cloud + rain), including Mie'
    units = '-'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_REF'
    longname = 'effective radius'
    units = 'microns'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_RFS'
    longname = 'effective radius standard deviation'
    units = 'microns'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_R'
    longname = 'mean radius'
    units = 'microns'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_RSD'
    longname = 'radius standard deviation'
    units = 'microns'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_QA'
    longname = 'liquid water on aerosols'
    units = 'g/kg'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_QC'
    longname = 'liquid water on cloud droplets'
    units = 'g/kg'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_QR'
    longname = 'liquid water on rain drops'
    units = 'g/kg'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_AR'
    longname = 'mean aerosol radius'
    units = 'nm'
    call add_to_namelist( count, microcount, name, longname, units, 0)

    name = 'lag_AM'
    longname = 'mean aerosol mass'
    units = 'ng'
    call add_to_namelist( count, microcount, name, longname, units, 0)

    IF ( switch_ice )  THEN

       name = 'lag_NIC'
       longname = 'ice concentration'
       units = 'mg-1'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_QI'
       longname = 'ice water'
       units = 'g/kg'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_NSC'
       longname = 'snow concentration'
       units = 'mg-1'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_QS'
       longname = 'snow water'
       units = 'g/kg'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_NGC'
       longname = 'graupel concentration'
       units = 'mg-1'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_QG'
       longname = 'graupel water'
       units = 'g/kg'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_NUC'
       longname = 'particle nucleation rate'
       units = '1/mg/s'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = "lag_IMM"
       longname = "immersion freezing rate"
       units = "1/mg/s"
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = "lag_NSUB"
       longname = "rate of particle conversion from sublimation"
       units = "1/mg/s"
       CALL add_to_namelist(count, microcount, name, longname, units, 0)

      ! name = "lag_HZN"
      ! longname = "deposition nucleation rate"
      ! units = "1/mg/s"
      ! CALL add_to_namelist(count, microcount, name, longname, units, 0 )

       name = "lag_HOM"
       longname = "Homogneous nucleation rate"
       units = "1/mg/s"
       CALL add_to_namelist(count, microcount, name, longname, units, 0 )

       name = 'lag_NML'
       longname = 'particle melting rate'
       units = '1/mg/s'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_DEP'
       longname = 'deposition rate'
       units = 'g/kg/s'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_SUB'
       longname = 'sublimation rate'
       units = 'g/kg/s'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_FRZ'
       longname = 'freezing rate'
       units = 'g/kg/s'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_MEL'
       longname = 'melting rate'
       units = 'g/kg/s'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_RIM'
       longname = 'riming rate'
       units = 'g/kg/s'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_iREF'
       longname = 'ice effective radius'
       units = 'microns'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_iRFS'
       longname = 'ice effective radius standard deviation'
       units = 'microns'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_iR'
       longname = 'ice mean radius'
       units = 'microns'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_iRSD'
       longname = 'ice radius standard deviation'
       units = 'microns'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_iAXS'
       longname = 'ice average a-axis length'
       units = 'microns'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_iCXS'
       longname = 'ice average c-axis length'
       units = 'microns'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_iPHI'
       longname = 'ice aspect ratio'
       units = '-'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_iRHO'
       longname = 'ice mass-weighted apparent density'
       units = 'kg/m3'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_iNmn'
       longname = 'number monomers'
       units = '-'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_iRH'
       longname = 'relative humidity ice'
       units = '-'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_iRMF'
       longname = 'rime mass fraction'
       units = '-'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_NINC'
       longname = 'ice nuclei concentration (part of lag_NAC, lag_NCC, lag_NRC)'
       units = 'mg-1'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = "lag_ITMX"
       longname = "ice particle mean cloud residence time"
       units = "s"
       CALL add_to_namelist(count, microcount, name, longname, units, 0)

       name = "lag_ITMN"
       longname = "ice particle max cloud residence time"
       units = "s"
       CALL add_to_namelist(count, microcount, name, longname, units, 0)

       name = "lag_TFRZ"
       longname = "mass-weighted particle freezing temperature"
       units = "K"
       CALL add_to_namelist(count, microcount, name, longname, units, 0)

       name = "lag_TFZC"
       longname = "mass-weighted cloud droplet freezing temperature"
       units = "K"
       CALL add_to_namelist(count, microcount, name, longname, units, 0)

       name = "lag_TFZR"
       longname = "mass-weighted rain freezing temperature"
       units = "K"
       CALL add_to_namelist(count, microcount, name, longname, units, 0)

       name = "lag_TFZI"
       longname = "mass-weighted cloud ice freezing temperature"
       units = "K"
       CALL add_to_namelist(count, microcount, name, longname, units, 0)

       name = "lag_TFZS"
       longname = "mass-weighted snow freezing temperature"
       units = "K"
       CALL add_to_namelist(count, microcount, name, longname, units, 0)

       name = "lag_TFZG"
       longname = "mass-weighted graupel freezing temperature"
       units = "K"
       CALL add_to_namelist(count, microcount, name, longname, units, 0)

       name = "lag_SSA"
       longname = "mean soot IN surface area per particle"
       units = "um2"
       call add_to_namelist(count, microcount, name, longname, units, 0)

       name = "lag_DSA"
       longname = "mean dust IN surface area per particle"
       units = "um2"
       call add_to_namelist(count, microcount, name, longname, units, 0)

       name = "lag_BSA"
       longname = "mean Snowmax IN surface area per particle"
       units = "um2"
       call add_to_namelist(count, microcount, name, longname, units, 0)

       name = "lag_PSA"
       longname = "mean pollen IN surface area per particle"
       units = "um2"
       call add_to_namelist(count, microcount, name, longname, units, 0)

    ENDIF

    IF ( switch_LEM_micro )  THEN

       name = 'lag_SPM'
       longname = 'perturbation supersaturation mean'
       units = '%'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_SPS'
       longname = 'perturbation supersaturation standard deviation'
       units = '%'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_TABS'
       longname = 'LEM absolute temperature'
       units = 'K'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

       name = 'lag_QV'
       longname = 'LEM water vapor mixing ratio'
       units = 'g/kg'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    ENDIF

    IF ( switch_particle_source )  THEN

       name = 'lag_NIN'
       longname = 'aerosol injection rate'
       units = '1/mg/s'
       CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    ENDIF

    name = 'lag_TMX'
    longname = 'maximum residence time'
    units = 's'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_TMN'
    longname = 'mean residence time'
    units = 's'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_LTMX'
    longname = 'liquid particle maximum residence time'
    units = 's'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_LTMN'
    longname = 'liquid particle mean residence time'
    units = 's'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    name = 'lag_ALB'
    longname = 'albedo (aerosol + cloud + rain), including Mie'
    units = '-'
    CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    IF (switch_collection) THEN

      name = "lag_COAL"
      longname = "coalescence rate"
      units = "g/kg/s"
      CALL add_to_namelist( count, microcount, name, longname, units, 0)

      name = "lag_BRK"
      longname = "breakup rate"
      units = "g/kg/s"
      CALL add_to_namelist( count, microcount, name, longname, units, 0)
      
      IF (switch_ice) THEN
         name = "lag_AGG"
         longname = "aggregation rate"
         units = "g/kg/s"
         CALL add_to_namelist( count, microcount, name, longname, units, 0)

         name = "lag_iBRK"
         longname = "ice collisional breakup rate"
         units = "1/kg/s"
         CALL add_to_namelist( count, microcount, name, longname, units, 0)

         name = "lag_DS"
         longname = "droplet shattering rate"
         units = "1/kg/s"
         CALL add_to_namelist( count, microcount, name, longname, units, 0)

         name = "lag_RS"
         longname = "rime-splintering rate"
         units = "1/kg/s"
         CALL add_to_namelist( count, microcount, name, longname, units, 0)
      ENDIF
   ENDIF

        name = 'QTFLUX'
        longname = 'Total ( resolved + SGS ) flux of total water'
        units = 'W/m2'
        CALL add_to_namelist( count, microcount, name, longname, units, 0 )

        name = 'QTFLUXS'
        longname = 'Subgrid flux of total water'
        units = 'W/m2'
        CALL add_to_namelist( count, microcount, name, longname, units, 0 )

    IF ( masterproc ) THEN
       WRITE(*,*) 'Added ', microcount, ' arrays to statistics for LAGRANGIAN microphysics'
    ENDIF

 END SUBROUTINE micro_hbuf_init

 SUBROUTINE micro_statistics()
!
!-- Collect microphysics history statistics (vertical profiles)
!-- Note that only the fields declared in micro_hbuf_init() are allowed to be
!-- collected. REQUIRED BY SAM.
    USE grid,              ONLY: dz, adz, nstat, nstep
    USE vars,              ONLY: rho
    USE hbuffer,           ONLY: hbuf_put
    USE params,            ONLY: lcond
#ifdef SAMTAK
    USE hbuf_1Dtimeseries, ONLY: ncondavg_1Dtimeseries, condavg_factor_1Dtimeseries,  &
                                 condavg_mask_1Dtimeseries, condavg_zcb_1Dtimeseries, &
                                 condavg_zct_1Dtimeseries, condavgname_1Dtimeseries,  &
                                 hbuf_1Dtimeseries_avg_put, hbuf_1Dtimeseries_put
#endif /*SAMTAK*/
    USE micro_vars,        ONLY: diss_LCM,                                         &
                                 ncc_prof, nrc_prof, cnd_prof, evp_prof, aut_prof, &
                                 acc_prof, nsd_prof, nac_prof, spm_prof, sps_prof, &
                                 sfx_prof, act_prof, dea_prof, std_SCLR_prof,      &
                                 T_prof, qv_prof, qa_prof, qc_prof, qr_prof,       &
                                 afx_prof, sgfx_prof, t_max_prof, t_mean_prof,     &
                                 t_count_prof, t_max_prof_liq, t_mean_prof_liq,    &
                                 t_count_prof_liq, t_max_prof_ice, t_mean_prof_ice,&
                                 t_count_prof_ice, rho_app_prof, iaavg_prof,       &
                                 icavg_prof, brk_prof, coal_prof, ds_prof, rs_prof,&
                                 icebrk_prof

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    REAL, DIMENSION(2)   :: tmp
    REAL                 :: factor_xy
    REAL, DIMENSION(nzm) :: tr0, tr2, diss, ncc
    INTEGER              :: k, n, ncond, i, j, kcb, kct, l, ierr

    CALL t_startf ('micro_statistics')

    factor_xy = 1.0 / FLOAT( nx * ny )

!
!-- SGS total water flux
    tr0(:) = ( mkwsb(1:nzm,1) * dz_LCM(:) + sgfx_prof / dV ) * lcond
    !-- Fine to leave dV here as an array, since we're dividing one array in z by another array in z, so it works correctly.
    CALL hbuf_put( 'QTFLUXS', tr0, factor_xy )
!
!-- Resolved + SGS total water flux
    tr0(:) = tr0(:) + ( mkwle(1:nzm,1) * dz_LCM(:) / dtn_LCM + afx_prof / dV ) * lcond
    CALL hbuf_put( 'QTFLUX', tr0, factor_xy )

    CALL hbuf_put( 'lag_NSD', nsd_prof, factor_xy )
    CALL hbuf_put( 'lag_NAC', nac_prof, factor_xy / 1.0E6 )
    CALL hbuf_put( 'lag_NCC', ncc_prof, factor_xy / 1.0E6 )
    CALL hbuf_put( 'lag_NRC', nrc_prof, factor_xy / 1.0E6 )
    CALL hbuf_put( 'lag_QA',  qa_prof/(dV*rho_air), factor_xy * 1.0E3 )
    CALL hbuf_put( 'lag_QC',  qc_prof/(dV*rho_air), factor_xy * 1.0E3 )
    CALL hbuf_put( 'lag_QR',  qr_prof/(dV*rho_air), factor_xy * 1.0E3 )
    CALL hbuf_put( 'lag_CND', cnd_prof/(dV*rho_air*dtn_LCM), 1.0E3 * factor_xy )
    CALL hbuf_put( 'lag_EVP', evp_prof/(dV*rho_air*dtn_LCM), 1.0E3 * factor_xy )
    CALL hbuf_put( 'lag_AUT', aut_prof/(dV*rho_air*dtn_LCM), 1.0E3 * factor_xy )
    CALL hbuf_put( 'lag_ACC', acc_prof/(dV*rho_air*dtn_LCM), 1.0E3 * factor_xy )

    CALL hbuf_put( 'lag_pQLF', sfx_prof/dV*lcond, factor_xy )
    CALL hbuf_put( 'lag_rQLF', afx_prof/dV*lcond, factor_xy )
    CALL hbuf_put( 'lag_sQLF', sgfx_prof/dV*lcond, factor_xy )

    CALL hbuf_put( 'lag_ACT', act_prof/(dV*rho_air*dtn_LCM), factor_xy / 1.0E6 )
    CALL hbuf_put( 'lag_DEA', dea_prof/(dV*rho_air*dtn_LCM), factor_xy / 1.0E6 )
    CALL hbuf_put( 'lag_TAU', tau_prof, factor_xy )
    CALL hbuf_put( 'lag_REF', REF_prof,  1.0E6 )
    CALL hbuf_put( 'lag_RFS', REFS_prof, 1.0E6 )
    CALL hbuf_put( 'lag_R',   R_prof,    1.0E6 )
    CALL hbuf_put( 'lag_RSD', RSD_prof,  1.0E6 )
    CALL hbuf_put( 'lag_AR', ar_prof, 1.0E6)
    call hbuf_put( 'lag_AM', am_prof/(dV*rho_air*dtn_LCM), factor_xy * 1.0E6)

    !Collection Quantities
    IF (switch_collection) THEN
      CALL hbuf_put( 'lag_COAL', coal_prof/(dV*rho_air*dtn_LCM), 1.0E3 * factor_xy)
      CALL hbuf_put( 'lag_BRK', brk_prof/(dV*rho_air*dtn_LCM), 1.0E3 * factor_xy)
    ENDIF

    IF ( switch_ice )  THEN
       CALL hbuf_put( 'lag_NIC',  nic_prof, factor_xy / 1.0E6 )
       CALL hbuf_put( 'lag_QI',   qi_prof/(dV*rho_air), factor_xy * 1.0E3 )
       CALL hbuf_put( 'lag_NSC',  nsc_prof, factor_xy / 1.0E6 )
       CALL hbuf_put( 'lag_QS',   qs_prof/(dV*rho_air), factor_xy * 1.0E3 )
       CALL hbuf_put( 'lag_NGC',  ngc_prof, factor_xy / 1.0E6 )
       CALL hbuf_put( 'lag_QG',   qg_prof/(dV*rho_air), factor_xy * 1.0E3 )

       CALL hbuf_put( 'lag_iRH',  rhice_prof, factor_xy )
       CALL hbuf_put( 'lag_NUC',  nuc_prof/(dV*rho_air*dtn_LCM), factor_xy / 1.0E6 )
       CALL hbuf_put( 'lag_IMM',  imm_prof/(dV*rho_air*dtn_LCM), factor_xy / 1.0E6 )
       !CALL hbuf_put( 'lag_HZN',  hzn_prof/(dV*rho_air*dtn_LCM), factor_xy / 1.0E6 )
       CALL hbuf_put( 'lag_HOM',  hom_prof/(dV*rho_air*dtn_LCM), factor_xy / 1.0E6 )
       CALL hbuf_put( 'lag_NML',  nml_prof/(dV*rho_air*dtn_LCM), factor_xy / 1.0E6 )
       CALL hbuf_put( 'lag_NSUB', nsub_prof/(dV*rho_air*dtn_LCM), factor_xy / 1.0E6)
       CALL hbuf_put( 'lag_DEP',  dep_prof/(dV*rho_air*dtn_LCM), 1.0E3 * factor_xy )
       CALL hbuf_put( 'lag_SUB',  sub_prof/(dV*rho_air*dtn_LCM), 1.0E3 * factor_xy )
       CALL hbuf_put( 'lag_FRZ',  frz_prof/(dV*rho_air*dtn_LCM), 1.0E3 * factor_xy )
       CALL hbuf_put( 'lag_MEL',  mel_prof/(dV*rho_air*dtn_LCM), 1.0E3 * factor_xy )
       CALL hbuf_put( 'lag_AGG',  agg_prof/(dV*rho_air*dtn_LCM), 1.0E3 * factor_xy )
       CALL hbuf_put( 'lag_iBRK',  icebrk_prof/(dV*rho_air*dtn_LCM), factor_xy )
       CALL hbuf_put( 'lag_DS',  ds_prof/(dV*rho_air*dtn_LCM), factor_xy )
       CALL hbuf_put( 'lag_RS',  rs_prof/(dV*rho_air*dtn_LCM), factor_xy )

       CALL hbuf_put( 'lag_RIM',  rim_prof/(dV*rho_air*dtn_LCM), 1.0E3 * factor_xy )
       CALL hbuf_put( 'lag_TFRZ', tfrz_prof, 1.0)
       CALL hbuf_put( 'lag_TFZC', tfrzc_prof, 1.0)
       CALL hbuf_put( 'lag_TFZR', tfrzr_prof, 1.0)
       CALL hbuf_put( 'lag_TFZI', tfrzi_prof, 1.0)
       CALL hbuf_put( 'lag_TFZS', tfrzs_prof, 1.0)
       CALL hbuf_put( 'lag_TFZG', tfrzg_prof, 1.0)
       CALL hbuf_put( 'lag_iRMF', rmf_prof, 1.0)

       CALL hbuf_put( 'lag_iREF', iREF_prof,  1.0E6 )
       CALL hbuf_put( 'lag_iRFS', iREFS_prof, 1.0E6 )
       CALL hbuf_put( 'lag_iR',   iR_prof,    1.0E6 )
       CALL hbuf_put( 'lag_iRSD', iRSD_prof,  1.0E6 )
       CALL hbuf_put( 'lag_iAXS', iaavg_prof, 1.0E6 )
       CALL hbuf_put( 'lag_iCXS', icavg_prof, 1.0E6 )
       CALL hbuf_put( 'lag_iPHI', iPHI_prof,  1.0 )
       CALL hbuf_put( 'lag_iRHO', rho_app_prof, 1.0 )
       CALL hbuf_put( 'lag_iNmn', Nmono_prof, 1.0 )
       CALL hbuf_put( 'lag_SSA', ssa_prof, 1.0E12)
       CALL hbuf_put( 'lag_DSA', dsa_prof, 1.0E12)
       CALL hbuf_put( 'lag_BSA', bsa_prof, 1.0E12)
       CALL hbuf_put( 'lag_PSA', psa_prof, 1.0E12)

       CALL hbuf_put( 'lag_NINC', ninc_prof, factor_xy / 1.0E6 )
    ENDIF

    IF ( switch_LEM_micro )  THEN
       CALL hbuf_put( 'lag_SPM',  spm_prof, 1.0E2 * factor_xy )
       CALL hbuf_put( 'lag_SPS',  sps_prof, 1.0E2 * factor_xy )
       CALL hbuf_put( 'lag_TABS', T_prof,   factor_xy )
       CALL hbuf_put( 'lag_QV',   qv_prof,  1.0E3 * factor_xy )
    ENDIF

    !
    !-- Particle lifetimes

    !
    !-- All particles
    CALL MPI_ALLREDUCE(t_max_prof,tr0,nzm,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    CALL hbuf_put( 'lag_TMX', tr0, 1.0)

    CALL MPI_ALLREDUCE(t_count_prof,tr0,nzm,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(t_mean_prof,tr2,nzm,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    DO k = 1, nzm
       IF ( tr0(k) .GT. 0.0 )  tr2(k) = tr2(k) / tr0(k)
    ENDDO
    CALL hbuf_put( 'lag_TMN', tr2, 1.0)

    !
    !-- Liquid particles only

    CALL MPI_ALLREDUCE(t_max_prof_liq,tr0,nzm,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    CALL hbuf_put( 'lag_LTMX', tr0, 1.0)

    CALL MPI_ALLREDUCE(t_count_prof_liq,tr0,nzm,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(t_mean_prof_liq,tr2,nzm,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    DO k = 1, nzm
       IF ( tr0(k) .GT. 0.0 )  tr2(k) = tr2(k) / tr0(k)
    ENDDO
    CALL hbuf_put( 'lag_LTMN', tr2, 1.0)

    !
    !-- Ice particles only
    IF (switch_ice) THEN
      CALL MPI_ALLREDUCE(t_max_prof_ice,tr0,nzm,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      CALL hbuf_put( 'lag_ITMX', tr0, 1.0)

      CALL MPI_ALLREDUCE(t_count_prof_ice,tr0,nzm,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(t_mean_prof_ice,tr2,nzm,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      
      DO k = 1, nzm
         IF ( tr0(k) .GT. 0.0 )  tr2(k) = tr2(k) / tr0(k)
      ENDDO
      CALL hbuf_put( 'lag_ITMN', tr2, 1.0)

      IF ( switch_particle_source )  THEN
         CALL hbuf_put( 'lag_NIN',  nin_prof/(dV*rho_air*dtn_LCM), factor_xy / 1.0E6 )
      ENDIF
    ENDIF

#ifdef SAMTAK
    DO ncond = 1, ncondavg_1Dtimeseries
!
!--    Mean cloud water droplet number
       tr0(1) = 0.0
       DO j = 1, ny
            DO i = 1, nx
                 IF ( condavg_mask_1Dtimeseries(i,j,ncond) > 0.0 ) THEN

                      kcb = condavg_zcb_1Dtimeseries(i,j,ncond)
                      kct = condavg_zct_1Dtimeseries(i,j,ncond)

          ncc(:)  = 0.0
                      DO k = kcb, kct
             ncc(k) = ncc(k) + micro_field(i,j,k,inc) + micro_field(i,j,k,inr)
                      ENDDO
                      tr0(1) = tr0(1) + SUM( ncc(kcb:kct) * adz(kcb:kct) ) / SUM( adz(kcb:kct) )

                  ENDIF
             ENDDO
        ENDDO
        CALL hbuf_1Dtimeseries_put( 'NCC'//TRIM(condavgname_1Dtimeseries(ncond)), tr0(1), 1.0 )

    ENDDO
#endif /*SAMTAK*/

    CALL t_stopf('micro_statistics')

 END SUBROUTINE micro_statistics

 END MODULE microphysics

 REAL FUNCTION Q_ext(r)
!
!-- Parameterization of Mie scattering based on Mitchell (2000, JAS)
    USE params,     ONLY: pi
    USE micro_vars, ONLY: lambda_micro
    IMPLICIT NONE

    REAL, INTENT(IN) :: r

    REAL :: D, x, rho, beta, a6, Q_ext_ADA, ra, rext, m, kmean, epsilon0, &
            kmax, C3, Q_edge, lambda

    REAL, SAVE :: nr, ni

    INTEGER :: n

    LOGICAL, SAVE :: initialized = .FALSE.
!
!-- Values from Hale and Querry (1973)
    REAL, DIMENSION(1:94), PARAMETER :: la_list =                              &
                           (/ 0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.350, &
                              0.375, 0.400, 0.425, 0.450, 0.475, 0.500, 0.525, &
                              0.550, 0.575, 0.600, 0.625, 0.650, 0.675, 0.700, &
                              0.725, 0.750, 0.775, 0.800, 0.825, 0.850, 0.875, &
                              0.900, 0.925, 0.950, 0.975, 1.000, 1.200, 1.400, &
                              1.600, 1.800, 2.000, 2.200, 2.400, 2.600, 2.650, &
                              2.700, 2.750, 2.800, 2.850, 2.900, 2.950, 3.000, &
                              3.050, 3.100, 3.150, 3.200, 3.250, 3.300, 3.350, &
                              3.400, 3.450, 3.500, 3.600, 3.700, 3.800, 3.900, &
                              4.000, 4.100, 4.200, 4.300, 4.400, 4.500, 4.600, &
                              4.700, 4.800, 4.900, 5.000, 5.100, 5.200, 5.300, &
                              5.400, 5.500, 5.600, 5.700, 5.800, 5.900, 6.000, &
                              6.100, 6.200, 6.300, 6.400, 6.500, 6.600, 6.700, &
                              6.800, 6.900, 7.000 /) * 1.0E-6

    REAL, DIMENSION(1:94), PARAMETER :: nr_list =                              &
                           (/ 1.396, 1.373, 1.362, 1.354, 1.349, 1.346, 1.343, &
                              1.341, 1.339, 1.338, 1.337, 1.336, 1.335, 1.334, &
                              1.333, 1.333, 1.332, 1.332, 1.331, 1.331, 1.331, &
                              1.330, 1.330, 1.330, 1.329, 1.329, 1.329, 1.328, &
                              1.328, 1.328, 1.327, 1.327, 1.327, 1.324, 1.321, &
                              1.317, 1.312, 1.306, 1.296, 1.279, 1.242, 1.219, &
                              1.188, 1.157, 1.142, 1.149, 1.201, 1.292, 1.371, &
                              1.426, 1.467, 1.483, 1.478, 1.467, 1.450, 1.432, &
                              1.420, 1.410, 1.400, 1.385, 1.374, 1.364, 1.357, &
                              1.351, 1.346, 1.342, 1.338, 1.334, 1.332, 1.330, &
                              1.330, 1.330, 1.328, 1.325, 1.322, 1.317, 1.312, &
                              1.305, 1.298, 1.289, 1.277, 1.262, 1.248, 1.265, &
                              1.319, 1.363, 1.357, 1.347, 1.339, 1.334, 1.329, &
                              1.324, 1.321, 1.317 /)

    REAL, DIMENSION(1:94), PARAMETER :: ni_list =                              &
            (/ 1.10E-7, 4.90E-8, 3.35E-8, 2.35E-8, 1.60E-8,  1.08E-8, 6.50E-9, &
               3.50E-9, 1.86E-9, 1.30E-9, 1.02E-9, 9.35E-10, 1.00E-9, 1.32E-9, &
               1.96E-9, 3.60E-9, 1.09E-8, 1.39E-8, 1.64E-8,  2.23E-8, 3.35E-8, &
               9.15E-8, 1.56E-7, 1.48E-7, 1.25E-7, 1.82E-7,  2.93E-7, 3.91E-7, &
               4.86E-7, 1.06E-6, 2.93E-6, 3.48E-6, 2.89E-6,  9.89E-6, 1.38E-4, &
               8.55E-5, 1.15E-4, 1.10E-3, 2.89E-4, 9.56E-4,  3.17E-3, 6.70E-5, &
               0.01900, 0.05900, 0.11500, 0.18500, 0.26800,  0.29800, 0.27200, &
               0.24000, 0.19200, 0.13500, 0.09240, 0.06100,  0.03680, 0.02610, &
               0.01950, 0.01320, 0.00940, 0.00515, 0.00360,  0.00340, 0.00380, &
               0.00460, 0.00562, 0.00688, 0.00845, 0.01030,  0.01340, 0.01470, &
               0.01570, 0.01500, 0.01370, 0.01240, 0.01110,  0.01010, 0.00980, &
               0.01030, 0.01160, 0.01420, 0.02030, 0.03300,  0.06220, 0.10700, &
               0.13100, 0.08800, 0.05700, 0.04490, 0.03920,  0.03560, 0.03370, &
               0.03270, 0.03220, 0.03200 /)

    lambda = lambda_micro

    IF ( .NOT. initialized )  THEN

       IF ( lambda .LE. la_list(1) )  THEN
          nr = nr_list(1)
          ni = ni_list(1)
       ELSEIF ( lambda .GE. la_list(94) )  THEN
          nr = nr_list(94)
          ni = ni_list(94)
       ELSE
          DO n = 1, 93
             IF ( lambda .GE. la_list(n)  .AND.  lambda .LE. la_list(n+1) )  THEN
                nr = nr_list(n) + ( nr_list(n+1) - nr_list(n) ) / ( la_list(n+1) - la_list(n) ) * ( lambda - la_list(n) )
                ni = ni_list(n) + ( ni_list(n+1) - ni_list(n) ) / ( la_list(n+1) - la_list(n) ) * ( lambda - la_list(n) )
                EXIT
             ENDIF
          ENDDO
       ENDIF

       initialized = .TRUE.

    ENDIF

    D    = 2.0 * r
    x    = 2.0 * pi * r / lambda
    rho  = 2.0 * x * ( nr - 1.0 )
    beta = ATAN( ni / ( nr - 1.0 ) )
    a6   = 1.0
!
!-- Anomalous diffraction approximation (ADA) based on van de Hulst (1981), see (14) in Mitchell (2000)
    Q_ext_ADA = 2.0 - 4.0 * EXP( -rho * TAN( beta ) ) * COS( beta ) / rho * SIN( rho - beta ) &
                - 4.0 * EXP( -rho * TAN( beta ) ) * ( COS( beta ) / rho )**2 * COS( rho - 2.0 * beta ) &
                + 4.0 * ( COS( beta ) / rho )**2 * COS( 2.0 * beta )
!
!-- Tunneling effects
    ra = 0.7393 * nr - 0.6069 ! (6)
    rext = ra / 2.0 ! (11)

    m = 0.5 ! dispersion gamma function, 0.5 has been found to be a good approcimation
    kmean = D / lambda
    epsilon0 = 0.25 + 0.6 * ( 1.0 - EXP( -8.0 * pi * ni / 3.0 ) )**2 ! (8)
    kmax = m / epsilon0

    C3 = rext * ( kmean**m * exp( -epsilon0 * kmean ) ) / ( kmax**m * EXP( -m ) ) ! (12)
!
!-- Edge effects
    Q_edge = a6 * ( 1.0 - exp( -0.06 * x ) ) * x**(-2.0/3.0) ! (17)
!
!-- Final approximation
    Q_ext = ( 1.0 + C3 ) * Q_ext_ADA + Q_edge ! (18)
    RETURN

 END FUNCTION Q_ext
