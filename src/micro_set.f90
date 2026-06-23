!TODO: adapt boundaries of u_LCM, u_LCM, u_LCM

  SUBROUTINE micro_set
!
!-- This subroutine has to be called in main.f90 when SAM starts.
!-- Throughout this subroutine, the specified variables are checked for
!-- inconsistencies.
!--
!-- 1. Set MPI_datatype = particle_array.
!-- 2. Obtain grid box information: actual distance of grid box center and side.
!-- 4. Allocate arrays
!-- 5. Terminate if an error is found
!-- 6. Setup output names, short names, and units. Check if nv01 == # of output scalar flags
!-- 7. Set flag_LCM_init = .FALSE. or read restart data
    USE micro_vars
    USE micro_restart,       ONLY: micro_restart_read
    USE micro_task_util_MPI, ONLY: micro_particle_array2mpi_datatype, micro_task_abort

    USE domain,        ONLY: YES3D, nx_gl, ny_gl, nsubdomains_x, nsubdomains_y
    USE grid,          ONLY: nx, ny, nzm, rank, RUN3D, dx, dy, dz, icycle, &
                             masterproc, ncycle, nelapse, nrestart_skip, nstat, &
                             nstep, nstop, nsubdomains, nz, zi, z, adz
    USE params,        ONLY: doprecip, docloud

    IMPLICIT NONE

    INTEGER :: i, j, k, ierr
    LOGICAL :: flag = .TRUE.
    LOGICAL :: flag_allocation = .TRUE.

    IF ( nzl .GT. nzm ) THEN
       IF (masterproc) PRINT*, 'LCM ERROR: Do not initialize particles higher than z(nzm).'
       flag = .FALSE.
    ENDIF

    IF ( switch_radiation_micro .AND. switch_analytical_micro ) THEN
       IF (masterproc) PRINT*, 'LCM ERROR: switch_radiation_micro = .TRUE. requires switch_analytical_micro = .FALSE.'
       flag = .FALSE.
    ENDIF

    IF (switch_LEM_micro .and. (maxval(adz) .gt. 1. .or. minval(adz) .lt. 1.)) then
      if (masterproc) print*, 'LEM cannot be run with vertical grid stretching'
      flag = .false.
    endif

    IF ( TRIM(init_aero) .NE. 'bin'  .AND.  TRIM(init_aero) .NE. 'random'  .AND.  TRIM(init_aero) .NE. 'USM' ) THEN
       IF (masterproc) PRINT*,'LCM ERROR: ', TRIM(init_aero), ' is not valid value for init_aero'
       flag = .FALSE.
    ENDIF

    IF ( .NOT. doprecip ) THEN
       IF (masterproc) PRINT*,'LCM ERROR: requires doprecip to be true'
       flag = .FALSE.
    ENDIF

    IF ( .NOT. docloud ) THEN
       IF (masterproc) PRINT*,'LCM ERROR: requires docloud to be true'
       flag = .FALSE.
    ENDIF

    IF ( switch_ice  .AND.  switch_LEM_micro ) THEN
       IF (masterproc) PRINT*,'LCM ERROR: ice does not work with LEM SGS model'
       flag = .FALSE.
    ENDIF

    IF ( switch_ice  .AND.  switch_radiation_micro ) THEN
       IF (masterproc) PRINT*,'LCM ERROR: radiation effects on ice need to be included'
       flag = .FALSE.
    ENDIF

    IF ( switch_ice  .AND.  TRIM(init_aero) .NE. 'random' ) THEN
       IF (masterproc) PRINT*,'LCM ERROR: ice requires init_aero = random'
       flag = .FALSE.
    ENDIF

    IF ( switch_ice  .AND.  ( frac_ice .LT. 0.0  .OR.  frac_ice .GT. 1.0 ) ) THEN
       IF (masterproc) PRINT*,'LCM ERROR: ice requires frac_ice between 0 and 1'
       flag = .FALSE.
    ENDIF

    IF ((.NOT. switch_ice) .AND. switch_particle_source .AND. (MAXVAL(p_source_in) .GT. 0 .OR. (MAXVAL(p_source_aero * n_source_aero * frac_in_source_particles) .NE. 0 ))) THEN
      PRINT*,'LCM ERROR: ice particle source present but ice microphysics disabled'
      STOP
    ENDIF
!
!-- Initialize Lagrangian particle information
    n_init  = nx_gl / nsubdomains_x * ny_gl / nsubdomains_y * np_xy * nzl ! # of particles for subdomain at initial time step
    n_total = nsubdomains * n_init ! # of particles for whole domain
    IF ( INT(nx_gl,KIND=8) / INT(nsubdomains_x,KIND=8) &
         * INT(ny_gl,KIND=8) / INT(nsubdomains_y,KIND=8) * INT(np_xy,KIND=8) *   &
         INT(nzl,KIND=8) .GT. 2147483647_8 ) THEN
     PRINT*,'LCM ERROR: Too many particles per subdomain (more than INTEGER(4) limit)'
     CALL micro_task_abort
    ENDIF
    IF ( INT(nsubdomains,KIND=8) * INT(nx_gl,KIND=8) / INT(nsubdomains_x,KIND=8) &
         * INT(ny_gl,KIND=8) / INT(nsubdomains_y,KIND=8) * INT(np_xy,KIND=8) *   &
         INT(nzl,KIND=8) .GT. 9223372036854775807_8 ) THEN
     PRINT*,'LCM ERROR: Too many particles (more than INTEGER(8) limit)'
     CALL micro_task_abort
    ENDIF
!
!-- Horiaontal width of whole domain
    xwidth = REAL( nx_gl ) * dx
    ywidth = REAL( ny_gl ) * dy
!
!-- Inverse of grid space
    idx = 1.0 / dx
    idy = 1.0 / dy
!
!-- Volume of a grid box
    DO k = 1, nz_LCM
       dV(k) = dx * dy * ( zi(k+1) - zi(k) )
       dz_LCM(k) = (zi(k+1)-zi(k))
    ENDDO
!
!-- Height at which particles will be deleted / reflected
    ztop    = zi(nz) - boundary_distance
    zbottom = 0.0 + boundary_distance
!
!-- Save LCM version of grid spacings
    dx_LCM = dx
    dy_LCM = dy
!
!-- Get global start id's of i and j
    CALL task_rank_to_index(rank,i_global_LCM,j_global_LCM)

    IF (masterproc) THEN
       PRINT*,''
       PRINT*,'==================================== Lagrangian Cloud Model ===================================='
       PRINT*,'  Number of particles per grid box:', INT(np_xy,KIND=2)
       PRINT*,'  Total number of Lagrangian particles:', n_total
       PRINT*,'  Particles have been initialized up to a height of:', REAL(z(nzl),KIND=4)
       PRINT*,'  LCM starts at time step:', INT(nstep_LCM_init,KIND=2)
       PRINT*,''
       PRINT*,'--------------------------------------- Particle motion ----------------------------------------'
       PRINT*,'  Number of resolved scale velocity iterations:', INT(n_uvw_iteration,KIND=2)
       PRINT*,'  Particles are allowed to sediment:', switch_sedimentation
       IF (switch_sedimentation) PRINT*,'  Sedimentation starts at timestep:', start_sedimentation
       PRINT*,''
       PRINT*,'----------------------------------------- SGS modeling -----------------------------------------'
       PRINT*,'  LCM SGS velocities are on:', switch_sgs_velocities
       PRINT*,'  LCM SGS mixing model is used:', switch_LEM_micro
       PRINT*,''
       PRINT*,'------------------------------------- Condensation growth --------------------------------------'
       IF ( switch_analytical_micro )  THEN
          PRINT*,'  Analytic solution is used to solve growth equation (switch_analytical_micro = .TRUE.).'
          PRINT*,'  Total aerosol/cloud droplet number concentration (in cm^-3):', REAL(n_aero/1.0E6,KIND=4)
       ELSE
          PRINT*,'  Numerical Euler method is used to solve growth equation (switch_analytical_micro = .FALSE.).'
          PRINT*,'  Total aerosol number concentration (in cm^-3):', REAL(n_aero/1.0E6,KIND=4)
          PRINT*,'  Geometric mean radius of aerosol size distribution (in µm):', REAL(rm_aero*1.0E6,KIND=4)
          PRINT*,'  Geometric stddev of aerosol size distribution (unitless):', REAL(sigma_aero,KIND=4)
          PRINT*,'  The aerosol distribution is initialized using the >',TRIM(init_aero),'< approach (init_aero).'
          IF ( switch_radiation_micro )  THEN
             PRINT*,'  Radiation effects on growth equation are considered (switch_radiation_micro =',switch_radiation_micro,').'
          ELSE
             PRINT*,'  Radiation effects on growth equation are NOT considered (switch_radiation_micro =',switch_radiation_micro,').'
          ENDIF
       ENDIF
       PRINT*,''
       PRINT*,'-------------------------------------- Collection growth ---------------------------------------'
       IF ( switch_collection )  THEN
          PRINT*,'  Collection growth is allowed (switch_collection =',switch_collection,').'
          PRINT*,'  Collection starts at timestep:', start_collection
       ELSE
          PRINT*,'  Collection growth is NOT allowed (switch_collection =',switch_collection,').'
       ENDIF
       PRINT*,''
       PRINT*,'----------------------------------------- Statistics -------------------------------------------'
       PRINT*,'  Cloud/rain water is separated at a radius of (in microns):', REAL(r_separate_micro*1.0E6,KIND=4)
       PRINT*,'  Droplet spectra are saved all ', spec_nsave, ' timesteps.'
       PRINT*,'  Saving droplet spectra starts after', spec_nsavestart, ' timesteps.'
    ENDIF
!
!-- STEP 1: Set MPI_datatype = particle_array
    IF ( nsubdomains > 1 )  CALL micro_particle_array2mpi_datatype
!
!-- STEP 2: Get grid position information
!--         SAM uses Arakawa C grid
!--
!--         SAM grid index orientation (s=scalar)
!--         - Horizontal:
!--          ---------------v(i=1,j=2)----------------
!--         |                                         |
!--         u(i=1,j=1)     w(i=1,j=1), s(i=1,j=1)     u(i=2,j=1)
!--         |                                         |
!--          ---------------v(i=1,j=1)----------------
!--
!--         - Vertical:
!--         w(k=2), zi(2) = gpzs(2) = z(0)+adz(1)*dz
!--           |
!--         u(k=1), v(k=1), s(k=1), gpzc(1) = z(1) = z(0)+0.5 * adz(1)*dz
!--           |
!--         w(k=1), zi(1) = gpzs(1) = 0
    nx_sb = MOD( rank, nsubdomains_x ) * nx
    ny_sb = ( rank / nsubdomains_x ) * ny
!
!-- Side
    DO i = 0, nx+1
       gpxs(i) = REAL(i-1+nx_sb) * dx  ! location of left side of grid box i, j, k
    ENDDO
    DO j = 0, ny+1
       gpys(j) = REAL(j-1+ny_sb) * dy  ! location of south side of grid box i, j, k
    ENDDO
    DO k = 1, nz
       gpzs(k) = zi(k)                 ! location of bottom of grid box i, j, k
    ENDDO
!
!-- Center
    DO i = -1, nx+1
       gpxc(i) = ( REAL(i+nx_sb) - 0.5 ) * dx  ! center's location of gridbox i, j, k in x-direction
    ENDDO
    DO j = -1, ny+1
       gpyc(j) = ( REAL(j+ny_sb) - 0.5 ) * dy  ! center's location of gridbox i, j, k in y-direction
    ENDDO
    DO k = 1, nzm
       gpzc(k) = z(k)                          ! center's location of gridbox i, j, k in z-direction
    ENDDO
!
!-- STEP 4:   Allocate arrays
    ALLOCATE( u_LCM( x1u:x2u, y1u:y2u, nzm ), &
              v_LCM( x1v:x2v, y1v:y2v, nzm ), &
              w_LCM( x1w:x2w, y1w:y2w, nzm ), &
              STAT=ierr )

    IF ( ierr /= 0 ) flag_allocation = .FALSE.

    ALLOCATE( tabs_LCM( 1:nx_LCM, 1:ny_LCM, nz_LCM ), &
              tke_LCM(  1:nx_LCM, 1:ny_LCM, nz_LCM ), &
              diss_LCM( 1:nx_LCM, 1:ny_LCM, nz_LCM ), &
              tk_LCM(   1:nx_LCM, 1:ny_LCM, nz_LCM ), &
              smix_LCM( 1:nx_LCM, 1:ny_LCM, nz_LCM ), &
              STAT=ierr )
    IF ( ierr /= 0 ) flag_allocation = .FALSE.

    IF ( .NOT.flag_allocation ) THEN
       PRINT*,'LCM ERROR: Failed to allocate arrays in micro_LCM_set on proc:', rank
       CALL micro_task_abort
    ENDIF
!
!-- STEP 5: Terminate if flag=.FALSE.
    IF ( .NOT.flag ) THEN
       IF (masterproc) PRINT*,'LCM ERROR: There is an error in the initialization. Abort.'
       CALL micro_task_abort
    ENDIF
!
!-- STEP 6: Setup output variable name list. Initialize history_length.
    history_length = 0
!
!-- STEP 7: Set the initialization flag false or read restart data.
!--         At this point, nstep is either 0 or last nstop
    IF ( nstep < nstep_LCM_init ) THEN
       flag_LCM_init = .FALSE.
    ELSE IF ( nstep >= nstep_LCM_init ) THEN ! read restart data
       CALL micro_restart_read
       flag_LCM_init = .TRUE.
    ENDIF

    IF (masterproc) THEN
       PRINT*,''
       PRINT*,'================================================================================================'
       PRINT*,'  The subroutine micro_LCM_set did not detect inconsistency.'
       PRINT*,'================================================================================================'
       PRINT*,''
    ENDIF

 END SUBROUTINE micro_set
