 SUBROUTINE micro_place_particle
!
!-- Places particles into domain
    USE micro_task_util_MPI, ONLY: micro_task_abort
    USE micro_random,        ONLY: randomu
    USE micro_vars,          ONLY: n_start, n_end, particle, iLEM, iUsgs,      &
                                   iVsgs, iWsgs, lcmp, n_init, n_particles,    &
                                   n_particles_alloc, np_xy, nx_sb, ny_sb, nzl,&
                                   pres_LCM, tabs_LCM, dx_LCM, dy_LCM, dz_LCM, &
                                   nx_LCM, ny_LCM, iU, iV, iW, iX, iY, iZ,     &
                                   n_aero, iilst, ijlst,                       &
                                   iklst, i_global_LCM, j_global_LCM, iID,     &
                                   niv0, nv0, nv1, n_ID, iWsedi, iSPR, itime,  &
                                   iTRACK, track_some_particles, n_particles,  &
                                   n_track_particles, inv_cdf, n_inv_cdf

    USE params,              ONLY: pi
    USE vars,                ONLY: qsatw
    USE grid,                ONLY: rank, zi, nsubdomains_x, nsubdomains_y
    USE domain,              ONLY: YES3D

    IMPLICIT NONE

    INTEGER    :: i, j, k, m, n, ierr
!
!-- Allocate Lagrangian particle array
    n_particles       = n_init
    IF ( n_init .EQ. 0 )  THEN
       n_particles_alloc = 1
    ELSE
       n_particles_alloc = n_init + MAX( n_init / INT(SQRT(REAL(n_init))) * 6 , 10 )
    ENDIF
    ALLOCATE( particle(1:n_particles_alloc), STAT=ierr )
    IF ( ierr /= 0 ) THEN
       PRINT*,'LCM ERROR: Failed to allocate arrays in micro_place_particle on proc:', rank
       CALL micro_task_abort
    ENDIF

    ! temporarily (until end of micro_place_particle) allocate lookup array of inverse aerosol cdf
    ALLOCATE(inv_cdf(4, 0:n_inv_cdf))

    n_start = 0
    n_end   = -1
!
!-- Place particles in model domain
    n  = 1
    DO k = 1, nzl   ! z-direction
       DO j = 1, ny_LCM   ! y-direction
          DO i = 1, nx_LCM   ! x-direction

             n_start(i,j,k) = n
             DO m = 1, np_xy   ! multiple particle on the same grid

                particle(n)%ivar(iLEM)  = m
                particle(n)%ivar(iilst) = i + i_global_LCM
                particle(n)%ivar(ijlst) = j + j_global_LCM
                particle(n)%ivar(iklst) = k
                particle(n)%ivar(iSPR)  = -999 ! sprayed particle mark initialization

                IF ( track_some_particles ) THEN
                  ! following: the condition to track particles, adjust to your preferences,
                  !            default is tracking evenly distributed thoughout the current particle array
                  IF ( MOD( n, INT(n_particles/MAX(1,n_track_particles/nsubdomains_x/nsubdomains_y)) ) == 0 ) THEN
                    !particle(n)%ivar(iTRACK) =  126 ! put marker for tracking particle
                    ! put marker for tracking particle, make it an ID, because iID seems to not work properly
                    particle(n)%ivar(iTRACK) = rank*1000000 + n/INT(n_particles/n_track_particles)
                  ELSE
                    particle(n)%ivar(iTRACK) = -999 ! default value, don't track particle
                  END IF
                END IF

                particle(n)%ivar(iID)   = n_ID
                n_ID = n_ID + 1_8

                particle(n)%var(itime) = 0.0_lcmp

                particle(n)%var(iX)    = REAL( ( i - 1 + nx_sb + randomu(rank+1) ) * dx_LCM, KIND=lcmp )
                IF ( YES3D .EQ. 1 )  THEN
                   particle(n)%var(iY) = REAL( ( j - 1 + ny_sb + randomu(rank+1) ) * dy_LCM, KIND=lcmp )
                ELSE
                   particle(n)%var(iY) = REAL( 0.5 * dy_LCM, KIND=lcmp )
                ENDIF
                !particle(n)%var(iZ)    = REAL( ( k - 1 + randomu(rank+1) ) * dz_LCM, KIND=lcmp )
                particle(n)%var(iZ) = REAL(zi(k) + randomu(rank+1)*dz_LCM(k), kind = lcmp) !Check this!

                particle(n)%var(iU)    = 0.0_lcmp
                particle(n)%var(iV)    = 0.0_lcmp
                particle(n)%var(iW)    = 0.0_lcmp

                particle(n)%var(iUsgs) = 0.0_lcmp
                particle(n)%var(iVsgs) = 0.0_lcmp
                particle(n)%var(iWsgs) = 0.0_lcmp

                particle(n)%var(iWsedi) = 0.0_lcmp  ! Sedimentation is always set to zero during initialization

                n = n + 1

             ENDDO
             n_end(i,j,k) = n-1
!
!--          Initilize further microphysical particle properties
             CALL micro_init_particles(i,j,k)
! !
! !--          Prepare calculation of resolved scale velocites
!              CALL micro_uvw_prep(i,j,k) ! not necessary here?!

          ENDDO
       ENDDO
    ENDDO

    ! deallocate inverse aerosol cdf to free up space
    DEALLOCATE(inv_cdf)

 END SUBROUTINE micro_place_particle

 SUBROUTINE micro_particle_source
   !-- Places particles into domain
       USE micro_task_util_MPI, ONLY: micro_task_abort
       USE micro_random,        ONLY: randomu
       USE micro_vars,          ONLY: n_start, n_end, particle, iLEM, iUsgs,      &
                                      iVsgs, iWsgs, lcmp, n_init, n_particles,    &
                                      n_particles_alloc, np_xy, nx_sb, ny_sb, nzl,&
                                      pres_LCM, tabs_LCM, dx_LCM, dy_LCM, dz_LCM, &
                                      nx_LCM, ny_LCM, iU, iV, iW, iX, iY, iZ,     &
                                      gpzs, n_aero, iilst, ijlst,                 &
                                      iklst, i_global_LCM, j_global_LCM, iID,     &
                                      niv0, nv0, nv1, n_ID, n_start_source_aero,  &
                                      particle_array, dtn_LCM, iA, iM, iN,        &
                                      in_startind, iTFZ,                          &
                                      ieta, iT, nz_LCM, iq, n_source_aero,        &
                                      nx_e_source_aero, nx_w_source_aero,         &
                                      ny_s_source_aero, ny_n_source_aero,         &
                                      nz_b_source_aero, nz_t_source_aero,         &
                                      n_source_in, type_source_in,                &
                                      r_source_aero, r_source_in,                 &
                                      nx_e_source_in, nx_w_source_in,             &
                                      ny_s_source_in, ny_n_source_in,             &
                                      nz_b_source_in, nz_t_source_in,             &
                                      rho_aero, rho_liq, dV, rho_air,             &
                                      p_source_aero, p_source_in,                 &
                                      nin_prof, sigma_source_aero,                &
                                      sigma_source_in, iWsedi,                    &
                                      iTYPE, switch_ice, iNmono, unit_source,     &
                                      switch_moving_source,                       &
                                      iSPR, vx_source_aero, vy_source_aero,       &
                                      vz_source_aero,                             &
                                      x_e_source_aero, x_w_source_aero,           &
                                      y_s_source_aero, y_n_source_aero,           &
                                      z_b_source_aero, z_t_source_aero,           &
                                      vx_source_in, vy_source_in, vz_source_in,   &
                                      x_e_source_in, x_w_source_in,               &
                                      y_s_source_in, y_n_source_in,               &
                                      z_b_source_in, z_t_source_in,               &
                                      type_source, salinity_source,               &
                                      n_start_source_in, frac_in_source_particles,&
                                      frac_in_source_weight, switch_phase_sources,&
                                      itime, n_end_source_aero, n_end_source_in,  &
                                      t_frz,  &
                                      iAero_portions, r_equi, search_sorted_array
       USE microphysics,        ONLY: micro_field, imv
   
       USE params,              ONLY: pi, fac_cond, cp, ggr
       USE vars,                ONLY: qsatw, t, gamaz
       USE grid,                ONLY: rank, nstep, nsubdomains, zi, z
       USE domain,              ONLY: YES3D, nx_gl, ny_gl
   
       IMPLICIT NONE
   
       INTEGER    :: i, j, k, m, n, ierr, n_particles_temp, n_alloc_source_aero,       &
                     n_alloc_source_in, n_alloc_source_total, pint_source_aero, pint_source_in, &
                     mm

       INTEGER(8) :: A_dummy
       REAL       :: SS_gb, dql, x_rand, r_aero, r_in, ierf_dummy, IERF, ERFI, sfc_in, tfz
   
       REAL :: temp_p_source_aero(4), temp_p_source_in(4)

      !  REAL, EXTERNAL :: t_frz
   
       TYPE(particle_array), DIMENSION(:), ALLOCATABLE :: particle_temp
   
       INTEGER :: n_start_temp(nx_LCM,ny_LCM,nz_LCM)
      !  REAL, DIMENSION(:,:), ALLOCATABLE :: zdiffarr_b, zdiffarr_t
       INTEGER :: n_end_temp(nx_LCM,ny_LCM,nz_LCM)
   !
   !-- Moving particle source
       DO m = 1, 4
   
   !-- Find indices of vertical grid at which to place particles
         ! allocate(zdiffarr_b(4, nz_LCM))
         ! allocate(zdiffarr_t(4, nz_LCM))
         ! do k = 1,nz_LCM
         !    zdiffarr_b(m, k) = abs(z(k)-z_b_source(m))
         !    zdiffarr_t(m, k) = abs(z(k)-z_t_source(m))
         ! end do
   
         ! deallocate(zdiffarr_b)
         ! deallocate(zdiffarr_t)
         temp_p_source_aero(m) = p_source_aero(m)
         temp_p_source_in(m) = p_source_in(m)
       END DO
   
       IF ((nstep .LT. n_start_source_aero .OR. (nstep .GT. n_end_source_aero .AND. n_end_source_aero .GT. 0)) .OR. (n_start_source_aero .LT. 0)) temp_p_source_aero = 0.0
       IF ((nstep .LT. n_start_source_in .OR. (nstep .GT. n_end_source_in .AND. n_end_source_in .GT. 0) .OR. (n_start_source_in .LT. 0)) .AND. switch_phase_sources .AND. switch_ice) temp_p_source_in = 0.0
   !
   !-- Timing controls on sources for CCN and IN
   !-- If the current timestep is less than the initial spraying timestep or greater than the final spraying timestep, 
   !-- then obviously don't spray for that species
   !-- Exceptions: If the initial spraying timestep is -1, never start spraying. If the final timestep is -1 and the initial
   !-- timestep is a positive integer, then spray until the end of the simulation. By default, both the starting and ending
   !-- timesteps are set to -1. This means that should a user wish to turn on a sprayer and have it run through the end of the
   !-- simulation, they only need to modify n_start_source_aero (or n_start_source_in), and not change n_end_source_aero (or n_end_source_in)
   !-- B. Ascher 08.05.2025

   !
   !-- IN type selection
   !-- Users can choose which species of IN superparticles emitted in source region should have. This is controlled by the type_source_in
   !-- variable. So, for example, if type_source_in is (/3, 1, 2, 4/), particles emitted in source region 1 will contain Snowmax IN,
   !-- particles emitted in region 2 will contain soot, those in region 3 will contain mineral dust, and those in region 4 will contain pollen.
   !-- B. Ascher 26.05.2025
       
       DO m = 1,4
   !
   !--   Exit the loop when temp_p_source(m) <= 0.0
         IF ( MAX(temp_p_source_aero(m), temp_p_source_in(m)) .LE. 0.0 )  EXIT
   !
   !--   Move vessle location according to their speed
         IF (switch_moving_source) THEN
            IF (temp_p_source_aero(m) .GT. 0.0) THEN
             !Move zonal boundaries
             x_e_source_aero(m) =  x_e_source_aero(m) + vx_source_aero(m) * dtn_LCM
             x_w_source_aero(m) =  x_w_source_aero(m) + vx_source_aero(m) * dtn_LCM
             !Apply zonal periodic boundary conditions
             x_e_source_aero(m) = MOD(x_e_source_aero(m), nx_gl * dx_LCM)
             x_w_source_aero(m) = MOD(x_w_source_aero(m), nx_gl * dx_LCM)
             !Move meridional boundaries
             y_n_source_aero(m) =  y_n_source_aero(m) + vy_source_aero(m) * dtn_LCM
             y_s_source_aero(m) =  y_s_source_aero(m) + vy_source_aero(m) * dtn_LCM
             !Apply meridional periodic boundary conditions
             IF (YES3D .EQ. 1) THEN
               y_n_source_aero(m) = MOD(y_n_source_aero(m), ny_gl * dy_LCM)
               y_s_source_aero(m) = MOD(y_s_source_aero(m), ny_gl * dy_LCM)
             ENDIF
             !Move vertical boundaries
             z_t_source_aero(m) =  z_t_source_aero(m) + vz_source_aero(m) * dtn_LCM
             z_b_source_aero(m) =  z_b_source_aero(m) + vz_source_aero(m) * dtn_LCM
            ENDIF
            IF (temp_p_source_in(m) .GT. 0.0) THEN !Move IN source if it exists
               !Move zonal boundaries
               x_e_source_in(m) =  x_e_source_in(m) + vx_source_in(m) * dtn_LCM
               x_w_source_in(m) =  x_w_source_in(m) + vx_source_in(m) * dtn_LCM
               !Apply zonal periodic boundary conditions
               x_e_source_in(m) = MOD(x_e_source_in(m), nx_gl * dx_LCM)
               x_w_source_in(m) = MOD(x_w_source_in(m), nx_gl * dx_LCM)
               !Move meridional boundaries
               y_n_source_in(m) =  y_n_source_in(m) + vy_source_in(m) * dtn_LCM
               y_s_source_in(m) =  y_s_source_in(m) + vy_source_in(m) * dtn_LCM
               !Apply meridional periodic boundary conditions
               IF (YES3D .EQ. 1) THEN
                  y_n_source_in(m) = MOD(y_n_source_in(m), ny_gl * dy_LCM)
                  y_s_source_in(m) = MOD(y_s_source_in(m), ny_gl * dy_LCM)
               ENDIF
               !Move vertical boundaries
               z_t_source_in(m) =  z_t_source_in(m) + vz_source_in(m) * dtn_LCM
               z_b_source_in(m) =  z_b_source_in(m) + vz_source_in(m) * dtn_LCM
            ENDIF
   !
   !--       Convert source position (x,y,z) to model indicies (i,j,k)
             !nz_b_source(m) = INT( z_b_source(m) / dz_LCM ) + 1
             !nz_t_source(m) = INT( z_t_source(m) / dz_LCM ) + 1
   
             IF (temp_p_source_aero(m) .GT. 0.0) THEN
               ny_s_source_aero(m) = INT( y_s_source_aero(m) / dy_LCM ) + 1
               ny_n_source_aero(m) = INT( y_n_source_aero(m) / dy_LCM ) + 1
               nx_w_source_aero(m) = INT( x_w_source_aero(m) / dx_LCM ) + 1
               nx_e_source_aero(m) = INT( x_e_source_aero(m) / dx_LCM ) + 1
               nz_b_source_aero(m) = search_sorted_array(gpzs, z_b_source_aero(m), nz_LCM+1)
               nz_t_source_aero(m) = search_sorted_array(gpzs, z_t_source_aero(m), nz_LCM+1)
             ENDIF
   
             IF (temp_p_source_in(m) .GT. 0.0) THEN
               ny_s_source_in(m) = INT( y_s_source_in(m) / dy_LCM ) + 1
               ny_n_source_in(m) = INT( y_n_source_in(m) / dy_LCM ) + 1
               nx_w_source_in(m) = INT( x_w_source_in(m) / dx_LCM ) + 1
               nx_e_source_in(m) = INT( x_e_source_in(m) / dx_LCM ) + 1
               nz_b_source_in(m) = search_sorted_array(gpzs, z_b_source_in(m), nz_LCM+1)
               nz_t_source_in(m) = search_sorted_array(gpzs, z_t_source_in(m), nz_LCM+1)
             ENDIF
         ENDIF
   !
   !--   Allocate temporary particle array
         IF (.NOT. switch_phase_sources) THEN
            IF ( n_particles .EQ. 0 )  THEN
               n_alloc_source_total = ( nz_t_source_aero(m) - nz_b_source_aero(m) + 1 ) * ( ny_n_source_aero(m) - ny_s_source_aero(m) + 1 ) * &
                                    ( nx_e_source_aero(m) - nx_w_source_aero(m) + 1 )
            ELSE
               n_alloc_source_total = MAX( n_particles, ( nz_t_source_aero(m) - nz_b_source_aero(m) + 1 ) * ( ny_n_source_aero(m) - ny_s_source_aero(m) + 1 ) &
                                 * ( nx_e_source_aero(m) - nx_w_source_aero(m) + 1 ) )
            ENDIF
         ELSE
            IF (temp_p_source_aero(m) .GT. 0.0) THEN
               n_alloc_source_aero = (nz_t_source_aero(m) - nz_b_source_aero(m) + 1 ) * (ny_n_source_aero(m) - ny_s_source_aero(m) + 1) * &
               (nx_e_source_aero(m) - nx_w_source_aero(m) + 1)
            ELSE
               n_alloc_source_aero = 0 !Need to reset this for every value of m
            ENDIF
            IF (temp_p_source_in(m) .GT. 0.0) THEN
               n_alloc_source_in = (nz_t_source_in(m) - nz_b_source_in(m) + 1 ) * (ny_n_source_in(m) - ny_s_source_in(m) + 1) * &
               (nx_e_source_in(m) - nx_w_source_in(m) + 1)
            ELSE
               n_alloc_source_in = 0 !Need to reset this for every value of m
            ENDIF
            IF (n_particles .EQ. 0) THEN
               n_alloc_source_total = n_alloc_source_aero + n_alloc_source_in
            ELSE
               n_alloc_source_total = MAX(n_particles, n_alloc_source_aero + n_alloc_source_in)
            ENDiF
         ENDIF
   
   
   
         ALLOCATE( particle_temp(1:2*n_alloc_source_total), STAT=ierr )
         IF ( ierr /= 0 ) THEN
            PRINT*,'LCM ERROR: Failed to allocate arrays in micro_particle_source on proc:', rank
            CALL micro_task_abort
         ENDIF
         n_start_temp     = 0
         n_end_temp       = -1
         n_particles_temp = 0

         nin_prof = 0.0
   !
   !--   Place particles in model domain
         DO k = 1, nz_LCM   ! z-direction
            DO j = 1, ny_LCM   ! y-direction
               DO i = 1, nx_LCM   ! x-direction
   
                  IF ( n_end(i,j,k) - n_start(i,j,k) .GE. 0 )  THEN
                     n_start_temp(i,j,k) = n_particles_temp + 1
                     n_end_temp(i,j,k)   = n_start_temp(i,j,k) + ( n_end(i,j,k) - n_start(i,j,k) )
                     mm = n_start(i,j,k)
                     DO n = n_start_temp(i,j,k), n_end_temp(i,j,k)
                        particle_temp(n)%var  = particle(mm)%var
                        particle_temp(n)%ivar = particle(mm)%ivar
                        mm = mm + 1
                     END DO
                     n_particles_temp = n_end_temp(i,j,k)
                  ENDIF
   
                  IF ( k .GE. nz_b_source_aero(m)  .AND.  k .LE. nz_t_source_aero(m)  .AND.                                &
                     j .GE. ny_s_source_aero(m) - j_global_LCM  .AND.  j .LE. ny_n_source_aero(m) - j_global_LCM  .AND.  &
                     i .GE. nx_w_source_aero(m) - i_global_LCM  .AND.  i .LE. nx_e_source_aero(m) - i_global_LCM .AND. &
                     temp_p_source_aero(m) .GT. 0.0)  THEN
   !
   !--               Within a source volume, check if you can add a new particle
                     IF ( TRIM(unit_source) .EQ. "1/m3/s" )  THEN
                        A_dummy = NINT( n_source_aero(m) * dV(k) * dtn_LCM / temp_p_source_aero(m), KIND=8 )
                     ELSEIF ( TRIM(unit_source) .EQ. "1/kg/s" )  THEN
                        A_dummy = NINT( n_source_aero(m) * dV(k) * rho_air(k) * dtn_LCM / temp_p_source_aero(m), KIND=8 )
                     ELSE
                        PRINT*,'LCM ERROR: unknown unit_source'
                        STOP
                     ENDIF
                     IF ( A_dummy .EQ. 0_8 )  THEN
                        PRINT*,'LCM ERROR: source produces particles with 0 weighting factor'
                        STOP
                     ENDIF
   !
   !--               Check how many particles will be released
                     x_rand = randomu(rank+1)
                     IF ( x_rand .LT. temp_p_source_aero(m) - FLOOR(temp_p_source_aero(m)) ) THEN
                        pint_source_aero = CEILING(temp_p_source_aero(m))
                     ELSE
                        pint_source_aero = FLOOR(temp_p_source_aero(m))
                     ENDIF

                     DO n = 1, pint_source_aero
   !
   !--                  New particle is added to the end of the list
                        IF ( n_end_temp(i,j,k) - n_start_temp(i,j,k) .GE. 0 )  THEN
                           !Really check whether we should be using n_end_temp and n_start_temp instead of n_end and n_start
                           n_end_temp(i,j,k)   = n_end_temp(i,j,k) + 1
                           n_particles_temp    = n_end_temp(i,j,k)
                        ELSE
                           n_start_temp(i,j,k) = n_particles_temp + 1
                           n_end_temp(i,j,k)   = n_start_temp(i,j,k)
                           n_particles_temp    = n_end_temp(i,j,k)
                        ENDIF
   !
   !--                  Mark sprayed particles
                        particle_temp(n_particles_temp)%ivar(iSPR)   = 100

                        particle_temp(n_particles_temp)%ivar(iLEM)  = MAXVAL( particle(n_start(i,j,k):n_end(i,j,k))%ivar(iLEM) ) + 1
                        particle_temp(n_particles_temp)%ivar(iilst) = i + i_global_LCM
                        particle_temp(n_particles_temp)%ivar(ijlst) = j + j_global_LCM
                        particle_temp(n_particles_temp)%ivar(iklst) = k

                        particle_temp(n_particles_temp)%ivar(iID)   = n_ID
                        n_ID = n_ID + 1_8

                        particle_temp(n_particles_temp)%var(itime) = 0.0_lcmp

                        particle_temp(n_particles_temp)%var(iX)    = REAL( ( i - 1 + nx_sb + randomu(rank+1) ) * dx_LCM, KIND=lcmp )
                        IF ( YES3D .EQ. 1 )  THEN
                           particle_temp(n_particles_temp)%var(iY) = REAL( ( j - 1 + ny_sb + randomu(rank+1) ) * dy_LCM, KIND=lcmp )
                        ELSE
                           particle_temp(n_particles_temp)%var(iY) = REAL( 0.5 * dy_LCM, KIND=lcmp )
                        ENDIF
                        particle_temp(n_particles_temp)%var(iZ)    = REAL(zi(k)+randomu(rank+1)*dz_LCM(k), kind = lcmp)
                        !particle_temp(n_particles_temp)%var(iZ)    = REAL( ( k - 1 + randomu(rank+1) ) * dz_LCM, KIND=lcmp )

                        particle_temp(n_particles_temp)%var(iU)    = 0.0_lcmp
                        particle_temp(n_particles_temp)%var(iV)    = 0.0_lcmp
                        particle_temp(n_particles_temp)%var(iW)    = 0.0_lcmp

                        particle_temp(n_particles_temp)%var(iUsgs) = 0.0_lcmp
                        particle_temp(n_particles_temp)%var(iVsgs) = 0.0_lcmp
                        particle_temp(n_particles_temp)%var(iWsgs) = 0.0_lcmp

                        particle_temp(n_particles_temp)%var(iWsedi) = 0.0_lcmp
   !
   !--                  Initialize particle properties and get radius of source particles, and adapt water/energy budget
                        SS_gb = ( micro_field(i,j,k,imv) - qsatw( tabs_LCM(i,j,k), pres_LCM(k) ) ) / qsatw( tabs_LCM(i,j,k), pres_LCM(k) )
                        SS_gb = MIN( SS_gb, -0.01 )

                        r_aero = -1.0
                        DO WHILE ( r_aero .LE. 0.001E-6 .OR. r_aero .GE. 10.0E-6 )
                           x_rand     = randomu(rank+1)
                           ierf_dummy = MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )
                           r_aero     = r_source_aero(m) * EXP( SQRT(2.0) * LOG(sigma_source_aero(m)) * ierf_dummy )
                        END DO

                        particle_temp(n_particles_temp)%ivar(iA)   = A_dummy
                        particle_temp(n_particles_temp)%var(iN)    = REAL( particle_temp(n_particles_temp)%ivar(iA) * 4.0 / 3.0 * pi * rho_aero(1) * r_aero**3, KIND=lcmp )
                        IF ( TRIM(type_source) .EQ. 'equilibrium' )  THEN
                           particle_temp(n_particles_temp)%var(iM)    = REAL( particle_temp(n_particles_temp)%ivar(iA) * 4.0 / 3.0 * pi * rho_liq * MAX( r_equi(SS_gb,tabs_LCM(i,j,k),r_aero), r_aero )**3, KIND=lcmp )
                        ELSEIF ( TRIM(type_source) .EQ. 'seawater' )  THEN
                           particle_temp(n_particles_temp)%var(iM)    = REAL( particle_temp(n_particles_temp)%ivar(iA) * 4.0 / 3.0 * pi * rho_liq * ( r_aero**3 / ( salinity_source * ( rho_liq * ( 1.0 + salinity_source ) ) / rho_aero(1) ) ), KIND=lcmp)
                        ELSE
                           PRINT*, 'LCM: type_source unknown'
                           STOP
                        ENDIF

                        ! assign source aerosol to first mode (arbitrary, but has to be something, could be made a user-specified thing)
                        particle_temp(n_particles_temp)%var(iAero_portions)    = 0.0
                        particle_temp(n_particles_temp)%var(iAero_portions(1)) = 1.0

                        particle_temp(n_particles_temp)%var(iT)   = REAL( tabs_LCM(i,j,k) + gamaz(k), KIND=lcmp ) !gamaz is better than z * ggr / cp because it avoids interpolation
                        particle_temp(n_particles_temp)%var(iq)   = REAL( micro_field(i,j,k,imv), KIND=lcmp )
                        particle_temp(n_particles_temp)%var(ieta) = REAL( micro_field(i,j,k,imv) - qsatw( tabs_LCM(i,j,k), pres_LCM(k) ), KIND=lcmp )

                        x_rand = randomu(rank+1)
                        IF ( x_rand .LT. frac_in_source_particles(m)  .AND.  switch_ice .AND. (.NOT. switch_phase_sources))  THEN
                           r_in = -1.0
                           DO WHILE ( r_in .LE. 0.001E-6 .OR. r_in .GE. 10.0E-4 )
                              x_rand     = randomu(rank+1)
                              ierf_dummy = MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )
                              r_in     = r_source_in(m) * EXP( SQRT(2.0) * LOG(sigma_source_in(m)) * ierf_dummy )
                           END DO
                           sfc_in = 4.0*pi*r_in**2
                           tfz = t_frz(type_source_in(m), sfc_in) !Surface area of each IN particle (assuming spherical shape)
                           particle_temp(n_particles_temp)%ivar(iA) = NINT(A_dummy * frac_in_source_weight(m)) !Weight the ice particles differently than the liquid ones
                           particle_temp(n_particles_temp)%var(iN) = particle_temp(n_particles_temp)%var(iN)*frac_in_source_weight(m)
                           !Adjust CCN mass to account for the different weighting
                           particle_temp(n_particles_temp)%var(iM) = particle_temp(n_particles_temp)%var(iM)*frac_in_source_weight(m)
                           !Ditto for water mass
                           particle_temp(n_particles_temp)%var(in_startind + type_source_in(m)) = REAL(A_dummy * frac_in_source_weight(m) * sfc_in)
                           !Total surface area of IN within the superdroplet
                           particle_temp(n_particles_temp)%var(iTFZ) = tfz
                           !Heterogeneous freezing temperature of the droplet
                            
                           particle_temp(n_particles_temp)%ivar(iTYPE) = 2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
                           nin_prof(k) = nin_prof(k) + NINT(A_dummy * frac_in_source_weight(m)) !Injection rate profile
                        ELSE
                           particle_temp(n_particles_temp)%var(in_startind+1:in_startind+4) = 0. !No IN surface area
                           particle_temp(n_particles_temp)%var(iTFZ) = 233.15 !-40 K; ensure only homogeneous freezing processes
                           particle_temp(n_particles_temp)%ivar(iTYPE) = 1 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
                           nin_prof(k) = nin_prof(k) + A_dummy !Injection rate profile
                        ENDIF
                        particle_temp(n_particles_temp)%ivar(iNmono) = 1

                        dql                    = particle_temp(n_particles_temp)%var(iM) / ( dV(k) * rho_air(k) )
                        micro_field(i,j,k,imv) = micro_field(i,j,k,imv) - dql
                        tabs_LCM(i,j,k)        = tabs_LCM(i,j,k) + fac_cond * dql
                        !-- Moved to after checking for IN. Previously would have overestimated temperature and vapor change
                        !-- as IN particles were added by sources - B.Ascher 26.05.2025
                     ENDDO
   
                  ENDIF
                  IF ( k .GE. nz_b_source_in(m)  .AND.  k .LE. nz_t_source_in(m)  .AND.                                &
                     j .GE. ny_s_source_in(m) - j_global_LCM  .AND.  j .LE. ny_n_source_in(m) - j_global_LCM  .AND.  &
                     i .GE. nx_w_source_in(m) - i_global_LCM  .AND.  i .LE. nx_e_source_in(m) - i_global_LCM .AND. &
                     temp_p_source_aero(m) .GT. 0.0)  THEN
         !
   !--               Within a source volume, check if you can add a new particle
                     IF ( TRIM(unit_source) .EQ. "1/m3/s" )  THEN
                        A_dummy = NINT( n_source_in(m) * dV(k) * dtn_LCM / temp_p_source_in(m), KIND=8 )
                     ELSEIF ( TRIM(unit_source) .EQ. "1/kg/s" )  THEN
                        A_dummy = NINT( n_source_in(m) * dV(k) * rho_air(k) * dtn_LCM / temp_p_source_in(m), KIND=8 )
                     ELSE
                        PRINT*,'LCM ERROR: unknown unit_source'
                        STOP
                     ENDIF
                     IF ( A_dummy .EQ. 0_8 )  THEN
                        PRINT*,'LCM ERROR: source produces particles with 0 weighting factor'
                        STOP
                     ENDIF
   !
   !--               Check how many particles will be released
                     x_rand = randomu(rank+1)
                     IF ( x_rand .LT. temp_p_source_in(m) - FLOOR(temp_p_source_in(m)) ) THEN
                        pint_source_in = CEILING(temp_p_source_in(m))
                     ELSE
                        pint_source_in = FLOOR(temp_p_source_in(m))
                     ENDIF

                     DO n = 1, pint_source_in
   !
   !--                  New particle is added to the end of the list
                        IF ( n_end_temp(i,j,k) - n_start_temp(i,j,k) .GE. 0 )  THEN
                           n_end_temp(i,j,k)   = n_end_temp(i,j,k) + 1
                           n_particles_temp    = n_end_temp(i,j,k)
                        ELSE
                           n_start_temp(i,j,k) = n_particles_temp + 1
                           n_end_temp(i,j,k)   = n_start_temp(i,j,k)
                           n_particles_temp    = n_end_temp(i,j,k)
                        ENDIF
   !
   !--                  Mark sprayed particles
                        particle_temp(n_particles_temp)%ivar(iSPR)   = 100

                        particle_temp(n_particles_temp)%ivar(iLEM)  = MAXVAL( particle(n_start(i,j,k):n_end(i,j,k))%ivar(iLEM) ) + 1
                        particle_temp(n_particles_temp)%ivar(iilst) = i + i_global_LCM
                        particle_temp(n_particles_temp)%ivar(ijlst) = j + j_global_LCM
                        particle_temp(n_particles_temp)%ivar(iklst) = k

                        particle_temp(n_particles_temp)%ivar(iID)   = n_ID
                        n_ID = n_ID + 1_8

                        particle_temp(n_particles_temp)%var(itime) = 0.0_lcmp

                        particle_temp(n_particles_temp)%var(iX)    = REAL( ( i - 1 + nx_sb + randomu(rank+1) ) * dx_LCM, KIND=lcmp )
                        IF ( YES3D .EQ. 1 )  THEN
                           particle_temp(n_particles_temp)%var(iY) = REAL( ( j - 1 + ny_sb + randomu(rank+1) ) * dy_LCM, KIND=lcmp )
                        ELSE
                           particle_temp(n_particles_temp)%var(iY) = REAL( 0.5 * dy_LCM, KIND=lcmp )
                        ENDIF
                        particle_temp(n_particles_temp)%var(iZ)    = REAL(zi(k)+randomu(rank+1)*dz_LCM(k), kind = lcmp)
                        !particle_temp(n_particles_temp)%var(iZ)    = REAL( ( k - 1 + randomu(rank+1) ) * dz_LCM, KIND=lcmp )

                        particle_temp(n_particles_temp)%var(iU)    = 0.0_lcmp
                        particle_temp(n_particles_temp)%var(iV)    = 0.0_lcmp
                        particle_temp(n_particles_temp)%var(iW)    = 0.0_lcmp

                        particle_temp(n_particles_temp)%var(iUsgs) = 0.0_lcmp
                        particle_temp(n_particles_temp)%var(iVsgs) = 0.0_lcmp
                        particle_temp(n_particles_temp)%var(iWsgs) = 0.0_lcmp

                        particle_temp(n_particles_temp)%var(iWsedi) = 0.0_lcmp
   !
   !--                  Initialize particle properties and get radius of source particles, and adapt water/energy budget
                        SS_gb = ( micro_field(i,j,k,imv) - qsatw( tabs_LCM(i,j,k), pres_LCM(k) ) ) / qsatw( tabs_LCM(i,j,k), pres_LCM(k) )
                        SS_gb = MIN( SS_gb, -0.01 )

                        r_aero = -1.0
                        DO WHILE ( r_aero .LE. 0.001E-6 .OR. r_aero .GE. 10.0E-6 )
                           x_rand     = randomu(rank+1)
                           ierf_dummy = MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )
                           r_aero     = r_source_aero(m) * EXP( SQRT(2.0) * LOG(sigma_source_aero(m)) * ierf_dummy )
                        END DO

                        r_in = -1.0
                        DO WHILE ( r_in .LE. 0.001E-6 .OR. r_in .GE. 100.0E-6 )
                           x_rand     = randomu(rank+1)
                           ierf_dummy = MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )
                           r_in     = r_source_in(m) * EXP( SQRT(2.0) * LOG(sigma_source_in(m)) * ierf_dummy )
                        END DO

                        !B. Ascher 06.05.2025: just use the same routine as for CCN sources, so that we don't have a monodispersive IN-containing droplet source

                        particle_temp(n_particles_temp)%ivar(iA)   = A_dummy
                        particle_temp(n_particles_temp)%var(iN)    = REAL( particle_temp(n_particles_temp)%ivar(iA) * 4.0 / 3.0 * pi * rho_aero(1) * r_aero**3, KIND=lcmp )
 
                        ! assign source aerosol to first mode (arbitrary, but has to be something, could be made a user-specified thing)
                        particle_temp(n_particles_temp)%var(iAero_portions)    = 0.0
                        particle_temp(n_particles_temp)%var(iAero_portions(1)) = 1.0

                        particle_temp(n_particles_temp)%var(iM) = REAL(particle_temp(n_particles_temp)%ivar(iA) * 4.0/3.0 * pi * rho_liq * MAX( r_equi(SS_gb,tabs_LCM(i,j,k),r_aero), r_aero )**3, KIND=lcmp )

                        dql                    = particle_temp(n_particles_temp)%var(iM) / ( dV(k) * rho_air(k) )
                        micro_field(i,j,k,imv) = micro_field(i,j,k,imv) - dql
                        tabs_LCM(i,j,k)        = tabs_LCM(i,j,k) + fac_cond * dql

                        particle_temp(n_particles_temp)%var(iT)   = REAL( tabs_LCM(i,j,k) + gamaz(k), KIND=lcmp ) !gamaz is better than z * ggr / cp because it avoids interpolation
                        particle_temp(n_particles_temp)%var(iq)   = REAL( micro_field(i,j,k,imv), KIND=lcmp )
                        particle_temp(n_particles_temp)%var(ieta) = REAL( micro_field(i,j,k,imv) - qsatw( tabs_LCM(i,j,k), pres_LCM(k) ), KIND=lcmp )
                        sfc_in = 4*pi*r_in**2 !Surface area of each IN particle (NOT the total IN surface area within the supedroplet)
                        tfz = t_frz(type_source_in(m), sfc_in)
                        particle_temp(n_particles_temp)%var(in_startind + type_source_in(m)) = REAL(A_dummy * sfc_in)
                        particle_temp(n_particles_temp)%var(iTFZ) = tfz
                        particle_temp(n_particles_temp)%ivar(ITYPE) = 2
                        particle_temp(n_particles_temp)%ivar(iNmono) = 1
                        nin_prof(k) = nin_prof(k) + A_dummy
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
   !
   !--   Store temporary particle array to main particle array
         IF ( n_particles_temp .GT. n_particles )  THEN
            DEALLOCATE( particle )
            n_particles       = n_particles_temp
            n_particles_alloc = n_particles_temp + MAX( n_particles_temp / INT(SQRT(REAL(n_particles_temp))) * 6 , 10 )
            ALLOCATE( particle(1:n_particles_alloc), STAT=ierr )
            IF ( ierr /= 0 ) THEN
               PRINT*,'LCM ERROR: Failed to allocate arrays in micro_particle_source on proc:', rank
               CALL micro_task_abort
            ENDIF
            DO n = 1, n_particles
               particle(n)%var  = particle_temp(n)%var
               particle(n)%ivar = particle_temp(n)%ivar
            ENDDO
            n_start = n_start_temp
            n_end   = n_end_temp
         ENDIF
         DEALLOCATE( particle_temp )
   
       ENDDO !end of multiple particle source
   
    END SUBROUTINE micro_particle_source

 SUBROUTINE micro_init_particles(i,j,k)
!
!-- Initializes aerosol sizes
    USE micro_task_util_MPI, ONLY: micro_task_abort
    USE micro_vars,   ONLY: particle, iA, iM, iN, iq, ieta,                    &
                            np_xy, dV, n_init, rho_liq, n_start, n_end, iT,    &
                            lcmp, rm_aero, sigma_aero, n_aero, rho_aero,       &
                            molecular_weight_aero, vanthoff_aero, shares_aero, &
                            sigma, molecular_weight_water, tabs_LCM,           &
                            pres_LCM, rho_air, init_aero, iZ, iTYPE,           &
                            switch_ice, frac_ice, n_ice, iNmono, unit_aero,    &
                            rm_in, sigma_in, iTFZ, in_startind, t_frz,         &
                            iAero_portions, get_bfactor, get_r_aero, r_equi,   &
                            inv_cdf, n_inv_cdf
    USE microphysics, ONLY: micro_field, imv
    USE micro_random, ONLY: randomu

    USE params,       ONLY: pi, diffelq, fac_cond, cp, ggr, rv
    USE grid,         ONLY: rank, masterproc
    USE vars,         ONLY: qsatw, t, gamaz

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: i, j, k
    INTEGER             :: n, nn, m, m_max, liqpart, intype, iMode, np_mode, np_done=0

    REAL                :: dlogr, r_l, r_r, r_aero, SS_gb, dql,        &
                           log_normal_pdf, log_normal_cdf, ERFI, x_rand,       &
                           ierf_dummy, IERF, gauss, N_dummy, r_r_m1, r_r_m2,   &
                           cdf_m1, cdf_m2, p_USM, f_ak, f_bk, f_ck, ak, bk, ck,&
                           p_random, f_m1, f_m2, icepart, afactor,    &
                           bfactor, koehler_requi, sumvols, r_in, sfc_in

    REAL,    SAVE       :: r_aero_min(4), r_aero_max(4)

    LOGICAL, SAVE       :: init_aero_finished = .FALSE.

   !  REAL, SAVE          :: inv_incdf(1:4,0:1000)
    REAL                :: in_typeprobs(4) !Probability of a particle with IN receiving a certain IN species.
                                           !B. Ascher 22.05.2025 - Consider making this allocatable 
                                           !to allow for an arbitrary number of IN species

    IF ( np_xy .EQ. 0 )  RETURN
    
    IF ( .NOT. init_aero_finished )  THEN
      IF (switch_ice) THEN
         in_typeprobs = [(sum(n_ice(1:m)), m = 1,4)]/sum(n_ice(:)) !B. Ascher 22.05.2025 - Same as above re arbitrary IN species
      ENDiF

       np_done = 0
       ! initialize aerosols for every mode separately (for the case of different compositions)
       DO iMode = 1, 4

          ! determine how many particles this mode gets in this grid box
          np_mode = FLOOR(np_xy * shares_aero(iMode) / SUM(shares_aero))
          ! add the remainder from rounding to the first mode
          IF (iMode==1) THEN
            n = 0
            DO nn = 1, 4
              n = n + FLOOR(np_xy * shares_aero(nn) / SUM(shares_aero))
            END DO
            np_mode = np_mode + np_xy - n
          END IF

          ! if no particles for this mode, cycle this loop iteration
          IF (np_mode==0) CYCLE

!
!--       Determine r_aero_min and r_aero_max using regula falsi
          r_aero_min(iMode) = 1.0E-9
          r_aero_max(iMode) = 1.0E-5

          ak   = r_aero_min(iMode)
          f_ak = log_normal_cdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),ak) / n_aero(iMode) - 5.0E-5

          bk   = r_aero_max(iMode)
          f_bk = log_normal_cdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),bk) / n_aero(iMode) - 5.0E-5

          DO

             ck   = ( ak * f_bk - bk * f_ak ) / ( f_bk - f_ak )
             f_ck = log_normal_cdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),ck) / n_aero(iMode) - 5.0E-5

             IF ( f_ak * f_ck .GT. 0.0 )  THEN
                ak   = ck
                f_ak = f_ck
             ELSE
                bk   = ck
                f_bk = f_ck
             ENDIF

             IF ( ABS( f_ck ) .LT. 1.0E-12 )  EXIT

          END DO
          r_aero_min(iMode) = ck

          ak   = r_aero_min(iMode)
          f_ak = log_normal_cdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),ak) / n_aero(iMode) - ( 1.0 - 5.0E-5 )

          bk   = r_aero_max(iMode)
          f_bk = log_normal_cdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),bk) / n_aero(iMode) - ( 1.0 - 5.0E-5 )

          DO

             ck   = ( ak * f_bk - bk * f_ak ) / ( f_bk - f_ak )
             f_ck = log_normal_cdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),ck) / n_aero(iMode) - ( 1.0 - 5.0E-5 )

             IF ( f_ak * f_ck .GT. 0.0 )  THEN
                ak   = ck
                f_ak = f_ck
             ELSE
                bk   = ck
                f_bk = f_ck
             ENDIF

             IF ( ABS( f_ck ) .LT. 1.0E-12 )  EXIT

          END DO
          r_aero_max(iMode) = ck

          IF ( TRIM(init_aero) .EQ. 'bin' )  THEN
!
!--          Aerosols are crated using a representation of a lognormal distribution by bins
             dlogr = ( LOG(r_aero_max(iMode)) - LOG(r_aero_min(iMode)) ) / np_mode

            DO n = np_done+1, np_done+np_mode !n_start(i,j,k), n_end(i,j,k)

               nn = n - np_done !n - n_start(i,j,k) + 1

               r_l    = EXP( LOG( r_aero_min(iMode) ) + (nn-1) * dlogr )
               r_r    = EXP( LOG( r_aero_min(iMode) ) + nn * dlogr )
               r_aero = SQRT( r_l * r_r )

               N_dummy = log_normal_pdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),r_aero)
               IF ( TRIM(unit_aero) .EQ. "1/m3" )  THEN
                  particle(n)%ivar(iA) = INT( N_dummy * r_aero * dlogr * dV(k), KIND=8 )
               ELSEIF ( TRIM(unit_aero) .EQ. "1/kg" )  THEN
                  particle(n)%ivar(iA) = INT( N_dummy * r_aero * dlogr * dV(k) * rho_air(k), KIND=8 )
               ELSE
                  PRINT*,'LCM ERROR: unknown unit_aero'
                  STOP
               ENDIF
               IF (switch_ice) particle(n)%ivar(iA) = particle(n)%ivar(iA) / (1.0 - frac_ice) ! Adjust multiplicity of CCN-containing particles based on IN particle fraction
               particle(n)%var(iN) = REAL(particle(n)%ivar(iA) * 4.0 / 3.0 * pi * r_aero**3 * rho_aero(iMode), KIND=lcmp)

               particle(n)%var(iAero_portions)        = 0.0
               particle(n)%var(iAero_portions(iMode)) = 1.0
               
               IF (switch_ice .and. (randomu(rank+1) .lt. frac_ice)) THEN
                  particle(n)%ivar(iTYPE) = 2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
                  x_rand = randomu(rank+1)
                  !First random variable; determine which species of IN the particle will get
                  intype = minloc(in_typeprobs - x_rand, mask = in_typeprobs>=x_rand, dim = 1)
                  N_dummy = n_ice(intype) * dV(k) / (np_xy * frac_ice)
                  !Each IN type has its own frequency of occurence; some IN should correspond to a greater
                  !superparticle multiplicity than others
                  IF ( TRIM(unit_aero) .EQ. "1/kg" )  THEN
                     N_dummy = N_dummy*rho_air(k)
                  ENDIF
                  particle(n)%var(iN) = N_dummy * 4.0 / 3.0 * pi * rho_aero(iMode) * r_aero**3
                  particle(n)%ivar(iA) = INT(N_dummy)
                  !Greatly reduce both aerosol mass and particle multiplicity to account for there being many fewer IN
                  !than CCN
                  !Use lognormal sampling to determine which size IN the droplet has
                  r_in = -1.0
                  DO WHILE ( r_in .LE. 0.001E-6 .OR. r_in .GE. 10.0E-4 )
                     x_rand     = randomu(rank+1)
                     ierf_dummy = MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )
                     r_in     = rm_in(intype) * EXP( SQRT(2.0) * LOG(sigma_in(intype)) * ierf_dummy )
                  END DO
                  sfc_in = 4.0*pi*r_in**2 !Area of each IN (assuming single spherical particle in each droplet)
                  !Now, determine the freezing temperature of the particle
                  particle(n)%var(in_startind + intype) = N_dummy * sfc_in !total surface area of IN in superdroplet
                  particle(n)%var(iTFZ) = t_frz(intype, sfc_in)
               ELSE
                  particle(n)%ivar(iTYPE) = 1
                  particle(n)%var(in_startind+1:in_startind+4) = 0. !0 corresponds to no IN
                  particle(n)%var(iTFZ) = 233.15 !This ensures that the droplet will freeze by homogeneous nucleation
               ENDIF
            END DO
!
!--          If a bin does not contain any aerosols, the assigned computational
!--          particle is reassigned to represent the bin with the highest number of
!--          aerosols
             DO n = np_done+1, np_done+np_mode !n_start(i,j,k), n_end(i,j,k)

                IF ( particle(n)%ivar(iA) .EQ. 0 )  THEN

                   m_max  = n
                   DO m = np_done+1, np_done+np_mode !n_start(i,j,k), n_end(i,j,k)

                      IF ( particle(m)%ivar(iA) .GT. particle(m_max)%ivar(iA)  .AND.  particle(m)%ivar(iA) .GE. 2 )  THEN
                         m_max  = m
                      ENDIF

                   ENDDO

                   IF ( particle(m_max)%ivar(iA) .EQ. 0 )  THEN
                      PRINT*, 'LCM ERROR: Initialization of aerosols not possible.'
                      STOP
                   ENDIF

                   particle(n)%ivar(iA)     = particle(m_max)%ivar(iA) / 2
                   particle(m_max)%ivar(iA) = particle(m_max)%ivar(iA) - particle(n)%ivar(iA)

                   particle(n)%var(iN)      = particle(m_max)%var(iN) / ( particle(m_max)%ivar(iA) + particle(n)%ivar(iA) ) * particle(n)%ivar(iA)
                   particle(m_max)%var(iN)  = particle(m_max)%var(iN) - particle(n)%var(iN)

                   particle(n)%var(iAero_portions) = particle(m_max)%var(iAero_portions)
                   particle(n)%ivar(iTYPE)         = particle(m_max)%ivar(iTYPE)

                   IF ( particle(n)%ivar(iA) .EQ. 0  .OR.  particle(m_max)%ivar(iA) .EQ. 0 )  THEN
                      PRINT*, 'LCM ERROR: Initialization of aerosols not possible.'
                      STOP
                   ENDIF

                ENDIF

          ENDDO
       ELSEIF ( TRIM(init_aero) .EQ. 'random' )  THEN
!
!--          Aerosol radii are chosen by a random generator following a lognormal distribution
             p_random = 1.0 / n_inv_cdf

             r_l               = r_aero_min(iMode)
             inv_cdf(iMode, 0) = r_aero_min(iMode)

             DO n = 1, n_inv_cdf-1
!
!--              Use regula falsi method to find the (n * p_random)-percentile of the distribution.
                 ak   = r_l
                 f_ak = log_normal_cdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),ak) / n_aero(iMode) - p_random * n

                 bk   = r_aero_max(iMode)
                 f_bk = log_normal_cdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),bk) / n_aero(iMode) - p_random * n

                 DO

                    ck   = ( ak * f_bk - bk * f_ak ) / ( f_bk - f_ak )
                    f_ck = log_normal_cdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),ck) / n_aero(iMode) - p_random * n

                    IF ( f_ak * f_ck .GT. 0.0 )  THEN
                       ak   = ck
                       f_ak = f_ck
                    ELSE
                       bk   = ck
                       f_bk = f_ck
                    ENDIF

                    IF ( ABS( f_ck ) .LT. 1.0E-12 )  EXIT

                 END DO

                 r_r        = ck
                 r_aero     = SQRT( r_l * r_r )
                 r_l        = r_r
                 inv_cdf(iMode, n) = r_aero

             ENDDO

             inv_cdf(iMode, n_inv_cdf) = r_aero_max(iMode)

       ELSEIF ( TRIM(init_aero) .EQ. 'USM' )  THEN
!
!--            Uniform sampling method by Shima et al. (2020), Eq. (85).
            p_USM = 1.0 / np_mode ! the first and last interval are ignored since they usually contain unrealistically small or large droplets

            r_l = r_aero_min(iMode)

            DO n = np_done+1, np_done+np_mode-1
               nn = n - np_done
!
!--               Use regula falsi method to find the (n * p_USM)-percentile of the distribution.
               ak   = r_l
               f_ak = log_normal_cdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),ak) / n_aero(iMode) - p_USM * nn

               bk   = r_aero_max(iMode)
               f_bk = log_normal_cdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),bk) / n_aero(iMode) - p_USM * nn

               DO

                  ck   = ( ak * f_bk - bk * f_ak ) / ( f_bk - f_ak )
                  f_ck = log_normal_cdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),ck) / n_aero(iMode) - p_USM * nn

                  IF ( f_ak * f_ck .GT. 0.0 )  THEN
                     ak   = ck
                     f_ak = f_ck
                  ELSE
                     bk   = ck
                     f_bk = f_ck
                  ENDIF

                  IF ( ABS( f_ck ) .LT. 1.0E-12 )  EXIT

               END DO

               r_r    = ck
               r_aero = SQRT( r_l * r_r )
               r_l    = r_r

               IF ( TRIM(unit_aero) .EQ. "1/m3" )  THEN
                  particle(n)%ivar(iA) = INT( dV(k) * n_aero(iMode) / np_mode, KIND=8 )
               ELSEIF ( TRIM(unit_aero) .EQ. "1/kg" )  THEN
                  particle(n)%ivar(iA) = INT( dV(k) * rho_air(k) * n_aero(iMode) / np_mode, KIND=8 )
               IF (switch_ice) particle(n)%ivar(iA) = particle(n)%ivar(iA) / (1.0 - frac_ice)
               ELSE
                  PRINT*,'LCM ERROR: unknown unit_aero'
                  STOP
               ENDIF
               particle(n)%var(iN)  = REAL( particle(n)%ivar(iA) * 4.0 / 3.0 * pi * rho_aero(iMode) * r_aero**3, KIND=lcmp )

               ! assign aerosol to mode number iMode
               particle(n)%var(iAero_portions)        = 0.0
               particle(n)%var(iAero_portions(iMode)) = 1.0

               IF (switch_ice .and. (randomu(rank+1) .lt. frac_ice)) THEN
                  particle(n)%ivar(iTYPE) = 2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
                  x_rand = randomu(rank+1)
                  !First random variable; determine which species of IN the particle will get
                  intype = minloc(in_typeprobs - x_rand, mask = in_typeprobs>=x_rand, dim = 1)
                  N_dummy = n_ice(intype) * dV(k) / (np_xy * frac_ice)
                  !Each IN type has its own frequency of occurence; some IN should correspond to a greater
                  !superparticle multiplicity than others
                  IF ( TRIM(unit_aero) .EQ. "1/kg" )  THEN
                     N_dummy = N_dummy*rho_air(k)
                  ENDIF
                  particle(n)%var(iN) = N_dummy * 4.0 / 3.0 * pi * rho_aero(iMode) * r_aero**3
                  particle(n)%ivar(iA) = INT(N_dummy)
                  !Greatly reduce both aerosol mass and particle multiplicity to account for there being many fewer IN
                  !than CCN
                  !Use lognormal sampling to determine which size IN the droplet has
                  r_in = -1.0
                  DO WHILE ( r_in .LE. 0.001E-6 .OR. r_in .GE. 10.0E-4 )
                     x_rand     = randomu(rank+1)
                     ierf_dummy = MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )
                     r_in     = rm_in(intype) * EXP( SQRT(2.0) * LOG(sigma_in(intype)) * ierf_dummy )
                  END DO
                  sfc_in = 4.0 * pi * r_in**2
                  !Now, determine the freezing temperature of the particle
                  particle(n)%var(in_startind+intype) = N_dummy * sfc_in !total surface area of IN in superdroplet
                  particle(n)%var(iTFZ) = t_frz(intype, sfc_in)
               ELSE
                  particle(n)%ivar(iTYPE) = 1
                  particle(n)%var(in_startind+1:in_startind+4) = 0. !No IN means no IN surface area
                  particle(n)%var(iTFZ) = 233.15 !This ensures that the droplet will freeze by homogeneous nucleation
               ENDIF

            ENDDO

!
!--       n = np_xy requires special treatment to avoid overestimation of the largest aerosols
            r_r    = r_aero_max(iMode)
            r_aero = SQRT( r_l * r_r )

            particle(np_done+np_mode)%ivar(iA) = INT( dV(k) * SUM(n_aero(:)) / (np_xy*(1.0 - frac_ice)), KIND=8 )
            particle(np_done+np_mode)%var(iN)  = REAL( particle(np_done+np_mode)%ivar(iA) * 4.0 / 3.0 * pi * rho_aero(iMode) * r_aero**3, KIND=lcmp )
            IF (switch_ice .and. (randomu(rank+1) .lt. frac_ice)) THEN
               particle(np_done+np_mode)%ivar(iTYPE) = 2 ! > 0: liquid; < 0: ice; |2|: with IN; |1|: without IN
               x_rand = randomu(rank+1)
               !First random variable; determine which species of IN the particle will get
               intype = minloc(in_typeprobs - x_rand, mask = in_typeprobs>=x_rand, dim = 1)
               N_dummy = n_ice(intype) * dV(k) / (np_xy * frac_ice)
               !Each IN type has its own frequency of occurence; some IN should correspond to a greater
               !superparticle multiplicity than others
               IF ( TRIM(unit_aero) .EQ. "1/kg" )  THEN
                  N_dummy = N_dummy*rho_air(k)
               ENDIF
               particle(np_done+np_mode)%var(iN) = N_dummy * 4.0 / 3.0 * pi * rho_aero(iMode) * r_aero**3
               particle(np_done+np_mode)%ivar(iA) = INT(N_dummy)
               !Greatly reduce both aerosol mass and particle multiplicity to account for there being many fewer IN
               !than CCN
               !Use lognormal sampling to determine which size IN the droplet has
               r_in = -1.0
               DO WHILE ( r_in .LE. 0.001E-6 .OR. r_in .GE. 10.0E-4 )
                  x_rand     = randomu(rank+1)
                  ierf_dummy = MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )
                  r_in     = rm_in(intype) * EXP( SQRT(2.0) * LOG(sigma_in(intype)) * ierf_dummy )
               END DO
               sfc_in = 4.0 * pi * r_in**2
               !Now, determine the freezing temperature of the particle
               particle(np_done+np_mode)%var(in_startind + intype) = N_dummy * sfc_in !total surface area of IN in superdroplet
               particle(np_done+np_mode)%var(iTFZ) = t_frz(intype, sfc_in)
            ELSE
               particle(np_done+np_mode)%ivar(iTYPE) = 1
               particle(np_done+np_mode)%var(in_startind+1:in_startind+4) = 0. !No IN means no IN surface area
               particle(np_done+np_mode)%var(iTFZ) = 233.15 !This ensures that the droplet will freeze by homogeneous nucleation
            ENDIF
         ENDIF
          ! track how many particles are done from previous modes already
          np_done = np_done + np_mode
      END DO
 

       init_aero_finished = .TRUE.

    ENDIF

    np_done = 0
    DO iMode = 1, 4
       ! determine how many particles this mode gets in this grid box
       np_mode = FLOOR(np_xy * shares_aero(iMode) / SUM(shares_aero))
       ! add the remainder from rounding to the first mode
       IF (iMode==1) THEN
         n = 0
         DO nn = 1, 4
           n = n + FLOOR(np_xy * shares_aero(nn) / SUM(shares_aero))
         END DO
         np_mode = np_mode + np_xy - n
       END IF

       ! if no particles for this mode, cycle this loop iteration
       IF (np_mode==0) CYCLE

       DO n = n_start(i,j,k)+np_done, n_start(i,j,k)+np_done+np_mode-1

          IF ( TRIM(init_aero) .EQ. 'USM'  .OR.   TRIM(init_aero) .EQ. 'bin' )  THEN
             ! particles in every grid box are exactly alike the particles in the first grid box
             m = n - n_start(i,j,k) + 1
             particle(n)%ivar(iA)            = particle(m)%ivar(iA)
             particle(n)%var(iN)             = particle(m)%var(iN)
             particle(n)%var(iAero_portions) = particle(m)%var(iAero_portions)
             particle(n)%ivar(iTYPE)         = particle(m)%ivar(iTYPE)
             IF (switch_ice) THEN
               particle(n)%var(in_startind+1:in_startind+4) = particle(m)%var(in_startind+1:in_startind+4)
            particle(n)%var(iTFZ) = particle(m)%var(iTFZ)
             ENDIF
          ENDIF

         IF ( TRIM(init_aero) .EQ. 'random'  )  THEN
            ! generate radius that is not out of the borders of the pdf
            ! (borders are determined by predetermined multiplicity of the particle, should not exceed pdf value)
            find_r: DO nn = 1, 1000 ! 1000 is arbitrary, just a couple tries
              x_rand = randomu(rank+1)
              r_aero = inv_cdf(iMode, MAX(MIN(NINT(x_rand*n_inv_cdf),n_inv_cdf),0))

              IF ( TRIM(unit_aero) .EQ. "1/m3" )  THEN
                 N_dummy = INT( dV(k) * n_aero(iMode) / (np_mode * (1.0 - frac_ice)), KIND=8 )
              ELSEIF ( TRIM(unit_aero) .EQ. "1/kg" )  THEN
                 N_dummy = INT( dV(k) * rho_air(k) * n_aero(iMode) / (np_mode * (1.0 - frac_ice)), KIND=8 )
              ENDIF

              ! test if N is acceptable for generated radius
              IF ( N_dummy < log_normal_pdf(n_aero(iMode),sigma_aero(iMode),rm_aero(iMode),r_aero) ) THEN
                 ! pdf-conform radius generated -> take it
                 EXIT find_r
              END IF
            END DO find_r

            IF ( TRIM(unit_aero) .EQ. "1/m3" )  THEN
               particle(n)%ivar(iA) = INT( dV(k) * n_aero(iMode) / (np_mode * (1.0 - frac_ice)), KIND=8 )
            ELSEIF ( TRIM(unit_aero) .EQ. "1/kg" )  THEN
               particle(n)%ivar(iA) = INT( dV(k) * n_aero(iMode) * rho_air(k) / (np_mode * (1.0 - frac_ice)), KIND=8 )
            ELSE
               PRINT*,'LCM ERROR: unknown unit_aero'
               STOP
            ENDIF


            IF (switch_ice) particle(n)%ivar(iA) = particle(n)%ivar(iA) / (1.0 - frac_ice) ! Adjust multiplicity of CCN-containing particles based on IN particle fraction
            particle(n)%var(iN)     = REAL( particle(n)%ivar(iA) * 4.0 / 3.0 * pi * rho_aero(iMode) * r_aero**3, KIND=lcmp )

            particle(n)%var(iAero_portions)        = 0.0
            particle(n)%var(iAero_portions(iMode)) = 1.0

            particle(n)%ivar(iTYPE) = 1


            IF (switch_ice .and. (randomu(rank+1) .lt. frac_ice)) THEN
               particle(n)%ivar(iTYPE) = 2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
               x_rand = randomu(rank+1)
               !First random variable; determine which species of IN the particle will get
               intype = minloc(in_typeprobs - x_rand, mask = in_typeprobs>=x_rand, dim = 1)
               N_dummy = n_ice(intype) * dV(k) / (np_xy * frac_ice)
               !Each IN type has its own frequency of occurence; some IN should correspond to a greater
               !superparticle multiplicity than others
               IF ( TRIM(unit_aero) .EQ. "1/kg" )  THEN
                  N_dummy = N_dummy*rho_air(k)
               ENDIF
               particle(n)%var(iN) = N_dummy * 4.0 / 3.0 * pi * rho_aero(iMode) * r_aero**3
               particle(n)%ivar(iA) = INT(N_dummy)
               !Greatly reduce both aerosol mass and particle multiplicity to account for there being many fewer IN
               !than CCN
               !Use lognormal sampling to determine which size IN the droplet has
               r_in = -1.0
               DO WHILE ( r_in .LE. 0.001E-6 .OR. r_in .GE. 10.0E-4 )
                  x_rand     = randomu(rank+1)
                  ierf_dummy = MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )
                  r_in     = rm_in(intype) * EXP( SQRT(2.0) * LOG(sigma_in(intype)) * ierf_dummy )
               END DO
               sfc_in = 4.0 * pi * r_in**2
               !Now, determine the freezing temperature of the particle
               particle(n)%var(in_startind + intype) = N_dummy * sfc_in !total surface area of IN in superdroplet
               particle(n)%var(iTFZ) = t_frz(intype, sfc_in)
            ELSE
               particle(n)%ivar(iTYPE) = 1
               particle(n)%var(in_startind+1:in_startind+4) = 0. !No IN means no IN surface area
               particle(n)%var(iTFZ) = 233.15 !This ensures that the droplet will freeze by homogeneous nucleation
            ENDIF
         ENDIF
      ENDDO
      np_done = np_done + np_mode
    ENDDO
!
!-- Get equilibrium radius of aerosols/droplets and adopt the fields of qv
!-- and T accordingly. Five iterations seem to be ok.
    DO nn = 1, 5

       SS_gb = ( micro_field(i,j,k,imv) - qsatw( tabs_LCM(i,j,k), pres_LCM(k) ) ) / qsatw( tabs_LCM(i,j,k), pres_LCM(k) )
       SS_gb = MIN( SS_gb, -0.01 )

       dql = 0.0
       DO n = n_start(i,j,k), n_end(i,j,k)
          r_aero = get_r_aero(n)
          dql    = dql - particle(n)%var(iM)
          particle(n)%var(iM)    = REAL( particle(n)%ivar(iA) * 4.0 / 3.0 * pi * rho_liq * MAX( r_equi(SS_gb,tabs_LCM(i,j,k),0.0,n), r_aero )**3, KIND=lcmp )
          dql    = dql + particle(n)%var(iM)

       ENDDO
       dql = dql / ( dV(k) * rho_air(k) )

       micro_field(i,j,k,imv) = micro_field(i,j,k,imv) - dql
       tabs_LCM(i,j,k)        = tabs_LCM(i,j,k) + fac_cond * dql

    ENDDO

    DO n = n_start(i,j,k), n_end(i,j,k)
       particle(n)%var(iT)      = REAL( tabs_LCM(i,j,k) + gamaz(k), KIND=lcmp ) !gamaz is better than z * ggr / cp because it avoids interpolation
       particle(n)%var(iq)      = REAL( micro_field(i,j,k,imv), KIND=lcmp )
       particle(n)%var(ieta)    = REAL( micro_field(i,j,k,imv) - qsatw( tabs_LCM(i,j,k), pres_LCM(k) ), KIND=lcmp )
       particle(n)%ivar(iNmono) = 1_8
    END DO

    CALL micro_user_init_particles(i,j,k)

 END SUBROUTINE micro_init_particles

 REAL FUNCTION t_frz (m, sfc)

   USE micro_random, ONLY : randomu
   USE grid, ONLY : rank

   IMPLICIT NONE

   REAL, INTENT(IN) :: sfc
   INTEGER, INTENT(IN) :: m
   REAL, DIMENSION(289) :: tempc, frzprob_dust, frzprob_soot, frzprob_snowmax, frzprob_pollen, &
                           ns_dust, ns_soot, ns_snowmax, ns_pollen
   REAL :: tfrz_dust, tfrz_soot, tfrz_snowmax, tfrz_pollen, tfrz, randnum

   INTEGER :: tcounter, loc14

   tempc(1:289) = [(-38. + (36./288.)*tcounter, tcounter = 0, 288)]
   loc14 = MINLOC(ABS(tempc+14.), dim = 1) !Need to split up the pollen parameterization at 14 C, due to a sharp bend in the data there
   !as 'Class B' ice-nucleating macromolecules become active

   ns_dust = exp(-0.517*MAX(-30., MIN(tempc, -14.)) + 8.781) 
   !INP site surface density following Eq. 5 of Ulrich et al. (2017)
   ns_soot = 7.463 * exp(-0.0101*MAX(-34., MIN(tempc, -18.))**2 - 0.8525*MAX(-34., MIN(tempc, -18.)) + 0.7667) 
   !INP site surface density following Eq. 6 of Ulrich et al. (2017)
   ns_snowmax = 4.53E11*(1-exp(-2E-10*exp(-MAX(-20., MIN(tempc, -2.))*2.34))) !INP site surface density, following Eq. 6 of Wex et al. (2015) and
   !assuming a bacterium diameter of 1.44 um (Wex et al. (2015) Fig. 9 where lambda = 3) and a bacterium dry density of 1350 kg/m^3
   ns_pollen(1:loc14) = EXP(14.701 + 6.0842*ATAN(-15.8840*MAX(MIN(tempc(1:loc14), -14.), -35.)) + 1.4573*LOG(-0.45627*MAX(MIN(tempc(1:loc14),-14.),-35.)))
   ns_pollen(loc14+1:289) = EXP(-17.864 - 6.9238*MAX(MIN(tempc(loc14+1:289), -2.), -14.) - 0.59183*MAX(MIN(tempc(loc14+1:289), -2.), -14.) - 1.7131E-2*MAX(MIN(tempc(loc14+1:289), -2.), -14.))
   !INP site surface density of birch (Betula pendata) pollen. Derived from measurements of Wieland et al. (2025), 
   !and converted from per mass to per surface area using birch pollen dimension estimates from Tatiana et al. (2014) 
   !and pollen grain masses from Brown and Irving (1973).
   frzprob_dust = 1.-exp(-sfc * ns_dust)
   frzprob_soot = 1.-exp(-sfc * ns_soot)
   frzprob_snowmax = 1.-exp(-sfc * ns_snowmax)
   frzprob_pollen = 1.-exp(-sfc * ns_pollen)

   randnum = randomu(rank+1)
   SELECT CASE (m)
   CASE (1)
       t_frz = tempc(MINLOC(ABS(frzprob_soot-randnum), dim = 1))
   CASE (2)
       t_frz = tempc(MINLOC(ABS(frzprob_dust-randnum), dim = 1))
   CASE (3)
       t_frz = tempc(MINLOC(ABS(frzprob_snowmax-randnum), dim = 1))
   CASE (4)
       t_frz = tempc(MINLOC(ABS(frzprob_pollen-randnum), dim = 1))
   CASE DEFAULT
       t_frz = -40.
   END SELECT
   t_frz = MAX(t_frz, -40.) !Make sure that regardless of IN surface area, t_frz doesn't go below -40 Celsius.
   !Shouldn't matter, because the homogeneous freezing routine should begin to dominate around -38, but want to be sure
   t_frz = t_frz + 273.15 !Convert from Celsius to Kelvin
 END FUNCTION t_frz

 REAL FUNCTION log_normal_pdf(n,s,m,x)

    USE params, ONLY: pi

    IMPLICIT NONE

    REAL, INTENT(IN) :: n, s, m, x
    INTEGER          :: nn

    log_normal_pdf = 0.0
    IF ( n .NE. 0.0  .AND.  s .NE. 0.0  .AND.  m .NE. 0.0 )  THEN
       log_normal_pdf = log_normal_pdf +                            &
                        n / ( SQRT( 2.0 * pi ) * LOG(s) * x ) *     &
                        EXP( -1.0 * LOG( x / m )**2 / ( 2.0 * LOG(s)**2 ) )
    ENDIF

 END FUNCTION log_normal_pdf

 REAL FUNCTION log_normal_cdf(n,s,m,x)

    IMPLICIT NONE

    REAL, INTENT(IN) :: n, s, m, x

    log_normal_cdf = 0.0
    IF ( n .NE. 0.0  .AND.  s .NE. 0.0  .AND.  m .NE. 0.0  .AND.  x .GT. 0.0 )  THEN
       log_normal_cdf = log_normal_cdf +                                    &
                        n * 0.5 * ( 1.0 + ERF( LOG( x / m ) / ( SQRT(2.0) * LOG(s) ) ) )
    ENDIF

 END FUNCTION log_normal_cdf

 REAL FUNCTION ERFI(x)
!
!-- Calculate the inverse error function by inverting a uniform approximation of
!-- the error function (Winitzki, International Conference on Computational
!-- Science and Its Applications (pp. 780-789). Springer.) The maximum error is
!-- less then 0.012 %.
    USE params, ONLY: pi

    IMPLICIT NONE

    REAL, INTENT(IN) :: x
    REAL, PARAMETER  :: a = 0.147

    erfi = 2.0 / ( pi * a ) + LOG( 1.0 - x**2 ) / 2.0
    erfi = SIGN( 1.0, x ) * SQRT( SQRT( erfi**2 - LOG( 1.0 - x**2 ) / a ) - erfi )

 END FUNCTION ERFI

 REAL FUNCTION IERF (arg)
!
!-- see also https://en.wikipedia.org/wiki/Error_function#Inverse_functions
    IMPLICIT NONE

    INTEGER,    PARAMETER :: QR_K = 16!selected_real_kind (32)
    REAL,       PARAMETER :: pi   = 3.141592653589793

    REAL,      INTENT(IN) :: arg
    REAL(KIND=QR_K)       :: arg_dp, IERF_dp, dummy_dp
    REAL(KIND=QR_K), SAVE :: coeff(0:5000)
    LOGICAL,         SAVE :: initialized = .FALSE.

    INTEGER :: mmm, nnn, ooo
!
!-- Precalculate coefficients for Maclaurin series for calculate ierf
    IF ( .NOT. initialized )  THEN
       coeff    = 0.0
       coeff(0) = 1.0
       DO mmm = 1, 5000
          DO nnn = 0, mmm - 1
              coeff(mmm) = coeff(mmm) + coeff(nnn) * coeff(mmm-1-nnn) / ( ( nnn + 1.0 ) * ( 2.0 * nnn + 1.0 ) )
          ENDDO
       ENDDO
       initialized = .TRUE.
    ENDIF

    IERF_dp = 0.0
    arg_dp  = arg * SQRT(pi) / 2.0
    DO mmm = 0, 5000
       dummy_dp = coeff(mmm) * arg_dp**( 2 * mmm + 1 ) / ( 2.0 * mmm + 1.0 )
       IERF_dp  = IERF_dp + dummy_dp
       IF ( dummy_dp / IERF_dp .LT. 1.0E-6 )  EXIT
    ENDDO

    IERF = IERF_dp

 END FUNCTION IERF

 REAL FUNCTION gauss( )
!
!-- Compute gauss random number with zero mean and unity standard deviation.
!-- Maximum absolute values are limited to 5.0.
    USE micro_random, ONLY: randomu
    USE grid,         ONLY: rank

    IMPLICIT NONE

    REAL :: x1, x2, y1, y2, w

    gauss = 10000.0

    DO WHILE ( ABS(gauss) .GT. 5.0 )
       w = 2.0
       DO WHILE ( w >= 1.0 )
          x1 = randomu(rank+1)
          x2 = randomu(rank+1)
          x1 = 2.0 * x1 - 1.0
          x2 = 2.0 * x2 - 1.0
          w  = x1 * x1 + x2 * x2
       END DO

       w = sqrt( (-2.0 * LOG( w ) ) / w )
       y1 = x1 * w
       y2 = x2 * w

       gauss = y1
    ENDDO

 END FUNCTION gauss
