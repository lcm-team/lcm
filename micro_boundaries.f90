SUBROUTINE micro_boundaries
!
!-- Delete paticles if they are above ztop or below zbottom
    USE micro_vars,          ONLY: iM, iA, iWsgs, rho_liq, particle,           &
                                   lcmp, particle, particle_array, iWsgs,      &
                                   n_particles, n_particles_alloc, zbottom,    &
                                   ztop, np_deleted, dx_LCM, dy_LCM, iN,       &
                                   rho_aero, iZ, gpzc, iU, iV, iW,             &
                                   iUsgs, iVsgs, iWsgs, iT, iq, ieta, iklst,   &
                                   iType, iMrime, iNMono,                      &
                                   idx, idy, nx_sb, ny_sb, nx_LCM, ny_LCM,     &
                                   iWsedi, iX, iY, nz_LCM, sigma_aero, rm_aero,   &
                                   tabs_LCM, pres_LCM, n_start, n_end,         &
                                   switch_reflective_walls_pichamber,          &
                                   switch_walls_pichamber, switch_steady_aerosol, &
                                   switch_steady_pichamber, dz_LCM,            &
                                   wall_x_e_pichamber, wall_x_w_pichamber,     &
                                   wall_y_n_pichamber, wall_y_s_pichamber,     &
                                   inject_bottom, inject_top,     &
                                   get_r_aero, get_r_particle, r_equi, gpzs,   &
                                   search_sorted_array
    USE microphysics,        ONLY: imv, micro_field
    USE grid,                ONLY: RUN3D, rank, z
    USE vars,                ONLY: total_water_prec, gamaz, qsatw
    USE params,              ONLY: pi
    USE micro_random,        ONLY: randomu

    IMPLICIT NONE

    INTEGER :: n, m, i, j, k, n_seeded
    REAL    :: radius, r_aero, ierf_dummy, x_rand, IERF, SS_gb
   !  REAL, DIMENSION(:,:), ALLOCATABLE :: zdiffarr

    TYPE(particle_array), ALLOCATABLE, DIMENSION(:) :: particle_temp

    CALL t_startf ('micro_boundaries')

    np_deleted = 0
    IF ( ANY( particle(1:n_particles)%var(iZ) .GT. ztop )     .OR.  &
         ANY( particle(1:n_particles)%var(iZ) .LT. zbottom )  .OR.  &
         switch_walls_pichamber )  THEN

       total_water_prec = total_water_prec + SUM( particle(1:n_particles)%var(iM) ) / ( dx_LCM * dy_LCM )

       n_particles_alloc = MAX( n_particles_alloc, n_particles )
       ALLOCATE( particle_temp(0:n_particles_alloc) )

       m = 1
       n_seeded = 1
       DO n = 1, n_particles

          radius = get_r_particle(n)

          IF ( particle(n)%var(iZ) .GE. zbottom  .AND. particle(n)%var(iZ) .LE. ztop )  THEN
!
!--          Nothing is done for particles within the model domain
             IF ( .NOT. switch_walls_pichamber )  THEN

                particle_temp(m)%ivar(:) = particle(n)%ivar(:)
                particle_temp(m)%var(:)  = particle(n)%var(:)
                m = m + 1
!
!--          Get interaction with walls of the pi chamber
             ELSEIF ( switch_walls_pichamber )  THEN

                IF ( particle(n)%var(iX) .GE. wall_x_w_pichamber  .AND.  &
                     particle(n)%var(iX) .LE. wall_x_e_pichamber  .AND.  &
                     particle(n)%var(iY) .GE. wall_y_s_pichamber  .AND.  &
                     particle(n)%var(iY) .LE. wall_y_n_pichamber )  THEN

                   particle_temp(m)%ivar(:) = particle(n)%ivar(:)
                   particle_temp(m)%var(:)  = particle(n)%var(:)
                   m = m + 1

                ELSE

                   IF ( .NOT. switch_reflective_walls_pichamber )  THEN ! particles get absorbed at wall

                      CYCLE

                   ELSE

                      particle_temp(m)%ivar(:) = particle(n)%ivar(:)
                      particle_temp(m)%var(:)  = particle(n)%var(:)

                      IF ( particle_temp(m)%var(iX) .LT. wall_x_w_pichamber )  THEN
                         particle_temp(m)%var(iX)    = REAL( wall_x_w_pichamber + ( wall_x_w_pichamber - particle_temp(m)%var(iX) ), KIND=lcmp )
                         particle_temp(m)%var(iUSGS) = -1.0_lcmp * particle_temp(m)%var(iUSGS)
                      ENDIF
                      IF ( particle_temp(m)%var(iX) .GT. wall_x_e_pichamber )  THEN
                         particle_temp(m)%var(iX)    = REAL( wall_x_e_pichamber + ( wall_x_e_pichamber - particle_temp(m)%var(iX) ), KIND=lcmp )
                         particle_temp(m)%var(iUSGS) = -1.0_lcmp * particle_temp(m)%var(iUSGS)
                      ENDIF
                      IF ( particle_temp(m)%var(iY) .LT. wall_y_s_pichamber )  THEN
                         particle_temp(m)%var(iY)    = REAL( wall_y_s_pichamber + ( wall_y_s_pichamber - particle_temp(m)%var(iY) ), KIND=lcmp )
                         particle_temp(m)%var(iVSGS) = -1.0_lcmp * particle_temp(m)%var(iVSGS)
                      ENDIF
                      IF ( particle_temp(m)%var(iY) .GT. wall_y_n_pichamber )  THEN
                         particle_temp(m)%var(iY)    = REAL( wall_y_n_pichamber + ( wall_y_n_pichamber - particle_temp(m)%var(iY) ), KIND=lcmp )
                         particle_temp(m)%var(iVSGS) = -1.0_lcmp * particle_temp(m)%var(iVSGS)
                      ENDIF

                      m = m + 1

                   ENDIF

                ENDIF

             ENDIF

          ELSEIF ( particle(n)%var(iZ) .LT. zbottom  .AND.  radius .GT. 1.0E-6 )  THEN
!
!--          Particles > 1.0 µm in radius are absorbed at bottom (2nd option, default) or
!--          are re-injected into the model domain at 50 % height, but without liquid water
!--          (1st option)
             IF ( .NOT. (switch_steady_pichamber .OR. switch_steady_aerosol))  THEN

                CYCLE

             ELSEIF ( switch_steady_pichamber .OR. switch_steady_aerosol )  THEN

                particle_temp(m)%ivar(:)    = particle(n)%ivar(:)
                particle_temp(m)%var(:)     = particle(n)%var(:)
                x_rand = randomu(rank+1)
                particle_temp(m)%var(iZ)    = REAL( inject_bottom + ( inject_top - inject_bottom ) * x_rand, KIND=lcmp )
                particle_temp(m)%var(iWsgs) = 0.0
!
!--             Get i, j, k indices of particle
                i = MAX( MIN( FLOOR( particle_temp(m)%var(iX) * idx ) + 1 - nx_sb, nx_LCM ), 1 )
                IF (RUN3D) THEN
                   j = MAX( MIN( FLOOR( particle_temp(m)%var(iY) * idy ) + 1 - ny_sb, ny_LCM ), 1 )
                ELSE
                   j = 1
                ENDIF
               !  allocate(zdiffarr(n_particles, nz_LCM))
               !  do k = 1, nz_LCM
               !    zdiffarr(m, k) = abs(particle_temp(m)%var(iZ)-z(k))
               !  end do
                k = search_sorted_array(gpzs, particle(n)%var(iZ), nz_LCM+1)
               !  k = minloc(zdiffarr(m))

                r_aero = get_r_aero(n)
                SS_gb  = ( micro_field(i,j,k,imv) - qsatw( tabs_LCM(i,j,k), pres_LCM(k) ) ) / qsatw( tabs_LCM(i,j,k), pres_LCM(k) )
                SS_gb  = MAX( -0.2, MIN( SS_gb, -0.01 ) )

                particle_temp(m)%var(iM) = REAL( particle_temp(m)%ivar(iA) * 4.0 / 3.0 * pi * rho_liq * MAX( r_equi(SS_gb,tabs_LCM(i,j,k),0.0,n), r_aero )**3, KIND=lcmp )
                particle_temp(m)%ivar(iTYPE) = ABS(particle(m)%ivar(iTYPE))
                particle_temp(m)%ivar(iNMono) = 1
                particle_temp(m)%var(iMrime) = 0.0

                m = m + 1

             ENDIF

          ELSE !reflection of aerosols
!
!--          Remaining particles are reflected at bottom or top
             IF ( particle(n)%var(iZ) .LT. zbottom )  THEN

                particle_temp(m)%ivar(:)    = particle(n)%ivar(:)
                particle_temp(m)%var(:)     = particle(n)%var(:)

                particle_temp(m)%var(iZ)    = REAL( zbottom + ( zbottom - particle_temp(m)%var(iZ) ), KIND=lcmp )
                particle_temp(m)%var(iWsgs) = ABS( particle_temp(m)%var(iWsgs) ) ! make sure that particle moves away from wall

                m = m + 1

             ENDIF

             IF ( particle(n)%var(iZ) .GT. ztop )  THEN

                particle_temp(m)%ivar(:)    = particle(n)%ivar(:)
                particle_temp(m)%var(:)     = particle(n)%var(:)

                particle_temp(m)%var(iZ)    = REAL( ztop + ( ztop - particle_temp(m)%var(iZ) ), KIND=lcmp )
                particle_temp(m)%var(iWsgs) = -1.0_lcmp * ABS( particle_temp(m)%var(iWsgs) ) ! make sure that particle moves away from wall

                m = m + 1

             ENDIF

          ENDIF

       ENDDO

       np_deleted  = n_particles - ( m - 1 )
       n_particles = m - 1

       DO n = 1, n_particles
          particle(n)%ivar(:) = particle_temp(n)%ivar(:)
          particle(n)%var(:)  = particle_temp(n)%var(:)
       ENDDO

       DEALLOCATE(particle_temp)

       total_water_prec = total_water_prec - SUM( particle(1:n_particles)%var(iM) ) / ( dx_LCM * dy_LCM )

    ENDIF

    CALL t_stopf ('micro_boundaries')

 END SUBROUTINE micro_boundaries
