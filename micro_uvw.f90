 SUBROUTINE micro_uvw_prep(i,j,k)
!
!-- Prepare resolved scale velocity calculation
    USE micro_vars, ONLY: particle, lcmp, n_start, n_end, u_LCM, v_LCM, w_LCM,    &
                          idx, idy, gpxs, gpys, gpzs, nz_LCM, iU, iV, iW, iX, iY, iZ

    USE grid,       ONLY: RUN3D

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: i, j, k
    INTEGER             :: n
    REAL                :: unext, vnext, wnext, alpha, beta, gamma

    CALL t_startf ('micro_uvw_prep')

    DO n = n_start(i,j,k), n_end(i,j,k)
       alpha = MAX( MIN( ( particle(n)%var(iX) - gpxs(i) ) * idx, 1.0 ), 0.0 )
       unext = u_LCM(i,j,k) * ( 1.0 - alpha ) + u_LCM(i+1,j,k) * alpha
       IF ( RUN3D ) THEN
          beta  = MAX( MIN( ( particle(n)%var(iY) - gpys(j) ) * idy, 1.0 ), 0.0 )
          vnext = v_LCM(i,j,k) * ( 1.0 - beta ) + v_LCM(i,j+1,k) * beta
       ELSE
          vnext = v_LCM(i,j,k)
       ENDIF
       gamma = MAX( MIN( ( particle(n)%var(iZ) - gpzs(k) ) / ( gpzs(k+1) - gpzs(k) ), 1.0 ), 0.0 )

       IF ( k .LE. nz_LCM - 1 )  THEN
          wnext = w_LCM(i,j,k) * ( 1.0 - gamma ) + w_LCM(i,j,k+1) * gamma
       ELSE
          wnext = w_LCM(i,j,k) * ( 1.0 - gamma ) ! assuming that w = 0 for k >= nz_LCM
       ENDIF

       particle(n)%var(iU) = REAL( unext, KIND=lcmp )
       particle(n)%var(iV) = REAL( vnext, KIND=lcmp )
       particle(n)%var(iW) = REAL( wnext, KIND=lcmp )
    ENDDO

    CALL t_stopf ('micro_uvw_prep')

 END SUBROUTINE micro_uvw_prep

 SUBROUTINE micro_uvw(i,j,k)
!
!-- This subroutine finds resolved velocity at t=t+1/2 with iterative method.
    USE micro_vars,   ONLY: particle, lcmp, n_start, n_end, dtn_LCM,           &
                            n_uvw_iteration, u_LCM, v_LCM, w_LCM,              &
                            gpzs, nz_LCM, iM, nz_LCM, nz_LCM,                  &
                            gpzs, iX, iY, iZ, iU, iV, iW, afx_prof, iTYPE
    USE microphysics, ONLY: micro_field, ipl, ipi
    USE micro_random, ONLY: num_gauss

    USE grid,         ONLY: RUN3D
    USE domain,       ONLY: YES3D

    IMPLICIT NONE

    REAL    :: alpha, beta, gamma, unow, unext, vnow, vnext, wnow, wnext, xp,  &
               yp, zp
    INTEGER :: i0, j0, k0, n, nn

    INTEGER, INTENT(IN) :: i, j, k

    CALL t_startf ('micro_uvw')
!
!-- Find u, v, w at t=t+1/2 by iterative method.
!-- If n_uvw_iteration = 0, particle position will be updated with a forward
!-- scheme, i.e., x(t+1) = x(t) + u(t) * dt using u(t) determined in micro_uvw_prep
    IF ( n_uvw_iteration .GE. 1 )  THEN

       DO n = n_start(i,j,k), n_end(i,j,k)

          unow = particle(n)%var(iU)
          vnow = particle(n)%var(iV)
          wnow = particle(n)%var(iW)

          DO nn = 1, n_uvw_iteration
!
!--          Guess new position
             xp = particle(n)%var(iX) + particle(n)%var(iU) * dtn_LCM
             yp = particle(n)%var(iY) + particle(n)%var(iV) * dtn_LCM * YES3D
             zp = particle(n)%var(iZ) + particle(n)%var(iW) * dtn_LCM
!
!--          Use "simple" interpolation as proposed by Grabowski et al. (2017, GMDD)
             CALL get_ijk(xp,yp,zp,i0,j0,k0,alpha,beta,gamma)

             unext = u_LCM(i0,j0,k0) * ( 1.0 - alpha ) + u_LCM(i0+1,j0,k0) * alpha
             IF ( RUN3D ) THEN
                vnext = v_LCM(i0,j0,k0) * ( 1.0 - beta )  + v_LCM(i0,j0+1,k0) * beta
             ELSE
                vnext = v_LCM(i0,j0,k0)
             ENDIF

             IF ( k0 .LE. nz_LCM - 1 )  THEN
                wnext = w_LCM(i0,j0,k0) * ( 1.0 - gamma ) + w_LCM(i0,j0,k0+1) * gamma
             ELSE
                wnext = w_LCM(i0,j0,k0) * ( 1.0 - gamma ) ! assuming that w = 0 for k >= nz_LCM
             ENDIF

             particle(n)%var(iU) = REAL( 0.5 * ( unow + unext ), KIND=lcmp )
             particle(n)%var(iV) = REAL( 0.5 * ( vnow + vnext ), KIND=lcmp )
             particle(n)%var(iW) = REAL( 0.5 * ( wnow + wnext ), KIND=lcmp )

          ENDDO

       ENDDO

    ENDIF
!
!-- Determine change of ql due to advection
    DO n = n_start(i,j,k), n_end(i,j,k)
       afx_prof(k) = afx_prof(k) - particle(n)%var(iM) * particle(n)%var(iW)
       IF (particle(n)%ivar(iTYPE) .gt. 0) THEN
         micro_field(i,j,k,ipl) = micro_field(i,j,k,ipl) - particle(n)%var(iM) * particle(n)%var(iW) !Negative because we're counting positive precipitation flux as going downwards
       ELSE
         micro_field(i,j,k,ipi) = micro_field(i,j,k,ipi) - particle(n)%var(iM) * particle(n)%var(iW)
       ENDIF
    END DO

    CALL t_stopf ('micro_uvw')

 END SUBROUTINE micro_uvw
