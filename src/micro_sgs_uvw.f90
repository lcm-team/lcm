 SUBROUTINE micro_sgs_uvw(i,j,k)
!
!-- This subroutine calculates random velocities for each particle if switch_sgs_velocities = .TRUE.
    USE micro_vars,   ONLY: particle, tke_LCM, tk_LCM, iUsgs, iVsgs, iWsgs,    &
                            lcmp, n_start, n_end, dtn_LCM, iZ, iM,     &
                            nz_LCM, gpzs, sgfx_prof, iTYPE
    USE micro_random, ONLY: num_gauss
    USE microphysics, ONLY: micro_field, ipi, ipl

    IMPLICIT NONE

    INTEGER, INTENT(in) :: i, j, k

    INTEGER :: n
    REAL    :: RL, sigma, SGS_pre, SGS_lmt

    CALL t_startf ('micro_sgs_uvw')
!
!-- Check if flow is turbulent to avoid divide by zero in the parameterization
!-- formula. Turbulent if tkes > 1.e-2 (velocity scale > 8 cm/s).
    IF( tke_LCM(i,j,k) < 1.0E-4 )  THEN

       particle(n_start(i,j,k):n_end(i,j,k))%var(iUsgs) = 0.0_lcmp   ! SGS u wind
       particle(n_start(i,j,k):n_end(i,j,k))%var(iVsgs) = 0.0_lcmp   ! SGS v wind
       particle(n_start(i,j,k):n_end(i,j,k))%var(iWsgs) = 0.0_lcmp   ! SGS w wind

    ELSE
!
!--    Lagrangian autocorrelation coefficient:
       RL = EXP( -dtn_LCM / ( tk_LCM(i,j,k) + 1.0E-20 ) * tke_LCM(i,j,k)  )
!
!--    SGS velocity standard deviation:
       sigma = SQRT( 2.0 / 3.0 * tke_LCM(i,j,k) )
!
!--    SGS u, v, w winds:
       DO n = n_start(i,j,k), n_end(i,j,k)

          SGS_pre = RL * particle(n)%var(iUsgs) + SQRT( 1.0 - RL**2 ) * num_gauss( sigma )
          SGS_lmt = SIGN( MIN( 0.4 * particle(n)%var(iZ), ABS( SGS_pre * dtn_LCM ) ) / dtn_LCM, SGS_pre )
          particle(n)%var(iUsgs) = REAL( SGS_lmt, KIND=lcmp )

          SGS_pre = RL * particle(n)%var(iVsgs) + SQRT( 1.0 - RL**2 ) * num_gauss( sigma )
          SGS_lmt = SIGN( MIN( 0.4 * particle(n)%var(iZ), ABS( SGS_pre * dtn_LCM ) ) / dtn_LCM, SGS_pre )
          particle(n)%var(iVsgs) = REAL( SGS_lmt, KIND=lcmp )

          SGS_pre = RL * particle(n)%var(iWsgs) + SQRT( 1.0 - RL**2 ) * num_gauss( sigma )
          SGS_lmt = SIGN( MIN( 0.4 * particle(n)%var(iZ), ABS( SGS_pre * dtn_LCM ) ) / dtn_LCM, SGS_pre )
          particle(n)%var(iWsgs) = REAL( SGS_lmt, KIND=lcmp )
!
!--       Determine change of ql due to SGS mixing
          sgfx_prof(k) = sgfx_prof(k) - particle(n)%var(iM) * particle(n)%var(iWsgs)
          IF (particle(n)%ivar(iTYPE) .gt. 0) THEN
            micro_field(i,j,k,ipl) = micro_field(i,j,k,ipl) - particle(n)%var(iM) * particle(n)%var(iWsgs) !Negative because we're counting positive precipitation flux as going downwards
          ELSE
            micro_field(i,j,k,ipi) = micro_field(i,j,k,ipi) - particle(n)%var(iM) * particle(n)%var(iWsgs)
          ENDIF

       ENDDO

    ENDIF

    CALL t_stopf ('micro_sgs_uvw')

 END SUBROUTINE micro_sgs_uvw
