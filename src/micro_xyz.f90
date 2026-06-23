 SUBROUTINE micro_xyz
!
!-- Update particle position
    USE micro_vars, ONLY: switch_sedimentation, start_sedimentation, lcmp,          &
                          particle, iUsgs, iVsgs, iWsgs, iWsedi, n_particles,     &
                          dtn_LCM, switch_sgs_velocities, iX, iY, iZ, iU, iV, iW

    USE domain,     ONLY: YES3D
    USE grid,       ONLY: nstep

    IMPLICIT NONE

    INTEGER :: n

    CALL t_startf ('micro_xyz')

    DO n = 1, n_particles
       particle(n)%var(iX) = REAL( particle(n)%var(iX) + particle(n)%var(iU) * dtn_LCM, KIND=lcmp )
       particle(n)%var(iY) = REAL( particle(n)%var(iY) + particle(n)%var(iV) * dtn_LCM * YES3D, KIND=lcmp )
       particle(n)%var(iZ) = REAL( particle(n)%var(iZ) + particle(n)%var(iW) * dtn_LCM, KIND=lcmp )
       IF ( switch_sgs_velocities )  THEN
          particle(n)%var(iX) = REAL( particle(n)%var(iX) + particle(n)%var(iUsgs) * dtn_LCM, KIND=lcmp )
          particle(n)%var(iY) = REAL( particle(n)%var(iY) + particle(n)%var(iVsgs) * dtn_LCM * YES3D, KIND=lcmp )
          particle(n)%var(iZ) = REAL( particle(n)%var(iZ) + particle(n)%var(iWsgs) * dtn_LCM, KIND=lcmp )
       ENDIF
       IF ( switch_sedimentation  .AND.  nstep .GE. start_sedimentation )  THEN
          particle(n)%var(iZ) = REAL( particle(n)%var(iZ) - particle(n)%var(iWsedi) * dtn_LCM, KIND=lcmp )
       ENDIF

    ENDDO

    CALL t_stopf ('micro_xyz')

 END SUBROUTINE micro_xyz
