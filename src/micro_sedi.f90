 SUBROUTINE micro_terminal_velocity(i,j,k)
!
!-- Calculates terminal velocity of droplets
    USE micro_vars,   ONLY: particle, rho_liq, n_start, n_end, lcmp, sigma,    &
                            rho_air, pres_LCM, iA, iM, iWsedi, gpzs, dtn_LCM,  &
                            dz_LCM, iZ, iAXS, iCXS, iTYPE, sfx_prof, iDNS, dV, &
                            rho_ice, get_r_particle
    USE microphysics, ONLY: micro_field, ipi, ipl

    USE params,       ONLY: pi, ggr, fac_cond, muelq
    USE vars,         ONLY: tabs0, t, precsfc, prec_xy

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: i, j, k
    INTEGER             :: n

    REAL                :: Bo, diameter, eta, l, Cac, C1, C2, C3, NDa, NRe, NP, &
                           XX, YY, phi, DD_cgs, alpha, beta, gamma, sigma_ice,  &
                           Bv, ws_cgs, aRE_KC, bRE_KC, SQRT_XX, KC_FACTOR,      &
                           muelq_cgs, g_cgs, rho_cgs, aRE, bRE, Av, aaxis,      &
                           caxis, A, V, A_CE, N_Re, N_Be, Psi, N_Re_turb,       &
                           C_drag, D, Wsedi
!
!-- Parameters for calculating the terminal velocity (Beard, 1976, J. Atmos. Sci.).
    REAL, PARAMETER, DIMENSION(0:6) :: b = (/ -0.318657e1,   0.992696,    &
                                              -0.153193e-2, -0.987059e-3, &
                                              -0.578878e-3,  0.855176e-4, &
                                              -0.327815e-5 /)
    REAL, PARAMETER, DIMENSION(0:5) :: c = (/ -0.500015e1,   0.523778e1,  &
                                              -0.204914e1,   0.475294,    &
                                              -0.542819e-1,  0.238449e-2 /)

    REAL, PARAMETER :: eta0 = 1.818e-5, l0 = 6.62e-8, p0 = 1013.25, &
                       T0 = 293.15, rho0 = 1.292509, rho_0 = 1.23, c1_ice = 0.0902, &
                       delta = 9.06

    CALL t_startf ('micro_terminal_velocity')

    DO n = n_start(i,j,k), n_end(i,j,k)

       IF ( particle(n)%ivar(iTYPE) .GT. 0 )  THEN  ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
!
!--       Terminal velocity (Beard, 1976, J. Atmos. Sci.)
          diameter = get_r_particle(n) * 2.0

          eta = rho_air(k) * eta0 / rho0
          l   = l0 * ( eta / eta0 ) * ( p0 / pres_LCM(k) ) * SQRT( tabs0(k) / T0 )
          Cac = 1.0 + 2.5 * l / diameter

          IF ( diameter .LE. 19.0e-6 )  THEN
             C1 = ( rho_liq - rho_air(k) ) * ggr / ( 18.0 * eta  )

             particle(n)%var(iWsedi) = REAL( C1 * Cac * diameter**2, KIND=lcmp )
          ELSEIF ( diameter .LE. 1070.0e-6 )  THEN
             C2  = 4.0 * rho_air(k) * ( rho_liq - rho_air(k) ) * ggr / ( 3.0 * eta**2 )
             NDa = C2 * diameter**3
             XX  = LOG(NDa)
             YY  = b(0) + b(1) * XX + b(2) * XX**2 + b(3) * XX**3 + b(4) * XX**4 + &
                   b(5) * XX**5 + b(6) * XX**6
             NRe = Cac * EXP(YY)

             particle(n)%var(iWsedi) = REAL( eta * NRe / ( rho_air(k) * diameter ), KIND=lcmp )
          ELSE
             C3  = 4.0 * ( rho_liq - rho_air(k) ) * ggr / ( 3.0 * sigma )
             Bo  = C3 * diameter**2
             NP  = sigma**3 * rho_air(k)**2 / ( eta**4 * ( rho_liq - rho_air(k) ) * ggr )
             XX  = LOG(Bo * NP**0.166666666)
             YY  = c(0) + c(1) * XX + c(2) * XX**2 + c(3) * XX**3 + c(4) * XX**4 + &
                   c(5) * XX**5
             NRe = NP**0.166666666 * EXP(YY)

             particle(n)%var(iWsedi) = REAL( eta * NRe / ( rho_air(k) * diameter ), KIND=lcmp )
          ENDIF

       ELSEIF ( particle(n)%ivar(iTYPE) .LT. 0 )  THEN  ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
! !
! !--       Ice crystal terminal velocities by Mitchell (1996, JAS), updated with Khvorostyanov and Curry (2002, JAS).
!           aaxis = particle(n)%var(iAXS)
!           caxis = particle(n)%var(iCXS)
!
!           phi  = caxis / aaxis
! !
! !--       Convert variables to CGS units
!           IF ( phi .GT. 1.0 )  THEN
!              DD_cgs = caxis * 2.0 * 100.0
!           ELSE
!              DD_cgs = aaxis * 2.0 * 100.0
!           ENDIF
!
!           muelq_cgs = muelq * 1000.0 / 100.0 ! dynamic viscosity ( kg/(m*s) to g/(cm*s) )
!           g_cgs     = ggr * 100.0 ! acceleration by gravity ( m/s^2 to cm/s^2 )
!           rho_cgs   = rho_air(k) * 1000.0 / 100.0**3 ! density of air ( kg/m^3 to g/cm^3 )
! !
! !--       Tables from Mitchell (1996, JAS) and Khvorostyanov and Curry (2002, JAS)
!           IF ( phi .EQ. 1.0 )  THEN ! spheres
!              alpha = 0.524
!              beta  = 3.0
!              gamma = 0.785
!              sigma_ice = 2.0
!
!           ELSEIF ( phi .LT. 1.0 )  THEN ! plates
!              IF ( DD_cgs .LT. 100.0E-4 )  THEN
!                 alpha = 0.00739
!                 beta  = 2.45
!                 gamma = 0.24
!                 sigma_ice = 1.85
!              ELSEIF ( DD_cgs .GE. 100.0E-4 )  THEN
!                 alpha = 0.00739
!                 beta  = 2.45
!                 gamma = 0.65
!                 sigma_ice = 2.0
!              ENDIF
!
!           ELSEIF ( phi .GT. 1.0 )  THEN ! columns
!              IF ( DD_cgs .LT. 100.0E-4 )  THEN
!                 alpha = 0.1677
!                 beta  = 2.91
!                 gamma = 0.684
!                 sigma_ice = 2.0
!              ELSEIF ( DD_cgs .LT. 300.0E-4 )  THEN
!                 alpha = 0.00166
!                 beta  = 1.91
!                 gamma = 0.0696
!                 sigma_ice = 1.50
!              ELSEIF ( DD_cgs .GE. 300.0E-4 )  THEN
!                 alpha = 0.000907
!                 beta  = 1.74
!                 gamma = 0.0512
!                 sigma_ice = 1.414
!              ENDIF
!           ENDIF
! !
! !--       Eq. (17) from Mitchell (1996, JAS), and some definitions
!           XX = 2.0 * alpha * g_cgs * rho_cgs / ( gamma * muelq_cgs**2 ) * DD_cgs**(beta-sigma_ice+2.0)
!           SQRT_XX = SQRT(XX)
!           KC_FACTOR = SQRT( 1.0 + c1_ice * SQRT_XX )
! !
! !--    Using Eqs. (11), (17), (18), (19), (20), (21) from Mitchell (1996, JAS)
!           ! IF ( XX .LE. 10.0 )  THEN
!           !    aRE = 0.04394
!           !    bRE = 0.970
!           ! ELSEIF ( XX .LE. 585.0 )  THEN
!           !    aRE = 0.06049
!           !    bRE = 0.831
!           ! ELSEIF ( XX .LE. 1.56E5 )  THEN
!           !    aRE = 0.2072
!           !    bRE = 0.638
!           ! ELSE
!           !    aRE = 1.0865
!           !    bRE = 0.499
!           ! ENDIF
! !
! !--       Eq. (22) from Mitchell (1996, JAS)
!           ! w_cgs = aRE * ( muelq_cgs / rho_cgs )**(1.0-2.0*bRE) *           &
!           !         ( 2.0 * alpha * g_cgs / ( rho_cgs * gamma ) )**bRE *     &
!           !         DD_cgs**( bRE * ( beta - sigma_ice + 2.0 ) - 1.0 )
! !
! !--       Eqs. (2.12) and (2.13) from Khvorostyanov and Curry (2002, JAS)
!           bRE_KC = 1.0 / 2.0 * c1_ice * SQRT_XX / ( KC_FACTOR - 1.0 ) / KC_FACTOR
!           aRE_KC = delta**2 / 4.0 * ( KC_FACTOR - 1.0 )**2 / XX**bRE_KC
! !
! !--       Eqs. (2.24) and (2.25) from Khvorostyanov and Curry (2002, JAS)
!           Av = aRE_KC * ( muelq_cgs / rho_cgs )**(1.0-2.0*bRE_KC) * &
!                ( 2.0 * alpha * g_cgs / ( rho_cgs * gamma ) )**bRE_KC
!           Bv = bRE_KC * ( beta - sigma_ice + 2.0 ) - 1.0
! !
! !--       Eqs. (2.23) and (2.26) from Khvorostyanov and Curry (2002, JAS)
!           ws_cgs = Av * DD_cgs**Bv * ( rho_0 / rho_air(k) )**bRE_KC
!
!           particle(n)%var(iWsedi) = REAL( ws_cgs / 100.0, KIND=lcmp )

!
!--       Terminal Velocity of ice particles, based on Khvorostyanov and Curry (2005), using
!--       mostly the definitions of Bradr and Seifert (2018)
          D    = 2.0 * MAX( particle(n)%var(iAXS), particle(n)%var(iCXS) ) ! maximum diameter
          A_CE = pi * particle(n)%var(iAXS) * MAX( particle(n)%var(iAXS), particle(n)%var(iCXS) ) ! are of circumscribed ellipse
          A    = A_CE * ( particle(n)%var(iDNS) / rho_ice )**EXP(-particle(n)%var(iCXS)/particle(n)%var(iAXS)) ! projected area , (5) in Shima et al. (2020)
          V    = 4.0 / 3.0 * pi * particle(n)%var(iAXS)**2 * particle(n)%var(iCXS) ! volume

          N_Be = 2.0 * ABS( particle(n)%var(iDNS) - rho_air(k) ) * V * ggr * D**2 / ( A * rho_air(k) * muelq**2 ) ! Best number, (6) in Bradr and Seifert (2018)
          N_Re = 8.0**2 / 4.0 * ( SQRT( 1.0 + 4.0 / 8.0**2 * SQRT( N_Be / 0.35 ) ) - 1.0 )**2 ! Reynolds number, (7) in Bradr and Seifert (2018)

          Psi       = ( 1.0 + ( N_Be / 2.8E6 )**2 ) / ( 1.0 + 1.6 * ( N_Be / 2.8E6 )**2 ) ! correction factor to determine turbulent Reynolds number, (8) in Bradr and Seifert (2018)
          N_Re_turb = N_Re * SQRT(Psi) ! turbulent Reynolds number
          C_drag    = 0.35 * ( 1.0 + 8.0 / SQRT( N_Re_turb ) )**2 / Psi ! drag coefficient, (9) in Bradr and Seifert (2018)

          Wsedi = SQRT( 2.0 * ggr * V * ABS( particle(n)%var(iDNS) - rho_air(k) ) / ( A * C_drag * rho_air(k) ) ) ! sedimentation velocity, (10) in Bradr and Seifert (2018)

          particle(n)%var(iWsedi) = REAL( Wsedi, KIND=lcmp )

       ENDIF
!
!--    Determine impact of sedimentation on liquid water static energy
       sfx_prof(k) = sfx_prof(k) + particle(n)%var(iM) * particle(n)%var(iWsedi)
       IF (particle(n)%ivar(iTYPE) .gt. 0) THEN
         micro_field(i,j,k,ipl) = micro_field(i,j,k,ipl) + particle(n)%var(iM) * particle(n)%var(iWsedi)
       ELSE
         micro_field(i,j,k,ipi) = micro_field(i,j,k,ipi) + particle(n)%var(iM) * particle(n)%var(iWsedi)
       ENDIF
       IF ( k .EQ. 1 )  THEN
           precsfc(i,j) = precsfc(i,j) + particle(n)%var(iM) * particle(n)%var(iWsedi) / ( dz_LCM(k) * dV(k) * rho_air(k) )
           prec_xy(i,j) = prec_xy(i,j) + particle(n)%var(iM) * particle(n)%var(iWsedi) / ( dz_LCM(k) * dV(k) * rho_air(k) )
       ENDIF

    ENDDO

    CALL t_stopf ('micro_terminal_velocity')

 END SUBROUTINE micro_terminal_velocity

!  SUBROUTINE micro_sedimentation(i,j,k)
!
!     USE micro_vars,   ONLY: n_particles, n_start, n_end, particle, gpzc, dtn_LCM
!     USE microphysics, ONLY: iM, iWsedi, idql_sed, micro_field
!     USE grid,         ONLY: nzm
!
!     IMPLICIT NONE
!
!     INTEGER, INTENT(IN) :: i, j, k
!     INTEGER             :: k_new, n
!
!     DO n = n_start(i,j,k), n_end(i,j,k)
! !
! !--    Move particle in z-direction by terminal velocity
!        particle(n)%var(3) = particle(n)%var(3) - particle(n)%var(iWsedi) * dtn_LCM
! !
! !--    Get new particle location to determine the change of ql by sedimentation
!        k_new = MINLOC( ABS( particle(n)%var(3) - gpzc ), DIM=1 ) !for non-unifrom grids
!        micro_field(i,j,k,idql_sed) = micro_field(i,j,k,idql_sed) - particle(n)%var(iM)
!        IF ( k_new .GE. 1 )  THEN
!           micro_field(i,j,k_new,idql_sed) = micro_field(i,j,k_new,idql_sed) + particle(n)%var(iM)
!        ENDIF
!
!     ENDDO
!
!  END SUBROUTINE micro_sedimentation
