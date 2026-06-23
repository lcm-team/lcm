 SUBROUTINE micro_collision(i,j,k)
!
!-- Calculates collision and coalescence of droplets
    USE micro_vars,   ONLY: n_start, n_end, particle, rho_liq, dV, diss_LCM,   &
                            lcmp, dtn_LCM, aut_prof, acc_prof,                 &
                            tke_LCM, switch_linear_sampling_collection,        &
                            r_separate_micro, IA, iM, iN, switch_breakup,      &
                            switch_icebreak, iWsedi, iTYPE, tabs_LCM,          &
                            icebreak_param, ds_param
    USE micro_random, ONLY: randomu
    USE microphysics, ONLY: irmf, irim
    USE micro_task_util_MPI, ONLY: micro_task_abort
   !  USE microphysics, ONLY: iA, iM, iN, iWsedi, iTYPE

    USE params,       ONLY: pi
    USE grid,         ONLY: rank

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  i, j, k
    INTEGER             ::  n, m

    REAL                ::  p_crit
    LOGICAL             ::  check_collection, check_shattering
    LOGICAL             ::  odd_n_particles
    INTEGER             ::  half_n_particles
    INTEGER             ::  tmp_n_ptcl, n_particles, mm
    INTEGER             ::  coalescence_count, breakup_count, agg_count, icebreak_count, &
                            nfrag_tot, nfrag_b, nfrag_s
    REAL                ::  smallest_radius

    INTEGER, DIMENSION(:), ALLOCATABLE :: ptcl_pair1, ptcl_pair2

    CALL t_startf ('micro_collision')
!
!-- To save computational resources, the largest droplet should be larger than
!-- 5.0 µm to cause any collisions.
    IF ( MAXVAL( particle(n_start(i,j,k):n_end(i,j,k))%var(iM) /               &
                 particle(n_start(i,j,k):n_end(i,j,k))%ivar(iA) )              &
         .LT. 5.0E-6**3 * 4.0 / 3.0 * pi * rho_liq )  THEN

       CALL t_stopf ('micro_collision')
       RETURN

    ENDIF

    n_particles = n_end(i,j,k) - n_start(i,j,k) + 1
!
!-- Linear Sampling method (LSM, Shima et al. 2009)
    IF ( switch_linear_sampling_collection )  THEN

       half_n_particles = FLOOR( n_particles / 2.0 )
!
!--    Special treatment when the number of particles in grid box is odd.
!--    Ignore the collision of smallest particle.
       IF ( MOD(n_particles,2) .NE. 0 ) THEN
          odd_n_particles = .TRUE.
          smallest_radius = MINVAL(particle(n_start(i,j,k):n_end(i,j,k))%var(iM)/ &
                                   particle(n_start(i,j,k):n_end(i,j,k))%ivar(iA))
          tmp_n_ptcl      = n_particles - 1
       ELSE
          odd_n_particles = .FALSE.
          tmp_n_ptcl      = n_particles
          smallest_radius = -999.0
       ENDIF

       ALLOCATE (ptcl_pair1(1:half_n_particles))
       ALLOCATE (ptcl_pair2(1:half_n_particles))
!
!--    Call permuation code to get 2 n_particle/2 particle index combination
       CALL permutation(odd_n_particles,n_start(i,j,k),n_end(i,j,k),tmp_n_ptcl,&
                        half_n_particles,ptcl_pair1,ptcl_pair2,smallest_radius)

       DO mm = 1, half_n_particles

          m = ptcl_pair2(mm)
          n = ptcl_pair1(mm)
!
!--       Coalescence/breakup among liquid particles
          IF ( particle(n)%ivar(iTYPE) .GT. 0  .AND.  particle(m)%ivar(iTYPE) .GT. 0 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
             CALL determine_collision(m,n,i,j,k,n_particles,coalescence_count,breakup_count)
             IF ( coalescence_count .GT. 0) THEN
                CALL update_coalescence(n,m,i,j,k,coalescence_count)
             ENDIF
             IF (breakup_count .GT. 0 .AND. switch_breakup) THEN
                CALL update_breakup(n,m,i,j,k)
             ENDIF
! !
!--       Aggregation/Breakup between ice particles
          ELSEIF ( particle(n)%ivar(iTYPE) .LT. 0  .AND.  particle(m)%ivar(iTYPE) .LT. 0 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
             CALL determine_icecoll(m,n,i,j,k,n_particles,agg_count,icebreak_count)
             IF (agg_count .gt. 0) THEN
               CALL update_aggregation(n,m,i,j,k)
             ENDIF

             IF (icebreak_count .gt. 0 .and. switch_icebreak) THEN
               IF (TRIM(icebreak_param) .EQ. "Phillips") THEN
                  CALL update_icebreak_phillips(n,m,i,j,k)
               ELSEIF (TRIM(icebreak_param) .EQ. "Sullivan") THEN
                  CALL update_icebreak_sullivan(n,m,i,j,k)
               ELSE
                  print*, "Ice Breakup Parameterization must be either 'Phillips' or 'Sullivan!'"
                  call micro_task_abort
               ENDIF
             ENDIF
!
!--       Riming among a frozen and a liquid particle
          ELSE

             CALL determine_riming(m,n,i,j,k,n_particles,check_collection,check_shattering,p_crit,nfrag_tot,nfrag_b,nfrag_s)

             IF ( check_collection ) THEN
                CALL update_riming(n,m,i,j,k,p_crit)
             ELSEIF ( check_shattering ) THEN
                CALL update_shattering(m,n,i,j,k,nfrag_tot,nfrag_b,nfrag_s)
             ENDIF

          ENDIF

       ENDDO

       DEALLOCATE(ptcl_pair1)
       DEALLOCATE(ptcl_pair2)

    ELSE
!
!--    Calculate the collectional growth following the "all-or-nothing" approach
       DO n = n_start(i,j,k), n_end(i,j,k)

          DO m = n+1, n_end(i,j,k)
!
!--          Coalescence/breakup among liquid particles
             IF ( particle(n)%ivar(iTYPE) .GT. 0  .AND.  particle(m)%ivar(iTYPE) .GT. 0 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN

                CALL determine_collision(m,n,i,j,k,n_particles,coalescence_count,breakup_count)

                IF ( coalescence_count .gt. 0 ) THEN
                   CALL update_coalescence(n,m,i,j,k,coalescence_count)
                ENDIF

                IF (breakup_count .gt. 0 .AND. switch_breakup) THEN
                   CALL update_breakup(n,m,i,j,k)
                ENDIF
!
!--          Aggregation/Breakup among ice particles
             ELSEIF ( particle(n)%ivar(iTYPE) .LT. 0  .AND.  particle(m)%ivar(iTYPE) .LT. 0 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
             CALL determine_icecoll(m,n,i,j,k,n_particles,agg_count,icebreak_count)
             IF (agg_count .gt. 0) THEN
               CALL update_aggregation(n,m,i,j,k)
             ENDIF

             IF (icebreak_count .gt. 0 .and. switch_icebreak) THEN
               IF (TRIM(icebreak_param) .EQ. "Phillips") THEN
                  CALL update_icebreak_phillips(n,m,i,j,k)
               ELSEIF (TRIM(icebreak_param) .EQ. "Sullivan") THEN
                  CALL update_icebreak_sullivan(n,m,i,j,k)
               ELSE
                  print*, "Ice Breakup Parameterization must be either 'Phillips' or 'Sullivan!'"
                  call micro_task_abort
               ENDIF
             ENDIF
!
!--          Riming among a frozen and a liquid particle
             ELSE

                CALL determine_riming(m,n,i,j,k,n_particles,check_collection,check_shattering,p_crit,nfrag_tot,nfrag_b,nfrag_s)

               IF ( check_collection ) THEN
                  CALL update_riming(n,m,i,j,k,p_crit)
               ELSEIF ( check_shattering ) THEN
                  CALL update_shattering(m,n,i,j,k,nfrag_tot,nfrag_b,nfrag_s)
               ENDIF

             ENDIF

          ENDDO

       ENDDO

    ENDIF

    CALL t_stopf ('micro_collision')

 END SUBROUTINE micro_collision

 REAL FUNCTION E_H80(r1,r2)
!
!-- Function that interpolates collision efficiencies from Hall (1980, JAS)
    IMPLICIT NONE

    INTEGER          ::  iq, ir, kr
    REAL             ::  pp, qq, rq
    REAL, INTENT(IN) ::  r1, r2

    REAL, DIMENSION(1:21),      SAVE ::  rat
    REAL, DIMENSION(1:15),      SAVE ::  r0
    REAL, DIMENSION(1:15,1:21), SAVE ::  ecoll

    r0  = (/   6.0,   8.0,  10.0, 15.0,  20.0,  25.0,  &
              30.0,  40.0,  50.0, 60.0,  70.0, 100.0,  &
             150.0, 200.0, 300.0 /)

    rat = (/ 0.00, 0.05, 0.10, 0.15, 0.20, 0.25,       &
             0.30, 0.35, 0.40, 0.45, 0.50, 0.55,       &
             0.60, 0.65, 0.70, 0.75, 0.80, 0.85,       &
             0.90, 0.95, 1.00 /)

    ecoll(:,1)  = (/ 0.001, 0.001, 0.001, 0.001, 0.001, &
                     0.001, 0.001, 0.001, 0.001, 0.001, &
                     0.001, 0.001, 0.001, 0.001, 0.001 /)
    ecoll(:,2)  = (/ 0.003, 0.003, 0.003, 0.004, 0.005, &
                     0.005, 0.005, 0.010, 0.100, 0.050, &
                     0.200, 0.500, 0.770, 0.870, 0.970 /)
    ecoll(:,3)  = (/ 0.007, 0.007, 0.007, 0.008, 0.009, &
                     0.010, 0.010, 0.070, 0.400, 0.430, &
                     0.580, 0.790, 0.930, 0.960, 1.000 /)
    ecoll(:,4)  = (/ 0.009, 0.009, 0.009, 0.012, 0.015, &
                     0.010, 0.020, 0.280, 0.600, 0.640, &
                     0.750, 0.910, 0.970, 0.980, 1.000 /)
    ecoll(:,5)  = (/ 0.014, 0.014, 0.014, 0.015, 0.016, &
                     0.030, 0.060, 0.500, 0.700, 0.770, &
                     0.840, 0.950, 0.970, 1.000, 1.000 /)
    ecoll(:,6)  = (/ 0.017, 0.017, 0.017, 0.020, 0.022, &
                     0.060, 0.100, 0.620, 0.780, 0.840, &
                     0.880, 0.950, 1.000, 1.000, 1.000 /)
    ecoll(:,7)  = (/ 0.030, 0.030, 0.024, 0.022, 0.032, &
                     0.062, 0.200, 0.680, 0.830, 0.870, &
                     0.900, 0.950, 1.000, 1.000, 1.000 /)
    ecoll(:,8)  = (/ 0.025, 0.025, 0.025, 0.036, 0.043, &
                     0.130, 0.270, 0.740, 0.860, 0.890, &
                     0.920, 1.000, 1.000, 1.000, 1.000 /)
    ecoll(:,9)  = (/ 0.027, 0.027, 0.027, 0.040, 0.052, &
                     0.200, 0.400, 0.780, 0.880, 0.900, &
                     0.940, 1.000, 1.000, 1.000, 1.000 /)
    ecoll(:,10) = (/ 0.030, 0.030, 0.030, 0.047, 0.064, &
                     0.250, 0.500, 0.800, 0.900, 0.910, &
                     0.950, 1.000, 1.000, 1.000, 1.000 /)
    ecoll(:,11) = (/ 0.040, 0.040, 0.033, 0.037, 0.068, &
                     0.240, 0.550, 0.800, 0.900, 0.910, &
                     0.950, 1.000, 1.000, 1.000, 1.000 /)
    ecoll(:,12) = (/ 0.035, 0.035, 0.035, 0.055, 0.079, &
                     0.290, 0.580, 0.800, 0.900, 0.910, &
                     0.950, 1.000, 1.000, 1.000, 1.000 /)
    ecoll(:,13) = (/ 0.037, 0.037, 0.037, 0.062, 0.082, &
                     0.290, 0.590, 0.780, 0.900, 0.910, &
                     0.950, 1.000, 1.000, 1.000, 1.000 /)
    ecoll(:,14) = (/ 0.037, 0.037, 0.037, 0.060, 0.080, &
                     0.290, 0.580, 0.770, 0.890, 0.910, &
                     0.950, 1.000, 1.000, 1.000, 1.000 /)
    ecoll(:,15) = (/ 0.037, 0.037, 0.037, 0.041, 0.075, &
                     0.250, 0.540, 0.760, 0.880, 0.920, &
                     0.950, 1.000, 1.000, 1.000, 1.000 /)
    ecoll(:,16) = (/ 0.037, 0.037, 0.037, 0.052, 0.067, &
                     0.250, 0.510, 0.770, 0.880, 0.930, &
                     0.970, 1.000, 1.000, 1.000, 1.000 /)
    ecoll(:,17) = (/ 0.037, 0.037, 0.037, 0.047, 0.057, &
                     0.250, 0.490, 0.770, 0.890, 0.950, &
                     1.000, 1.000, 1.000, 1.000, 1.000 /)
    ecoll(:,18) = (/ 0.036, 0.036, 0.036, 0.042, 0.048, &
                     0.230, 0.470, 0.780, 0.920, 1.000, &
                     1.020, 1.020, 1.020, 1.020, 1.020 /)
    ecoll(:,19) = (/ 0.040, 0.040, 0.035, 0.033, 0.040, &
                     0.112, 0.450, 0.790, 1.010, 1.030, &
                     1.040, 1.040, 1.040, 1.040, 1.040 /)
    ecoll(:,20) = (/ 0.033, 0.033, 0.033, 0.033, 0.033, &
                     0.119, 0.470, 0.950, 1.300, 1.700, &
                     2.300, 2.300, 2.300, 2.300, 2.300 /)
    ecoll(:,21) = (/ 0.027, 0.027, 0.027, 0.027, 0.027, &
                     0.125, 0.520, 1.400, 2.300, 3.000, &
                     4.000, 4.000, 4.000, 4.000, 4.000 /)

!
!-- Calculate the radius class index of particles with respect to array r0.
!-- Radius has to be in microns
    DO  kr = 1, 15
       IF ( MAX( r1, r2 ) * 1.0E6 < r0(kr) )  THEN
          ir = kr
          EXIT
       ENDIF
    ENDDO
    IF ( MAX( r1, r2 ) * 1.0E6 >= r0(15) )  ir = 16
!
!-- Two-dimensional linear interpolation of the collision efficiency.
!-- Radius has to be in microns
    rq = MIN( r1 / r2, r2 / r1 )
    iq = INT( rq * 20 ) + 1
    iq = MAX( iq , 2)

    IF ( ir < 16 )  THEN
       IF ( ir >= 2 )  THEN
          pp = ( ( MAX( r1, r2 ) * 1.0E6 ) - r0(ir-1) ) / &
               ( r0(ir) - r0(ir-1) )
          qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )

          E_H80 = ( 1.0 - pp ) * ( 1.0 - qq )                &
                  * ecoll(ir-1,iq-1)                         &
                  + pp * ( 1.0 - qq ) * ecoll(ir,iq-1)       &
                  + qq * ( 1.0 - pp ) * ecoll(ir-1,iq)       &
                  + pp * qq * ecoll(ir,iq)
       ELSE
          qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )

          E_H80 = ( 1.0 - qq ) * ecoll(1,iq-1) + qq * ecoll(1,iq)
       ENDIF
    ELSE
       qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )

       E_H80 = MIN( ( 1.0 - qq ) * ecoll(15,iq-1) + qq * ecoll(15,iq), 1.0 )
    ENDIF

    IF ( E_H80 < 1.0E-20 )  E_H80 = 0.0

 END FUNCTION E_H80

 REAL FUNCTION gck(r1,r2,ws1,ws2,epsilon,tke)

    USE params, ONLY : pi, muelq, ggr

    IMPLICIT NONE

    REAL ::  ao, ao_gr, bbb, be, b1, b2, ccc, c1, c1_gr, c2, d1, d2, eta, e1,  &
             e2, fao_gr, fr, grfin, lambda, lambda_re, lf, rc, rrp, sst, tauk, &
             tl, t2, tt, t1, urms, vk, vrms1xy, vrms2xy, v1, v1v2xy, v1xysq,   &
             v2, v2xysq, wrfin, wrgrav2, wrtur2xy, xx, yy, z, st1, tau1, st2,  &
             tau2, vis_kin

    REAL, INTENT(IN) ::  r1, r2, ws1, ws2, epsilon, tke
    REAL             ::  phi_w, zhi !other functions to be called

    REAL, PARAMETER  :: rho_dummy = 1.2

!    urms      = 2.02 * ( epsilon / 0.04 )**( 1.0 / 3.0 ) !from Riechelmann et al. (2012), assumes a lengthscale L = 374 m, tke ~ ( L * epsilon )**(2/3) 
    urms      = SQRT( 2.0 / 3.0 * tke ) ! relates root mean square velocity to sgs TKE.

    vis_kin   = muelq / rho_dummy  ! get kinematic viscosity

    lambda    = urms * SQRT( 15.0 * vis_kin / epsilon )     ! Taylor microscale
    lambda_re = urms**2 * SQRT( 15.0 / epsilon / vis_kin )  ! Taylor microscale Reynolds number
    tl        = urms**2 / epsilon
    lf        = 0.5 * urms**3 / epsilon
    tauk      = SQRT( vis_kin / epsilon )
    eta       = ( vis_kin**3 / epsilon )**0.25
    vk        = eta / tauk

    ao = ( 11.0 + 7.0 * lambda_re ) / ( 205.0 + lambda_re )
    tt = SQRT( 2.0 * lambda_re / ( SQRT( 15.0 ) * ao ) ) * tauk

    tau1 = ws1 / ggr   ! inertial time scale
    st1  = tau1 / tauk  ! Stokes number
    tau2 = ws2 / ggr   ! inertial time scale
    st2  = tau2 / tauk  ! Stokes number

!
!-- Calculate average radial relative velocity at contact (wrfin)
    z   = tt / tl
    be  = SQRT( 2.0 ) * lambda / lf
    bbb = SQRT( 1.0 - 2.0 * be**2 )
    d1  = ( 1.0 + bbb ) / ( 2.0 * bbb )
    e1  = lf * ( 1.0 + bbb ) * 0.5
    d2  = ( 1.0 - bbb ) * 0.5 / bbb
    e2  = lf * ( 1.0 - bbb ) * 0.5
    ccc = SQRT( 1.0 - 2.0 * z**2 )
    b1  = ( 1.0 + ccc ) * 0.5 / ccc
    c1  = tl * ( 1.0 + ccc ) * 0.5
    b2  = ( 1.0 - ccc ) * 0.5 / ccc
    c2  = tl * ( 1.0 - ccc ) * 0.5

    v1  = ws1
    t1  = tau1
    v2  = ws2
    t2  = tau2
    rrp = r1 + r2

    v1xysq  = b1 * d1 * phi_w(c1,e1,v1,t1) - b1 * d2 * phi_w(c1,e2,v1,t1) &
              - b2 * d1 * phi_w(c2,e1,v1,t1) + b2 * d2 * phi_w(c2,e2,v1,t1)
    v1xysq  = v1xysq * urms**2 / t1
              vrms1xy = SQRT( v1xysq )

    v2xysq  = b1 * d1 * phi_w(c1,e1,v2,t2) - b1 * d2 * phi_w(c1,e2,v2,t2) &
              - b2 * d1 * phi_w(c2,e1,v2,t2) + b2 * d2 * phi_w(c2,e2,v2,t2)
    v2xysq  = v2xysq * urms**2 / t2
    vrms2xy = SQRT( v2xysq )

    IF ( ws1 >= ws2 )  THEN
       v1 = ws1
       t1 = tau1
       v2 = ws2
       t2 = tau2
    ELSE
       v1 = ws2
       t1 = tau2
       v2 = ws1
       t2 = tau1
    ENDIF

    v1v2xy   =  b1 * d1 * zhi(c1,e1,v1,t1,v2,t2) - &
                b1 * d2 * zhi(c1,e2,v1,t1,v2,t2) - &
                b2 * d1 * zhi(c2,e1,v1,t1,v2,t2) + &
                b2 * d2* zhi(c2,e2,v1,t1,v2,t2)
    fr       = d1 * EXP( -rrp / e1 ) - d2 * EXP( -rrp / e2 )
    v1v2xy   = v1v2xy * fr * urms**2 / tau1 / tau2
    wrtur2xy = vrms1xy**2 + vrms2xy**2 - 2.0 * v1v2xy
    IF ( wrtur2xy < 0.0 )  wrtur2xy = 0.0
    wrgrav2  = pi / 8.0 * ( ws2 - ws1 )**2
    wrfin    = SQRT( ( 2.0 / pi ) * ( wrtur2xy + wrgrav2) )

!
!-- Calculate radial distribution function (grfin)
    IF ( st1 > st1 )  THEN
       sst = st1
    ELSE
       sst = st2
    ENDIF

    xx = -0.1988 * sst**4 + 1.5275 * sst**3 - 4.2942 *       &
         sst**2 + 5.3406 * sst
    IF ( xx < 0.0 )  xx = 0.0
    yy = 0.1886 * EXP( 20.306 / lambda_re )

    c1_gr  =  xx / ( ggr / vk * tauk )**yy

    ao_gr  = ao + ( pi / 8.0 ) * ( ggr / vk * tauk )**2
    fao_gr = 20.115 * SQRT( ao_gr / lambda_re )
    rc     = SQRT( fao_gr * ABS( st2 - st1 ) ) * eta

    grfin  = ( ( eta**2 + rc**2 ) / ( rrp**2 + rc**2) )**( c1_gr*0.5 )
    IF ( grfin < 1.0 )  grfin = 1.0

!
!-- Calculate general collection kernel (without the consideration of
!-- collection efficiencies)
    gck = 2.0 * pi * rrp**2 * wrfin * grfin

 END FUNCTION gck

 REAL FUNCTION phi_w( a, b, vsett, tau0 )
!
!-- Function used in the Ayala et al. (2008, NJP) analytical model for turbulent
!-- effects on the collision kernel
    IMPLICIT NONE

    REAL, INTENT(IN) ::  a, b, tau0, vsett
    REAL             ::  aa1

    aa1 = 1.0 / tau0 + 1.0 / a + vsett / b
    phi_w = 1.0 / aa1  - 0.5 * vsett / b / aa1**2

 END FUNCTION phi_w

 REAL FUNCTION zhi( a, b, vsett1, tau1, vsett2, tau2 )
!
!-- Function used in the Ayala et al. (2008, NJP) analytical model for turbulent
!-- effects on the collision kernel
    IMPLICIT NONE

    REAL, INTENT(IN) ::  a, b, tau1, tau2, vsett1, vsett2
    REAL             ::  aa1, aa2, aa3, aa4, aa5, aa6

    aa1 = vsett2 / b - 1.0 / tau2 - 1.0 / a
    aa2 = vsett1 / b + 1.0 / tau1 + 1.0 / a
    aa3 = ( vsett1 - vsett2 ) / b + 1.0 / tau1 + 1.0 / tau2
    aa4 = ( vsett2 / b )**2 - ( 1.0 / tau2 + 1.0 / a )**2
    aa5 = vsett2 / b + 1.0 / tau2 + 1.0 / a
    aa6 = 1.0 / tau1 - 1.0 / a + ( 1.0 / tau2 + 1.0 / a) * vsett1 / vsett2

    zhi = (1.0 / aa1 - 1.0 / aa2 ) * ( vsett1 - vsett2 ) * 0.5 /    &
          b / aa3**2 + ( 4.0 / aa4 - 1.0 / aa5**2 - 1.0 / aa1**2 )  &
          * vsett2 * 0.5 / b /aa6 + ( 2.0 * ( b / aa2 - b / aa1 ) - &
          vsett1 / aa2**2 + vsett2 / aa1**2 ) * 0.5 / b / aa3

 END FUNCTION zhi

 REAL FUNCTION E_turb(r1,r2,epsilon)
!
!-- Function to get turbulent collision efficiencies (Wang and Grabowski, 2009, QJRMS)
    IMPLICIT NONE

    INTEGER :: iq, ir, kr, kq

    INTEGER :: ira1, ira2

    REAL ::  pp, qq, rq, y1, y2, y3, r1_micro, r2_micro

    REAL, DIMENSION(1:11) ::  rat           !<
    REAL, DIMENSION(1:7)  ::  r0            !<

    REAL, DIMENSION(1:7,1:11) ::  ecoll_100, ecoll_400 !<

    REAL, INTENT(IN) ::  r1, r2, epsilon

!
!-- Initial assignment of constants
    r0  = (/  10.0, 20.0, 30.0, 40.0, 50.0, 60.0,  &
             100.0 /)

    rat = (/ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, &
             0.7, 0.8, 0.9, 1.0 /)
!
!-- Tabulated turbulent enhancement factor at 100 cm**2/s**3
    ecoll_100(:,1)  = (/  1.74,   1.74,   1.773, 1.49,  &
                          1.207,  1.207,  1.0 /)
    ecoll_100(:,2)  = (/  1.46,   1.46,   1.421, 1.245, &
                          1.069,  1.069,  1.0 /)
    ecoll_100(:,3)  = (/  1.32,   1.32,   1.245, 1.123, &
                          1.000,  1.000,  1.0 /)
    ecoll_100(:,4)  = (/  1.250,  1.250,  1.148, 1.087, &
                          1.025,  1.025,  1.0 /)
    ecoll_100(:,5)  = (/  1.186,  1.186,  1.066, 1.060, &
                          1.056,  1.056,  1.0 /)
    ecoll_100(:,6)  = (/  1.045,  1.045,  1.000, 1.014, &
                          1.028,  1.028,  1.0 /)
    ecoll_100(:,7)  = (/  1.070,  1.070,  1.030, 1.038, &
                          1.046,  1.046,  1.0 /)
    ecoll_100(:,8)  = (/  1.000,  1.000,  1.054, 1.042, &
                          1.029,  1.029,  1.0 /)
    ecoll_100(:,9)  = (/  1.223,  1.223,  1.117, 1.069, &
                          1.021,  1.021,  1.0 /)
    ecoll_100(:,10) = (/  1.570,  1.570,  1.244, 1.166, &
                          1.088,  1.088,  1.0 /)
    ecoll_100(:,11) = (/ 20.3,   20.3,   14.6,   8.61,  &
                          2.60,   2.60,   1.0 /)
!
!-- Tabulated turbulent enhancement factor at 400 cm**2/s**3
    ecoll_400(:,1)  = (/  4.976,  4.976,  3.593,  2.519, &
                          1.445,  1.445,  1.0 /)
    ecoll_400(:,2)  = (/  2.984,  2.984,  2.181,  1.691, &
                          1.201,  1.201,  1.0 /)
    ecoll_400(:,3)  = (/  1.988,  1.988,  1.475,  1.313, &
                          1.150,  1.150,  1.0 /)
    ecoll_400(:,4)  = (/  1.490,  1.490,  1.187,  1.156, &
                          1.126,  1.126,  1.0 /)
    ecoll_400(:,5)  = (/  1.249,  1.249,  1.088,  1.090, &
                          1.092,  1.092,  1.0 /)
    ecoll_400(:,6)  = (/  1.139,  1.139,  1.130,  1.091, &
                          1.051,  1.051,  1.0 /)
    ecoll_400(:,7)  = (/  1.220,  1.220,  1.190,  1.138, &
                          1.086,  1.086,  1.0 /)
    ecoll_400(:,8)  = (/  1.325,  1.325,  1.267,  1.165, &
                          1.063,  1.063,  1.0 /)
    ecoll_400(:,9)  = (/  1.716,  1.716,  1.345,  1.223, &
                          1.100,  1.100,  1.0 /)
    ecoll_400(:,10) = (/  3.788,  3.788,  1.501,  1.311, &
                          1.120,  1.120,  1.0 /)
    ecoll_400(:,11) = (/ 36.52,  36.52,  19.16,  22.80,  &
                          26.0,   26.0,    1.0 /)
!
!-- Calculate the radius class index of particles with respect to array r0
!-- The droplet radius has to be given in microns.
    r1_micro = r1 * 1.0E6
    DO  kr = 1, 7
       IF ( r1_micro < r0(kr) )  THEN
          ira1 = kr
          EXIT
       ENDIF
    ENDDO
    IF ( r1_micro >= r0(7) )  ira1 = 8
    r2_micro = r2 * 1.0E6
    DO  kr = 1, 7
       IF ( r2_micro < r0(kr) )  THEN
          ira2 = kr
          EXIT
       ENDIF
    ENDDO
    IF ( r2_micro >= r0(7) )  ira2 = 8

!
!-- Two-dimensional linear interpolation of the turbulent enhancement factor.
!-- The droplet radius has to be given in microns.
    ir = MAX( ira1, ira2 )
    rq = MIN( r1 / r2, r2 / r1 )

    DO  kq = 2, 11
       IF ( rq <= rat(kq) )  THEN
          iq = kq
          EXIT
       ENDIF
    ENDDO

    y1 = 1.0  ! turbulent enhancement factor at 0 m**2/s**3

    IF ( ir < 8 )  THEN
       IF ( ir >= 2 )  THEN
          pp = ( MAX( r1, r2 ) * 1.0E6 -  &
               r0(ir-1) ) / ( r0(ir) - r0(ir-1) )
          qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
          y2 = ( 1.0 - pp ) * ( 1.0 - qq ) * ecoll_100(ir-1,iq-1) + &
               pp * ( 1.0 - qq ) * ecoll_100(ir,iq-1)        + &
               qq * ( 1.0 - pp ) * ecoll_100(ir-1,iq)        + &
               pp * qq           * ecoll_100(ir,iq)
          y3 = ( 1.0-pp ) * ( 1.0 - qq ) * ecoll_400(ir-1,iq-1)      + &
               pp * ( 1.0 - qq ) * ecoll_400(ir,iq-1)        + &
               qq * ( 1.0 - pp ) * ecoll_400(ir-1,iq)        + &
               pp * qq           * ecoll_400(ir,iq)
       ELSE
          qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
          y2 = ( 1.0 - qq ) * ecoll_100(1,iq-1) + qq * ecoll_100(1,iq)
          y3 = ( 1.0 - qq ) * ecoll_400(1,iq-1) + qq * ecoll_400(1,iq)
       ENDIF
    ELSE
       qq = ( rq - rat(iq-1) ) / ( rat(iq) - rat(iq-1) )
       y2 = ( 1.0 - qq ) * ecoll_100(7,iq-1) + qq * ecoll_100(7,iq)
       y3 = ( 1.0 - qq ) * ecoll_400(7,iq-1) + qq * ecoll_400(7,iq)
    ENDIF
!
!-- Linear interpolation of turbulent enhancement factor
    IF ( epsilon <= 0.01 )  THEN
       E_turb = ( epsilon - 0.01 ) / ( 0.0  - 0.01 ) * y1 &
                + ( epsilon - 0.0  ) / ( 0.01 - 0.0  ) * y2
    ELSEIF ( epsilon <= 0.06 )  THEN
       E_turb = ( epsilon - 0.04 ) / ( 0.01 - 0.04 ) * y2 &
                + ( epsilon - 0.01 ) / ( 0.04 - 0.01 ) * y3
    ELSE
       E_turb = ( 0.06 - 0.04 ) / ( 0.01 - 0.04 ) * y2 &
                + ( 0.06 - 0.01 ) / ( 0.04 - 0.01 ) * y3
    ENDIF

    IF ( E_turb < 1.0 )  E_turb = 1.0

 END FUNCTION E_turb

 SUBROUTINE permutation( odd_switch, n_start, n_end, n_ptcl, &
                         half_n_ptcl, ptcl_pair1, ptcl_pair2, smallest_radius )

    USE microphysics, ONLY: micro_field
    USE micro_vars,   ONLY: particle, iA, iM, iN, iWsedi
    USE grid,         ONLY: rank
    USE micro_random, ONLY: randomu

    IMPLICIT NONE

    REAL    :: mm, smallest_radius
    INTEGER :: half_n_ptcl, temp, k, m, n_ptcl, i
    INTEGER :: n_start,n_end
    LOGICAL :: odd_switch

    INTEGER, DIMENSION(n_ptcl)      :: tmp_ptcl_arr
    INTEGER, DIMENSION(half_n_ptcl) :: ptcl_pair1, ptcl_pair2
!
!-- Store particle indices to array. If there is an odd number of particles per
!-- grid box, replace the smallest particle index with the last particle index,
    k = 1
    DO i = n_start, n_ptcl+n_start-1
       tmp_ptcl_arr(k) = i
       IF ( odd_switch )  THEN
          IF( particle(i)%var(iM) / particle(i)%ivar(iA) .EQ. smallest_radius ) THEN
             tmp_ptcl_arr(k) = n_end
          ENDIF
       ENDIF
       k = k + 1
    ENDDO
!
!-- Using a random generator, the order of perticles are permutated.
    DO i = 2, n_ptcl
       mm = randomu(rank+1)
       m  = 1 + FLOOR( (n_ptcl-1) * mm )
       IF ( m .NE. i )  THEN
          temp            = tmp_ptcl_arr(m)
          tmp_ptcl_arr(m) = tmp_ptcl_arr(i)
          tmp_ptcl_arr(i) = temp
       ENDIF
    ENDDO
!
!-- The particles are divided into two arrays of the same length.
    ptcl_pair1 = tmp_ptcl_arr(:half_n_ptcl)
    ptcl_pair2 = tmp_ptcl_arr(half_n_ptcl+1:)

 END SUBROUTINE permutation

 SUBROUTINE determine_collision(m,n,i,j,k,n_particles,coalescence_count,breakup_count)

    USE micro_vars,   ONLY: n_start, n_end, particle, rho_liq, dV, diss_LCM,   &
                            lcmp, dtn_LCM, r_separate_micro, tke_LCM, iX, iY, iZ,    &
                            switch_linear_sampling_collection, iWsedi, iM, iA, sigma, &
                            switch_breakup, get_r_particle
    USE micro_random, ONLY: randomu
    USE microphysics, ONLY: micro_field
    USE params,       ONLY: pi
    USE grid,         ONLY: rank, time, nrestart, caseid
    USE vars,         ONLY: esatw, qsatw

    IMPLICIT NONE

    REAL                 ::  xn, xm, xsn, xsm, KK_coll, KK_coal, wsn, wsm, rn, rm, E_H80, &
                             x_separate, xsn_seeded, xsm_seeded, radius, p_crit, p_crit_lim, &
                             CKE, St, Sc, DelSSig, W1, W2, ET, Ecoal, x_rand, & !Parameters for determining likelihood of coalescence vs breakup
                             rmin, rmax, Ds, Dl, Vs, Vl, ecoal_s, &
                             ecoal_ll, sf_dmin, sf_dmax
    REAL, EXTERNAL       ::  ecoal_seifert
    INTEGER              ::  half_n_particles, pcount
    INTEGER, INTENT(IN)  ::  m, n, i, j, k, n_particles
    INTEGER, INTENT(OUT) :: coalescence_count, breakup_count
   !  REAL,    INTENT(OUT) ::  p_crit !No longer need to output this, since we've split it up into coalescence_count and breakup_count
    !B. Ascher 03.06.2025
!
!-- Collision is only calculated if one weighting factor is larger than 1 to avoid
!-- zero weighting factors
    coalescence_count = 0
    breakup_count = 0
    IF ( MAX( particle(n)%ivar(iA), particle(m)%ivar(iA) ) .GT. 1_8 ) THEN

!-- Compute mass of superdroplets m and n
      xn = REAL(particle(n)%var(iM)/particle(n)%ivar(iA))
      xm = REAL(particle(m)%var(iM)/particle(m)%ivar(iA))
!--    Probability that the necessary collisions take place
       wsn = particle(n)%var(iWsedi)
       wsm = particle(m)%var(iWsedi)
!
!--    Radii of droplets n and m
       rn = get_r_particle(n)
       rm = get_r_particle(m)
       Ds = 2.0*MIN(rn, rm)
       Dl = 2.0*MAX(rn, rm)
       IF (xn .LE. xm) THEN
         Vs = wsn
         Vl = wsm
       ELSE
         Vs = wsm
         Vl = wsn
       ENDIF
!
!--    Get collision kernel using coalescence efficiencies from Hall (1980)
       KK_coll = pi * ( rn + rm )**2 * E_H80(rn,rm) * ABS( wsn - wsm )
!
!--    Get collection kernel using coalescence efficiencies from Hall (1980) and
!--    turbulence enhancement by Ayala et al. (2008) and Grabowski and Wang (2009)
       ! IF ( epsilon .LT. 1.0E-10 )  THEN
       !    KK_coal = pi * ( rn + rm )**2 * E_H80(rn,rm) * ABS( wsn - wsm )
       ! ELSE
       !    KK_coal = E_H80(rn,rm) * gck(rn,rm,wsn,wsm,diss_LCM(i,j,k),tke_LCM(i,j,k)) * E_turb(rn,rm,diss_LCM(i,j,k))
       ! ENDIF
!
!--    Get collection kernel based on Long (1974)
       ! IF ( MAX(rn,rm) .LE. 50.0E-6 )  THEN
       !    KK_coal = 9.44E9 * ( xn**2 + xm**2 ) / rho_liq**2 * 1.0E6 !1.0E6 necessary for conversion to m3/s
       ! ELSE
       !    KK_coal = 5.78E3 * ( xn + xm ) / rho_liq
       ! ENDIF
!
!--    Probability of collision (not yet coalescence!)
       p_crit = MAX( particle(n)%ivar(iA), particle(m)%ivar(iA) ) / dV(k) * KK_coll * dtn_LCM
!
!--    Collision probability needs to be upscaled if the LSM method is used.
       IF ( switch_linear_sampling_collection )  THEN
          half_n_particles = FLOOR( n_particles / 2.0 )
          p_crit           = p_crit * n_particles * ( n_particles - 1.0 ) / ( half_n_particles * 2.0 )
       ENDIF

       IF ( p_crit .GT. randomu(rank+1) )  THEN

         !Now, need to figure out the probability that colliding droplets coalesce or breakup. Much of this breakup parameterization
         !is based on McFarquhar (2004), itself based on the studies of Low and List (1982) (referred to as LL82). However,
         !LL82 only applies for larger droplets. Therefore, we'll use Seifert et al.(2005)'s parameteriation for droplets under 300um,
         !then smoothly vary into the LL82 coalescence efficiency between 300um and 600um, above which only the LL82 values
         !will be used

         IF (Ds .lt. 0.1E-3) THEN
            Ecoal = 1.0
         ELSEIF (Ds .lt. 0.6E-3) THEN

            ecoal_s = ecoal_seifert(Ds, Dl, Vs, Vl)
            CALL coll_energy(Ds, Dl, Vs, Vl, CKE, St, Sc, W1, W2, DelSsig)

            ! CKE = (2. * rho_liq * pi / 3.) * (rn**3*rm**3 / (rn**3 + rm**3))*(wsn - wsm)**2 !Collision kinetic energy
            ! !From LL82 Part I Eq. 3.3. Their equation uses 1/12 rather than 2 because they use diameter rather than radius
            ! St = pi * sigma * ((2.0*rn)**2 + (2.0*rm)**2) !Total surface energy of the two separate droplets
            ! !Eq. 4.2 in LL82 Part I
            ! Sc = pi * sigma * ((2.0*rn)**3 + (2.0*rm)**3)**(2./3.)
            ! !Eq. 4.3 in LL82 Part I
            ! DelSsig = St-Sc
            ! !Eq. 4.1 in LL82 Part I
            Et = CKE+DelSSig
            !Eq 4.4 in LL82 Part I
            IF (Et .lt. 5E-6) THEN
               ecoal_ll = 0.778*(1+min(rn, rm)/max(rn, rm))**(-2)*EXP(-2.61E6*sigma*Et**2/Sc)
            ELSE
               ecoal_ll = 0.0
            ENDIF

            Ecoal = 0.5*((1+ COS((MAX(Ds, 0.3E-3)-0.3E-3)*pi/0.3E-3))*ecoal_s + (1-COS((MAX(Ds,0.3E-3)-0.3E-3)*pi/(0.3E-3)))*ecoal_ll)
         ELSE
            CALL coll_energy(Ds, Dl, Vs, Vl, CKE, St, Sc, W1, W2, DelSsig)

            Et = CKE+DelSSig
            !Eq 4.4 in LL82 Part I
            IF (Et .lt. 5E-6) THEN
               Ecoal = 0.778*(1+min(rn, rm)/max(rn, rm))**(-2)*EXP(-2.61E6*sigma*Et**2/Sc)
            ELSE
               Ecoal = 0.0
            ENDIF
         ENDIF
         IF ( p_crit .LE. 1.0 )  THEN
!
!--          If a single collision takes place, p_crit must be 1
            p_crit = 1.0
         ELSE
!
!--          Multiple collisions
            p_crit = MAX( NINT( p_crit ), 1 )
!
!--          Calculate a limiting p_crit that ensures that no more particles are collected than
!--          available ( MAX(A,B) - p * MIN(A,B) >= 1 <=> p <= ( MAX(A,B) - 1 ) / MIN(A,B) )
            p_crit_lim = ( MAX( particle(n)%ivar(iA), particle(m)%ivar(iA) ) - 1.0 ) &
                           / MIN( particle(n)%ivar(iA), particle(m)%ivar(iA) )
            p_crit_lim = MAX( FLOOR( p_crit_lim ), 1 )
!
!--          p_crit must be the minimum of the limiter and the original value
            p_crit = MIN( p_crit_lim, p_crit )
         ENDIF

         DO pcount = 1, NINT(p_crit)
            x_rand = randomu(rank+1)
            IF (.not. (switch_breakup) .or. x_rand .lt. Ecoal) THEN
               coalescence_count = coalescence_count + 1
            ELSE
               breakup_count = breakup_count + 1
            ENDIF
         END DO

      ENDIF

    ENDIF
    
 END SUBROUTINE determine_collision

 SUBROUTINE update_breakup(n,m,i,j,k)
   USE micro_vars, ONLY: rho_liq, particle, sigma, lcmp, &
                         iM, iA, iN, iWsedi, breakup_conserve_mode, &
                         n_start, n_end, in_startind, iTFZ, iTYPE,  &
#ifdef SAMCHEM
                         iAero_portions, naqspc, iaqspc, brk_prof
#else
                         iAero_portions, brk_prof
#endif /*SAMCHEM*/
   USE micro_random, ONLY: randomu
   USE params, ONLY: pi
   USE grid, ONLY: rank
   USE microphysics, ONLY: micro_field, ibreak
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, m, i, j, k
   INTEGER             :: incounter
   INTEGER(KIND=8)     :: s_ind, l_ind, collnum, p, MinN, Mults, Multf
   REAL :: Dl, Ds, CKE, W1, W2, St, Sc, DelSsig, &
           Rf, Rs, Rd, Ff, Fs, Fd, Vs, Vl,       &
           Ms0, Ml0, Ms1, Ml1, Mf, mn, mm,  &
           cs, cl, cf, cm, Df, Dp, Minddiff,           &
           newnum, newmass, newccn
#ifdef SAMCHEM
   REAL, DIMENSION(naqspc) :: aqs, aql, aqf, aqm, newspc
#endif /*SAMCHEM*/
   REAL, DIMENSION(4) :: in_m, in_s, in_l, in_f, newin !IN mass arrays
   REAL, DIMENSION(4) :: ps, pl, pf, pm          !CCN proportion arrays
   REAL, EXTERNAL :: t_frz


   !Figure out the smaller weighting factor
   collnum = MIN(particle(n)%ivar(iA), particle(m)%ivar(iA))
   
   !Figure out which index corresponds to the larger particle and which to the smaller particle
   !For now, only allow a single breakup event per timestep. Might change this later, but for now I worry that particle
   !properties change too much by the first breakup to make a second one valid. - B.Ascher 11.06.2025
   mn = particle(n)%var(iM)/REAL(particle(n)%ivar(iA))
   mm = particle(m)%var(iM)/REAL(particle(m)%ivar(iA))

   IF (mn .gt. mm) THEN
      s_ind = m
      l_ind = n
      Ms0 = mm
      Ml0 = mn
   ELSE
      s_ind = n
      l_ind = m
      Ms0 = mn
      Ml0 = mm
   ENDIF
   Ds = ((6./(pi*rho_liq)*Ms0)**(1./3.)) !small particle diameter
   Dl = ((6./(pi*rho_liq)*Ml0)**(1./3.)) !large particle diameter
   Vs = particle(s_ind)%var(iWsedi)                                       !small particle fall speed
   Vl = particle(l_ind)%var(iWsedi)                                       !large particle fall speed
   cs = particle(s_ind)%var(iN)/particle(s_ind)%ivar(iA)                  !small particle aerosol mass
   cl = particle(l_ind)%var(iN)/particle(l_ind)%ivar(iA)                  !large particle aerosol mass
   in_s = particle(s_ind)%var(in_startind+1:in_startind+4)/particle(s_ind)%ivar(iA) !small particle IN surface area
   in_l = particle(l_ind)%var(in_startind+1:in_startind+4)/particle(l_ind)%ivar(iA) !large particle IN surface area
   ps = particle(s_ind)%var(iAero_portions) !small particle aerosol mass fractions
   pl = particle(l_ind)%var(iAero_portions) !large particle aerosol mass fractions
#ifdef SAMCHEM
   aqs = particle(s_ind)%var(iaqspc)/particle(s_ind)%ivar(iA) !small particle aqueous molecules
   aql = particle(l_ind)%var(iaqspc)/particle(l_ind)%ivar(iA) !large particle aqueous molecules
#endif /*SAMCHEM*/

   CALL coll_energy(Ds, Dl, Vs, Vl, CKE, St, Sc, W1, W2, DelSsig)
   CALL breakup_fractions(Ds, Dl, CKE, St, Sc, W1, W2, DelSsig, Ff, Fs, Fd, Rf, Rs, Rd)
   CALL breakup_mcfarquhar(Ds, Dl, Vs, Vl, CKE, St, Sc, W1, W2, DelSsig, Ff, Fs, Fd, Rf, Rs, Rd, Mults, Multf, Ms1, Ml1, Mf)
   pm = (ps*cs + pl*cl)/SUM(ps*cs + pl*cl) !Aerosol portions of the merged droplet - proportional to the aerosol mass in each colliding droplet
   cm = cs + cl
   in_m = in_s + in_l
   cs = Ms1/(Ms0 + Ml0)*cm
   cl = Ml1/(Ms0 + Ml0)*cm
   cf = Mf/(Ms0 + Ml0)*cm
   in_s = Ms1/(Ms0 + Ml0)*in_m
   in_l = Ml1/(Ms0 + Ml0)*in_m
   in_f = Mf/(Ms0 + Ml0)*in_m
   Df = (6/(pi*rho_liq)*Mf)**(1./3.)
#ifdef SAMCHEM
   aqm = aqs + aql
   aqs = Ms1/(Ms0 + Ml0)*aqm
   aql = Ml1/(Ms0 + Ml0)*aqm
   aqf = Mf/(Ms0 + Ml0)*aqm
#endif /*SAMCHEM*/

   !breakup_conserve_flag determines whether we proceed along an approach preserving droplet sizes as much as possible,
   !or whether we follow the De Jong et al. (2023) approach of averaging fragment sizes, leaving the SD with the larger weighting
   !factor alone, and assigning all droplets within the SD with the smaller weighting factor to the average fragment mass
   IF (breakup_conserve_mode .EQ. 1) THEN !First option; try to preserve large and small droplet masses as much as possible
      MinN = -1 !Index of particle to merge fragments with
      particle(l_ind)%var(iM) = REAL((Ml0*(particle(l_ind)%ivar(iA)-collnum))+ Ml1*collnum, KIND=lcmp)
      particle(l_ind)%var(iAero_portions) = REAL(pl * (particle(l_ind)%var(iN)/particle(l_ind)%ivar(iA))*(particle(l_ind)%ivar(iA) - collnum) + pm * cl * collnum, KIND=lcmp)
      particle(l_ind)%var(iAero_portions) = particle(l_ind)%var(iAero_portions)/SUM(particle(l_ind)%var(iAero_portions))
#ifdef SAMCHEM
      particle(l_ind)%var(iaqspc) = REAL((particle(l_ind)%var(iaqspc)/particle(l_ind)%ivar(iA)) * (particle(l_ind)%ivar(iA) - collnum) + aql * collnum, KIND=lcmp)
#endif /*SAMCHEM*/
      particle(l_ind)%var(iN) = REAL((particle(l_ind)%var(iN)/particle(l_ind)%ivar(iA))*(particle(l_ind)%ivar(iA) - collnum) &
                              + cl*collnum, KIND=lcmp)
      particle(l_ind)%var(in_startind+1:in_startind+4) = REAL(&
      (particle(l_ind)%var(in_startind+1:in_startind+4)/particle(l_ind)%ivar(iA))*(particle(l_ind)%ivar(iA) - collnum) &
      + in_l*collnum, KIND=lcmp)
      
      !Here, set the small particle's mass to be the pre-collision mass times the number of remaining uncollided droplets plus
      !the post-collision mass times the number of collided small droplets

      IF (Mults .gt. 0) THEN !Filament breakup, where the small droplet survives
         Minddiff = ABS(Df - (6/(pi*rho_liq)*Ms1)**(1./3.)) !Difference between small droplet size and fragment size
         IF (Multf .gt. 0) THEN
            DO p = n_start(i,j,k), n_end(i,j,k)
               Dp = (6/(pi*rho_liq)*particle(p)%var(iM)/particle(p)%ivar(iA))**(1./3.)
               IF (ABS(Dp - Df) .lt. Minddiff .AND. (ABS(Dp-Df) .LT. 0.1E-3) .AND. (particle(p)%ivar(iTYPE) .gt. 0)) THEN !Condition: merging particle must have a diameter within
                                                                                 !0.1 mm of the fragment and be the best-matching diameter in the gridbox
                  Minddiff = ABS(Dp-Df)
                  MinN = p
               ENDIF
            END DO
         ENDIF

         IF (MinN .gt. 0) THEN !Successfully assigned the fragments to another particle within the gridbox
            particle(MinN)%ivar(iA) = INT(particle(MinN)%ivar(iA) + Multf*collnum, KIND=8)
            particle(MinN)%var(iM) = REAL(particle(MinN)%var(iM) + Mf*Multf*collnum, KIND=lcmp)
            particle(MinN)%var(iAero_portions) = REAL(particle(MinN)%var(iAero_portions)*particle(MinN)%var(iN) + pm*cf*Multf*collnum, KIND=lcmp)
            particle(MinN)%var(iAero_portions) = REAL(particle(MinN)%var(iAero_portions)/SUM(particle(MinN)%var(iAero_Portions)), KIND=lcmp)
            particle(MinN)%var(iN) = REAL(particle(MinN)%var(iN) + cf*Multf*collnum, KIND=lcmp)
            particle(MinN)%var(in_startind+1:in_startind+4) = REAL(particle(MinN)%var(in_startind+1:in_startind+4) &
                                                            + in_f*Multf*collnum, KIND=lcmp)
#ifdef SAMCHEM
            particle(MinN)%var(iaqspc) = REAL(particle(MinN)%var(iaqspc) + aqf*Multf*collnum, KIND=lcmp)
#endif /*SAMCHEM*/
            particle(s_ind)%var(iM) = REAL((Ms0*(particle(s_ind)%ivar(iA) - collnum)) + Ms1*collnum, KIND=lcmp)
            particle(s_ind)%var(iAero_portions) = REAL(particle(s_ind)%var(iAero_portions)*(particle(s_ind)%var(iN)/particle(s_ind)%ivar(iA) * (particle(s_ind)%ivar(iA) - collnum)) &
                                                + pm*cs*collnum, KIND=lcmp)
            particle(s_ind)%var(iAero_portions) = REAL(particle(s_ind)%var(iAero_portions)/SUM(particle(s_ind)%var(iAero_portions)), KIND=lcmp)
            particle(s_ind)%var(iN) = REAL((particle(s_ind)%var(iN)/particle(s_ind)%ivar(iA))*(particle(s_ind)%ivar(iA) - collnum) &
                                    + cs*collnum, KIND=lcmp)
            particle(s_ind)%var(in_startind+1:in_startind+4) = REAL(&
            (particle(s_ind)%var(in_startind+1:in_startind+4)/particle(s_ind)%ivar(iA))*(particle(s_ind)%ivar(iA) - collnum) &
            +in_s*collnum, KIND=lcmp)
#ifdef SAMCHEM
            particle(s_ind)%var(iaqspc) = REAL((particle(s_ind)%var(iaqspc)/particle(s_ind)%ivar(iA))*(particle(s_ind)%ivar(iA) - collnum) + aqs*collnum, KIND=lcmp)
#endif /*SAMCHEM*/
         ELSE !In case no appropriate merging particle within the gridbox was found, merge with the small drop remnant
            particle(s_ind)%var(iM) = REAL((Ms0*(particle(s_ind)%ivar(iA) - collnum)) + Ms1*collnum + &
                                       Mf*Multf*collnum, KIND=lcmp)
            particle(s_ind)%var(iAero_portions) = REAL(particle(s_ind)%var(iAero_portions)*(particle(s_ind)%var(iN)/particle(s_ind)%ivar(iA))*(particle(s_ind)%ivar(iA) - collnum) &
                                                  + pm*collnum*(cs + cf*Multf), KIND=lcmp)
            particle(s_ind)%var(iAero_portions) = REAL(particle(s_ind)%var(iAero_portions)/SUM(particle(s_ind)%var(iAero_portions)), KIND=lcmp)
            particle(s_ind)%var(iN) = REAL((particle(s_ind)%var(iN)/particle(s_ind)%ivar(iA))*(particle(s_ind)%ivar(iA) - collnum) &
                                    + cs*collnum + cf*Multf*collnum, KIND=lcmp)
            particle(s_ind)%var(in_startind+1:in_startind+4) = REAL(&
            (particle(s_ind)%var(in_startind+1:in_startind+4)/particle(s_ind)%ivar(iA))*(particle(s_ind)%ivar(iA) - collnum) &
            +in_s*collnum + in_f*Multf*collnum, KIND=lcmp)
#ifdef SAMCHEM
            particle(s_ind)%var(iaqspc) = REAL((particle(s_ind)%var(iaqspc)/particle(s_ind)%ivar(iA))*(particle(s_ind)%ivar(iA) - collnum) + collnum*(aqs + aqf*Multf), KIND=lcmp)
#endif /*SAMCHEM*/
            particle(s_ind)%ivar(iA) = INT(particle(s_ind)%ivar(iA) + Multf*collnum, KIND=8)
         ENDIF
      ELSE !Sheet or disk breakup, where the small droplet is destroyed and turns into another fragment
         
         particle(s_ind)%var(iM) = REAL(Ms0*(particle(s_ind)%ivar(iA) - collnum) + Mf*Multf*collnum, KIND=lcmp)
         particle(s_ind)%var(iAero_portions) = REAL(ps*(particle(s_ind)%var(iN)/particle(s_ind)%ivar(iA))*(particle(s_ind)%ivar(iA) - collnum) &
                                    + pm*cf*Multf*collnum, KIND=lcmp)
         particle(s_ind)%var(iAero_portions) = REAL(particle(s_ind)%var(iAero_portions)/SUM(particle(s_ind)%var(iAero_portions)), KIND=lcmp)
         particle(s_ind)%var(iN) = REAL((particle(s_ind)%var(iN)/particle(s_ind)%ivar(iA))*(particle(s_ind)%ivar(iA) - collnum) &
                                    + cf*Multf*collnum, KIND=lcmp)
         particle(s_ind)%var(in_startind+1:in_startind+4) = REAL(&
         (particle(s_ind)%var(in_startind+1:in_startind+4)/particle(s_ind)%ivar(iA))*(particle(s_ind)%ivar(iA) - collnum) &
         + in_f*Multf*collnum, KIND=lcmp)
#ifdef SAMCHEM
         particle(s_ind)%var(iaqspc) = REAL((particle(s_ind)%var(iaqspc)/particle(s_ind)%ivar(iA))*(particle(s_ind)%ivar(iA) - collnum) + aqf*multf*collnum, KIND=lcmp)
#endif /*SAMCHEM*/
         particle(s_ind)%ivar(iA) = INT(particle(s_ind)%ivar(iA) + (Multf-1)*collnum, KIND=8) !-1 because the small droplet was destroyed in the collision
      ENDIF
         
   ELSE !Here, follow the De Jong et al. (2023) approach. Say that all of the smaller-weighting factor SD's droplets have fragmented
      !And reduce the weighting factor of the initially larger-weighting-factor droplet. Incorporate collision-involved droplets in the smaller-weighting factor droplet
      IF (particle(n)%ivar(iA) .LT. particle(m)%ivar(iA)) THEN
         IF (Mults .gt. 0) THEN !Proxy for a filament breakup, since that's the only type of collision in which the small droplet survives
            particle(n)%ivar(iA) = INT(collnum*(Multf+2), KIND=8)
            particle(n)%var(iM) = REAL(collnum*(Ms1 + Ml1 + Multf*Mf), KIND=lcmp)
            particle(n)%var(iN) = REAL(collnum*(cs + cl + Multf*cf), KIND=lcmp)
            particle(n)%var(in_startind+1:in_startind+4) = REAL(collnum*(in_s + in_l + Multf*in_f), KIND=lcmp)
            particle(n)%var(iAero_portions) = REAL(pm, KIND=lcmp)
#ifdef SAMCHEM
            particle(n)%var(iaqspc) = REAL(collnum*(aqs + aql + Multf*aqf), KIND=lcmp)
#endif /*SAMCHEM*/
         ELSE !Sheet or disk breakup
            particle(n)%ivar(iA) = INT(collnum*(Multf+1), KIND=8)
            particle(n)%var(iM) = REAL(collnum*(Ml1 + Multf*Mf), KIND=lcmp)
            particle(n)%var(iN) = REAL(collnum*(cl + Multf*cf), KIND=lcmp)
            particle(n)%var(in_startind+1:in_startind+4) = REAL(collnum*(in_l + Multf*in_f), KIND=lcmp)
            particle(n)%var(iAero_portions) = REAL(pm, KIND=lcmp)
#ifdef SAMCHEM
            particle(n)%var(iaqspc) = REAL(collnum*(aql + Multf*aqf), KIND=lcmp)
#endif /*SAMCHEM*/
         ENDIF
         
         particle(m)%var(iM) = REAL(mm * (particle(m)%ivar(iA) - collnum), KIND=lcmp)
         particle(m)%var(iN) = REAL(particle(m)%var(iN)/particle(m)%ivar(iA)*(particle(m)%ivar(iA) - collnum), KIND=lcmp)
         particle(m)%var(in_startind+1:in_startind+4) = REAL(particle(m)%var(in_startind+1:in_startind+4)/particle(m)%ivar(iA) * &
                                                      (particle(m)%ivar(iA) - collnum), KIND=lcmp)
#ifdef SAMCHEM
         particle(m)%var(iaqspc) = REAL(particle(m)%var(iaqspc)/particle(m)%ivar(iA) * (particle(m)%ivar(iA) - collnum), KIND=lcmp)
#endif /*SAMCHEM*/
         particle(m)%ivar(iA) = INT(particle(m)%ivar(iA) - collnum, KIND=8)
      ELSEIF (particle(m)%ivar(iA) .LT. particle(n)%ivar(iA)) THEN
         IF (Mults .gt. 0) THEN !Filament breakup, where the small droplet survives
            particle(m)%ivar(iA) = INT(collnum*(Multf+2), KIND=8)
            particle(m)%var(iM) = REAL(collnum*(Ms1 + Ml1 + Multf*Mf), KIND=lcmp)
            particle(m)%var(iN) = REAL(collnum*(cs + cl + Multf*cf), KIND=lcmp)
            particle(m)%var(in_startind+1:in_startind+4) = REAL(collnum*(in_s + in_l + Multf*in_f), KIND=lcmp)
            particle(m)%var(iAero_portions) = REAL(pm, KIND=lcmp)
#ifdef SAMCHEM
            particle(m)%var(iaqspc) = REAL(collnum*(aqs + aql + Multf*aqf), KIND=lcmp)
#endif /*SAMCHEM*/
         ELSE !Sheet or disk breakup
            particle(m)%ivar(iA) = INT(collnum*(Multf+1), KIND=8)
            particle(m)%var(iM) = REAL(collnum*(Ml1 + Multf*Mf), KIND=lcmp)
            particle(m)%var(iN) = REAL(collnum*(cl + Multf*cf), KIND=lcmp)
            particle(m)%var(in_startind+1:in_startind+4) = REAL(collnum*(in_l + Multf*in_f), KIND=lcmp)
            particle(m)%var(iAero_portions) = REAL(pm, KIND=lcmp)
#ifdef SAMCHEM
            particle(m)%var(iaqspc) = REAL(collnum*(aql + Multf*aqf), KIND=lcmp)
#endif /*SAMCHEM*/
         ENDIF
         particle(n)%var(iM) = REAL(mn * (particle(n)%ivar(iA) - collnum), KIND=lcmp)
         particle(n)%var(iN) = REAL(particle(n)%var(iN)/particle(n)%ivar(iA)*(particle(n)%ivar(iA) - collnum), KIND=lcmp)
         particle(n)%var(in_startind+1:in_startind+4) = REAL(particle(n)%var(in_startind+1:in_startind+4)/particle(n)%ivar(iA) * &
                                                      (particle(n)%ivar(iA) - collnum), lcmp)
         particle(n)%ivar(iA) = INT(particle(n)%ivar(iA) - collnum, KIND=8)
#ifdef SAMCHEM
         particle(n)%var(iaqspc) = REAL(particle(n)%var(iaqspc)/particle(n)%ivar(iA) * (particle(n)%ivar(iA) - collnum), KIND=lcmp)
#endif /*SAMCHEM*/
      ELSE !Particles have equal weighting factor; in this case, spread the fragmentation equally among them
        !No droplets left unaffected - therefore, new mass is simply the mass of the new large, small, and fragment droplets
         IF (Mults > 0)THEN !Filament collision, where the small droplet survives
               newmass = collnum*(Ms1 + Ml1 + Mf*Multf)
               newccn = collnum*(cs + cl + Multf*cf)
               newin = collnum*(in_s + in_l + Multf*in_f)
               newnum = collnum*(Multf+2)
#ifdef SAMCHEM
               newspc = collnum*(aqs + aql + Multf*aqf)
#endif /*SAMCHEM*/
         ELSE
               newmass = collnum*(Ml1 + Mf*Multf)
               newccn = collnum*(cl + Multf*cf)
               newin = collnum*(in_l + Multf*in_f)
               newnum = collnum*(Multf+1)
#ifdef SAMCHEM
               newspc = collnum*(aql + Multf*aqf)
#endif /*SAMCHEM*/
         ENDIF
         particle(n)%ivar(iA) = floor(newnum / 2., KIND=8)
         particle(m)%ivar(iA) = INT(newnum - particle(n)%ivar(iA), KIND=8)
         particle(n)%var(iM) = REAL(newmass/newnum * particle(n)%ivar(iA), KIND=lcmp)
         particle(m)%var(iM) = REAL(newmass/newnum * particle(m)%ivar(iA), KINd=lcmp)
         particle(n)%var(iN) = REAL(newccn/newnum*particle(n)%ivar(iA), KIND=lcmp)
         particle(m)%var(iN) = REAL(newccn/newnum*particle(m)%ivar(iA), KIND=lcmp)
         particle(n)%var(in_startind+1:in_startind+4) = REAL(newin/newnum*particle(n)%ivar(iA), KIND=lcmp)
         particle(m)%var(in_startind+1:in_startind+4) = REAL(newin/newnum*particle(m)%ivar(iA), KIND=lcmp)
         particle(n)%var(iAero_portions) = REAL(pm, KIND=lcmp)
         particle(m)%var(iAero_portions) = REAL(pm, KIND=lcmp)
#ifdef SAMCHEM
         particle(n)%var(iaqspc) = REAL(newspc/newnum*particle(n)%ivar(iA), KIND=lcmp)
         particle(m)%var(iaqspc) = REAL(newspc/newnum*particle(m)%ivar(iA), KIND=lcmp)
#endif /*SAMCHEM*/
      ENDIF
   ENDIF

   !Update IN and freezing temperature information
   IF (sum(particle(n)%var(in_startind+1:in_startind+4)) .GT. 0.0) THEN
      particle(n)%ivar(iTYPE) = 2
      DO incounter = 1, 4
         particle(n)%var(iTFZ) = MAX(particle(n)%var(iTFZ), t_frz(incounter, particle(n)%var(in_startind+incounter)/particle(n)%ivar(iA)))
      END DO
   ENDIF
   IF (sum(particle(m)%var(in_startind+1:in_startind+4)) .GT. 0.0) THEN
      particle(m)%ivar(iTYPE) = 2
      DO incounter = 1, 4
         particle(m)%var(iTFZ) = MAX(particle(m)%var(iTFZ), t_frz(incounter, particle(m)%var(in_startind+incounter)/particle(m)%ivar(iA)))
      END DO
   ENDIF

   IF ((breakup_conserve_mode .eq. 1) .and. (MinN .gt. 0)) THEN
      IF (sum(particle(MinN)%var(in_startind+1:in_startind+4)) .GT. 0.0) THEN
         particle(MinN)%ivar(iTYPE) = 2
         DO incounter = 1, 4
            particle(MinN)%var(iTFZ) = MAX(particle(MinN)%var(iTFZ), t_frz(incounter, particle(MinN)%var(in_startind+incounter)/particle(MinN)%ivar(iA)))
         END DO
      ENDIF
   ENDIF

   micro_field(i,j,k,ibreak) = micro_field(i,j,k,ibreak) + collnum*(mn + mm)
   brk_prof(k) = brk_prof(k) + collnum*(mn + mm)

 END SUBROUTINE update_breakup

 SUBROUTINE update_coalescence(n,m,i,j,k,coalescence_count)

    USE micro_vars,   ONLY: n_start, n_end, particle, rho_liq, dV, diss_LCM,   &
                            lcmp, dtn_LCM, aut_prof, acc_prof,                 &
                            r_separate_micro, iA, iM, iN,                      &
                            iTYPE, in_startind, iTFZ, t_frz,                   &
                            coal_prof, iAero_portions

    USE micro_random, ONLY: randomu
    USE microphysics, ONLY: micro_field, icoal
    USE params,       ONLY: pi
    USE grid,         ONLY: rank
    USE vars,         ONLY: esatw, qsatw
#ifdef SAMCHEM
    USE microphysics, ONLY: iaqspc
    USE domain,       ONLY: naqspc
#endif /*SAMCHEM*/

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n, m, i, j, k
    INTEGER             :: intype
    INTEGER, INTENT(IN) :: coalescence_count
    REAL                :: xn, xm, xsn, xsm

    REAL, DIMENSION(4)  :: xim, xin
    REAL                :: x_separate, sfc_in, frztemp

#ifdef SAMCHEM
    REAL                :: nConc(naqspc), mConc(naqspc)
#endif /*SAMCHEM*/
!
!-- calculate separation mass
    x_separate = 4.0 / 3.0 * pi * rho_liq * r_separate_micro**3
!
!-- Compute mass of superdroplets n and m
    xn = particle(n)%var(iM) / particle(n)%ivar(iA)
    xm = particle(m)%var(iM) / particle(m)%ivar(iA)
!
!-- Compute aerosol mass of superdroplet n and m
    xsn = particle(n)%var(iN) / particle(n)%ivar(iA)
    xsm = particle(m)%var(iN) / particle(m)%ivar(iA)
!
!-- Compute IN surface area of superdroplets n and m
    xin = particle(n)%var(in_startind+1:in_startind+4) / particle(n)%ivar(iA)
    xim = particle(m)%var(in_startind+1:in_startind+4) / particle(m)%ivar(iA)

#ifdef SAMCHEM
!
!-- Compute number of atoms per represented droplet
!-- (analogous calculation to xm, xn, xsn, xsm for aqueous chemical species)
    IF (naqspc>0) THEN
      nConc = particle(n)%var(iaqspc) / particle(n)%ivar(iA)
      mConc = particle(m)%var(iaqspc) / particle(m)%ivar(iA)
    END IF
#endif /*SAMCHEM*/
!
!-- Each droplet of the super-droplet with the smaller weighting
!-- factor collects one droplet of the super-droplet with the larger
!-- weighting factor

!-- Currently, the way IN is handled in this routine is quite complicated; It would be simpler if we could
!-- keep track of IN surface area as a len-4 array within each particle. Ask Fabian about this
!-- B.Ascher 26.05.2025
    IF ( particle(n)%ivar(iA) .LT. particle(m)%ivar(iA) )  THEN

       particle(n)%var(iM)  = REAL( particle(n)%var(iM) + particle(n)%ivar(iA) * xm * coalescence_count, KIND=lcmp )
       particle(m)%var(iM)  = REAL( particle(m)%var(iM) - particle(n)%ivar(iA) * xm * coalescence_count, KIND=lcmp )

       ! this must be done before the new iN is calculated
       ! aerosol mode portions change for colliding particles (only)
       particle(n)%var(iAero_portions)  = REAL(   particle(n)%var(iAero_portions)*particle(n)%var(iN) &
                                        &       + particle(m)%var(iAero_portions)*particle(n)%ivar(iA) * xsm * coalescence_count &
                                        & , KIND=lcmp )
       particle(n)%var(iAero_portions)  = particle(n)%var(iAero_portions) / SUM(particle(n)%var(iAero_portions))

       particle(n)%var(iN)  = REAL( particle(n)%var(iN) + particle(n)%ivar(iA) * xsm * coalescence_count, KIND=lcmp )
       particle(m)%var(iN)  = REAL( particle(m)%var(iN) - particle(n)%ivar(iA) * xsm * coalescence_count, KIND=lcmp )
         
       particle(n)%var(in_startind+1:in_startind+4) = &
       REAL( particle(n)%var(in_startind+1:in_startind+4) + particle(n)%ivar(iA) * xim * coalescence_count, KIND=lcmp)
       particle(m)%var(in_startind+1:in_startind+4)  = &
       REAL( particle(m)%var(in_startind+1:in_startind+4) - particle(n)%ivar(iA) * xim * coalescence_count, KIND=lcmp)

       particle(m)%ivar(iA) = particle(m)%ivar(iA) - particle(n)%ivar(iA) * coalescence_count
#ifdef SAMCHEM
       IF (naqspc>0) THEN
         particle(n)%var(iaqspc) = REAL( particle(n)%var(iaqspc) + particle(n)%ivar(iA) * mConc * coalescence_count, KIND=lcmp )
         particle(m)%var(iaqspc) = REAL( particle(m)%var(iaqspc) - particle(n)%ivar(iA) * mConc * coalescence_count, KIND=lcmp )
       END IF
#endif /*SAMCHEM*/
      !
      !--    If true, particle n has or will get IN due to this interaction
      IF (particle(m)%ivar(iTYPE) .EQ. 2 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
         particle(n)%ivar(iTYPE) = 2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
         DO intype = 1,4
            sfc_in = particle(n)%var(in_startind+intype) / REAL(particle(n)%ivar(iA))
            particle(n)%var(iTFZ) = MAX(particle(n)%var(iTFZ), t_frz(intype, sfc_in))
         END DO
      ENDIF
!
!--    Analysis to determine autoconversion and accretion rates
       IF ( xm .LT. x_separate  .AND.  xn .LT. x_separate  .AND.  xm * coalescence_count + xn .GE. x_separate )  THEN
          aut_prof(k) = aut_prof(k) + ( xm * coalescence_count + xn ) * particle(n)%ivar(iA)
       ENDIF
       IF ( xm .GE. x_separate .AND. xn .LT. x_separate )  THEN
          acc_prof(k) = acc_prof(k) + xn * particle(n)%ivar(iA)
       ENDIF
       IF ( xm .LT. x_separate .AND. xn .GE. x_separate )  THEN
          acc_prof(k) = acc_prof(k) + xm * coalescence_count * particle(n)%ivar(iA)
       ENDIF

       micro_field(i,j,k,icoal) = micro_field(i,j,k,icoal) + ( xm * coalescence_count + xn ) * particle(n)%ivar(iA)
       coal_prof(k) = coal_prof(k) + ( xm * coalescence_count + xn ) * particle(n)%ivar(iA)

    ELSEIF ( particle(m)%ivar(iA) .LT. particle(n)%ivar(iA) )  THEN

       particle(m)%var(iM) = REAL( particle(m)%var(iM) + particle(m)%ivar(iA) * xn * coalescence_count, KIND=lcmp )
       particle(n)%var(iM) = REAL( particle(n)%var(iM) - particle(m)%ivar(iA) * xn * coalescence_count, KIND=lcmp )

       ! aerosol mode portions change for colliding particles (only)
       ! this must be done before the new iN is calculated
       particle(m)%var(iAero_portions)  = REAL(   particle(m)%var(iAero_portions)*particle(m)%var(iN) &
                                        &       + particle(n)%var(iAero_portions)*particle(m)%ivar(iA) * xsn * coalescence_count &
                                        & , KIND=lcmp )
       particle(m)%var(iAero_portions)  = particle(m)%var(iAero_portions) / SUM(particle(m)%var(iAero_portions))

       particle(m)%var(iN) = REAL( particle(m)%var(iN) + particle(m)%ivar(iA) * xsn * coalescence_count, KIND=lcmp )
       particle(n)%var(iN) = REAL( particle(n)%var(iN) - particle(m)%ivar(iA) * xsn * coalescence_count, KIND=lcmp )

       particle(m)%var(in_startind+1:in_startind+4) = &
       REAL (particle(m)%var(in_startind+1:in_startind+4) + particle(m)%ivar(iA) * xin * coalescence_count, KIND=lcmp)
       particle(n)%var(in_startind+1:in_startind+4) = &
       REAL (particle(n)%var(in_startind+1:in_startind+4) - particle(m)%ivar(iA) * xin * coalescence_count, KIND=lcmp)

#ifdef SAMCHEM
       IF (naqspc>0) THEN
         particle(m)%var(iaqspc) = REAL( particle(m)%var(iaqspc) + particle(m)%ivar(iA) * nConc * coalescence_count, KIND=lcmp )
         particle(n)%var(iaqspc) = REAL( particle(n)%var(iaqspc) - particle(m)%ivar(iA) * nConc * coalescence_count, KIND=lcmp )
       ENDIF
#endif /*SAMCHEM*/

       particle(n)%ivar(iA) = particle(n)%ivar(iA) - particle(m)%ivar(iA) * coalescence_count
      !
      !--    If true, particle m has or will get IN due to this interaction
       IF (particle(n)%ivar(iTYPE) .EQ. 2 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
         particle(m)%ivar(iTYPE) = 2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
         DO intype = 1,4
            sfc_in = particle(m)%var(in_startind+intype) / REAL(particle(m)%ivar(iA))
            particle(m)%var(iTFZ) = MAX(particle(m)%var(iTFZ), t_frz(intype, sfc_in))
         END DO
      ENDIF
!
!--    Analysis to determine autoconversion and accretion rates
       IF ( xm .LT. x_separate  .AND.  xn .LT. x_separate  .AND.  xm + xn * coalescence_count .GE. x_separate )  THEN
          aut_prof(k) = aut_prof(k) + ( xm + xn * coalescence_count ) * particle(m)%ivar(iA)
       ENDIF
       IF ( xn .GE. x_separate .AND. xm .LT. x_separate )  THEN
          acc_prof(k) = acc_prof(k) + xm * particle(m)%ivar(iA)
       ENDIF
       IF ( xn .LT. x_separate .AND. xm .GE. x_separate )  THEN
          acc_prof(k) = acc_prof(k) + xn * coalescence_count * particle(m)%ivar(iA)
       ENDIF

       micro_field(i,j,k,icoal) = micro_field(i,j,k,icoal) + ( xn * coalescence_count + xm ) * particle(m)%ivar(iA)
       coal_prof(k) = coal_prof(k) + ( xn * coalescence_count + xm ) * particle(m)%ivar(iA)
!
!-- If the weighting factors are equal or the droplets of the same superdroplet collide,
!-- the resultant droplet and aerosol mass is distributed evenly on the two interacting
!-- computational particles.
!--
!-- Multiple collisions cannot be considered here as their collection kernel changes
!-- significantly after the first collection since all droplets will have the same
!-- size, i.e., no differential sedimentation and hence no collsions in the second interaction
    ELSEIF ( particle(n)%ivar(iA) .EQ. particle(m)%ivar(iA) )  THEN

       particle(n)%ivar(iA) = FLOOR( particle(n)%ivar(iA) / 2.0, KIND=8 )
       particle(m)%ivar(iA) = particle(m)%ivar(iA) - particle(n)%ivar(iA)

       particle(n)%var(iM) = REAL(( xn + xm ) * particle(n)%ivar(iA), KIND=lcmp)
       particle(m)%var(iM) = REAL(( xn + xm ) * particle(m)%ivar(iA), KIND=lcmp)

       ! aerosol mode portions change for colliding particles (both super-droplets in this case)
       particle(m)%var(iAero_portions)  = (xsn*particle(n)%var(iAero_portions)+xsm*particle(m)%var(iAero_portions))
       particle(m)%var(iAero_portions)  = particle(m)%var(iAero_portions) / SUM(particle(m)%var(iAero_portions))
       particle(n)%var(iAero_portions)  = particle(m)%var(iAero_portions)

       particle(n)%var(iN) = REAL(( xsn + xsm ) * particle(n)%ivar(iA), KIND=lcmp)
       particle(m)%var(iN) = REAL(( xsn + xsm ) * particle(m)%ivar(iA), KIND=lcmp)

       particle(n)%var(in_startind+1:in_startind+4) = REAL((xin + xim) * particle(n)%ivar(iA), KIND=lcmp)
       particle(m)%var(in_startind+1:in_startind+4) = REAL((xin + xim) * particle(m)%ivar(iA), KIND=lcmp)

#ifdef SAMCHEM
       IF (naqspc>0) THEN
         particle(n)%var(iaqspc) = ( nConc + mConc ) * particle(n)%ivar(iA)
         particle(m)%var(iaqspc) = ( nConc + mConc ) * particle(m)%ivar(iA)
       ENDIF
#endif /*SAMCHEM*/

      !
      !--    If true, particle n and m will get IN due to this interaction
         IF ( particle(n)%ivar(iTYPE) .EQ. 2  .OR.  particle(m)%ivar(iTYPE) .EQ. 2 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
            particle(n)%ivar(iTYPE) = 2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
            particle(m)%ivar(iTYPE) = 2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
            DO intype = 1,4
               sfc_in = particle(n)%var(in_startind+intype) / REAL(particle(n)%ivar(iA))
               particle(n)%var(iTFZ) = MAX(particle(n)%var(iTFZ), t_frz(intype, sfc_in))
               particle(m)%var(iTFZ) = MAX(particle(m)%var(iTFZ), t_frz(intype, sfc_in))
            ENDDO
         ENDIF
!
!--    Analysis to determine autoconversion and accretion rates
       IF ( xm .LT. x_separate  .AND.  xn .LT. x_separate  .AND.  xm + xn .GE. x_separate )  THEN
          IF ( n .NE. m )  THEN
             aut_prof(k) = aut_prof(k) + ( xm + xn ) * ( particle(n)%ivar(iA) + particle(m)%ivar(iA) )
          ELSE ! internal collisions (n = m)
             aut_prof(k) = aut_prof(k) + ( xm + xn ) * particle(n)%ivar(iA)
          ENDIF
       ENDIF
       IF ( ( xm .GE. x_separate .AND. xn .LT. x_separate )  .OR.  ( xm .LT. x_separate .AND. xn .GE. x_separate ) )  THEN
          acc_prof(k) = acc_prof(k) + MIN( xm, xn ) * ( particle(n)%ivar(iA) + particle(m)%ivar(iA) )
       ENDIF

       micro_field(i,j,k,icoal) = micro_field(i,j,k,icoal) + ( xm + xn ) * ( particle(n)%ivar(iA) + particle(m)%ivar(iA) )
       coal_prof(k) = coal_prof(k) + ( xm + xn ) * ( particle(n)%ivar(iA) + particle(m)%ivar(iA) )

    ENDIF

 END SUBROUTINE update_coalescence

 SUBROUTINE determine_icecoll(m,n,i,j,k,n_particles,agg_count,icebreak_count)

    USE micro_vars,   ONLY: n_start, n_end, particle, rho_liq, dV, diss_LCM,   &
                            lcmp, dtn_LCM, switch_linear_sampling_collection,  &
                            iA, iM, iWsedi, iAXS, iCXS, iDNS, rho_ice, Imrime, &
                            switch_icebreak
    USE micro_random, ONLY: randomu
    USE microphysics, ONLY: micro_field
    USE params,       ONLY: pi
    USE grid,         ONLY: rank
    USE vars,         ONLY: esatw, qsatw

    IMPLICIT NONE

    REAL                 ::  An, Am, KK_agg, wsn, wsm, p_crit_lim, p_crit, &
                             rimefrac_n, rimefrac_m, KK_coll
    INTEGER              ::  half_n_particles
    INTEGER, INTENT(IN)  ::  m, n, i, j, k, n_particles
    INTEGER, INTENT(OUT) ::  agg_count, icebreak_count
    LOGICAL              ::  breakup_check1, breakup_check2

    REAL, PARAMETER      :: E_agg = 0.1 ! Morrison and Grabowski (2010) !Think about replacing this with temperature-dependent efficiency from Connolly et al. (2012)
!
    agg_count = 0
    icebreak_count = 0
!-- Collision is only calculated if one weighting factor is larger than 1 to avoid
!-- zero weighting factors
    IF ( MAX( particle(n)%ivar(iA), particle(m)%ivar(iA) ) .GT. 1_8 ) THEN
!
!--    Compute projected area of superdroplet n and m (Eqn. (5) and (6) in Shima et al. 2020)
       An = pi * MAX( particle(n)%var(iCXS), particle(n)%var(iAXS) )**2 * ( particle(n)%var(iDNS) / rho_ice )**EXP(-particle(n)%var(iCXS)/particle(n)%var(iAXS))
       Am = pi * MAX( particle(m)%var(iCXS), particle(m)%var(iAXS) )**2 * ( particle(m)%var(iDNS) / rho_ice )**EXP(-particle(m)%var(iCXS)/particle(m)%var(iAXS))
!
!--    Terminal velocity of superdroplet an and m
       wsn = particle(n)%var(iWsedi)
       wsm = particle(m)%var(iWsedi)
!
!--    Aggregation kernel (Eqn. (59) in Shima et al. (2020))
       KK_coll = ( SQRT(An) + SQRT(Am) )**2 * ABS( wsn - wsm )
!
!--    Probability of collision
       p_crit = MAX( particle(n)%ivar(iA), particle(m)%ivar(iA) ) / dV(k) * KK_coll * dtn_LCM
!
!--    Collision probability needs to be upscaled if the LSM method is used.
       IF ( switch_linear_sampling_collection )  THEN
          half_n_particles = FLOOR( n_particles / 2.0 )
          p_crit           = p_crit * n_particles * ( n_particles - 1.0 ) / ( half_n_particles * 2.0 )
       ENDIF

!--    Determine whether one of the particles is either graupel/hail or a large frozen droplet. If yes, treat this as a breakup event. If not, treat it as an aggregation event.
!--    In both Phillips et al., (2017) and the TAKAHrim scheme of Georgakaki et al., (2020), only collisions involving at least one rimed particle can lead to SIP
       rimefrac_n = particle(n)%var(iMrime)/particle(n)%var(iM)
       rimefrac_m = particle(m)%var(iMrime)/particle(m)%var(iM)
       IF (p_crit .gt. randomu(rank+1)) THEN
         breakup_check1 = ((MAX(rimefrac_n, rimefrac_m) .gt. 0.5) .and. (MAX(MAX(particle(n)%var(iCXS), particle(n)%var(iAXS)), MAX(particle(m)%var(iCXS), particle(m)%var(iAXS))) .ge. 50.0E-6)) !Rime check
         breakup_check2 = ((MAX(particle(n)%var(iAXS), particle(n)%var(iCXS)) .gt. 50.0E-6 .and. (particle(n)%var(iDNS) .eq. rho_ice)) .or. &
         (MAX(particle(m)%var(iAXS), particle(m)%var(iCXS)) .gt. 50.0E-6 .and. (particle(m)%var(iDNS) .eq. rho_ice))) !Large frozen drop check
         IF ((breakup_check1 .or. breakup_check2).and. switch_icebreak) THEN !Breakup
            icebreak_count = 1
         ELSE
            agg_count = MIN(NINT(E_agg * p_crit), 1)
         ENDIF
      ENDIF
   ENDIF

 END SUBROUTINE determine_icecoll

 SUBROUTINE update_icebreak_phillips(n,m,i,j,k)
   USE micro_vars,    ONLY: particle, icebrk_prof, iA, iM, iN, iMrime, &
                            iAXS, iCXS, iDNS, iTFZ, rho_ice, rho_liq, dV, t_frz, &
                            in_startind, IWsedi, tabs_LCM, n_start, n_end, iTYPE
   USE micro_random,  ONLY: randomu
   USE microphysics,  ONLY: micro_field, iicebreak
   USE params,        ONLY: pi
   use grid,          ONLY: rank, time

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: i,j,k,n,m
   INTEGER(KIND=8)     :: j_int, k_int, MinN, p, Nmax,            &
                          frag_ind, collnum, Nfrag, counter=0
   REAL                :: C, Ds, Dl, Vs, Vl, rimefrac_s, rimefrac_l, &
                          dd, deq, alpha, Am, a0, K0,     &
                          colldiam, T0, gamma, psi, zeta,      &
                          fragmass, efactor
   REAL                   mindiff, fragradius, phi,              &
                          oldaaxis, remaaxis, oldcaxis, remcaxis, &
                          oldmass, remmass, oldvolume, remvolume, &
                          oldrho, remrho, newaaxis, newcaxis,     &
                          newmass, newvolume, newrho, oldrime,    &
                          remrime, newrime, fragrime, fragphi,    &
                          ERFI, x_rand, ierf_dummy, realnfrag
   INTEGER             :: colltype !For now, unused as parameterization only considers "Mode 1" collisions

   REAL(KIND=8)        :: mm, mn

   !-- Ice collisional breakup as parameterized by Phillips et al. (2017) with modifications by Grzegorcyk et al. (2023) for graupel-dendrite collisions

   !Ds and Dl are the diameters of the smaller and larger, respectively
   !Vs and Vl are the fall speeds of the smaller and larger ice particles, respectively
   !rimefrac_s and rimefrac_l are the rime mass fractions of the smaller and larger ice particles, respectively
   !mm, mn are the masses per particle of SDs m and n, respectively
   !dd is the maximum dimension of the smaller particle
   !colldiam is similar to dd, but forced to be within the bounds of each collision class
   !deq is the spherical equivalent diameter of the smaller particle
   !alpha is the spherical equivalent surface area of the smaller particle
   !C is the kinetic energy required to break each fragment
   !Am is the surface density of breakable fragments
   !a0 is a threshold surface density
   !K0 is the collisional kinetic energy (Eq. 1)
   !T is the current grid point temperature
   !T0 is a threshold temperature (in K)
   !Nmax is the maximum number of fragments each collision type can produce
   !psi is a coefficient on the C term to correct for sublimation effects in observations
   !zeta is the fraction of the smaller particle's mass contained in each fragment
   !gamma is an exponential factor relating the ratio of collision kinetic energy to C to the number of fragments
   !realnfrag is the number of fragments produced by each collision (can be much less than 1)
   !Nfrag is the actual number of fragments per collision used in update of particle prroperties. It is rounded to the nearest integer
   !fragmass is the mass of each fragment
   !frag_ind is the index of the fragmenting particle
   !efactor is the term inside the exponential in Eq. 13 (separated here to reduce risk of typos)

   collnum = MIN(particle(n)%ivar(iA), particle(m)%ivar(iA))

   mn = particle(n)%var(iM)/REAL(particle(n)%ivar(iA), KIND=8)
   mm = particle(m)%var(iM)/REAL(particle(m)%ivar(iA), KIND=8)

   IF (particle(n)%var(iAXS) .gt. particle(m)%var(iAXS)) THEN
      j_int = n !Larger particle
      k_int = m !smaller
   ELSE
      j_int = m
      k_int = n
   ENDIF

   vs = particle(k_int)%var(iWsedi)
   vl = particle(j_int)%var(iWsedi)

   rimefrac_s = particle(k_int)%var(iMrime)/particle(k_int)%var(iM)
   rimefrac_l = particle(j_int)%var(iMrime)/particle(j_int)%var(iM)

   dd = 2.0*particle(k_int)%var(iAXS) !Currently, this parameterization is only specified for planar crystals, not columnar ones
   deq = 2*(particle(k_int)%var(iAXS)**2 * particle(k_int)%var(iCXS))**(1./3.)
   alpha = pi*deq**2
   K0 = 0.5*((mn*mm)/(mn+mm))*(vs - vl)**2 !Though we're mixing up m,n, and s and l, it doesn't matter in this equation

!-- Now, figure out which case of collision we need to simulate

   IF ((rimefrac_s .ge. 0.5) .and. (rimefrac_l .ge. 0.5) .and. (dd .lt. 5.0E-3)) THEN !Collision between graupel (0.5mm to 0.5cm) and hail
      colldiam = MAX(dd, 5.0E-4) !If dd is less than 0.5mm, set it to 0.5mm
      a0 = (3.78E4)*(1. + (0.0079/(dd**1.5)))
      T0 = -15.+273.16
      !psi = 3.5E-3 !Psi is included in the Phillips et al. (2017) parameterization in an attempt to account for boosting of fragment number concentrations by sublimation.
      !However, Grzegorczyk et al. (2023) clearly showed that sublimation is not necessary for large fragment numbers, and do not use a psi factor in their parameterizations.
      psi = 1.0
      C = (6.3E6)*psi
      Am = a0/3. + MAX(2.*a0/3. - (a0/9.)*ABS(tabs_lcm(i,j,k) - T0), 0.)
      gamma = 0.3
      zeta = 1.0E-3
      fragphi = 1.0 !Aspect ratio of fragments
      Nmax = 100
      IF ((tabs_LCM(i,j,k) .gt. 256.16) .and. (tabs_LCM(i,j,k) .lt. 260.16)) THEN !Use Grzegorczyk et al. (2023) parameterization around -15C
         zeta = 0.0
         C = 9.7E5
         Am = 1.6E6
         gamma = 0.55
         x_rand     = randomu(rank+1)
         ierf_dummy = MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )
         fragradius = 0.5E-6*EXP(2.0*SQRT(2.0)*0.69 * ierf_dummy + 4.86) !Exponential fit parameters to Eq. 4 in Grzegorczyk et al. (2023) for graupel with dendrites
         fragphi = 0.7
         Nmax = 500
      ENDIF
      frag_ind = k_int
   ELSEIF ((rimefrac_s .ge. 0.5) .and. (rimefrac_l .ge. 0.5) .and. (dd .ge. 5.0E-3)) THEN !Collision between hail particles (>0.5 cm)
      a0 = 4.35E5
      T0 = -15+273.16
      C = 3.31E5
      Am = a0/3. + MAX(2.*a0/3. - (a0/9.)*ABS(tabs_lcm(i,j,k) - T0), 0.)
      gamma = 0.54
      fragphi = 1.0
      zeta = 1.0E-6
      Nmax = 1000
      frag_ind = k_int
   ELSEIF ((rimefrac_s .lt. 0.5) .and. (particle(k_int)%var(iDNS) .le. 400)) THEN !Collision between dendrites and graupel/hail with smaller particle being the dendrite
      colldiam = MIN(MAX(dd, 5.0E-4),5.0E-3) !Force dd to be between 0.5mm and 0.5cm
      ! gamma = 0.5 - 0.25*rimefrac_s
      gamma = 0.78
      !psi = 3.5E-3 !Psi is included in the Phillips et al. (2017) parameterization in an attempt to account for boosting of fragment number concentrations by sublimation.
      !However, Grzegorczyk et al. (2023) clearly showed that sublimation is not necessary for large fragment numbers, and do not use a psi factor in their parameterizations.
      psi = 1.0
      ! C = (3.09E6)*psi
      C = 5.8E8
      ! Am = (1.41E6)*(1.+100.*(rimefrac_s**2))*(1 + (3.98E-5)/(colldiam**1.5))
      Am = 5.0E6
      ! zeta = 1.0E-3 !Not used by Grzegorczyk et al. (2023); size is instead parameterized with a lognormal distribution.
      x_rand     = randomu(rank+1)
      ierf_dummy = MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )
      fragradius = 0.5E-6*EXP(2.0*SQRT(2.0)*1.18 * ierf_dummy + 6.38) !Exponential fit parameters to Eq. 4 in Grzegorczyk et al. (2023) for graupel with dendrites
      fragphi = 0.5 !Aspect ratio of fragments
      zeta = 0.0
      Nmax = 500 !Hundreds of fragments were seen in Grzegorzyk et al. (2025)'s experiments, making a 100-fragment limit seem unrealistcally low
      frag_ind = k_int
   ELSEIF ((rimefrac_l .lt. 0.5) .and. (particle(j_int)%var(iDNS) .le. 400)) THEN !Collision between dendrites and graupel/hail with larger particle being the dendrite
      dd = 2.0*particle(j_int)%var(iAXS)
      colldiam = MIN(MAX(dd, 5.0E-4), 5.0E-3)
      ! gamma = 0.5 - 0.25*rimefrac_l
      gamma = 0.78
      ! psi = 3.5E-3
      psi = 1.0
      ! C = (3.09E6)*psi
      C = 5.8E8
      ! Am = (1.41E6)*(1.+100.*(rimefrac_l**2))*(1 + (3.98E-5)/(colldiam**1.5))
      Am = 5.0E6
      ! zeta = 1.0E-3 !Not used by Grzegorczyk et al. (2025); size is instead parameterized with a lognormal distribution.
      x_rand     = randomu(rank+1)
      ierf_dummy = MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )
      fragradius = 0.5E-6*EXP(2.0*SQRT(2.0)*1.18 * ierf_dummy + 6.38) !Exponential fit parameters to Eq. 4 in Grzegorczyk et al. (2023) for graupel with dendrites
      fragphi = 0.5 !Aspect ratio of fragments
      zeta = 0.0
      Nmax = 500
      frag_ind = j_int
   ELSEIF ((rimefrac_s .lt. 0.5) .and. (particle(k_int)%var(iDNS) .ge. 400)) THEN !Collision between plates and graupel/hail with the smaller particle being the plate
      colldiam = MIN(MAX(dd, 5.0E-4), 5.0E-3)
      gamma = 0.50 - 0.25*rimefrac_s
      ! psi = 3.5E-3
      psi = 1.0
      C = (7.08E6)*psi
      Am = (1.58E7)*(1.0 + 100.*(rimefrac_s**2))*(1.0 + (1.33E-4)/(dd**1.5))
      fragphi = 1.0
      zeta = 1.0E-3
      Nmax = 100
      frag_ind = k_int
   ELSEIF ((rimefrac_l .lt. 0.5) .and. (particle(j_int)%var(iDNS) .ge. 400)) THEN !Collision between plates and graupel/hail with the larger particle being the plate
      dd = particle(j_int)%var(iAXS)
      gamma = 0.50 - 0.25*rimefrac_l
      ! psi = 3.5E-3
      psi = 1.0
      C = (7.08E6)*psi
      Am = (1.58E7)*(1.0 + 100.*(rimefrac_l**2))*(1.0 + (1.33E-4)/(dd**1.5))
      fragphi = 1.0
      zeta = 1.0E-3
      Nmax = 100
      frag_ind = j_int
   ENDIF

   efactor = (C*K0/(alpha*Am))**gamma
   realnfrag = MIN(alpha*Am*(1.0 - EXP(-1*efactor)), REAL(Nmax, KIND=8))

   Nfrag = FLOOR(realnfrag)
   x_rand = randomu(rank+1)
   IF (x_rand .lt. (realnfrag - Nfrag)) Nfrag = Nfrag + 1

   IF (Nfrag .gt. 0) THEN !Only do particle updates if we actually produced fragments from the collision. Otherwise, treat it as though the particles "bounced off" each other
   
!-- Now, try to assign the fragments to an existing ice particle within the gridbox, similar to our approach with liquid droplet collisional breakup
      oldaaxis = particle(frag_ind)%var(iAXS) !Old a-axis length of the fragile crystal
      oldcaxis = particle(frag_ind)%var(iCXS) !Old c-axis length of the fragile crystal
      oldvolume = 4./3. * pi * oldaaxis**2*oldcaxis !Old per-crystal volume of the fragile crystal
      oldmass = particle(frag_ind)%var(iM)/particle(frag_ind)%ivar(iA) !Old per-crystal mass of the fragile crystal
      oldrho = particle(frag_ind)%var(iDNS) !Old effective density of the fragile crystal
      oldrime = particle(frag_ind)%var(iMrime)/particle(frag_ind)%ivar(iA) !Old per-crystal rime mass of the fragile crystal

      IF (fragphi .lt. 0.9) THEN !Check whether we're using Grzegorczyk et al. (2023)'s parameterization with a lognormal size distribution
         fragmass = rho_ice * 4.0/3.0 * pi * fragradius**3 * fragphi
         !Now, need to make sure that we're not making fragments which have more mass than the initial ice crystal!
         DO WHILE(Nfrag*fragmass .gt. oldmass)
            counter = counter+1
            newmass = fragmass*(0.9**counter)
            fragradius = fragradius*(newmass/fragmass)**(1./3.)
            fragmass = newmass
         END DO
      ELSE
         fragmass = zeta * particle(frag_ind)%var(iM)/particle(frag_ind)%ivar(iA) !Mass per fragment
         !Figure out what to do with particle aspect ratio, caxis, axis. My first instinct would be to say that each fragment, being a spike, should be somewhat columnar. 
         !But for now, let's do isormophic for simplicity (revisit this with some better data on typical shapes of fragments - 14.10.2025)
         fragradius = (fragmass * 3./(4. * pi * rho_ice))**(1./3.) !Assume fragments are solid ice
      ENDIF
      fragrime = oldrime * (fragmass/oldmass) !Mass of rime per fragment

      remmass = oldmass - fragmass*Nfrag !Per-crystal mass of remnant fragile crystals (those which have been impacted)
      remvolume = oldvolume - Nfrag*4./3.*pi*(fragradius**3 * fragphi) !Per-crystal volume of remnant fragile crystals (those which have been impacted)
      remrho = remmass/remvolume !Effective density of remnant fragile crystals (those which have been impacted)
      remrime = oldrime - fragrime*Nfrag !Per-crystal rime mass of remnant fragile crystals (those which have been impacted)

      MinN = -1
      Mindiff = ABS(fragmass - particle(frag_ind)%var(iM)/particle(frag_ind)%ivar(iA))
      DO p = n_start(i,j,k), n_end(i,j,k)
         IF ((ABS(particle(p)%var(iM)/particle(p)%ivar(iA) - fragmass) .lt. Mindiff) .and. ((particle(p)%var(iM)/particle(p)%ivar(iA))/fragmass .le. 2.0) .and. &
         ((particle(p)%var(iM)/particle(p)%ivar(iA))/fragmass .ge. 0.5) .and. (particle(p)%ivar(iTYPE) .lt. 0)) THEN
            MinN = p
            Mindiff = ABS(particle(p)%var(iM)/particle(p)%ivar(iA) - fragmass)
         ENDIF
      END DO

      IF (MinN .gt. 0) THEN !Successfully assigned fragments to another particle in the gridbox
         !Adjust aspect ratio based on volume-weighted aspect ratios of pre-existing ice crystal and new fragments
         phi = particle(MinN)%var(iCXS)/particle(MinN)%var(iAXS)
         phi = (phi*(particle(MinN)%var(iCXS)*particle(MinN)%var(iAXS)**2)*(particle(MinN)%ivar(iA)) + fragphi**2*fragradius**3*Nfrag*collnum)/((particle(MinN)%var(iCXS)*particle(MinN)%var(iAXS)**2)*(particle(MinN)%ivar(iA)) + fragphi*fragradius**3*Nfrag*collnum)
         !fragradius**3*fragphi is the fragment's volume; the additional factor of fragphi in the numerator is to weight the fragment's aspect ratio by the volume
         particle(MinN)%var(iDNS) = (particle(MinN)%var(iDNS)*particle(MinN)%var(iM) + rho_ice*fragmass*Nfrag*collnum)/(particle(MinN)%var(iM)+fragmass*Nfrag*collnum)
         particle(MinN)%var(iM) = particle(MinN)%var(iM) + fragmass*Nfrag*collnum
         particle(MinN)%ivar(iA) = particle(MinN)%ivar(iA) + Nfrag*collnum
         particle(MinN)%var(iCXS) = (3.*particle(MinN)%var(iM)*phi**2/(4.*pi*particle(MinN)%var(iDNS)*particle(MinN)%ivar(iA)))**(1./3.)
         particle(MinN)%var(iAXS) = particle(MinN)%var(iCXS)/phi
         particle(MinN)%var(iMrime) = particle(MinN)%var(iMrime) + fragrime*Nfrag*collnum
         !Now, update properties of the old fragile crystal

         particle(frag_ind)%var(iM) = MAX(particle(frag_ind)%var(iM) - fragmass*Nfrag*collnum, 4.0/3.0*pi*rho_ice*(0.1E-6)**3)

         newrho = (oldrho*oldmass*(particle(frag_ind)%ivar(iA) - collnum) + remrho*remmass*collnum)/(particle(frag_ind)%var(iM))
         particle(frag_ind)%var(iDNS) = newrho
         newmass = particle(frag_ind)%var(iM)/particle(frag_ind)%ivar(iA)
         phi = oldcaxis/oldaaxis
         particle(frag_ind)%var(iCXS) = (3.*newmass*phi**2/(4.*pi*newrho))**(1./3.)
         particle(frag_ind)%var(iAXS) = particle(frag_ind)%var(iCXS)/phi
         particle(frag_ind)%var(iMrime) = oldrime*(particle(frag_ind)%ivar(iA) - collnum) + remrime*collnum
         particle(frag_ind)%ivar(iA) = particle(frag_ind)%ivar(iA)
         !Keep ice crystal aspect ratio the same. Revisit this assumption - should they be preferentially removed from the major axis? -B.Ascher 14.10.2025

         !Leave aerosol mass, IN surface area, temperature, specific humidity unchanged. Don't really need to worry about velocity, since it'll be updated in the next timestep anyway
      ELSE !If unsuccessful, treat the fragments as remaining with the shattered particle
         phi = oldcaxis/oldaaxis
         phi = (phi*oldvolume*(particle(frag_ind)%ivar(iA)-collnum) + phi*remvolume*collnum+ 1.0*4./3.*pi*fragradius**3*Nfrag*collnum)/(oldvolume*(particle(frag_ind)%ivar(iA)-collnum) + remvolume*collnum + 4./3.*pi*fragradius**3*Nfrag*collnum)
         particle(frag_ind)%var(iDNS) = (oldrho*oldmass*(particle(frag_ind)%ivar(iA) - collnum) + remrho*remmass*collnum + rho_ice*Nfrag*fragmass*collnum)/(oldmass*(particle(frag_ind)%ivar(iA) - collnum) + remmass*collnum + fragmass*Nfrag*collnum)
         newrho = particle(frag_ind)%var(iDNS)
         particle(frag_ind)%ivar(iA) = particle(frag_ind)%ivar(iA) + Nfrag*collnum
         newmass = particle(frag_ind)%var(iM)/particle(frag_ind)%ivar(iA)
         particle(frag_ind)%var(iCXS) = (3.*newmass*phi**2/(4.*pi*newrho))**(1./3.)
         particle(frag_ind)%var(iAXS) = particle(frag_ind)%var(iCXS)/phi
      ENDIF
      micro_field(i,j,k,iicebreak) = micro_field(i,j,k,iicebreak) + collnum*Nfrag
      icebrk_prof(k) = icebrk_prof(k) + collnum*Nfrag

   ENDIF

 END SUBROUTINE update_icebreak_phillips

 SUBROUTINE update_icebreak_sullivan(n,m,i,j,k)
   USE micro_vars,    ONLY: particle, dtn_LCM, icebrk_prof, iA, iM, iN, iTYPE, iMrime, &
                            iAXS, iCXS, iDNS, iTFZ, rho_ice, rho_liq, dV, t_frz,       &
                            in_startind, IWsedi, tabs_LCM, n_start, n_end, iAero_portions
   USE micro_random,  ONLY: randomu
   USE microphysics,  ONLY: micro_field, iicebreak
   USE params,        ONLY: pi
   use grid,          ONLY: rank

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: i,j,k,n,m
   INTEGER(KIND=8)     :: j_int, k_int, MinN, p, frag_ind, collnum
   REAL                :: xn, xm, rimefrac_s, rimefrac_l
   REAL                   mindiff, fragradius, phi,              &
                          oldaaxis, oldcaxis, oldmass, oldvolume, &
                          oldrho, newaaxis, newcaxis, newmass,    &
                          newvolume, newrho, oldrime, newrime,    &
                          oldaeromass, fragmass, fragrime, nfrag, &
                          dfragile
   REAL, PARAMETER :: d0 = 0.02, Tmin = 252 !D0 from eq. B14 in Sotiropoulou et al. (2021), Tmin from Table 1 in Sullivan et al. (2018)
   REAL, DIMENSION(4)  :: oldaeroports, oldins

   collnum = MIN(particle(n)%ivar(iA), particle(m)%ivar(iA))

   xn = particle(n)%var(iM)/particle(n)%ivar(iA)
   xm = particle(m)%var(iM)/particle(n)%ivar(iA)

   IF (particle(n)%var(iAXS) .gt. particle(m)%var(iAXS)) THEN
      j_int = n !Larger particle
      k_int = m !smaller
   ELSE
      j_int = m
      k_int = n
   ENDIF

   rimefrac_s = particle(k_int)%var(iMrime)/particle(k_int)%var(iM)
   rimefrac_l = particle(j_int)%var(iMrime)/particle(j_int)%var(iM)

!-- Though not specified in Sullivan et al. (2018), we'll assign fragility the same way as in Phillips et al. (2018) - the 'unrimed' particle, or the smaller particle if both are heavily rimed
!-- Calculate number of fragments produced by ice-ice collisions, following the parameterization of Sullivan et al. (2018) modified by Sotiropoulou et al. (2021)

   IF (rimefrac_s .le. 0.5) THEN
      frag_ind = k_int
   ELSEIF ((rimefrac_s .ge. 0.5) .and. (rimefrac_l .le. 0.5)) THEN
      frag_ind = j_int
   ENDIF
   dfragile = MAX(particle(frag_ind)%var(iAXS), particle(frag_ind)%var(iCXS)) * 2.0
   Nfrag = 280 * (MAX(tabs_LCM(i,j,k), Tmin) - Tmin)**(1.2)*EXP((-1*(MAX(tabs_LCM(i,j,k), Tmin)-Tmin)/5.)) * dfragile / d0 !Temperature-dependent fragment number, with no fragmentation below 252 K.
   
   !In Sullivan et al. (2018), there is no information about the size of fragments produced by the breakups. Therefore, assume here that the fragile ice crystal completely fragments
!-- The coefficient in front would normally be 280, as in Sullivan et al. (2018), but Sotiropoulou et al. (2021) found that the parameterization is more accurate when scaled by particle size.

   fragmass = 1./Nfrag * particle(frag_ind)%var(iM)/particle(frag_ind)%ivar(iA)
   fragradius = (3. * fragmass / (4. * pi * rho_ice))**(1./3.)

   oldaaxis = particle(frag_ind)%var(iAXS) !Old a-axis length of the fragile crystal
   oldcaxis = particle(frag_ind)%var(iCXS) !Old c-axis length of the fragile crystal
   oldvolume = 4./3. * pi * oldaaxis**2*oldcaxis !Old per-crystal volume of the fragile crystal
   oldmass = particle(frag_ind)%var(iM)/particle(frag_ind)%ivar(iA) !Old per-crystal mass of the fragile crystal
   oldrho = particle(frag_ind)%var(iDNS) !Old effective density of the fragile crystal
   oldrime = particle(frag_ind)%var(iMrime)/particle(frag_ind)%ivar(iA) !Old per-crystal rime mass of the fragile crystal
   oldaeroports = particle(frag_ind)%var(iAero_portions)
   oldaeromass = particle(frag_ind)%var(iN) / particle(frag_ind)%ivar(iA)
   oldins = particle(frag_ind)%var(in_startind+1:in_startind+4) / particle(frag_ind)%ivar(iA)

   MinN = -1
   Mindiff = ABS(fragmass - particle(frag_ind)%var(iM)/particle(frag_ind)%ivar(iA))
   DO p = n_start(i,j,k), n_end(i,j,k)
      IF ((ABS(particle(p)%var(iM)/particle(p)%ivar(iA) - fragmass) .lt. Mindiff) .and. ((particle(p)%var(iM)/particle(p)%ivar(iA))/fragmass .le. 2.0) .and. &
      ((particle(p)%var(iM)/particle(p)%ivar(iA))/fragmass .ge. 0.5) .and. (particle(p)%ivar(iTYPE) .lt. 0)) THEN
         MinN = p
         Mindiff = ABS(particle(p)%var(iM)/particle(p)%ivar(iA) - fragmass)
      ENDIF
   END DO

   IF ((MinN .gt. 0) .and. (particle(frag_ind)%ivar(iA) .gt. collnum)) THEN !Successfully assigned fragments to another particle in the gridbox. Only do this if there are unfragmented superparticles left
      
      !Adjust aspect ratio based on volume-weighted aspect ratios of pre-existing ice crystal and new fragments
      phi = particle(MinN)%var(iCXS)/particle(MinN)%var(iAXS)
      phi = (phi*(particle(MinN)%var(iCXS)*particle(MinN)%var(iAXS)**2)*(particle(MinN)%ivar(iA)) + 1.0*fragradius**3*Nfrag*collnum)/((particle(MinN)%var(iCXS)*particle(MinN)%var(iAXS)**2)*(particle(MinN)%ivar(iA)) + fragradius**3*Nfrag*collnum)
      particle(MinN)%var(iDNS) = (particle(MinN)%var(iDNS)*particle(MinN)%var(iM) + rho_ice*oldmass*collnum)/(particle(MinN)%var(iM)+oldmass*collnum)
      particle(MinN)%var(iAero_portions) = (particle(MinN)%var(iAero_portions)*particle(MinN)%var(iN) + oldaeroports*oldaeromass*collnum)/(particle(MinN)%var(iN) + oldaeromass*collnum)
      particle(MinN)%var(iN) = particle(MinN)%var(iN) + oldaeromass*collnum
      particle(MinN)%var(in_startind+1:in_startind+4) = particle(MinN)%var(in_startind+1:in_startind+4) + oldins * collnum
      particle(MinN)%var(iM) = particle(MinN)%var(iM) + oldmass*collnum
      particle(MinN)%ivar(iA) = particle(MinN)%ivar(iA) + FLOOR(Nfrag*collnum, KIND=8)
      particle(MinN)%var(iCXS) = (3.*particle(MinN)%var(iM)*phi**2/(4.*pi*particle(MinN)%var(iDNS)*particle(MinN)%ivar(iA)))**(1./3.)
      particle(MinN)%var(iAXS) = particle(MinN)%var(iCXS)/phi
      particle(MinN)%var(iMrime) = particle(MinN)%var(iMrime) + oldrime*collnum
      IF (SUM(particle(MinN)%var(in_startind+1:in_startind+4)) .gt. 0.0) THEN
         particle(MinN)%ivar(iTYPE) = -2
         particle(MinN)%var(iTFZ) = MAX(particle(MinN)%var(iTFZ), particle(frag_ind)%var(iTFZ))
      ENDIF


      !Now, update properties of unfragmented fragile crystals
      particle(frag_ind)%var(iN) = particle(frag_ind)%var(iN) - oldaeromass*collnum
      particle(frag_ind)%var(in_startind+1:in_startind+4) = particle(frag_ind)%var(in_startind+1:in_startind+4) - oldins*collnum
      particle(frag_ind)%var(iMrime) = particle(frag_ind)%var(iMrime) - oldrime*collnum
      particle(frag_ind)%var(iM) = particle(frag_ind)%var(iM) - oldmass*collnum
      particle(frag_ind)%ivar(iA) = particle(frag_ind)%ivar(iA) - collnum
   ELSE !Either no similar particle successfully found in gridbox, or the more fragile particle is also less numerous (all particles fragment)
      phi = oldcaxis/oldaaxis
      particle(frag_ind)%ivar(iA) = particle(frag_ind)%ivar(iA) + FLOOR(Nfrag*collnum, KIND=8)
      newmass = particle(frag_ind)%var(iM)/particle(frag_ind)%ivar(iA)
      particle(frag_ind)%var(iCXS) = (3.*newmass*phi**2/(4.*pi*particle(n)%var(iDNS)))**(1./3.)
      particle(frag_ind)%var(iAXS) = particle(frag_ind)%var(iCXS)/phi
   ENDIF

   micro_field(i,j,k,iicebreak) = micro_field(i,j,k,iicebreak) + collnum*Nfrag
   icebrk_prof(k) = icebrk_prof(k) + collnum*Nfrag

END SUBROUTINE update_icebreak_sullivan

 SUBROUTINE update_aggregation(n,m,i,j,k)

    USE micro_vars,   ONLY: n_start, n_end, particle, rho_liq, dV, diss_LCM,   &
                            lcmp, dtn_LCM, agg_prof, r_separate_micro, iMrime,   &
                            iNmono, iA, iM, iN, iTYPE, iAXS, iCXS, iDNS, &
                            rho_ice, iTFZ, in_startind, t_frz
    USE micro_random, ONLY: randomu
    USE microphysics, ONLY: micro_field, iagg
    USE params,       ONLY: pi
    USE grid,         ONLY: rank

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n, m, i, j, k
    INTEGER             :: j_int, k_int, intype

    REAL :: rho_apparent, rho_max, rho_min, rho_crit, V_j, V_k, V_max, x_j, &
            x_k, xrm, xrn, aaxis, caxis, xn, xm, xsn, xsm, sfc_in
    REAL, DIMENSION(4) :: xim, xin

    REAL, PARAMETER :: rho_ice_crit = 10.0 ! in kg/m3, critical ice density
!
!-- Check for the larger ice crystal such that j > k
    IF ( MAX(particle(n)%var(iAXS),particle(n)%var(iCXS)) .GT. &
         MAX(particle(m)%var(iAXS),particle(m)%var(iCXS)) )  THEN
       j_int = n !larger
       k_int = m !smaller
    ELSE
       j_int = m !larger
       k_int = n !smaller
    ENDIF
    x_j = particle(j_int)%var(iM) / particle(j_int)%ivar(iA)
    x_k = particle(k_int)%var(iM) / particle(k_int)%ivar(iA)
!
!-- Volumina of particles (before (26) in Shima et al. (2020))
    V_j = 4.0 / 3.0 * pi * particle(j_int)%var(iAXS)**2 * particle(j_int)%var(iCXS)
    V_k = 4.0 / 3.0 * pi * particle(k_int)%var(iAXS)**2 * particle(k_int)%var(iCXS)
!
!-- Estimate maximum density of aggregate
    rho_max  = ( x_j + x_k ) / ( V_j + V_k )
!
!-- Update axes for plates (c < a) first, then columns (c > a). It is assumed
!-- that the maximum dimension will be preserved during aggregation.
    IF ( particle(j_int)%var(iCXS) / particle(j_int)%var(iAXS) .LE. 1.0 )  THEN
!
!--    Maximum volume of aggregate
       V_max = 4.0 / 3.0 * pi * particle(j_int)%var(iAXS)**2 * &
               ( particle(j_int)%var(iCXS) + MIN( particle(k_int)%var(iAXS), particle(k_int)%var(iCXS) ) )
!
!--    Minimum density of aggregate, (61) in Shima et al. (2020)
       rho_min = ( x_j + x_k ) / V_max
!
!--    New apparent density of aggregate, (62) in Shima et al. (2020)
       rho_apparent = ( ( rho_ice - rho_max ) * rho_max + &
                        ( rho_max - rho_crit ) * rho_min  &
                      ) / ( rho_ice - rho_crit )
!
!--    A-axis of aggregate (kept)
       aaxis = particle(j_int)%var(iAXS)
!
!--    C-axis of aggregate, (63) in Shima et al. (2020)
       caxis = ( x_j + x_k ) / ( 4.0 / 3.0 * pi * rho_apparent * aaxis**2 )

    ELSE !column
!
!--    Maximum volume of aggregate
       V_max = 4.0 / 3.0 * pi * particle(j_int)%var(iCXS) * &
               ( particle(j_int)%var(iAXS) + MIN( particle(k_int)%var(iAXS), particle(k_int)%var(iCXS) ) ) * &
               MAX( particle(j_int)%var(iAXS), particle(k_int)%var(iAXS), particle(k_int)%var(iCXS) )
!
!--    Minimum density of aggregate, (61) in Shima et al. (2020)
       rho_min = ( x_j + x_k ) / V_max
!
!--    New apparent density if aggregate, (62) in Shima et al. (2020)
       rho_apparent = ( ( rho_ice - rho_max ) * rho_max + &
                        ( rho_max - rho_crit ) * rho_min  &
                      ) / ( rho_ice - rho_crit )
!
!--    C-axis of aggregate (kept)
       caxis = particle(j_int)%var(iCXS)
!
!--    A-axis of aggregate, (70) in Shima et al. (2020)
       aaxis = SQRT( ( x_j + x_k ) / ( 4.0 / 3.0 * pi * rho_apparent * caxis ) )

    ENDIF
!
!-- Compute mass of superdroplets n and m
    xn = particle(n)%var(iM) / particle(n)%ivar(iA)
    xm = particle(m)%var(iM) / particle(m)%ivar(iA)
!
!-- Compute aerosol mass of superdroplet n and m
    xsn = particle(n)%var(iN) / particle(n)%ivar(iA)
    xsm = particle(m)%var(iN) / particle(m)%ivar(iA)
!
!-- Compute IN surface area of superdroplets n and m
    xin = particle(n)%var(in_startind+1:in_startind+4) / particle(n)%ivar(iA)
    xim = particle(m)%var(in_startind+1:in_startind+4) / particle(m)%ivar(iA)
!
!-- Compute rimed water mass of superdroplet n and m
    xrn = particle(n)%var(iMrime) / particle(n)%ivar(iA)
    xrm = particle(m)%var(iMrime) / particle(m)%ivar(iA)
!
!-- Each droplet of the super-droplet with the smaller weighting
!-- factor collects one droplet of the super-droplet with the larger
!-- weighting factor
    IF ( particle(n)%ivar(iA) .LT. particle(m)%ivar(iA) )  THEN

       particle(n)%var(iM)  = REAL( particle(n)%var(iM) + particle(n)%ivar(iA) * xm, KIND=lcmp )
       particle(m)%var(iM)  = REAL( particle(m)%var(iM) - particle(n)%ivar(iA) * xm, KIND=lcmp )

       particle(n)%var(iN)  = REAL( particle(n)%var(iN) + particle(n)%ivar(iA) * xsm, KIND=lcmp )
       particle(m)%var(iN)  = REAL( particle(m)%var(iN) - particle(n)%ivar(iA) * xsm, KIND=lcmp )
      
       particle(n)%var(in_startind+1:in_startind+4) = &
       REAL( particle(n)%var(in_startind+1:in_startind+4) + particle(n)%ivar(iA) * xim, KIND=lcmp)
       particle(m)%var(in_startind+1:in_startind+4)  = &
       REAL( particle(m)%var(in_startind+1:in_startind+4) - particle(n)%ivar(iA) * xim, KIND=lcmp)

       particle(m)%ivar(iA) = particle(m)%ivar(iA) - particle(n)%ivar(iA)
      !
      !--    If true, particle n has or will get IN due to this interaction
      IF (particle(m)%ivar(iTYPE) .EQ. -2 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
         particle(n)%ivar(iTYPE) = -2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
         DO intype = 1,4
            sfc_in = particle(n)%var(in_startind+intype) / REAL(particle(n)%ivar(iA))
            particle(n)%var(iTFZ) = MAX(particle(n)%var(iTFZ), t_frz(intype, sfc_in))
         END DO
      ENDIF
!
!--    Update of monomer number
       particle(n)%ivar(iNmono) = particle(n)%ivar(iNmono) + particle(m)%ivar(iNmono)
!
!--    Update of rimed water mass
       particle(n)%var(iMrime)  = REAL( particle(n)%var(iMrime) + particle(n)%ivar(iA) * xrm, KIND=lcmp )
       particle(m)%var(iMrime)  = REAL( particle(m)%var(iMrime) - particle(n)%ivar(iA) * xrm, KIND=lcmp )
!
!--    Update of axes and densities
       particle(n)%var(iDNS)  = REAL( rho_apparent, KIND=lcmp )
       particle(n)%var(iAXS)  = REAL( aaxis, KIND=lcmp )
       particle(n)%var(iCXS)  = REAL( caxis, KIND=lcmp )
!
!--    Aggregation profile
       IF ( particle(n)%ivar(iNmono) .GT. 10 )  THEN ! new aggregate is snow ( number monomers > 10 )
          IF ( particle(n)%ivar(iNmono) - particle(m)%ivar(iNmono) .LE. 10  .AND.  particle(m)%ivar(iNmono) .LE. 10 )  THEN ! n and m were ice before aggregation
             agg_prof(k) = agg_prof(k) + ( xm + xn ) * particle(n)%ivar(iA)
             micro_field(i,j,k,iagg) = micro_field(i,j,k,iagg) + (xm + xn) * particle(n)%ivar(iA)
          ENDIF
          IF ( particle(m)%ivar(iNmono) .LE. 10 )  THEN ! m was ice before aggregation
             agg_prof(k) = agg_prof(k) + xm * particle(n)%ivar(iA)
             micro_field(i,j,k,iagg) = micro_field(i,j,k,iagg) + xm * particle(n)%ivar(iA)
          ENDIF
          IF ( particle(n)%ivar(iNmono) - particle(m)%ivar(iNmono) .LE. 10 )  THEN ! n was ice before aggregation
             agg_prof(k) = agg_prof(k) + xn * particle(n)%ivar(iA)
             micro_field(i,j,k,iagg) = micro_field(i,j,k,iagg) + xn * particle(n)%ivar(iA)
          ENDIF
       ENDIF

    ELSEIF ( particle(m)%ivar(iA) .LT. particle(n)%ivar(iA) )  THEN

       particle(m)%var(iM) = REAL( particle(m)%var(iM) + particle(m)%ivar(iA) * xn, KIND=lcmp )
       particle(n)%var(iM) = REAL( particle(n)%var(iM) - particle(m)%ivar(iA) * xn, KIND=lcmp )

       particle(m)%var(iN) = REAL( particle(m)%var(iN) + particle(m)%ivar(iA) * xsn, KIND=lcmp )
       particle(n)%var(iN) = REAL( particle(n)%var(iN) - particle(m)%ivar(iA) * xsn, KIND=lcmp )
      
       particle(m)%var(in_startind+1:in_startind+4) = &
       REAL (particle(m)%var(in_startind+1:in_startind+4) + particle(m)%ivar(iA) * xin, KIND=lcmp)
       particle(n)%var(in_startind+1:in_startind+4) = &
       REAL (particle(n)%var(in_startind+1:in_startind+4) - particle(m)%ivar(iA) * xin, KIND=lcmp)

       particle(n)%ivar(iA) = particle(n)%ivar(iA) - particle(m)%ivar(iA)
      !
      !--    If true, particle m has or will get IN due to this interaction
       IF (particle(n)%ivar(iTYPE) .EQ. -2 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
         particle(m)%ivar(iTYPE) = -2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
         DO intype = 1,4
            sfc_in = particle(m)%var(in_startind+intype) / REAL(particle(m)%ivar(iA))
            particle(m)%var(iTFZ) = MAX(particle(m)%var(iTFZ), t_frz(intype, sfc_in))
         END DO
      ENDIF
!
!--    Update of monomer number
       particle(m)%ivar(iNmono) = particle(m)%ivar(iNmono) + particle(n)%ivar(iNmono)
!
!--    Update of rimed water mass
       particle(m)%var(iMrime)  = REAL( particle(m)%var(iMrime) + particle(m)%ivar(iA) * xrn, KIND=lcmp )
       particle(n)%var(iMrime)  = REAL( particle(n)%var(iMrime) - particle(m)%ivar(iA) * xrn, KIND=lcmp )
!
!--    Update of axes and densities
       particle(m)%var(iDNS)  = REAL( rho_apparent, KIND=lcmp )
       particle(m)%var(iAXS)  = REAL( aaxis, KIND=lcmp )
       particle(m)%var(iCXS)  = REAL( caxis, KIND=lcmp )
!
!--    Aggregation profile
       IF ( particle(m)%ivar(iNmono) .GT. 10 )  THEN ! new aggregate is snow ( number monomers > 10 )
          IF ( particle(m)%ivar(iNmono) - particle(n)%ivar(iNmono) .LE. 10  .AND.  particle(n)%ivar(iNmono) .LE. 10 )  THEN ! n and m were ice before aggregation
             agg_prof(k) = agg_prof(k) + ( xm + xn ) * particle(m)%ivar(iA)
             micro_field(i,j,k,iagg) = micro_field(i,j,k,iagg) + ( xm + xn ) * particle(m)%ivar(iA)
          ENDIF
          IF ( particle(n)%ivar(iNmono) .LE. 10 )  THEN ! n was ice before aggregation
             agg_prof(k) = agg_prof(k) + xn * particle(m)%ivar(iA)
             micro_field(i,j,k,iagg) = micro_field(i,j,k,iagg) + xn * particle(m)%ivar(iA)
          ENDIF
          IF ( particle(m)%ivar(iNmono) - particle(n)%ivar(iNmono) .LE. 10 )  THEN ! m was ice before aggregation
             agg_prof(k) = agg_prof(k) + xm * particle(m)%ivar(iA)
             micro_field(i,j,k,iagg) = micro_field(i,j,k,iagg) + xm * particle(m)%ivar(iA)
          ENDIF
       ENDIF
!
!-- If the weighting factors are equal or the droplets of the same superdroplet collide,
!-- the resultant droplet and aerosol mass is distributed evenly on the two interacting
!-- computational particles.
    ELSEIF ( particle(n)%ivar(iA) .EQ. particle(m)%ivar(iA) )  THEN

       particle(n)%ivar(iA) = FLOOR( particle(n)%ivar(iA) / 2.0, KIND=8 )
       particle(m)%ivar(iA) = particle(m)%ivar(iA) - particle(n)%ivar(iA)

       particle(n)%var(iM) = REAL(( xn + xm ) * particle(n)%ivar(iA), KIND=lcmp)
       particle(m)%var(iM) = REAL(( xn + xm ) * particle(m)%ivar(iA), KIND=lcmp)

       particle(n)%var(iN) = REAL(( xsn + xsm ) * particle(n)%ivar(iA), KIND=lcmp)
       particle(m)%var(iN) = REAL(( xsn + xsm ) * particle(m)%ivar(iA), KIND=lcmp)
       
       particle(n)%var(in_startind+1:in_startind+4) = REAL((xin + xim) * particle(n)%ivar(iA), KIND=lcmp)
       particle(m)%var(in_startind+1:in_startind+4) = REAL((xin + xim) * particle(m)%ivar(iA), KIND=lcmp)
      !
      !--    If true, particle n and m will get IN due to this interaction
         IF ( particle(n)%ivar(iTYPE) .EQ. -2  .OR.  particle(m)%ivar(iTYPE) .EQ. -2 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
            particle(n)%ivar(iTYPE) = -2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
            particle(m)%ivar(iTYPE) = -2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
            DO intype = 1,4
               sfc_in = particle(n)%var(in_startind+intype) / REAL(particle(n)%ivar(iA))
               particle(n)%var(iTFZ) = MAX(particle(n)%var(iTFZ), t_frz(intype, sfc_in))
               particle(m)%var(iTFZ) = MAX(particle(m)%var(iTFZ), t_frz(intype, sfc_in))
            ENDDO
         ENDIF
!
!--    Aggregation profile (calculation needs to be placed here because of changes in n and m particle)
       IF ( particle(n)%ivar(iNmono) + particle(m)%ivar(iNmono) .GT. 10 )  THEN ! new aggregate is snow ( number monomers > 10 )
          IF ( particle(m)%ivar(iNmono) .LE. 10  .AND.  particle(n)%ivar(iNmono) .LE. 10 )  THEN ! n and m were ice before aggregation
             agg_prof(k) = agg_prof(k) + ( xm + xn ) * particle(m)%ivar(iA)
             micro_field(i,j,k,iagg) = micro_field(i,j,k,iagg) + ( xm + xn ) * particle(m)%ivar(iA)
          ENDIF
          IF ( particle(n)%ivar(iNmono) .LE. 10  .AND.  particle(m)%ivar(iNmono) .GT. 10 )  THEN ! n was ice before aggregation
             agg_prof(k) = agg_prof(k) + xn * particle(m)%ivar(iA)
             micro_field(i,j,k,iagg) = micro_field(i,j,k,iagg) + xn * particle(m)%ivar(iA)
          ENDIF
          IF ( particle(m)%ivar(iNmono) .LE. 10  .AND.  particle(n)%ivar(iNmono) .GT. 10 )  THEN ! m was ice before aggregation
             agg_prof(k) = agg_prof(k) + xm * particle(m)%ivar(iA)
             micro_field(i,j,k,iagg) = micro_field(i,j,k,iagg) + xm * particle(m)%ivar(iA)
          ENDIF
       ENDIF
!
!--    Update of monomer number
       particle(n)%ivar(iNmono) = particle(n)%ivar(iNmono) + particle(m)%ivar(iNmono)
       particle(m)%ivar(iNmono) = particle(n)%ivar(iNmono)
!
!--    Update of rimed water mass
       particle(m)%var(iMrime)  = REAL(( xrn + xrm ) * particle(n)%ivar(iA), KIND=lcmp)
       particle(n)%var(iMrime)  = REAL(( xrn + xrm ) * particle(m)%ivar(iA), KIND=lcmp)
!
!--    Update of axes and densities
       particle(m)%var(iDNS)  = REAL( rho_apparent, KIND=lcmp )
       particle(n)%var(iDNS)  = REAL( rho_apparent, KIND=lcmp )
       particle(m)%var(iAXS)  = REAL( aaxis, KIND=lcmp )
       particle(n)%var(iAXS)  = REAL( aaxis, KIND=lcmp )
       particle(m)%var(iCXS)  = REAL( caxis, KIND=lcmp )
       particle(n)%var(iCXS)  = REAL( caxis, KIND=lcmp )

    ENDIF

 END SUBROUTINE update_aggregation

 SUBROUTINE determine_riming(m,n,i,j,k,n_particles,check_collection,check_shattering,p_crit,nfrag_tot,nfrag_b,nfrag_s)

    USE micro_vars,   ONLY: n_start, n_end, particle, rho_liq, dV, diss_LCM,   &
                            lcmp, dtn_LCM, switch_linear_sampling_collection,  &
                            iA, iM, iN, iWsedi, iAXS, iCXS, iDNS, rho_ice,     &
                            iTYPE, iMrime, iNmono, rho_air, ds_param,          &
                            switch_ds, ds_prof
    USE micro_random, ONLY: randomu
    USE microphysics, ONLY: micro_field, ilastt, iss, isLES, imv, tabs_LCM,    &
                            pres_LCM, ids
    USE params,       ONLY: pi, muelq, ggr, fac_fus
    USE grid,         ONLY: rank
    USE vars,         ONLY: esatw, qsatw
    USE micro_ds

    IMPLICIT NONE

    REAL                 ::  An, Am, KK_agg, wsn, wsm, p_crit_lim, eta, St,    &
                             A_j, A_rime, Ace_j, E_rime, KK_rime, phi_j, r_k,  &
                             v_j, v_k, lambda, N_Re, N_St, r_j, C_SC, pp,      &
                             E_BG74, N_Fr, N_Re_columns, E_EM17_plates,        &
                             E_EM17_columns, Nfrag, dqf, tinymass, mindiff, phi
    INTEGER              ::  half_n_particles, j_int, k_int, mm
    INTEGER, INTENT(IN)  ::  m, n, i, j, k, n_particles
    LOGICAL, INTENT(OUT) ::  check_collection, check_shattering
    REAL,    INTENT(OUT) ::  p_crit
    INTEGER, INTENT(OUT) ::  nfrag_b, nfrag_s, nfrag_tot
    INTEGER              ::  p, MinN

    check_collection = .FALSE.
    check_shattering = .FALSE.
    nfrag_b = 0
    nfrag_s = 0
    nfrag_tot = 0
!
!-- Collision is only calculated if one weighting factor is larger than 1 to avoid
!-- zero weighting factors
    IF ( MAX( particle(n)%ivar(iA), particle(m)%ivar(iA) ) .GT. 1_8 ) THEN

       IF ( particle(m)%ivar(iTYPE) .LT. 0 )  THEN  ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
          j_int = m !ice
          k_int = n !liquid
       ELSE
          j_int = n
          k_int = m
       ENDIF
!
!--    Fall velocities
       v_k = particle(k_int)%var(iWsedi)
       v_j = particle(j_int)%var(iWsedi)
!
!--    Droplet and ice radius
       r_k   = ( particle(k_int)%var(iM) / particle(k_int)%ivar(iA) / ( 4.0 / 3.0 * pi * rho_liq ) )**0.33333333333
       r_j   = ( particle(j_int)%var(iAXS)**2 * particle(j_int)%var(iCXS))**0.33333333
!
!--    Aspect ratio of ice particle
       phi_j = particle(j_int)%var(iCXS) / particle(j_int)%var(iAXS)
!
!--    Area of circumscribed ellipse of ice particle
       Ace_j = pi * particle(j_int)%var(iAXS) * MAX( particle(j_int)%var(iAXS), particle(j_int)%var(iCXS) )
!
!--    True area of ice particle, (5) and (6) in Shima et al. (2020)
       A_j   = Ace_j * ( particle(j_int)%var(iDNS) / rho_ice )**EXP(-phi_j)
!
!--    Geometric cross section area, (38) in Shima et al. (2020)
       A_rime = pi * ( particle(j_int)%var(iAXS) + r_k ) * ( MAX( particle(j_int)%var(iAXS), particle(j_int)%var(iCXS) ) + r_k ) &
                - ( Ace_j - A_j )
!
!--    Determine riming efficiency, (39) - (41) in Shima et al. (2020)
       IF ( v_j .LT. v_k )  THEN
!
!--       Ratio of ice particle radius to droplet radius
          pp = r_j / r_k
!
!--       Reynolds number of the droplet
          N_Re = v_k * 2.0 * r_k / ( muelq / rho_air(k) )
!
!--       Cunningham slip correction factor for ice particle
          lambda = 6.62E-6 * 1.0E-2 ! mean free path, conversion from cm to m
          C_SC   = 1.0 + 1.26 * lambda / r_j
!
!--       Stokes impaction parameter
          N_St   = pp**2 * particle(j_int)%var(iDNS) * N_Re * C_SC / ( 9.0 * rho_air(k) )
!
!--       Get riming efficiency
          E_rime = E_BG74(pp,N_Re,N_St)

       ELSE
!
!--       Ratio of droplet radius to ice radius
          pp = r_k / r_j
!
!--       Reynolds number of the ice particle
          N_Re = v_j * 2.0 * r_j / ( muelq / rho_air(k) )
!
!--       Mixed Froude number
          N_Fr = ( v_j - v_k ) * v_k / ( ggr * MAX( particle(j_int)%var(iAXS), particle(j_int)%var(iCXS) ) )
!
!--       Reynolds number of ice columns
          N_Re_columns = v_j * 2.0 * particle(j_int)%var(iAXS) / ( muelq / rho_air(k) )
!
!--       Get riming efficiency
          IF ( phi_j .LE. 1.0 )  THEN
             E_rime = phi_j * E_BG74(pp,N_Re,N_Fr) + ( 1.0 - phi_j ) * E_EM17_plates(N_Re,N_Fr)
          ELSE
             E_rime = ( 1.0 / phi_j ) * E_BG74(pp,N_Re,N_Fr) + ( 1.0 - ( 1.0 / phi_j ) ) * E_EM17_columns(N_Re_columns,N_Fr)
          ENDIF

       ENDIF
!
!--    Cober and List (1993) riming efficiency
       ! N_St     = 2.0 * rho_liq * ABS( v_j - v_k ) * r_k**2 / ( 9.0 * ( muelq / rho_air(k) ) * MAX( particle(j_int)%var(iAXS), particle(j_int)%var(iCXS) ) )
       ! E_rime = MAX( 0.0, 0.55 * LOG10( 2.51 * N_St ) )
!
!--    Riming kernel, (37) in Shima et al. (2020)
       KK_rime = E_rime * A_rime * ABS( v_j - v_k )
!
!--    Probability of aggregation
       p_crit = MAX( particle(n)%ivar(iA), particle(m)%ivar(iA) ) / dV(k) * KK_rime * dtn_LCM
!
!--    Collision probability needs to be upscaled if the LSM method is used.
       IF ( switch_linear_sampling_collection )  THEN
          half_n_particles = FLOOR( n_particles / 2.0 )
          p_crit           = p_crit * n_particles * ( n_particles - 1.0 ) / ( half_n_particles * 2.0 )
       ENDIF
!
!--    Multiple collisions are not possible for riming!
       IF ( p_crit .GT. randomu(rank+1) )  THEN

!--    First, need to check for droplet shattering (which cannot occur simultaneously with riming)
          IF (switch_ds .and. (TRIM(ds_param) .EQ. "Phillips") .and. (tabs_LCM(i,j,k) .ge. 248.16) .and. (tabs_LCM(i,j,k) .le. 270.16)) THEN
             CALL ds_phillips(i,j,k,k_int,tabs_LCM(i,j,k), nfrag_tot, nfrag_b)
             nfrag_s = nfrag_tot - nfrag_b
          ELSEIF (switch_ds .and. (TRIM(ds_param) .EQ. "Lachapelle") .and. (tabs_LCM(i,j,k) .ge. 248.16) .and. (tabs_LCM(i,j,k) .le. 270.16)) THEN
            CALL ds_lachapelle(i,j,k,k_int,tabs_LCM(i,j,k),nfrag_tot)
          ELSEIF (switch_ds .and. (TRIM(ds_param) .EQ. "Lauber") .and. (tabs_LCM(i,j,k) .ge. 248.16) .and. (tabs_LCM(i,j,k) .le. 270.16)) THEN
            CALL ds_lauber(i,j,k,k_int,tabs_LCM(i,j,k),nfrag_tot)
          ENDIF
          IF (nfrag_tot .gt. 0) THEN !Freezing fragmentation occurred - no riming
             check_collection = .FALSE.
             check_shattering = .TRUE.
             p_crit = 0.0
          ELSE
            check_collection = .TRUE.  ! collections will take place
            check_shattering = .FALSE.
            p_crit = 1.0
          ENDIF
       ENDIF
    ENDIF
   END SUBROUTINE determine_riming
 SUBROUTINE update_shattering(m,n,i,j,k,nfrag_tot,nfrag_b,nfrag_s)
   USE micro_vars,   ONLY: n_start, n_end, particle, rho_liq, dV, lcmp, dtn_LCM, &
                            iA, iM, iN, iWsedi, iAXS, iCXS, iDNS, rho_ice,       &
                            iTYPE, iMrime, iNmono, rho_air, ds_param,            &
                            switch_ds, ds_prof
    USE micro_random, ONLY: randomu
    USE microphysics, ONLY: micro_field, ilastt, iss, isLES, imv, tabs_LCM,      &
                            pres_LCM, ids
    USE params,       ONLY: pi, muelq, ggr, fac_fus
    USE grid,         ONLY: rank
    USE vars,         ONLY: esatw, qsatw
    USE micro_ds

    IMPLICIT NONE

    REAL                 ::  Nfrag, dqf, tinymass, remicemass, fragradius, &
                             mindiff, phi, totfrzmass
    INTEGER              ::  j_int, k_int, mm, collnum
    INTEGER, INTENT(IN)  ::  m, n, i, j, k, nfrag_tot, nfrag_b, nfrag_s
    INTEGER              ::  p, MinN
   collnum = MIN(particle(n)%ivar(iA), particle(m)%ivar(iA))
   IF (particle(n)%ivar(iTYPE) .lt. 0) THEN
      j_int = n !ice
      k_int = m !liquid
   ELSE
      j_int = m !ice
      k_int = n !liquid
   ENDIF
   totfrzmass = 0

   IF ((nfrag_s .gt. 0) .and. (nfrag_b .gt. 0)) THEN !Fragmentation occurred with small fragments
      MinN = -1
      tinymass = pi/6.*rho_ice*(1E-6)**3
      mindiff = 1 !kg. Minimum difference in mass (per particle) between fragments and ice particles in gridbox
      !Search for other ice particles in the grid box to assign the small fragments to
      DO p = n_start(i,j,k), n_end(i,j,k) !Yes, this is a double loop, but it should only trigger seldomly
         IF ((ABS(particle(p)%var(iM)/particle(p)%ivar(iA) - tinymass) .lt. mindiff) .and. (particle(p)%ivar(iTYPE) .lt. 0) .and. &
         (particle(p)%var(iM)/particle(p)%ivar(iA)/tinymass .lt. 5) .and. (particle(p)%var(iM)/particle(p)%ivar(iA)/tinymass .gt. 0.5)) THEN !Merge with another particle in gridbox
            mindiff = ABS(particle(p)%var(iM)/particle(p)%ivar(iA) - tinymass)
            MinN = p
         ENDIF
      END DO
      IF (MinN .gt. 0) THEN !Found an appropriate particle within the gridbox which the fragments can merge to
         ! totfrzmass = totfrzmass + tinymass*nfrag_s*collnum
         phi = particle(MinN)%var(iCXS)/particle(MinN)%var(iAXS) !Set aspect ratio of fragments to the same as that of the merging particle, for simplicity
         particle(MinN)%var(iM) = particle(MinN)%var(iM) + tinymass*nfrag_s*collnum
         particle(MinN)%var(iDNS) = rho_ice !Set all particles as solid ice. This should almost certainly be the case anyway for small ice
         particle(MinN)%ivar(iA) = particle(MinN)%ivar(iA) + nfrag_s*collnum
         particle(MinN)%var(iAXS) = (3*particle(MinN)%var(iM)/(4*pi*rho_ice*particle(MinN)%ivar(iA)*phi))**(1./3.)
         particle(MinN)%var(iCXS) = particle(MinN)%var(iAXS)*phi
         !For simplicity, we'll assume that these ice fragments are pure water vapor, without any trace gases, CCN, or IN, so they don't affect the quantities of these species in the merged particle
         particle(k_int)%var(iM) = particle(k_int)%var(iM) - tinymass*nfrag_s*collnum
         !Now, need to transfer the number of water droplets which actually collided (and froze and shattered)
         !into the ice phase
         MinN = -1
         remicemass = particle(k_int)%var(iM) / (particle(k_int)%ivar(iA) * (nfrag_b+1))
         fragradius = (3.0 * remicemass / (4.0 * pi * rho_ice))**(0.3333333)
         mindiff = 1
         DO p = n_start(i,j,k), n_end(i,j,k)
            IF ((ABS(particle(p)%var(iM)/particle(p)%ivar(iA) - remicemass) .lt. mindiff) .and. (particle(p)%ivar(iTYPE) .lt. 0) .and. &
            (particle(p)%var(iM)/particle(p)%ivar(iA)/remicemass .lt. 5) .and. (particle(p)%var(iM)/particle(p)%ivar(iA)/remicemass .gt. 0.5)) THEN !Merge with another particle in gridbox
               mindiff = ABS(particle(p)%var(iM)/particle(p)%ivar(iA) - remicemass)
               MinN = p
            ENDIF
         END DO
         IF (MinN .gt. 0 .and. (collnum .lt. particle(j_int)%ivar(iA))) THEN !Found a matching particle for large fragments
            ! totfrzmass = totfrzmass + remicemass * (nfrag_b + 1) * collnum
            phi = particle(MinN)%var(iCXS)/particle(MinN)%var(iAXS)
            phi = (phi*(particle(MinN)%var(iCXS)*particle(MinN)%var(iAXS)**2)*(particle(MinN)%ivar(iA)) + 1.0*fragradius**3*(nfrag_b+1)*collnum)/((particle(MinN)%var(iCXS)*particle(MinN)%var(iAXS)**2)*(particle(MinN)%ivar(iA)) + 1.0*fragradius**3*(nfrag_b+1)*collnum)
            !Assuming the fragments have an aspect ratio of 1
            particle(MinN)%var(iDNS) = (particle(MinN)%var(iDNS)*particle(MinN)%var(iM) + rho_ice*remicemass*(nfrag_b+1)*collnum)/(particle(MinN)%var(iM)+remicemass*(nfrag_b+1)*collnum)
            particle(MinN)%var(iM) = particle(MinN)%var(iM) + remicemass*(nfrag_b+1)*collnum
            particle(MinN)%ivar(iA) = particle(MinN)%ivar(iA) + (nfrag_b+1)*collnum
            particle(MinN)%var(iCXS) = (3.*particle(MinN)%var(iM)*phi**2/(4.*pi*particle(MinN)%var(iDNS)*particle(MinN)%ivar(iA)))**(1./3.)
            particle(MinN)%var(iAXS) = particle(MinN)%var(iCXS)/phi
            particle(k_int)%var(iM) = particle(k_int)%var(iM) - remicemass * (nfrag_b+1) * collnum
            particle(k_int)%ivar(iA) = particle(k_int)%ivar(iA) - collnum
         ELSEIF (MinN .lt. 0 .and. (collnum .lt. particle(j_int)%ivar(iA))) THEN !No matching particle found - merge large ice fragments with colliding ice particle j
            phi = particle(j_int)%var(iCXS)/particle(j_int)%var(iAXS)
            phi = (phi*(particle(j_int)%var(iCXS)*particle(j_int)%var(iAXS)**2)*(particle(j_int)%ivar(iA)) + 1.0*fragradius**3*(nfrag_b+1)*collnum)/((particle(j_int)%var(iCXS)*particle(j_int)%var(iAXS)**2)*(particle(j_int)%ivar(iA)) + 1.0*fragradius**3*(nfrag_b+1)*collnum)
            !Assuming the fragments have an aspect ratio of 1
            particle(j_int)%var(iDNS) = (particle(j_int)%var(iDNS)*particle(j_int)%var(iM) + rho_ice*remicemass*(nfrag_b+1)*collnum)/(particle(j_int)%var(iM)+remicemass*(nfrag_b+1)*collnum)
            particle(j_int)%var(iM) = particle(j_int)%var(iM) + remicemass*(nfrag_b+1)*collnum
            particle(j_int)%ivar(iA) = particle(j_int)%ivar(iA) + (nfrag_b+1)*collnum
            particle(j_int)%var(iCXS) = (3.*particle(j_int)%var(iM)*phi**2/(4.*pi*particle(j_int)%var(iDNS)*particle(j_int)%ivar(iA)))**(1./3.)
            particle(j_int)%var(iAXS) = particle(j_int)%var(iCXS)/phi
            particle(k_int)%var(iM) = particle(k_int)%var(iM) - remicemass * (nfrag_b+1) * collnum
            particle(k_int)%ivar(iA) = particle(k_int)%ivar(iA) - collnum
         ELSE !Particle k has a lower multiplicity than particle j - assign all fragments to particle k
            particle(k_int)%var(iDNS) = rho_ice
            particle(k_int)%var(iAXS) = (3 * particle(k_int)%var((iM)) / (4* pi * particle(k_int)%ivar(iA) * rho_ice))**(0.33333333)
            particle(k_int)%var(iCXS) = particle(k_int)%var(iAXS)
            particle(k_int)%ivar(iA) = particle(k_int)%ivar(iA) * (nfrag_b+1)
            particle(k_int)%ivar(iTYPE) = -1*ABS(particle(k_int)%ivar(iTYPE))
         ENDIF
         !This isn't (nfrag_b+1)*collnum since we're just removing the number of particles which were affected by the collision
      ELSEIF (MinN .lt. 0 .and. (collnum .lt. particle(k_int)%ivar(iA))) THEN !No matching particle found - merge all ice fragments with colliding ice particle j
         remicemass = particle(k_int)%var(iM) / (particle(k_int)%ivar(iA) * (nfrag_tot+1))
         fragradius = (3.0 * remicemass / (4.0 * pi * rho_ice))**(0.3333333)
         phi = particle(j_int)%var(iCXS)/particle(j_int)%var(iAXS)
         phi = (phi*(particle(j_int)%var(iCXS)*particle(j_int)%var(iAXS)**2)*(particle(j_int)%ivar(iA)) + 1.0*fragradius**3*(nfrag_tot+1)*collnum)/((particle(j_int)%var(iCXS)*particle(j_int)%var(iAXS)**2)*(particle(j_int)%ivar(iA)) + 1.0*fragradius**3*(nfrag_tot+1)*collnum)
         !Assuming the fragments have an aspect ratio of 1
         particle(j_int)%var(iDNS) = (particle(j_int)%var(iDNS)*particle(j_int)%var(iM) + rho_ice*remicemass*(nfrag_tot+1)*collnum)/(particle(j_int)%var(iM)+remicemass*(nfrag_tot+1)*collnum)
         particle(j_int)%var(iM) = particle(j_int)%var(iM) + remicemass*(nfrag_tot+1)*collnum
         particle(j_int)%ivar(iA) = particle(j_int)%ivar(iA) + (nfrag_tot+1)*collnum
         particle(j_int)%var(iCXS) = (3.*particle(j_int)%var(iM)*phi**2/(4.*pi*particle(j_int)%var(iDNS)*particle(j_int)%ivar(iA)))**(1./3.)
         particle(j_int)%var(iAXS) = particle(j_int)%var(iCXS)/phi
         particle(k_int)%var(iM) = particle(k_int)%var(iM) - remicemass * (nfrag_tot+1) * collnum
         particle(k_int)%ivar(iA) = particle(k_int)%ivar(iA) - collnum
      ELSE !Particle k has a lower multiplicity than particle j - assign all fragments to particle k
         particle(k_int)%var(iDNS) = rho_ice
         particle(k_int)%var(iAXS) = (3 * particle(k_int)%var((iM)) / (4* pi * particle(k_int)%ivar(iA) * rho_ice))**(0.33333333)
         particle(k_int)%var(iCXS) = particle(k_int)%var(iAXS)
         particle(k_int)%ivar(iA) = particle(k_int)%ivar(iA) * (nfrag_tot+1)
         particle(k_int)%ivar(iTYPE) = -1*ABS(particle(k_int)%ivar(iTYPE))
      ENDIF
      
   ELSE !fragmentation, but without one size category of fragments
      MinN = -1
      remicemass = particle(k_int)%var(iM) / (particle(k_int)%ivar(iA) * (nfrag_tot+1))
      fragradius = (3.0 * remicemass / (4.0 * pi * rho_ice))**(0.3333333)
      mindiff = 1
      DO p = n_start(i,j,k), n_end(i,j,k)
         IF ((ABS(particle(p)%var(iM)/particle(p)%ivar(iA) - remicemass) .lt. mindiff) .and. (particle(p)%ivar(iTYPE) .lt. 0) .and. &
         (particle(p)%var(iM)/particle(p)%ivar(iA)/remicemass .lt. 5) .and. (particle(p)%var(iM)/particle(p)%ivar(iA)/remicemass .gt. 0.5)) THEN !Merge with another particle in gridbox
            mindiff = ABS(particle(p)%var(iM)/particle(p)%ivar(iA) - remicemass)
            MinN = p
         ENDIF
      END DO
      IF (MinN .gt. 0 .and. (collnum .lt. particle(k_int)%ivar(iA))) THEN !Found a matching particle for fragments
         ! totfrzmass = totfrzmass + remicemass * (nfrag_b + 1) * collnum
         phi = particle(MinN)%var(iCXS)/particle(MinN)%var(iAXS)
         phi = (phi*(particle(MinN)%var(iCXS)*particle(MinN)%var(iAXS)**2)*(particle(MinN)%ivar(iA)) + 1.0*fragradius**3*(nfrag_tot + 1)*collnum)/((particle(MinN)%var(iCXS)*particle(MinN)%var(iAXS)**2)*(particle(MinN)%ivar(iA)) + 1.0*fragradius**3*(nfrag_tot + 1)*collnum)
         !Assuming the fragments have an aspect ratio of 1
         particle(MinN)%var(iDNS) = (particle(MinN)%var(iDNS)*particle(MinN)%var(iM) + rho_ice*remicemass*(nfrag_tot + 1)*collnum)/(particle(MinN)%var(iM)+remicemass*(nfrag_tot + 1)*collnum)
         particle(MinN)%var(iM) = particle(MinN)%var(iM) + remicemass*(nfrag_tot + 1)*collnum
         particle(MinN)%ivar(iA) = particle(MinN)%ivar(iA) + (nfrag_tot + 1)*collnum
         particle(MinN)%var(iCXS) = (3.*particle(MinN)%var(iM)*phi**2/(4.*pi*particle(MinN)%var(iDNS)*particle(MinN)%ivar(iA)))**(1./3.)
         particle(MinN)%var(iAXS) = particle(MinN)%var(iCXS)/phi
         particle(k_int)%var(iM) = particle(k_int)%var(iM) - remicemass * (nfrag_tot+1) * collnum
         particle(k_int)%ivar(iA) = particle(k_int)%ivar(iA) - collnum
      ELSEIF (MinN .lt. 0 .and. (collnum .lt. particle(k_int)%ivar(iA))) THEN !No matching particle found - merge large ice fragments with colliding ice particle j
         phi = particle(j_int)%var(iCXS)/particle(j_int)%var(iAXS)
         phi = (phi*(particle(j_int)%var(iCXS)*particle(j_int)%var(iAXS)**2)*(particle(j_int)%ivar(iA)) + 1.0*fragradius**3*(nfrag_tot + 1)*collnum)/((particle(j_int)%var(iCXS)*particle(j_int)%var(iAXS)**2)*(particle(j_int)%ivar(iA)) + 1.0*fragradius**3*(nfrag_tot + 1)*collnum)
         !Assuming the fragments have an aspect ratio of 1
         particle(j_int)%var(iDNS) = (particle(j_int)%var(iDNS)*particle(j_int)%var(iM) + rho_ice*remicemass*(nfrag_tot + 1)*collnum)/(particle(j_int)%var(iM)+remicemass*(nfrag_tot + 1)*collnum)
         particle(j_int)%var(iM) = particle(j_int)%var(iM) + remicemass*(nfrag_tot + 1)*collnum
         particle(j_int)%ivar(iA) = particle(j_int)%ivar(iA) + (nfrag_tot + 1)*collnum
         particle(j_int)%var(iCXS) = (3.*particle(j_int)%var(iM)*phi**2/(4.*pi*particle(j_int)%var(iDNS)*particle(j_int)%ivar(iA)))**(1./3.)
         particle(j_int)%var(iAXS) = particle(j_int)%var(iCXS)/phi
         particle(k_int)%var(iM) = particle(k_int)%var(iM) - remicemass * (nfrag_tot+1) * collnum
         particle(k_int)%ivar(iA) = particle(k_int)%ivar(iA) - collnum
      ELSE !Particle k has lower multiplicity than particle j - assign all fragments to particle k
         particle(k_int)%var(iDNS) = rho_ice
         particle(k_int)%var(iAXS) = (3 * particle(k_int)%var((iM)) / (4* pi * particle(k_int)%ivar(iA) * rho_ice))**(0.33333333)
         particle(k_int)%var(iCXS) = particle(k_int)%var(iAXS)
         particle(k_int)%ivar(iA) = particle(k_int)%ivar(iA) * (nfrag_tot+1)
         particle(k_int)%ivar(iTYPE) = -1*ABS(particle(k_int)%ivar(iTYPE))
      ENDIF
   ENDIF
   ds_prof(k) = ds_prof(k) + (nfrag_tot+1)*collnum
   micro_field(i,j,k,ids) = micro_field(i,j,k,ids) + (nfrag_tot + 1)*collnum
   tabs_LCM(i,j,k)             = tabs_LCM(i,j,k) + ( fac_fus * (remicemass * (nfrag_tot + 1)*collnum) ) / ( dV(k) * rho_air(k) )
   micro_field(i,j,k,iss)      = micro_field(i,j,k,imv) - qsatw( tabs_LCM(i,j,k), pres_LCM(k) ) !absolute supersaturation
   micro_field(i,j,k,ilastT)   = tabs_LCM(i,j,k) ! temperature at last timestep
   micro_field(i,j,k,isLES)    = micro_field(i,j,k,iss) / qsatw( tabs_LCM(i,j,k), pres_LCM(k) ) !relative supersaturation

 END SUBROUTINE update_shattering

 SUBROUTINE update_riming(n,m,i,j,k,p_crit)

    USE micro_vars,   ONLY: n_start, n_end, particle, rho_liq, dV, diss_LCM,   &
                            lcmp, dtn_LCM, frz_prof, r_separate_micro, iMrime, &
                            iNmono, iA, iM, iN, iTYPE, iAXS, iCXS, iDNS,  &
                            rho_ice, rim_prof, tabs_LCM, rho_air, pres_lcm,    &
                            iWsedi, iTFZ, in_startind, t_frz
    USE micro_random, ONLY: randomu
    USE microphysics, ONLY: micro_field, imv, isLES, iss, irim, ilastT
    USE params,       ONLY: pi, fac_fus, muelq, therco, diffelq, &
                            variable_therco, variable_diffelq
    USE grid,         ONLY: rank

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n, m, i, j, k
    INTEGER             :: j_int, k_int, mm, intype
    REAL, INTENT(IN)    :: p_crit

    REAL :: rho_apparent, rho_max, rho_min, rho_crit, V_j, V_k, V_max, x_j,    &
            x_k, xrm, xrn, aaxis, caxis, xn, xm, xsn, xsm, dqf,                &
            phi_j, r_k, rho_rime, qsatw, N_Re, N_St, r_j, T_surface,           &
            f_v_impact, T_surface_local, v_impact, rho_rime_HP85, sfc_in
    REAL, DIMENSION(4) :: xim, xin
!
!-- Check which particles are ice and liquid
    IF ( particle(m)%ivar(iTYPE) .LT. 0 )  THEN  ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
       j_int = m !ice
       k_int = n !liquid
    ELSE
       j_int = n
       k_int = m
    ENDIF
!
!-- Droplet and ice radius
    r_k   = ( particle(k_int)%var(iM) / particle(k_int)%ivar(iA) / ( 4.0 / 3.0 * pi * rho_liq ) )**0.33333333333
    r_j   = ( particle(j_int)%var(iAXS)**2 * particle(j_int)%var(iCXS))**0.33333333
!
!-- Droplet and ice radius
    x_k = particle(k_int)%var(iM) / particle(k_int)%ivar(iA)
    x_j = particle(j_int)%var(iM) / particle(j_int)%ivar(iA)
!
!-- Aspect ratio of ice particle
    phi_j = particle(j_int)%var(iCXS) / particle(j_int)%var(iAXS)
!
!-- Volumes of particles (before (26) in Shima et al. (2020))
    V_j = 4.0 / 3.0 * pi * particle(j_int)%var(iAXS)**2 * particle(j_int)%var(iCXS)
!
!-- Update for ice smaller than liquid particle
    IF ( MAX( particle(j_int)%var(iCXS), particle(j_int)%var(iAXS) ) .LT. r_k )  THEN
!
!--    New apparent density of rimed particle, (42) in Shima et al. (2020)
       rho_apparent = rho_ice
!
!--    Axes of rimed particle, (43) in Shima et al. (2020)
       aaxis = ( ( x_k + x_j ) / ( 4.0 * pi / 3.0 * rho_ice ) )**0.33333333
       caxis = aaxis

    ELSE

       N_Re = particle(j_int)%var(iWsedi) * 2.0 * r_j / ( muelq / rho_air(k) )
       N_St = ( r_k / r_j )**2 * rho_liq * N_Re / ( 9.0 * rho_air(k) )

       v_impact        = ABS(particle(j_int)%var(iWsedi)-particle(k_int)%var(iWsedi))*f_v_impact(N_Re,N_St)
       IF (variable_diffelq .and. variable_therco) THEN
         T_surface_local = T_surface(tabs_LCM(i,j,k),micro_field(i,j,k,imv),pres_LCM(k),therco(i,j,k),diffelq(i,j,k))
       ELSEIF (variable_diffelq) THEN
         T_surface_local = T_surface(tabs_LCM(i,j,k),micro_field(i,j,k,imv),pres_LCM(k),therco(1,1,1),diffelq(i,j,k))
       ELSEIF (variable_therco) THEN
         T_surface_local = T_surface(tabs_LCM(i,j,k),micro_field(i,j,k,imv),pres_LCM(k),therco(i,j,k),diffelq(1,1,1))
       ELSE
         T_surface_local = T_surface(tabs_LCM(i,j,k),micro_field(i,j,k,imv),pres_LCM(k),therco(1,1,1),diffelq(1,1,1))
       ENDIF
       rho_rime        = rho_rime_HP85(T_surface_local,v_impact,r_k)

       IF ( phi_j .LE. 0.8  .OR.  ( phi_j .GT. 1.0  .AND.  phi_j .LE. 1.25 ) )  THEN
!
!--       New apparent density of rimed particle, (49) in Shima et al. (2020)
          rho_apparent = ( x_j + x_k ) / ( V_j + x_k / rho_rime )
!
!--       A axis of rimed particle, (108) in Shima et al. (2020)
          aaxis = MAX( particle(j_int)%var(iAXS), r_k * ( rho_liq / rho_rime )**0.33333333 )
!
!--       C axis of rimed particle, (51) in Shima et al. (2020)
          caxis = ( V_j + x_k / rho_rime ) / ( 4.0 * pi / 3.0 * particle(j_int)%var(iAXS)**2 )

       ELSE
!
!--       New apparent density of rimed particle, (52) in Shima et al. (2020)
          rho_apparent = ( x_j + x_k ) / ( V_j + x_k / rho_rime )
!
!--       C axis of rimed particle, (109) in Shima et al. (2020)
          caxis = MAX( particle(j_int)%var(iCXS), r_k * ( rho_liq / rho_rime )**0.33333333 )
!
!--       A axis of rimed particle, (53) in Shima et al. (2020)
          aaxis = SQRT( ( V_j + x_k / rho_rime ) / ( 4.0 * pi / 3.0 * particle(j_int)%var(iCXS) ) )

       ENDIF

    ENDIF
!
!-- Compute mass of superdroplets n and m
    xn = particle(n)%var(iM) / particle(n)%ivar(iA)
    xm = particle(m)%var(iM) / particle(m)%ivar(iA)
!
!-- Compute aerosol mass of superdroplet n and m
    xsn = particle(n)%var(iN) / particle(n)%ivar(iA)
    xsm = particle(m)%var(iN) / particle(m)%ivar(iA)
!
!-- Compute IN surface area of superdroplets n and m
    xin = particle(n)%var(in_startind+1:in_startind+4) / particle(n)%ivar(iA)
    xim = particle(m)%var(in_startind+1:in_startind+4) / particle(m)%ivar(iA)
!
!-- Compute rimed water mass of superdroplet n and m
    xrn = particle(n)%var(iMrime) / particle(n)%ivar(iA)
    xrm = particle(m)%var(iMrime) / particle(m)%ivar(iA)
!
!-- If n = j, ice particle (n) collects liquid droplet (m): number monomers equals number monomers of n; rime mass is sum of rime mass of n plus water mass of m
!-- If m = j, liquid droplet (n) collects ice particle (m): number monomers equals number monomers of m; rime mass is sum of water mass of n plus rime mass of m
    IF ( particle(n)%ivar(iA) .LT. particle(m)%ivar(iA) )  THEN

       particle(n)%var(iM)  = REAL( particle(n)%var(iM) + particle(n)%ivar(iA) * xm, KIND=lcmp )
       particle(m)%var(iM)  = REAL( particle(m)%var(iM) - particle(n)%ivar(iA) * xm, KIND=lcmp )

       particle(n)%var(iN)  = REAL( particle(n)%var(iN) + particle(n)%ivar(iA) * xsm, KIND=lcmp )
       particle(m)%var(iN)  = REAL( particle(m)%var(iN) - particle(n)%ivar(iA) * xsm, KIND=lcmp )
      
       particle(n)%var(in_startind+1:in_startind+4) = &
       REAL( particle(n)%var(in_startind+1:in_startind+4) + particle(n)%ivar(iA) * xim * p_crit, KIND=lcmp)
       particle(m)%var(in_startind+1:in_startind+4)  = &
       REAL( particle(m)%var(in_startind+1:in_startind+4) - particle(n)%ivar(iA) * xim * p_crit, KIND=lcmp)

       particle(m)%ivar(iA) = particle(m)%ivar(iA) - particle(n)%ivar(iA) * p_crit
      !
      !--    If true, particle n has or will get IN due to this interaction
      IF (ABS(particle(m)%ivar(iTYPE)) .EQ. 2 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
         particle(n)%ivar(iTYPE) = -2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
         DO intype = 1,4
            sfc_in = particle(n)%var(in_startind+intype) / REAL(particle(n)%ivar(iA))
            particle(n)%var(iTFZ) = MAX(particle(n)%var(iTFZ), t_frz(intype, sfc_in))
         END DO
      ELSE
         particle(n)%ivar(iTYPE) = -1*ABS(particle(n)%ivar(iTYPE))
      ENDIF
!--    Update of monomer number
       IF ( n .EQ. j_int )  THEN
          particle(n)%ivar(iNmono) = particle(n)%ivar(iNmono)
       ELSE
          particle(n)%ivar(iNmono) = particle(m)%ivar(iNmono)
       ENDIF
!
!--    Update of rimed water mass
       IF ( n .EQ. j_int )  THEN
          particle(n)%var(iMrime)  = REAL( particle(n)%var(iMrime) + particle(n)%ivar(iA) * xm, KIND=lcmp ) != particle(n)%ivar(iA) * ( xrn + xm )
          dqf = particle(n)%ivar(iA) * xm
       ELSE
          particle(n)%var(iMrime)  = REAL( particle(n)%ivar(iA) * xn + particle(n)%ivar(iA) * xrm, KIND=lcmp ) != particle(n)%ivar(iA) * ( xn + xrm )
          particle(m)%var(iMrime)  = REAL( particle(m)%var(iMrime) - particle(n)%ivar(iA) * xrm, KIND=lcmp )
          dqf = particle(n)%ivar(iA) * xn
       ENDIF
!
!--    Update of axes and densities
       particle(n)%var(iDNS)  = REAL( rho_apparent, KIND=lcmp )
       particle(n)%var(iAXS)  = REAL( aaxis, KIND=lcmp )
       particle(n)%var(iCXS)  = REAL( caxis, KIND=lcmp )

    ELSEIF ( particle(m)%ivar(iA) .LT. particle(n)%ivar(iA) )  THEN

       particle(m)%var(iM) = REAL( particle(m)%var(iM) + particle(m)%ivar(iA) * xn, KIND=lcmp )
       particle(n)%var(iM) = REAL( particle(n)%var(iM) - particle(m)%ivar(iA) * xn, KIND=lcmp )

       particle(m)%var(iN) = REAL( particle(m)%var(iN) + particle(m)%ivar(iA) * xsn, KIND=lcmp )
       particle(n)%var(iN) = REAL( particle(n)%var(iN) - particle(m)%ivar(iA) * xsn, KIND=lcmp )
      
       particle(m)%var(in_startind+1:in_startind+4) = &
       REAL (particle(m)%var(in_startind+1:in_startind+4) + particle(m)%ivar(iA) * xin * p_crit, KIND=lcmp)
       particle(n)%var(in_startind+1:in_startind+4) = &
       REAL (particle(n)%var(in_startind+1:in_startind+4) - particle(m)%ivar(iA) * xin * p_crit, KIND=lcmp)

       particle(n)%ivar(iA) = particle(n)%ivar(iA) - particle(m)%ivar(iA) * p_crit
      !
      !--    If true, particle m has or will get IN due to this interaction
       IF (ABS(particle(n)%ivar(iTYPE)) .EQ. 2 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
         particle(m)%ivar(iTYPE) = -2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
         DO intype = 1,4
            sfc_in = particle(m)%var(in_startind+intype) / REAL(particle(m)%ivar(iA))
            particle(m)%var(iTFZ) = MAX(particle(m)%var(iTFZ), t_frz(intype, sfc_in))
         END DO
      ELSE
         particle(m)%ivar(iTYPE) = -1*ABS(particle(m)%ivar(iTYPE))
      ENDIF
!
!--    Update of monomer number
       IF ( m .EQ. j_int )  THEN
          particle(m)%ivar(iNmono) = particle(m)%ivar(iNmono)
       ELSE
          particle(m)%ivar(iNmono) = particle(n)%ivar(iNmono)
       ENDIF
!
!--    Update of rimed water mass
       IF ( m .EQ. j_int )  THEN
          particle(m)%var(iMrime)  = REAL( particle(m)%var(iMrime) + particle(m)%ivar(iA) * xn, KIND=lcmp ) != particle(m)%ivar(iA) * ( xrm + xn )
          dqf = particle(m)%ivar(iA) * xn
       ELSE
          particle(m)%var(iMrime)  = REAL( particle(m)%ivar(iA) * xm + particle(m)%ivar(iA) * xrn, KIND=lcmp ) != particle(m)%ivar(iA) * ( xm + xrn )
          particle(n)%var(iMrime)  = REAL( particle(n)%var(iMrime) - particle(m)%ivar(iA) * xrn, KIND=lcmp )
          dqf = particle(m)%ivar(iA) * xm
       ENDIF
!
!--    Update of axes and densities
       particle(m)%var(iDNS)  = REAL( rho_apparent, KIND=lcmp )
       particle(m)%var(iAXS)  = REAL( aaxis, KIND=lcmp )
       particle(m)%var(iCXS)  = REAL( caxis, KIND=lcmp )
!
!-- If the weighting factors are equal or the droplets of the same superdroplet collide,
!-- the resultant droplet and aerosol mass is distributed evenly on the two interacting
!-- computational particles.
    ELSEIF ( particle(n)%ivar(iA) .EQ. particle(m)%ivar(iA) )  THEN

       particle(n)%ivar(iA) = FLOOR( particle(n)%ivar(iA) / 2.0, KIND=8 )
       particle(m)%ivar(iA) = particle(m)%ivar(iA) - particle(n)%ivar(iA)

       particle(n)%var(iM) = ( xn + xm ) * particle(n)%ivar(iA)
       particle(m)%var(iM) = ( xn + xm ) * particle(m)%ivar(iA)

       particle(n)%var(iN) = ( xsn + xsm ) * particle(n)%ivar(iA)
       particle(m)%var(iN) = ( xsn + xsm ) * particle(m)%ivar(iA)
       
       particle(n)%var(in_startind+1:in_startind+4) = (xin + xim) * particle(n)%ivar(iA)
       particle(m)%var(in_startind+1:in_startind+4) = (xin + xim) * particle(m)%ivar(iA)
      !
      !--    If true, particle n and m will get IN due to this interaction
         IF ( ABS(particle(n)%ivar(iTYPE)) .EQ. 2  .OR.  ABS(particle(m)%ivar(iTYPE)) .EQ. 2 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
            particle(n)%ivar(iTYPE) = -2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
            particle(m)%ivar(iTYPE) = -2 ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
            DO intype = 1,4
               sfc_in = particle(n)%var(in_startind+intype) / REAL(particle(n)%ivar(iA))
               particle(n)%var(iTFZ) = MAX(particle(n)%var(iTFZ), t_frz(intype, sfc_in))
               particle(m)%var(iTFZ) = MAX(particle(m)%var(iTFZ), t_frz(intype, sfc_in))
            ENDDO
         ELSE
         !Both particles homogeneously nuclated
            particle(n)%ivar(iTYPE) = -1
            particle(m)%ivar(iTYPE) = -1
         ENDIF
!
!--    Update of monomer number
       IF ( n .EQ. j_int )  THEN
          particle(n)%ivar(iNmono) = particle(n)%ivar(iNmono)
          particle(m)%ivar(iNmono) = particle(n)%ivar(iNmono)
       ELSE
          particle(n)%ivar(iNmono) = particle(m)%ivar(iNmono)
          particle(m)%ivar(iNmono) = particle(m)%ivar(iNmono)
       ENDIF
!
!--    Update of rimed water mass
       IF ( n .EQ. j_int )  THEN
          particle(m)%var(iMrime)  = ( xrn + xm ) * particle(n)%ivar(iA)
          particle(n)%var(iMrime)  = ( xrn + xm ) * particle(m)%ivar(iA)
          dqf = xm * ( particle(m)%ivar(iA) + particle(n)%ivar(iA) )
       ELSE
          particle(m)%var(iMrime)  = ( xrm + xn ) * particle(n)%ivar(iA)
          particle(n)%var(iMrime)  = ( xrm + xn ) * particle(m)%ivar(iA)
          dqf = xn * ( particle(m)%ivar(iA) + particle(n)%ivar(iA) )
       ENDIF
!
!--    Update of axes and densities
       particle(m)%var(iDNS)  = REAL( rho_apparent, KIND=lcmp )
       particle(n)%var(iDNS)  = REAL( rho_apparent, KIND=lcmp )
       particle(m)%var(iAXS)  = REAL( aaxis, KIND=lcmp )
       particle(n)%var(iAXS)  = REAL( aaxis, KIND=lcmp )
       particle(m)%var(iCXS)  = REAL( caxis, KIND=lcmp )
       particle(n)%var(iCXS)  = REAL( caxis, KIND=lcmp )

    ENDIF
!
!--    Riming and freezing profile
       rim_prof(k) = rim_prof(k) + dqf
       frz_prof(k) = frz_prof(k) + dqf
       micro_field(i,j,k,irim) = micro_field(i,j,k,irim) + dqf
!
!-- Compute change of temperature due to freezing
    tabs_LCM(i,j,k)             = tabs_LCM(i,j,k) + ( fac_fus * dqf ) / ( dV(k) * rho_air(k) )
    micro_field(i,j,k,iss)      = micro_field(i,j,k,imv) - qsatw( tabs_LCM(i,j,k), pres_LCM(k) ) !absolute supersaturation
    micro_field(i,j,k,ilastT)   = tabs_LCM(i,j,k) ! temperature at last timestep

    micro_field(i,j,k,isLES)    = micro_field(i,j,k,iss) / qsatw( tabs_LCM(i,j,k), pres_LCM(k) ) !relative supersaturation

 END SUBROUTINE update_riming

 REAL FUNCTION T_surface(T_int,q_int,p_int,therco,diffelq)
!
!-- Estimate surface temperature of rimed ice crystal (no gas kinetic effects)
    USE params, ONLY: lsub, rv, rgas

    IMPLICIT NONE

    REAL, INTENT(IN)  :: T_int, q_int, p_int,therco,diffelq
    REAL              :: esati_pre, esati, e_a, supersat_ice, C_ice_pre, drho
!
!-- Get ambient supersaturation with respect to ice
    esati_pre    = esati( T_int )
    e_a          = q_int * p_int / ( q_int + rgas / rv )
    supersat_ice = e_a / esati_pre - 1.0
!
!-- Prefactor in diffusional growth equation
    C_ice_pre  = 1.0 / ( rv * T_int / ( 100.0 * esati_pre * diffelq ) +  &
                         ( lsub / ( rv * T_int ) - 1.0 ) * lsub /        &
                         ( therco * T_int )                              &
                      )
!
!-- Excess water vapor density abover surface (without gas kinetic effects)
    drho = C_ice_pre * supersat_ice / diffelq ! because dmdt = 4 * pi * C_pre * diff_coeff * f_vent * drho = 4 * pi * C_pre * G_pre * f_vent * supersat_ice <==> drho = G_pre * supersat_ice / diff_coeff
!
!-- Surface temperature
    T_surface = T_int + lsub * diffelq / therco * drho
    RETURN

 END FUNCTION T_surface

 REAL FUNCTION f_v_impact(N_Re, N_St)
!
!-- Calculate impact velocity of liquid droplets on ice crystal (riming) following Rasmussen and Heymsfield (1985, JAS)
    IMPLICIT NONE

    REAL, INTENT(IN) :: N_Re, N_St
    REAL             :: f_v_impact_300, f_v_impact_100,        &
                        f_v_impact_30, f_v_impact_10, WW


    WW = LOG10(N_St)

    IF ( N_Re .GT. 300.0 )  THEN

       f_v_impact = 0.356  + 0.4738 * WW - 0.1233  * WW**2 - 0.1618 * WW**3 + 0.08087 * WW**4
       IF ( N_St .GT. 10.0 )  f_v_impact = 0.63
       IF ( N_St .LT. 0.1 )   f_v_impact = 0.0

    ELSEIF ( N_Re .GT. 100.0  .AND.  N_Re .LE. 300.0 )  THEN

       f_v_impact_300 = 0.356  + 0.4738 * WW - 0.1233  * WW**2 - 0.1618 * WW**3 + 0.08087 * WW**4
       IF ( N_St .GT. 10.0 )  f_v_impact_300 = 0.63
       IF ( N_St .LT. 0.1 )   f_v_impact_300 = 0.0

       f_v_impact_100 = 0.3272 + 0.4907 * WW - 0.09452 * WW**2 - 0.1906 * WW**3 + 0.07105 * WW**4
       IF ( N_St .GT. 10.0 )  f_v_impact_100 = 0.61
       IF ( N_St .LT. 0.1 )   f_v_impact_100 = 0.0

       f_v_impact = f_v_impact_100 + ( f_v_impact_300 - f_v_impact_100 ) / 200.0  * ( N_Re - 100.0 )

    ELSEIF ( N_Re .GT. 30.0  .AND.  N_Re .LE. 100.0 )  THEN

       f_v_impact_100 = 0.3272 + 0.4907 * WW - 0.09452 * WW**2 - 0.1906 * WW**3 + 0.07105 * WW**4
       IF ( N_St .GT. 10.0 )  f_v_impact_100 = 0.61
       IF ( N_St .LT. 0.1 )   f_v_impact_100 = 0.0

       f_v_impact_30  = 0.2927 + 0.5085 * WW - 0.03453 * WW**2 - 0.2184 * WW**3 + 0.03595 * WW**4
       IF ( N_St .GT. 10.0 )  f_v_impact_30 = 0.59
       IF ( N_St .LT. 0.1 )   f_v_impact_30 = 0.0

       f_v_impact = f_v_impact_30 + ( f_v_impact_100 - f_v_impact_30 ) / 70.0  * ( N_Re - 30.0 )

    ELSEIF ( N_Re .GT. 10.0  .AND.  N_Re .LE. 30.0 )  THEN

       f_v_impact_30 = 0.3272 + 0.4907 * WW - 0.09452 * WW**2 - 0.1906 * WW**3 + 0.07105 * WW**4
       IF ( N_St .GT. 10.0 )  f_v_impact_30 = 0.59
       IF ( N_St .LT. 0.1 )   f_v_impact_30 = 0.0

       f_v_impact_10 = 0.1701 + 0.7246 * WW + 0.2257  * WW**2 - 1.13   * WW**3 + 0.5756  * WW**4
       IF ( N_St .GT. 10.0 )  f_v_impact_10 = 0.57
       IF ( N_St .LT. 0.4 )   f_v_impact_10 = 0.0

       f_v_impact = f_v_impact_10 + ( f_v_impact_30 - f_v_impact_10 ) / 20.0  * ( N_Re - 10.0 )

    ELSE ! ELSEIF ( N_Re .LE. 10.0 )  THEN

       f_v_impact = 0.1701 + 0.7246 * WW + 0.2257  * WW**2 - 1.13   * WW**3 + 0.5756  * WW**4
       IF ( N_St .GT. 10.0 )  f_v_impact = 0.57
       IF ( N_St .LT. 0.4 )   f_v_impact = 0.0

    ENDIF

    f_v_impact = MAX( f_v_impact, 0.0 )

    RETURN

 END FUNCTION f_v_impact

 REAL FUNCTION rho_rime_HP85(T_surface_int,v_impact_int,radius_int)
!
!-- Rime density following Heymsfield and Pflaum (1985, JAS), Eqns. (13), (15), and (16)
!-- T_parcel, q_parcel, p_env: must be provided
    IMPLICIT NONE

    REAL, INTENT(IN) :: T_surface_int, v_impact_int, radius_int
    REAL             :: AA
!
!-- Parameter describes impact of rimed particles on density: the larger the droplets, the higher the density; the higher the impact velocity, the higher the density
    AA  = ( -1.0 ) * ( radius_int * 1.0E6 ) * v_impact_int / ( T_surface_int - 273.15 )   ! in [ µm * m/s / °C ] (13)

    IF ( T_surface_int - 273.15 .GT. -5.0  .OR.  AA .GT. 1.60 )  THEN
       AA = MIN( AA, 3.5 )  ! limiter as suggested by Shima et al. (2020) to avoid too small densities for temperatures close to zero, also supported by measurements presented in Pflaum and Pruppacher (1979)
       rho_rime_HP85 = EXP( -0.03115 - 1.7030 * AA + 0.9116 * AA**2 - 0.1224 * AA**3 ) ! (16)
    ELSE
       rho_rime_HP85 = 0.30 * AA**0.44  ! (15)  -> see Fig. 5
    ENDIF
    rho_rime_HP85 = rho_rime_HP85 * 1000.0 ! from g/cm^3 to kg/m^3
    rho_rime_HP85 = MAX( MIN( rho_rime_HP85, 910.0 ), 100.0 ) ! limit values as suggested by Heymsfield and Pflaum (1985, JAS)

    RETURN

 END FUNCTION rho_rime_HP85

 REAL FUNCTION E_BG74(pp,N_Re,N_St)
!
!-- Collision efficiencies by Beard and Grover (1974, JAS)
    USE params, ONLY: pi

    IMPLICIT NONE

    REAL, INTENT(IN) :: pp, N_Re, N_St
    REAL             :: FF, GG, N_St0, ZZ, HH, y_c0
!
!-- Characteristic Stokes number
    FF    = LOG(N_Re)
    GG    = -0.1007 - 0.358 * FF + 0.0261 * FF**2
    N_St0 = EXP(GG)
!
!-- Collision efficiency
    ZZ     = LOG(N_St/N_St0)
    HH     = 0.1465 + 1.302 * ZZ - 0.607 * ZZ**2 + 0.293 * ZZ**3
    y_c0   = MAX( 2.0 / pi * ATAN(HH), 0.0 )
    E_BG74 = ( y_c0 + pp )**2 / ( 1.0 + pp )**2

    RETURN

 END FUNCTION E_BG74

 REAL FUNCTION E_EM17_plates(N_Re_in,N_Fr)
!
!-- Riming efficiencies for plates by Erfani and Mitchell (2017, ACP)
    IMPLICIT NONE

    REAL, INTENT(IN) :: N_Re_in, N_Fr
    REAL             :: N_Fr_crit, E_EM17, N_Fr_thres, N_Re
!
!-- Bound Reynolds number
    N_Re = MIN( MAX( N_Re_in, 2.0 ), 120.0 )
!
!-- Critical mixed Froude number
    IF ( N_Re .GT. 1.0  .AND.  N_Re .LE. 10.0 )  THEN
       N_Fr_crit = 1.250 * N_Re**(-0.350)
    ELSEIF ( N_Re .GT. 10.0  .AND.  N_Re .LE. 40.0 )  THEN
       N_Fr_crit = 1.072 * N_Re**(-0.301)
    ELSEIF ( N_Re .GT. 40.0  .AND.  N_Re .LE. 120.0 )  THEN
       N_Fr_crit = 0.356 * N_Re**(-0.003)
    ENDIF
!
!-- Threshold mixed Froude number
    N_Fr_thres = -5.07E-10 * N_Re**5 + 1.73E-7 * N_Re**4 - 2.17E-5 * N_Re**3 + 0.0013 * N_Re**2 - 0.037 * N_Re + 0.8355
!
!-- Collision efficiency
    IF ( N_Fr .GE. 0.01  .AND.  N_Fr .LE. 0.35  .AND.  N_Re .GT. 2.0  .AND.  N_Re .LE. 120.0 )  THEN
       E_EM17 = ( 0.787 * N_Fr**0.988 ) * ( 0.263 * LOG( N_Re ) - 0.264 )
    ELSEIF ( N_Fr .GE. 0.35  .AND.  N_Fr .LE. N_Fr_thres  .AND.  N_Re .GT. 2.0  .AND.  N_Re .LE. 120.0 )  THEN
       E_EM17 = ( 0.7475 * LOG10( N_Fr ) + 0.620 ) * ( 0.263 * LOG( N_Re ) - 0.264 )
    ELSEIF ( N_Fr .GT. N_Fr_thres  .AND.  N_Fr .LT. 35.0  .AND.  N_Re .GE. 1.0  .AND.  N_Re .LE. 120.0 )  THEN
       E_EM17 = SQRT( MAX( 0.0, 1.0 - 1.0 / 5.0 * ( LOG10( N_Fr / N_Fr_crit ) - SQRT(5.0) )**2 ) )
    ELSE
       E_EM17 = 0.0
    ENDIF

    E_EM17_plates = E_EM17

    RETURN

 END FUNCTION E_EM17_plates

 REAL FUNCTION E_EM17_columns(N_Re_in, N_Fr)
!
!-- Riming efficiencies for columns by Erfani and Mitchell (2017, ACP)
    IMPLICIT NONE

    REAL, INTENT(IN) :: N_Re_in, N_Fr
    REAL             :: N_Fr_crit, N_Fr_thres, rrr, E_EM17, N_Re
!
!-- Bound Reynolds number
    N_Re = MIN( MAX( N_Re_in, 0.2 ), 20.0 )
!
!-- Critical mixed Froude number
    IF ( N_Re .GE. 0.2  .AND.  N_Re .LE. 1.7 )  THEN
       N_Fr_crit = 0.7797 * N_Re**(-0.009)
    ELSEIF ( N_Re .GT. 1.7  .AND.  N_Re .LE. 20.0 )  THEN
       N_Fr_crit = 1.0916 * N_Re**(-0.635)
    ENDIF
!
!-- Threshold mixed Froude number
    IF ( N_Re .GE. 0.2  .AND.  N_Re .LE. 2.0 )  THEN
       N_Fr_thres = 0.0251 * N_Re**2 - 0.0144 * N_Re + 0.811
    ELSEIF ( N_Re .GT. 2.0  .AND.  N_Re .LE. 20.0 )  THEN
       N_Fr_thres = -0.0003 * N_Re**3 + 0.0124 * N_Re**2 - 0.1634 * N_Re + 1.0075
    ENDIF
!
!-- r parameter related to the major radius of the ellipse (typo in Eqn. (19) -> switch lines)
    IF ( N_Re .GE. 0.2  .AND.  N_Re .LE. 1.7 )  THEN
       rrr = 0.7422 * N_Re**0.2111
    ELSEIF ( N_Re .GT. 1.7  .AND.  N_Re .LE. 20.0 )  THEN
       rrr = 0.8025 * N_Re**0.0604
    ENDIF
!
!-- Collision efficiency
    IF ( N_Fr .GE. 0.01  .AND.  N_Fr .LE. N_Fr_thres  .AND.  N_Re .GE. 0.2  .AND.  N_Re .LE. 3.0 )  THEN
       E_EM17 = ( 0.787 * N_Fr**0.988 ) * ( -0.0121 * N_Re**2 + 0.1297 * N_Re + 0.0598 )
    ELSEIF ( N_Fr .GE. 0.01  .AND.  N_Fr .LE. N_Fr_thres  .AND.  N_Re .GT. 3.0  .AND.  N_Re .LE. 20.0 )  THEN
       E_EM17 = ( 0.787 * N_Fr**0.988 ) * ( -0.0005 * N_Re**2 + 0.1028 * N_Re + 0.0359 )
    ELSEIF ( N_Fr .GT. N_Fr_thres  .AND.  N_Fr .LT. 20.0  .AND.  N_Re .GE. 0.2  .AND.  N_Re .LE. 20.0 )  THEN
       E_EM17 = rrr * SQRT( MAX( 0.0, 1.0 - 1.0 / 3.5 * ( LOG10( N_Fr / N_Fr_crit ) - SQRT(3.5) )**2 ) )
    ELSE
       E_EM17 = 0.0
    ENDIF

    E_EM17_columns = E_EM17

    RETURN

 END FUNCTION E_EM17_columns

 SUBROUTINE coll_energy(Ds, Dl, Vs, Vl, CKE, St, Sc, W1, W2, DelSsig)
   !Calculate several variables related to particle collision which are used by the breakup calculation routine
   !The definition and calculation of these variables comes from Low and List, 1982: 
   !Collision, Coalescence and Breakup of Raindrops. Parts I and II (referred to below as LL82)
   USE params, ONLY: pi
   USE micro_vars, ONLY: rho_liq
   IMPLICIT NONE
   REAL, INTENT(IN) :: Ds, Dl, Vs, Vl
   REAL, INTENT(OUT) :: CKE, St, Sc, W1, W2, DelSsig
   REAL              :: sigma = 0.0728 !Water surface tension, J/m^2
   !Ds is the small particle diameter, in m
   !Dl is the large particle diameter, in m
   !Vs is the small particle terminal velocity, in m/s
   !Vl is the large particle terminal velocity, in m/s
   !CKE is the collision kinetic energy, in J - Eq. 3.1 in LL82 Part I
   !St is the total surface energy of the two particles - Eq. 4.2 in LL82 Part I
   !Sc is the coalesced surface energy of the two particles - Eq. 4.3 in LL82 Part I
   !W1 is the Weber number, the ratio of collision kinetic energy to the coalesced surface energy (unitless) - Eq. 4.33 in LL82 Part II
   !W2 is an alternative Weber number, instead being the ratio of collision kinetic energy to 
   !the total surface energy (of the two separate particles) - Eq. 4.42 in LL82 Part II
   !DelSsig is the difference between the Total surface energy and the Coalesced surface energy - Eq. 4.1 in LL82 Part I
   !Sigma is the water surface tension constant of 7.82 * 10**(-2) J/m^2 - LL82 Part I
   !Consider using a temperature-dependent sigma here instead - B. Ascher 30.05.2025

   CKE = (pi*rho_liq / 12) * (Ds**3 * Dl**3) / (Ds**3 + Dl**3)* (Vl - Vs)**2
   St = pi * sigma * (Ds**2 + Dl**2)
   Sc = pi * sigma * (Ds**3 + Dl**3)**(2./3.)
   DelSsig = St - Sc
   W1 = CKE / Sc
   W2 = CKE / St

 END SUBROUTINE coll_energy

 REAL FUNCTION ecoal_seifert(Ds, Dl, Vs, Vl)
 !Calculate the coalescence efficiency according to Eq. (1) of Seifert et al. (2005)
 USE params, ONLY: pi
 REAL, INTENT(IN) :: Ds, Dl, Vs, Vl
 REAL, PARAMETER  :: sigma = 76.1-0.155*(293.15-273.15)
 REAL             :: N_w, U, gg, p, xx, rs
   rs = 50.*Ds
   p  = Ds/Dl
   U = 100.*abs(Vs-Vl)
   N_w = rs * U**2 / sigma
   gg  = 2**(3./2.)/(6.*pi) * p**4 * (1.+ p) / ((1.+p**2)*(1.+p**3))
   xx  = N_w**(0.5) * gg
   ecoal_seifert = max(0.767 - 10.14 * xx,0.0)
   RETURN
 END FUnCTION ecoal_seifert


 SUBROUTINE breakup_fractions(Ds, Dl, CKE, St, Sc, W1, W2, DelSsig, Ff, Fs, Fd, Rf, Rs, Rd)
 !Calculate which fraction of breakups are filament, sheet, or disk breakups as defined in Low and List, 1982 Part I
 !Also calculate the number of fragments released by each type of breakup
   IMPLICIT NONE
   REAL, INTENT(IN) :: Ds, Dl, CKE, St, Sc, W1, W2, DelSsig
   REAL, INTENT(OUT) :: Ff, Fs, Fd, Rf, Rs, Rd
   REAL :: Dlcm, Dscm, Ff1, Ff2, Ds0, astar, bstar

   !Fraction of filament breakups - Eqs 4.40, 4.41 in LL82 Part I
   IF (CKE .GE. 0.893E-6) THEN
      Rf = 1.11E-4 * CKE**(-0.654)
   ELSE
      Rf = 1.0
   ENDIF

   !Fraction of sheet breakups - Eqs 4.43, 4.44 in LL82 Part I
   IF (W2 .ge. 0.86) THEN
      Rs = 0.685 * (1 - EXP(-1.63 * (W2 - 0.86)))
   ELSE
      Rs = 0
   ENDIF

   !Fraction of disk breakups - Eqs 4.45, 4.46 in LL82 Part I
   Rd = 1.0 - MIN(Rf + Rs, 1.0)

   !Convert droplet diameters to centimeters for use in fragment number calculations
   Dlcm = Dl * 100.
   Dscm = Ds * 100.

   Ff1 = MAX((-2.25E4*(Dlcm - 0.403)**2 - 37.9)*(Dscm)**2.5 + 9.67*(Dlcm + 0.170)**2 + 4.95, 2.0) !LL82 Part I Eq. 3.3
   !Number of fragments from a filament collision if the small droplet size is above some critical value
   astar = 1.02E4
   bstar = 2.83
   Ff2 = astar * Dscm**bstar + 2 !LL82 Part I Eq. 3.4
   !Number of fragments from a filament collision if the small droplet size is below some critical value
   Ds0 = ((Ff1-2.0)/astar)**(1/bstar) !Find the critical small droplet size - LL82 Part I Eq. 3.5

   IF (Dscm .ge. Ds0) THEN
      Ff = Ff1
   ELSE
      Ff = Ff2
   ENDIF

   Fs = MAX(5.0*ERF((St - 2.53E-6)/1.85E-6)+6, 2.0) !Number of fragments from a sheet collision - LL82 Part I Eq. 3.7

   Fd = MAX(297.5 + 23.7*LOG(CKE), 2.0) !Number of fragments from a disk collision - LL82 Part I Eq. 3.9

END SUBROUTINE breakup_fractions

SUBROUTINE breakup_mcfarquhar(Ds0, Dl0, Vs, Vl, CKE, St, Sc, W1, W2, DelSsig, Ff, Fs, Fd, Rf, Rs, Rd, Mults, Multf, Ms1, Ml1, Mf)
   !Calculate weighting factors and masses of the small, large, and new (if necessary) particles based on McFarquhar (2004)
   !Referred to below as M2004
   USE micro_vars, ONLY: rho_liq, particle
   USE micro_random, ONLY: randomu
   use params, ONLY: pi
   USE grid, ONLY: rank
   IMPLICIT NONE
   ! INTEGER, INTENT(IN) :: sind, lind !sind are small, large particle indices, respectively. May be unnecessary in this function
   REAL, INTENT(IN)    :: Ds0, Dl0, Vs, Vl, CKE, St, Sc, W1, W2, DelSsig, Ff, Fs, Fd, Rf, Rs, Rd
   INTEGER(KIND=8), INTENT(OUT):: Mults, Multf
   REAL, INTENT(OUT)   :: Ms1, Ml1, Mf
   REAL                :: Dscm, Dlcm, Ds, Dl, Df,          &
                          Multl, Dscm_new, Dlcm_new, Dfcm, &
                          mulnf, sigmalnf, hlnf,           &
                          mulns, sigmalns, hlns,           &
                          mulnd, sigmalnd, hlnd,           &
                          musf, sigmasf, hsf,              &
                          muls, sigmals, hls,              &
                          muld, sigmald, hld,              &
                          Dcoal3
   REAL                :: x_rand
   REAL, EXTERNAL      :: ERFI
   INTEGER             :: nancheckcount

   !First, figure out which type of collision we're dealing with
      Ds = Ds0
      Dl = Dl0

      Dscm = Ds * 100. !Droplet diameter needs to be in centimeters, because M2004 uses cm
      Dlcm = Dl * 100. !Ditto

      ! CALL coll_energy(Ds, Dl, Vs, Vl, CKE, St, Sc, W1, W2, DelSsig)
      ! CALL breakup_fractions(Ds, Dl, CKE, St, Sc, W1, W2, DelSsig, Ff, Fs, Fd, Rf, Rs, Rd)
      nancheckcount = 0
      IF (randomu(rank+1) .LE. Rf) THEN !Do filament collisions
         mulnf = -3.884 + 1.15*(1 - EXP(-6.60E6*DelSsig)) !M2004 Eq. 17
         sigmalnf = 0.478 - 0.105*(W1 -1.546)**2 !M2004 Eq. 18
         sigmalnf = MAX(sigmalnf, 0.186)
         !Actually, we don't need to calculate Hln, do we, since we already have Ff, Fs, Fd?
         musf = Dscm !M2004 Eq. 19
         Hsf = 54.1*W1**(-0.851) + 10.4 !M2004 Eq. 20
         sigmasf = 1.0/(sqrt(2.*pi)*Hsf) !M2004 Eq. 11

         Dlcm_new = -1.0
         DO WHILE(Dlcm_new .LE. 0.0)
            x_rand = randomu(rank+1)
            Dfcm = EXP((sqrt(2.0)*sigmalnf)*MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )+mulnf)
            nancheckcount = nancheckcount + 1
            IF (nancheckcount .GT. 10) Dfcm = Dfcm * 0.9**(nancheckcount-10)
            !After 10 iterations, begin artificially reducing fragment size so that we have a nonzero large droplet size
            x_rand = randomu(rank+1)
            Dscm_new = (sqrt(2.0)*sigmasf)*MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )+musf
            Dcoal3 = (Dscm**3 + Dlcm**3)
            Dlcm_new = (Dcoal3 - MIN(Dfcm**3*(INT(Ff-2.0)) + Dscm_new**3, Dcoal3))**(1.0/3.0)
         END DO

         Mults = 1
         Multl = 1 !Do we even need to keep track of this? Maybe not. Multl is always 1
         Multf = INT(Ff-2.0)



         !At the end of the filament collison, have 3 superdroplets
         !1. The only slightly changed large droplets - same number as were involved in collision
         !2. The slightly changed small droplets - same number as were involved in collision
         !3. The fragments - some number of droplets probably between 2 and 10
         !So, if we want to return the large droplets to their original superdroplet, there shouldn't be a change in multiplicity
         !compared to before the collision. We'll just need to adjust the mass (mean radius) of those droplets slightly. Ditto
         !for the small droplets
         !For the filaments, though, we should see if we can find particles within the gridbox with similar masses
         !However, Fabian noted that, while this approach is probably better for droplet size distributions, it's bad for
         !Aerosol, and we might be better served following De Jong et al. (2023)'s approach, where they combined all fragments
         !together into a single superparticle with an average radius. This approach overestimates the fragment sizes, but is better
         !for aerosol conservation. Might be best for now to set a user-controlled flag which determines which method to use.
         !'B.Ascher 02.06.2025

      ELSEIF (randomu(rank+1) .LE. Rs) THEN !Do sheet collisions
         mulns = -3.65 + 0.51*W1
         sigmalns = 0.279 + 1.91E4*CKE
         !sigmals = -0.00612 + 8.65E-4*(Dlcm/Dscm) + 3.84E4*DelSsig !This actually doesn't matter, since we're
         !Going to just consider the large droplet mass as "whatever is left over"
         Dlcm_new = -1.0
         DO WHILE(Dlcm_new .LE. 0.0)
            x_rand = randomu(rank+1)
            Dfcm = EXP((sqrt(2.0)*sigmalns)*MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )+mulns)
            nancheckcount = nancheckcount + 1
            IF (nancheckcount .GT. 10) Dfcm = Dfcm * 0.9**(nancheckcount-10)
            !After 10 iterations, begin artificially reducing fragment size so that we have a nonzero large droplet size
            Dcoal3 = (Dscm**3 + Dlcm**3)
            Dlcm_new = (Dcoal3 - MIN(MAX(INT(Fs-1.0), 1)*Dfcm**3, Dcoal3))**(1./3.)
            Dscm_new = 0.0 !No new small droplets; since we need to output their size, just set it to zero
         END DO

         Mults = 0 !Small droplet completely destroyed
         Multf = MAX(INT(Fs-1.0), 1) !Number of fragments is one large drop plus a bunch of small droplets
         !Need to have at least one 'small' fragment to prevent divide by zero errors
         Multl = 1
      ELSE !Do disk collisions
         mulnd = -2.905 + 1.224 * (Dscm / Dlcm)
         sigmalnd = 0.0049 + 2.131 * Dlcm - 2.986 * Dscm
         Dlcm_new = -1.0 !Just consider the large droplet mass as "whatever is left over" from the collision
         DO WHILE(Dlcm_new .LE. 0.0)
            x_rand = randomu(rank+1)
            Dfcm = EXP((sqrt(2.0)*sigmalnd)*MIN( MAX( ERFI( 2.0 * x_rand - 1.0 ), -2.0 ), 2.0 )+mulnd)
            nancheckcount = nancheckcount + 1
            IF (nancheckcount .GT. 10) Dfcm = Dfcm * 0.9**(nancheckcount-10)
            !After 10 iterations, begin artificially reducing fragment size so that we have a nonzero large droplet size
            Dcoal3 = (Dscm**3 + Dlcm**3)
            Dlcm_new = (Dcoal3 - MIN(MAX(INT(Fs-1.0), 1)*Dfcm**3, Dcoal3))**(1./3.)
            Dscm_new = 0.0 !No new small droplets; since we need to output their size, just set it to zero
         END DO
         
         Mults = 0 !Small droplet completely destroyed by collision
         Multf = MAX(INT(Fs-1.0), 1)!Number of fragments is one large drop plus a bunch of small droplets
         !Need to have at least one 'small' fragment to prevent divide by zero errors
         Multl = 1
      ENDIF

      Ds = Dscm_new / 100.
      Dl = Dlcm_new / 100.
      Df = Dfcm / 100.

      Ms1 = pi/6. * rho_liq * (Ds)**3
      Ml1 = pi/6. * rho_liq * (Dl)**3
      Mf = pi/6. * rho_liq * (Df)**3

   END SUBROUTINE breakup_mcfarquhar












