#ifdef SAMCHEM
 SUBROUTINE micro_condensation(i,j,k, chem_new_radii)
#else
 SUBROUTINE micro_condensation(i,j,k)
#endif /*SAMCHEM*/
!
!-- Calculates condensational growth
    USE micro_vars,   ONLY: n_start, n_end, particle, rho_liq, dV, iT, iq, iN, &
                            iTFZ, switch_LEM_micro, lcmp, rho_aero,            &
                            molecular_weight_water, vanthoff_aero, sigma,      &
                            pres_LCM, rho_air, dtn_LCM,                        &
                            cnd_prof, evp_prof,                                &
                            switch_analytical_micro, act_prof, dea_prof,       &
                            ieta, iZ, r_act, qnmin_LCM, tabs_LCM,              &
                            switch_koehler_micro, switch_FFC_micro,            &
                            k_ice, A_ice, B_ice, iAXS, N_Avogadro,             &
                            iCXS, iDNS, iTYPE, rho_ice, imm_prof, hzn_prof,    &
                            hom_prof, nuc_prof, nml_prof, switch_habit_ice,    &
                            sub_prof, dep_prof, frz_prof,                      &
                            mel_prof, nsub_prof, molecular_weight_aero,        &
                            switch_radiation_micro, dz_LCM, switch_ice,        &
                            start_ice, iMrime, iA, iM, iWsedi, iID, iNmono,    &
                            nv, boltzk, planckh, av, bv, dw, srfthick,         &
                            habitparam, get_bfactor,   &
                            get_r_aero, switch_ds, ds_param, ds_prof, np_xy,   &
                            in_startind, n_ice, frac_ice, unit_aero, iAero_portions, &
                            n_aero, nstep_cut_supersat, supersat_threshold
    USE micro_random, ONLY: randomu
    USE microphysics, ONLY: micro_field, isLES, imv, iss, iDEA, iACT,inuc, ids,   &
                            iimm, ihzn, ihom, ifrz, imel, icon, ievp, idep, isub, &
#ifdef SAMCHEM
                            ilastT, ilastp, ipH, idelR, iSbyO3t, iSbyO3m, iSbyH2O2t, iSbyH2O2m
#else
                            ilastT, ilastp
#endif /*SAMCHEM*/
    USE params,       ONLY: pi, rv, lcond, cp, fac_cond, ggr, &
                            muelq, lsub, fac_sub, fac_fus, &
                            diffelq, therco, variable_diffelq, variable_therco
    USE vars,         ONLY: esatw, qsatw, dtqsatw, esati, qsati, dtqsati
    USE grid,         ONLY: rank, nstep, z, masterproc
    USE rad,          ONLY: qrad

#ifdef SAMCHEM
    USE chemistry,    ONLY: chem_proc, flag_chem_init
    USE params,       ONLY: dochemistry
#endif /*SAMCHEM*/

    USE micro_ds
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: i, j, k

    INTEGER             :: n
    INTEGER             :: nliq, nice, m
    REAL                :: venteff, new_radius, eta_dyn, eta_lst, eta_new_tmp, &
                           eta_int_tmp, eta_new, eta_mean, tau, fsw, eta_int,  &
                           forcing, qsatw_pre, dql, D_pre, esatw_pre,          &
                           radiation, dV_n, dql_n, qrad_liq, k_liq, k_tot,     &
                           aaxis, caxis, dm, A_ice_pre, A_liq_pre, C_ice_pre,  &
                           mass, p_nucleate, C_liq_pre, rho_apparent, dsd_ice, &
                           dsd_liq, t_nucleate, w_sedi, capacitance, dqi_n,    &
                           dqi, dqf, dqf_n, esati_pre, qsati_pre, pre_ice,     &
                           pre_liq, ql_local, supersat_dummy, therco_local, diffelq_local
                           
    INTEGER             :: iequiv, pren_nuc
    REAL                :: slopelg,aresid,cresid,uresid,nup1,nup2,ndwn1,ndwn2, &
                           mfp, vapvel, scrit_a, scrit_c, &
                           gtp1, gtp2, gtp_stand, S_a, S_c, rhoa, alphat, lt, &
                           requiv, alpha_c, alpha_a
    !-- Parameters for Diskice formulation (Harrington et al. (2019))
                           
    real                :: jv, jt, homofrzprob

    !nv is 3.35*10**(-28), m**-3, the number of water molecules per cubic meter in liquid water
    !boltzk is 1.38*10**(-23), J/k, Boltzmann's constant; this is already in non-si units in the radiation part of the code, but it's not much extra memory to put here
    !planckh is 6.626*10**(-34), J*s, Planck's constant; this is already in non-si units in the radiation part of the code, but it's not much extra memory to put here
    !av is -2.323*10**(-18), J, fit parameter for homogeneous freezing equation from Kuhn et. al. (2011)
    !bv is -1.075*10**(-20), J/kg, fit parameter for homogeneous freezing equation from Kuhn et. al. (2011)
    !dw is -2.13*10**(-20), J, some change in free energy as the droplet freezes
    !jv is the volume-specific homogeneous freezing rate for a single drop; can be interpreted as the number of freezing events expected
    !to happen per cubic meter of water in the droplet per second
    !jt is the total homogeneous freezing rate for a single drop; can be interpreted as the number of freezing events expected to happen
    !within the drop per second. This will increase as the cube of droplet radius
    !s is the water droplet's "surface thickness" - the distance in which surface freezing can happen
    !homofrzprob is the probability of the liquid water drop freezing homogeneously in a given timestep

    REAL                :: afactor, bfactor, alpha, A_droplet, N_FFC

    REAL, PARAMETER     :: alpha_mordy = 0.036 ! Mordy (1959)
    REAL, PARAMETER     :: molecular_weight_FFC = 243.0E-3 ! molecular weight of FFCs (kg/mol)
    REAL, PARAMETER     :: f_FFC = 0.05 ! mass fraction of FFC
    REAL, PARAMETER     :: min_mass_aero = 4.0 / 3.0 * pi * 1800.0 * (0.1E-9)**3 ! mass of minmum aerosol (for heterogeneous ice nucleation)
    REAL, PARAMETER     :: mslope = 10 !Unitless coefficient determining the response of alpha_a and alpha_c to supersat_ice/scrit. 10 comes from Harrington et al. (2019)



    REAL, DIMENSION(n_start(i,j,k):n_end(i,j,k)) :: f_vent_mass_liq, radius,   &
                                                    r_aero, Q_a,               &
                                                    particle_tabs, r0,         &
                                                    supersat_liq,              &
                                                    supersat_ice,              &
                                                    ssi,                       & !Almost identical to supersat_ice, but based on vapor pressure rather than mixing ratio
                                                    f_vent_mass_ice,           &
                                                    capacitance_ice,           &
                                                    f_vent_axes_ice

#ifdef SAMCHEM
    REAL, DIMENSION(n_start(i,j,k):n_end(i,j,k)) :: chem_new_radii
    LOGICAL, DIMENSION(n_start(i,j,k):n_end(i,j,k)) :: did_chem
#endif /*SAMCHEM*/

    INTEGER, DIMENSION(:), allocatable :: liqind, iceind
    INTEGER :: liqpos, icepos
!-- Variables for droplet shattering
    INTEGER :: nfrag_tot=0, nfrag_b=0, nfrag_s=0, MinN=-1, p=0
    REAL    :: phi, tinymass, mindiff

    nliq = count(particle(n_start(i,j,k):n_end(i,j,k))%ivar(iTYPE) .GT. 0)
    nice = count(particle(n_start(i,j,k):n_end(i,j,k))%ivar(iTYPE) .LT. 0)
    liqpos = 1; icepos = 1
    allocate(liqind(nliq))
    allocate(iceind(nice))
    do n = n_start(i,j,k),n_end(i,j,k)
      if (particle(n)%ivar(iTYPE) .gt. 0) THEN
         liqind(liqpos) = n
         liqpos = liqpos + 1
      elseif ( particle(n)%ivar(iTYPE) .lt. 0) THEN
         iceind(icepos) = n
         icepos = icepos + 1
      end if
    end do
    liqpos = 1; icepos = 1
    CALL t_startf ('micro_condensation')
!
!-- Initialize some variables
    dql = 0.0
    dqi = 0.0
    dqf = 0.0
    IF (trim(habitparam) .eq. "DISKICE") THEN
      call getscrit(tabs_LCM(i,j,k), scrit_a, scrit_c)
    ENDIF
!
!-- Parameters determining the diffusion of water vapor and the conduction of heat (C_pre)
!-- psychrometric correction (A_pre), factor for radiation consideration (D_pre), curvature (afactor) and
!-- solute effect (bfactor) for Köhler theory
    IF (variable_therco) THEN
      therco_local = therco(i,j,k)
    ELSE
      therco_local = therco(1,1,1)
    ENDIF
    IF (variable_diffelq) THEN
      diffelq_local = diffelq(i,j,k)
    ELSE
      diffelq_local = diffelq(1,1,1)
    ENDIF
    qsatw_pre  = qsatw( tabs_LCM(i,j,k), pres_LCM(k) )
    esatw_pre  = esatw( tabs_LCM(i,j,k) )
    C_liq_pre  = 1.0 / ( rv * tabs_LCM(i,j,k) / ( 100.0 * esatw_pre * diffelq_local ) +  &
                         ( lcond / ( rv * tabs_LCM(i,j,k) ) - 1.0 ) * lcond /      &
                         ( therco_local * tabs_LCM(i,j,k) )                              &
                       )
    A_liq_pre  = 1.0 + dtqsatw( tabs_LCM(i,j,k), pres_LCM(k) ) * lcond / cp
    D_pre      = 1.0 / ( therco_local * tabs_LCM(i,j,k) ) * ( lcond / ( rv * tabs_LCM(i,j,k) ) - 1.0 )

    qsati_pre  = qsati( tabs_LCM(i,j,k), pres_LCM(k) )
    esati_pre  = esati( tabs_LCM(i,j,k) )
    C_ice_pre  = 1.0 / ( rv * tabs_LCM(i,j,k) / ( 100.0 * esati_pre * diffelq_local ) +  &
                         ( lsub / ( rv * tabs_LCM(i,j,k) ) - 1.0 ) * lsub /        &
                         ( therco_local * tabs_LCM(i,j,k) )                              &
                      )
    A_ice_pre  = 1.0 + dtqsati( tabs_LCM(i,j,k), pres_LCM(k) ) * lsub  / cp

    IF ( switch_koehler_micro )  THEN
       afactor = 2.0 * sigma / ( rho_liq * rv * tabs_LCM(i,j,k) )
       ! bfactor is calculated later depending on particle composition
    ELSE
       afactor = 0.0
       bfactor = 0.0
    ENDIF

    ql_local = 0.0 ! local liquid water mixing ratio
    DO n = n_start(i,j,k), n_end(i,j,k)

       r_aero(n)        = get_r_aero(n)
       particle_tabs(n) = particle(n)%var(iT) - particle(n)%var(iZ) * ggr / cp ! only for switch_LEM_micro

       IF ( particle(n)%ivar(iTYPE) .GT. 0 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
          radius(n)          = ( particle(n)%var(iM) / ( 4.0 / 3.0 * pi * particle(n)%ivar(iA) * rho_liq ) )**0.33333333333
          !radius(n)          = ( (particle(n)%var(iM)/rho_liq + particle(n)%var(iN)/rho_aero) / ( 4.0 / 3.0 * pi * particle(n)%ivar(iA) ) )**0.33333333333
          f_vent_mass_liq(n) = venteff(radius(n),particle(n)%var(iWsedi),rho_air(k))
          ql_local           = ql_local + particle(n)%var(iM) / ( dV(k) * rho_air(k) )
       ELSEIF ( particle(n)%ivar(iTYPE) .LT. 0 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
          aaxis              = particle(n)%var(iAXS)
          caxis              = particle(n)%var(iCXS)
          w_sedi             = particle(n)%var(iWsedi)
          capacitance_ice(n) = capacitance(aaxis,caxis)
          CALL ice_ventilation(f_vent_mass_ice(n),f_vent_axes_ice(n),aaxis,caxis,w_sedi,capacitance_ice(n),rho_air(k),diffelq_local)
       ENDIF

    ENDDO

    IF ( switch_radiation_micro )  THEN
       radiation = 0.0
       DO n = n_start(i,j,k), n_end(i,j,k)
!
!--       Calculate absoption coefficient Q_a using equation (8) from Roach (1976, QJRMS), and heating rate for each droplet
          Q_a(n) = 1.18 * ( 1.0 - EXP( -0.28E6 * radius(n) ) )
!
!--       Integral droplet surface
          radiation = radiation + 4.0 * pi * radius(n)**2 * Q_a(n) * particle(n)%ivar(iA) / dV(k)

       ENDDO
!
!--    Determine liquid contribution to qrad, (10) in Klinger et al. (2019)
       k_tot    = 1.0 ! to be adapted
       k_liq    = 1.0 ! 0.090361 * SUM( particle(n_start(i,j,k):n_end(i,j,k))%var(iM) ) / dV(k) * dz_LCM * 1000.0 ! based on CAM model (https://www2.cesm.ucar.edu/models/atm-cam/docs/description/node23.html#SECTION00895000000000000000)
       qrad_liq = qrad(i,j,k) * k_liq / k_tot
!
!--    Heating rate ( qrad(i,j,k) ) is converted from [K/s] to [W/m^3]. Then, it is distributed among all droplets considering their radius and absorption coefficient, (12) in Klinger et al. (2019)
       radiation = qrad_liq * cp * rho_air(k) / radiation
!
!--    If there is no water, set radiation to zero
       IF ( SUM( particle(n_start(i,j,k):n_end(i,j,k))%var(iM) ) / dV(k) / rho_air(k) .LT. 1.0E-5 ) THEN
          radiation = 0.0
       ENDIF

    ELSE
       Q_a       = 0.0
       radiation = 0.0
    ENDIF
!
!-- Lengthscale for gas-kinetic effects (see Kogan, 1991)
    IF ( switch_FFC_micro )  THEN
!
!--    Film-forming compounds (Feingold and Chuang, 2002)
       DO n = n_start(i,j,k), n_end(i,j,k)
          A_droplet = 4.0 * pi * radius(n)**2
          N_FFC     = f_FFC * ( particle(n)%var(iN) / particle(n)%ivar(iA) ) / molecular_weight_FFC * N_Avogadro
          If ( A_droplet / N_FFC .LT. 22.0E-20 )  THEN
             alpha = 1.0E-5
          ELSE
             alpha = alpha_mordy
          ENDIF
          r0(n) = diffelq_local / alpha * SQRT( 2.0 * pi / ( rv * tabs_LCM(i,j,k) ) ) / &
               ( 1.0 + diffelq_local * lcond**2 * 100.0 * esatw_pre / ( therco_local * rv**2 * tabs_LCM(i,j,k)**3 ) )
       ENDDO
    ELSE
       alpha = alpha_mordy
       r0(:) = diffelq_local / alpha * SQRT( 2.0 * pi / ( rv * tabs_LCM(i,j,k) ) ) / &
               ( 1.0 + diffelq_local * lcond**2 * 100.0 * esatw_pre / ( therco_local * rv**2 * tabs_LCM(i,j,k)**3 ) )
    ! r0 = 1.86E-6 ! use a constant value from Grabowski et al. (2017, GMDD), valid for about 10°C
    ENDIF
!
!-- Determine average supersaturation during timestep
    IF ( .NOT. switch_LEM_micro )  THEN

       dsd_ice = 0.0
       dsd_liq = 0.0
       pre_ice = A_ice_pre * 4.0 * pi * C_ice_pre / ( qsati_pre * dV(k) )
       pre_liq = A_liq_pre * 4.0 * pi * C_liq_pre / ( qsatw_pre * dV(k) )
       DO n = n_start(i,j,k), n_end(i,j,k)
          IF ( particle(n)%ivar(iTYPE) .GT. 0 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
             dsd_liq = dsd_liq + particle(n)%ivar(iA) * radius(n)**2 / ( radius(n) + r0(n) ) * f_vent_mass_liq(n)
          ELSE
             dsd_ice = dsd_ice + particle(n)%ivar(iA) * capacitance_ice(n) * f_vent_mass_ice(n)
          ENDIF
       ENDDO
       tau     = 1.0 / ( pre_liq * dsd_liq + pre_ice * dsd_ice )
       eta_dyn = micro_field(i,j,k,imv) - qsatw( tabs_LCM(i,j,k)          , pres_LCM(k) ) !Get absolute supersaturation after LES dynamics (eta_dyn)
       eta_lst = micro_field(i,j,k,imv) - qsatw( micro_field(i,j,k,ilastT), micro_field(i,j,k,ilastp) ) !Get the supersaturation of the last timestep (eta_lst) -> advected quantity
       forcing = ( eta_dyn - eta_lst ) / dtn_LCM - ( qsati_pre - qsatw_pre ) * pre_ice * dsd_ice !the first term represents the change of the absolute supersaturation in the LES (dynamics and radiation), the second term considers differences between ice/liquid saturation pressures
!
!--    Get mean supersaturation during timestep using the analytical/numerical method
!--    by Tzivion, Feingold, and Levin (1989, JAS). First step, get the new and integral
!--    supersaturation:
       eta_new_tmp = eta_new( eta_lst, forcing, tau, dtn_LCM )
       eta_int_tmp = eta_int( eta_lst, forcing, tau, dtn_LCM ) ! this is the supersaturation that gets depleted by the hydrometers during the timesteps. That's why it is 0 for dt = 0!!!
!
!--    Bound the integral supersaturation by its values at the beginning and end of the timestep
       fsw = MAX( MIN( 1.0, MAX( ( MIN( eta_lst, eta_new_tmp ) * dtn_LCM - eta_int_tmp ) * 1.0E35, 0.0 ) ), &
                  MIN( 1.0, MAX( ( eta_int_tmp - MAX( eta_lst, eta_new_tmp ) * dtn_LCM ) * 1.0E35, 0.0 ) ) )

       eta_mean = ( eta_int_tmp * ( 1.0 - fsw ) + fsw * ( eta_lst + eta_new_tmp ) * dtn_LCM * 0.5 ) / dtn_LCM

       supersat_liq(:) = eta_mean / qsatw_pre

       IF (nstep <= nstep_cut_supersat) THEN
         DO n = n_start(i,j,k), n_end(i,j,k)
           IF (supersat_liq(n)>supersat_threshold) supersat_liq(n) = supersat_threshold
         ENDDO
       END IF

       IF ( switch_ice )  THEN
          supersat_ice(:) = ( eta_mean + qsatw_pre - qsati_pre ) / qsati_pre
          ssi(:) = ((eta_mean + qsatw_pre) * pres_LCM(k) / (0.622 + (eta_mean + qsatw_pre)) - esati_pre)/esati_pre
       ENDIF


    ELSE

       DO n = n_start(i,j,k), n_end(i,j,k)

          dV_n    = dV(k) / ( n_end(i,j,k) - n_start(i,j,k) + 1 )

          dsd_ice = 0.0
          dsd_liq = 0.0
          pre_ice = A_ice_pre * 4.0 * pi * C_ice_pre / ( qsati_pre * dV_n )
          pre_liq = A_liq_pre * 4.0 * pi * C_liq_pre / ( qsatw_pre * dV_n )
          IF ( particle(n)%ivar(iTYPE) .GT. 0 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
             dsd_liq = particle(n)%ivar(iA) * radius(n)**2 / ( radius(n) + r0(n) ) * f_vent_mass_liq(n)
          ELSE
             dsd_ice = particle(n)%ivar(iA) * capacitance_ice(n) * f_vent_mass_ice(n)
          ENDIF
          tau     = 1.0 / ( pre_liq * dsd_liq + pre_ice * dsd_ice )
          eta_dyn = particle(n)%var(iq) - qsatw( particle_tabs(n), pres_LCM(k) )
          eta_lst = particle(n)%var(ieta)
          forcing = ( eta_dyn - eta_lst ) / dtn_LCM - ( qsati( particle_tabs(n), pres_LCM(k) ) - qsatw( particle_tabs(n), pres_LCM(k) ) ) * pre_ice * dsd_ice  !the first term represents the change of the absolute supersaturation in the LES (dynamics and radiation), the second term combines effects on droplet growth by droplet curvature, aerosol mass, and radiation effects
!
!--       Get mean supersaturation during timestep using the analytical/numerical method
!--       by Tzivion, Feingold, and Levin (1989, JAS). First step, get the new and integral
!--       supersaturation:
          eta_new_tmp = eta_new( eta_lst, forcing, tau, dtn_LCM )
          eta_int_tmp = eta_int( eta_lst, forcing, tau, dtn_LCM ) ! this is the supersaturation that gets depleted by the hydrometers during the timesteps. That's why it is 0 for dt = 0!!!
!
!--       Bound the integral supersaturation by its values at the beginning and end of the timestep
          fsw = MAX( MIN( 1.0, MAX( ( MIN( eta_lst, eta_new_tmp ) * dtn_LCM - eta_int_tmp ) * 1.0E35, 0.0 ) ), &
                     MIN( 1.0, MAX( ( eta_int_tmp - MAX( eta_lst, eta_new_tmp ) * dtn_LCM ) * 1.0E35, 0.0 ) ) )

          eta_mean = ( eta_int_tmp * ( 1.0 - fsw ) + fsw * ( eta_lst + eta_new_tmp ) * dtn_LCM * 0.5 ) / dtn_LCM

          supersat_liq(n) = eta_mean / qsatw( particle_tabs(n), pres_LCM(k) )
          IF ( switch_ice )  THEN
             supersat_ice(n) = ( eta_mean + qsatw( particle_tabs(n), pres_LCM(k) ) - qsati( particle_tabs(n), pres_LCM(k) ) ) / qsati( particle_tabs(n), pres_LCM(k) )
             ssi(n) = ((eta_mean + qsatw(particle_tabs(n), pres_LCM(k))) * pres_LCM(k) / (0.622 + (eta_mean + qsatw(particle_tabs(n), pres_LCM(k)))) - esati(particle_tabs(n)))/esati(particle_tabs(n))
          ENDIF

       ENDDO
    ENDIF
    IF (.not. switch_LEM_micro) THEN
       call diffparams(tabs_LCM(i,j,k), pres_LCM(k), therco_local, diffelq_local, esati_pre, ssi(n_start(i,j,k)), mfp, &
                vapvel, scrit_a, scrit_c, gtp1, gtp2, gtp_stand, S_a, S_c, rhoa, lt)
    ENDIF
    call fittingparams(mslope,slopelg,aresid,cresid,uresid,nup1,nup2,ndwn1,ndwn2)

#ifdef SAMCHEM
    ! call chemistry procedure and get new radii from it instead radius_liquid_euler
    did_chem = .FALSE.
    micro_field(i,j,k,ipH)   = 0.0
    micro_field(i,j,k,idelR) = 0.0
    IF ( dochemistry .AND. flag_chem_init .AND. k>10 ) THEN
      ! NOTE: doesn't work if LEM is on!
      CALL chem_proc( i,j,k, &
                    & dtn_LCM,                                                          &
                    & micro_field(i,j,k,imv),                                           &
                    & micro_field(i,j,k,ilastT),                                        &
                    & (tabs_LCM(i,j,k) - micro_field(i,j,k,ilastT))/dtn_LCM,            &
                    & micro_field(i,j,k,ilastp)*100.0,                                  &
                    & (pres_LCM(k)     - micro_field(i,j,k,ilastp))*100.0/dtn_LCM,      &
                    & chem_new_radii,                                                   &
                    & did_chem,                                                         &
                    & f_vent_mass_liq,                                                  &
                    & micro_field(i,j,k,ipH),                                           &
                    & micro_field(i,j,k,idelR),                                         &
                    & micro_field(i,j,k,iSbyO3t),                                       &
                    & micro_field(i,j,k,iSbyO3m),                                       &
                    & micro_field(i,j,k,iSbyH2O2t),                                     &
                    & micro_field(i,j,k,iSbyH2O2m)                                      )
    END IF
#endif /*SAMCHEM*/

!
!-- Now, start the integration of the radius growth equation
    DO m = 1, nliq
   !  DO n = n_start(i,j,k), n_end(i,j,k) !This will restrict the loop to only liquid particles based on the partitioning above. No need for checking below
      n = liqind(m)
      ! IF (particle(n)%ivar(ITYPE) .gt. 0) THEN
!
#ifdef SAMCHEM
         IF ((.NOT. dochemistry) .OR. (.NOT. flag_chem_init) .OR. (.NOT. did_chem(n))) THEN
#endif /*SAMCHEM*/
           IF ( .NOT. switch_analytical_micro )  THEN
              new_radius = radius(n)
              ! calculate bfactor (solute term, depends on aerosol composition)
              IF ( switch_koehler_micro ) bfactor = get_bfactor(n)
              CALL radius_liquid_euler(new_radius,dtn_LCM, &
                                       r0(n),C_liq_pre,supersat_liq(n),f_vent_mass_liq(n), &
                                       afactor,bfactor,r_aero(n),D_pre,Q_a(n)*radiation)
           ELSE
              new_radius = MAX( ( radius(n) + r0(n) )**2 + 2.0 * dtn_LCM * C_liq_pre / rho_liq *  &
                                  supersat_liq(n) * f_vent_mass_liq(n), ( r_aero(n) + r0(n) )**2 )
              new_radius = SQRT( new_radius ) - r0(n)
           ENDIF
#ifdef SAMCHEM
         ELSE
           new_radius = chem_new_radii(n)
         ENDIF
#endif /*SAMCHEM*/

         particle(n)%var(iM) = REAL( 4.0 / 3.0 * pi * rho_liq * new_radius**3, lcmp ) * particle(n)%ivar(iA)

         dm = 4.0 / 3.0 * pi * rho_liq * ( new_radius**3 - radius(n)**3 ) 
         dql_n = dm * particle(n)%ivar(iA)

         dql   = dql + dql_n
         IF ( switch_LEM_micro )  THEN
            dV_n   = dV(k) / ( n_end(i,j,k) - n_start(i,j,k) + 1.0 )
            dql_n  = dql_n / ( dV_n * rho_air(k) )
            particle(n)%var(iT) = particle(n)%var(iT) + dql_n * fac_cond
            particle(n)%var(iq) = particle(n)%var(iq) - dql_n
         ENDIF
   !
   !--       Some domain-averaged analysis
         IF ( new_radius .GT. radius(n) )  THEN
            cnd_prof(k) = cnd_prof(k) + dm * particle(n)%ivar(iA) !Changed from dql_n because that can be overwritten if using LEM
            micro_field(i,j,k,icon) = micro_field(i,j,k,icon) + dm * particle(n)%ivar(iA)
         ELSE
            evp_prof(k) = evp_prof(k) + dm * particle(n)%ivar(iA)
            micro_field(i,j,k,ievp) = micro_field(i,j,k,ievp) + dm * particle(n)%ivar(iA)
         ENDIF

         IF ( new_radius .GT. r_act  .AND.  radius(n) .LE. r_act )  THEN
         act_prof(k) = act_prof(k) + particle(n)%ivar(iA)
         micro_field(i,j,k,iACT) = micro_field(i,j,k,iACT) + particle(n)%ivar(iA)
         ENDIF
         IF ( new_radius .LE. r_act  .AND.  radius(n) .GT. r_act )  THEN
         dea_prof(k) = dea_prof(k) + particle(n)%ivar(iA)
         micro_field(i,j,k,iDEA) = micro_field(i,j,k,iDEA) + particle(n)%ivar(iA)
         ENDIF
      ! ENDIF
   end do

   if (trim(habitparam) .eq. "Chen") then

      do m = 1, nice
      ! do n = n_start(i,j,k), n_end(i,j,k)
         ! IF (particle(n)%ivar(ITYPE) .lt. 0) THEN
            n = iceind(m)
            mass         = REAL(particle(n)%var(iM) / particle(n)%ivar(iA))
            aaxis        = particle(n)%var(iAXS)
            caxis        = particle(n)%var(iCXS)
            rho_apparent = particle(n)%var(iDNS)
            w_sedi       = particle(n)%var(iWsedi)

            dm = mass
            CALL mass_ice_euler(mass,dtn_LCM,C_ice_pre,capacitance_ice(n),f_vent_mass_ice(n),supersat_ice(n))
            dm = mass - dm

            particle(n)%var(iM)   = mass * particle(n)%ivar(iA)
            IF (dm .lt. 0) THEN
               particle(n)%var(iMrime) = mass/(mass - dm) * particle(n)%var(iMrime) !B. Ascher 04.11.2025 - Sublimate rimed mass as well; previously, had errors with nonsensical rime fractions!
            ENDIF
            IF (particle(n)%var(iMrime)/particle(n)%var(iM) .gt. 1.0) THEN
               dm = dm + (particle(n)%var(iM) - particle(n)%var(iMrime))
               particle(n)%var(iMrime) = particle(n)%var(iM)
            ENDIF

            IF ( switch_habit_ice )  THEN
               supersat_dummy = MIN( supersat_ice(n), esatw_pre / esati_pre - 1.0 ) ! limit ice supersaturation to ice supersaturation at water saturated conditions, Miller and Young (1979)
               CALL acd_update(aaxis,caxis,rho_apparent,mass,dm,C_ice_pre,supersat_dummy,tabs_LCM(i,j,k),diffelq_local,f_vent_axes_ice(n))
               particle(n)%var(iAXS) = aaxis
               particle(n)%var(iCXS) = caxis
               particle(n)%var(iDNS) = rho_apparent
            ELSE
               particle(n)%var(iAXS) = ( mass / ( 4.0 / 3.0 * pi * rho_ice ) )**0.333333333
               particle(n)%var(iCXS) = ( mass / ( 4.0 / 3.0 * pi * rho_ice ) )**0.333333333
               particle(n)%var(iDNS) = rho_ice
            ENDIF


            IF ( dm .GT. 0.0 )  THEN
               dep_prof(k) = dep_prof(k) + particle(n)%ivar(iA) * dm
               micro_field(i,j,k,idep) = micro_field(i,j,k,idep) + particle(n)%ivar(iA) * dm
            ELSE
               sub_prof(k) = sub_prof(k) + particle(n)%ivar(iA) * dm
               micro_field(i,j,k,isub) = micro_field(i,j,k,isub) + particle(n)%ivar(iA) * dm
            ENDIF

            dqi_n = dm * particle(n)%ivar(iA)
            dqi   = dqi + dqi_n
            IF ( switch_LEM_micro )  THEN
               dV_n   = dV(k) / ( n_end(i,j,k) - n_start(i,j,k) + 1.0 )
               dqi_n  = dqi_n / ( dV_n * rho_air(k) )
               particle(n)%var(iT) = particle(n)%var(iT) + dqi_n * fac_sub
               particle(n)%var(iq) = particle(n)%var(iq) - dqi_n
            ENDIF
         ! ENDIF
      ENDDO
   ELSEIF (trim(habitparam) .eq. "DISKICE") then
      IF (.not. switch_habit_ice) THEN
         iequiv = 1
      ENDIF
      ! DO n = n_start(i,j,k), n_end(i,j,k)
      Do m = 1, nice
         ! IF (particle(n)%ivar(iTYPE) .lt. 0) THEN
            n = iceind(m)
            mass         = REAL(particle(n)%var(iM) / particle(n)%ivar(iA))
            aaxis        = particle(n)%var(iAXS)
            caxis        = particle(n)%var(iCXS)
            rho_apparent = particle(n)%var(iDNS)
            w_sedi       = particle(n)%var(iWsedi)

            dm = mass
            IF (switch_LEM_micro) THEN
               call diffparams(tabs_LCM(i,j,k), pres_LCM(k), therco_local, diffelq_local, esati_pre, ssi(n), mfp, &
               vapvel, scrit_a, scrit_c, gtp1, gtp2, gtp_stand, S_a, S_c, rhoa, lt)
            ENDIF
            call ice_growth_diskice(ssi(n),esati_pre,tabs_LCM(i,j,k),pres_LCM(k), therco_local, diffelq_local, &
               f_vent_mass_ice(n), f_vent_axes_ice(n), mfp, vapvel, scrit_a, scrit_c, gtp1, gtp2, gtp_stand, S_a, S_c, rhoa, lt, &
               mass,aaxis,caxis,rho_apparent,capacitance_ice(n),mslope,slopelg,aresid,cresid,uresid,nup1,nup2, &
               ndwn1,ndwn2,iequiv,requiv,alpha_a,alpha_c)
            dm = mass - dm
            particle(n)%var(iM) = mass * particle(n)%ivar(iA)
            IF (dm .lt. 0) THEN
               particle(n)%var(iMrime) = mass/(mass - dm) * particle(n)%var(iMrime) !B. Ascher 04.11.2025 - Sublimate rimed mass as well; previously, had errors with nonsensical rime fractions!
            ENDIF
            IF (particle(n)%var(iMrime)/particle(n)%var(iM) .gt. 1.0) THEN
               dm = dm + (particle(n)%var(iM) - particle(n)%var(iMrime))
               particle(n)%var(iMrime) = particle(n)%var(iM)
            ENDIF
            particle(n)%var(iAXS) = aaxis
            particle(n)%var(iCXS) = caxis
            particle(n)%var(iDNS) = rho_apparent

            IF ( dm .GT. 0.0 )  THEN
               dep_prof(k) = dep_prof(k) + particle(n)%ivar(iA) * dm
               micro_field(i,j,k,idep) = micro_field(i,j,k,idep) + particle(n)%ivar(iA) * dm
            ELSE
               sub_prof(k) = sub_prof(k) + particle(n)%ivar(iA) * dm
               micro_field(i,j,k,isub) = micro_field(i,j,k,isub) + particle(n)%ivar(iA) * dm
            ENDIF

            dqi_n = dm * particle(n)%ivar(iA)
            dqi   = dqi + dqi_n
            IF ( switch_LEM_micro )  THEN
               dV_n   = dV(k) / ( n_end(i,j,k) - n_start(i,j,k) + 1.0 )
               dqi_n  = dqi_n / ( dV_n * rho_air(k) )
               particle(n)%var(iT) = particle(n)%var(iT) + dqi_n * fac_sub
               particle(n)%var(iq) = particle(n)%var(iq) - dqi_n
            ENDIF
         ! ENDIF
      END DO
   ENDIF

   DO n = n_start(i,j,k), n_end(i,j,k)

       !Ice nucleation paramterizations
      IF(nstep .GE. start_ice .and. particle(n)%ivar(ITYPE) .gt. 0.0 .and. switch_ice) THEN
         !Heterogeneous nucleation - immersion freezing for liquid particles with ice nuclei
         IF ((particle(n)%ivar(iTYPE) .EQ. 2) .and. (particle(n)%var(iM) .gt. 0.0) .and. (tabs_LCM(i,j,k) .le. particle(n)%var(iTFZ)) .and. supersat_liq(n) .GT. 0.0) THEN !immersion freezing when ambient temperature < particle freezing temp
            pren_nuc = particle(n)%ivar(iA)
            IF (switch_ds .and. (TRIM(ds_param) .EQ. "Phillips") .and. (tabs_LCM(i,j,k) .ge. 248.16) .and. (tabs_LCM(i,j,k) .le. 270.16)) THEN
               CALL ds_phillips(i,j,k, n, tabs_LCM(i,j,k), nfrag_tot, nfrag_b)
               nfrag_s = nfrag_tot - nfrag_b
               IF ((nfrag_s .gt. 0) .and. (nfrag_b .gt. 0)) THEN !Fragmentation occurred with small fragments
                  mindiff = ABS(tinymass - particle(n)%var(iM)/particle(n)%ivar(iA))
                  MinN = -1
                  tinymass = pi/6.*rho_ice*(1E-6)**3
                  !Search for other ice particles in the grid box to assign the small fragments to
                  DO p = n_start(i,j,k), n_end(i,j,k) !Yes, this is a double loop, but it should only trigger seldomly
                        IF ((ABS(particle(p)%var(iM)/particle(p)%ivar(iA) - tinymass) .lt. mindiff) .and. (particle(p)%ivar(iTYPE) .lt. 0) .and. &
                        (particle(p)%var(iM)/particle(p)%ivar(iA)/tinymass .lt. 5) .and. (particle(p)%var(iM)/particle(p)%ivar(iA)/tinymass .gt. 0.5)) THEN !Merge with another particle in gridbox
                           mindiff = ABS(particle(p)%var(iM)/particle(p)%ivar(iA) - tinymass)
                           MinN = p
                        ENDIF
                  END DO
                  IF (MinN .gt. 0) THEN !Found an appropriate particle within the gridbox which the fragments can merge to
                        phi = particle(MinN)%var(iCXS)/particle(MinN)%var(iAXS) !Set aspect ratio of fragments to the same as that of the merging particle, for simplicity
                        particle(MinN)%var(iM) = particle(MinN)%var(iM) + tinymass*nfrag_s*particle(n)%ivar(iA)
                        particle(MinN)%var(iDNS) = rho_ice !Set all particles as solid ice. This should almost certainly be the case anyway for small ice
                        particle(MinN)%ivar(iA) = particle(MinN)%ivar(iA) + nfrag_s*particle(n)%ivar(iA)
                        particle(MinN)%var(iAXS) = (3*particle(MinN)%var(iM)/(4*pi*rho_ice*particle(MinN)%ivar(iA)*phi))**(1./3.)
                        particle(MinN)%var(iCXS) = particle(MinN)%var(iAXS)*phi
                        !For simplicity, we'll assume that these ice fragments are pure water vapor, without any trace gases, CCN, or IN, so they don't affect the quantities of these species in the merged particle

                        particle(n)%var(iM) = particle(n)%var(iM) - tinymass*nfrag_s*particle(n)%ivar(iA)
                        particle(n)%ivar(iA) = particle(n)%ivar(iA) * (nfrag_b + 1)
                        !Don't need to change anything else; that'll be taken care of in the main micro_condensation routine

                  ELSE !No appropriate other particle found - assign all fragments to this particle
                        particle(n)%ivar(iA) = particle(n)%ivar(iA) * (nfrag_b + nfrag_s + 1)
                        !Don't need to change anything else; this'll be taken care of in the main micro_condensation routine
                  ENDIF
               ELSEIF ((nfrag_b + nfrag_s) .gt. 0) THEN !Fragmentation, but without one size category of fragments - assign all fragments to the particle
                  particle(n)%ivar(iA) = particle(n)%ivar(iA) * (nfrag_b + nfrag_s + 1)
               ENDIF
            ELSEIF (switch_ds .and. (TRIM(ds_param) .EQ. "Lachapelle") .and. (tabs_LCM(i,j,k) .ge. 248.16) .and. (tabs_LCM(i,j,k) .le. 270.16)) THEN
               CALL ds_lachapelle(i,j,k, n, tabs_LCM(i,j,k), nfrag_tot) !Lachapelle's approach does not distinguish between large and small fragments
               IF (nfrag_tot .gt. 0) THEN
                     particle(n)%ivar(iA) = particle(n)%ivar(iA)*(1 + nfrag_tot)
               ENDIF
            ELSEIF (switch_ds .and. (TRIM(ds_param) .EQ. "Lauber") .and. (tabs_LCM(i,j,k) .ge. 248.16) .and. (tabs_LCM(i,j,k) .le. 270.16)) THEN
               CALL ds_lauber(i,j,k, n, tabs_LCM(i,j,k), nfrag_tot) !Lachapelle's approach does not distinguish between large and small fragments
               IF (nfrag_tot .gt. 0) THEN
                     particle(n)%ivar(iA) = particle(n)%ivar(iA)*(1 + nfrag_tot)
               ENDIF
            ENDIF

            IF (nfrag_tot .gt. 0) THEN
               ds_prof(k) = ds_prof(k) + pren_nuc * nfrag_tot
               micro_field(i,j,k,ids) = micro_field(i,j,k,ids) + pren_nuc * nfrag_tot
            ENDIF

            particle(n)%ivar(iTYPE) = -2 !We'll keep using the convention of -2 being ice with IN, -1 being ice without IN
                                          !since we can use that to keep track of homogeneously vs heterogeneously nucleated ice particles
            particle(n)%var(iAXS) = MAX((3.0*particle(n)%var(iM) / (4.*particle(n)%ivar(iA) * rho_ice * pi))**0.33333333, 0.1E-6)
            particle(n)%var(iCXS) = particle(n)%var(IAXS) !Set an initial spherical ice particle which will then grow with habits
            particle(n)%var(iDNS) = rho_ice !Likewise, set initial density to that of a pure ice crystal
            particle(n)%ivar(iNmono) = 1
            particle(n)%var(iMrime)  = 0.0
            dqf_n = particle(n)%var(iM) !Mass of the superdroplet which is freezing 
            dqf = dqf + dqf_n !Add this to the total mass of water freezing in this timestep
            IF (switch_LEM_micro) THEN
               dV_n = dV(k) / (n_end(i,j,k) - n_start(i,j,k) + 1.0)
               dqf_n = dqf_n / (dV_n * rho_air(k))
               particle(n)%var(iT) = particle(n)%var(iT) + dqf_n * fac_fus
            ENDIF
            nuc_prof(k) = nuc_prof(k) + particle(n)%ivar(iA)
            imm_prof(k) = imm_prof(k) + particle(n)%ivar(iA)
            frz_prof(k) = frz_prof(k) + particle(n)%var(iM)

            micro_field(i,j,k,inuc) = micro_field(i,j,k,inuc) + particle(n)%ivar(iA)
            micro_field(i,j,k,iimm) = micro_field(i,j,k,iimm) + particle(n)%ivar(iA)
            micro_field(i,j,k,ifrz) = micro_field(i,j,k,ifrz) + particle(n)%var(iM)
         ENDIF


!
!--   Homogeneous nucleation rate, with parameterization from T. Kuhn et. al. (2011)
         IF ((particle(n)%ivar(iType) .gt. 0) .and. (tabs_LCM(i,j,k) .lt. 239.15) .and. (particle(n)%var(iM) .gt. 0.0) .and. (supersat_liq(n) .gt. -0.2)) THEN
            !B. Ascher 12.11.2025 - Adjust homogeneous freezing routine to require water mass and at least 80% liquid supersaturation
            jv = nv*boltzk*tabs_LCM(i,j,k)/planckh*exp((-1.*av+bv*tabs_LCM(i,j,k))/(boltzk*tabs_LCM(i,j,k)))
            !Volume-specific homogeneous nucleation rate; #/(m^3 * s)
            jt = jv*4.1888*radius(n)**3*(1.+3.*(srfthick/radius(n))*exp(-1.*dw/(boltzk*tabs_LCM(i,j,k))))
            !Total homogeneous freezing rate for a single drop; #/s
            homofrzprob = 1-exp(-dtn_LCM*jt)
            
            IF (homofrzprob .gt. randomu(rank+1)) THEN
               particle(n)%var(iAXS)   = radius(n) * ( rho_liq / rho_ice )**0.333333333
               particle(n)%var(iCXS)   = radius(n) * ( rho_liq / rho_ice )**0.333333333
               particle(n)%var(iDNS)   = rho_ice
               particle(n)%ivar(iTYPE) = -1*ABS(particle(n)%ivar(iTYPE)) ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
               particle(n)%ivar(iNmono) = 1
               particle(n)%var(iMrime) = 0.0
               dqf_n = particle(n)%var(iM)
               dqf   = dqf + dqf_n
               IF ( switch_LEM_micro )  THEN
                  dV_n   = dV(k) / ( n_end(i,j,k) - n_start(i,j,k) + 1.0 )
                  dqf_n  = dqf_n / ( dV_n * rho_air(k) )
                  particle(n)%var(iT) = particle(n)%var(iT) + dqf_n * fac_fus
               ENDIF
               nuc_prof(k) = nuc_prof(k) + particle(n)%ivar(iA)
               hom_prof(k) = hom_prof(k) + particle(n)%ivar(iA)
               frz_prof(k) = frz_prof(k) + particle(n)%var(iM)

               micro_field(i,j,k,inuc) = micro_field(i,j,k,inuc) + particle(n)%ivar(iA)
               micro_field(i,j,k,ihom) = micro_field(i,j,k,ihom) + particle(n)%ivar(iA)
               micro_field(i,j,k,ifrz) = micro_field(i,j,k,ifrz) + particle(n)%var(iM)
            ENDIF

         ENDIF
      ENDIF
!
!--    Determine melting: instantanous melting
       IF ( particle(n)%ivar(iTYPE) .LT. 0  .AND.  tabs_LCM(i,j,k) .GT. 273.15 )  THEN ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
          particle(n)%ivar(iTYPE)  = ABS(particle(n)%ivar(iTYPE)) ! > 0: liquid, < 0: ice; |2|: with IN, |1|: without IN
          particle(n)%ivar(iNmono) = 1
          particle(n)%var(iMrime)  = 0.0

          dqf_n = -particle(n)%var(iM)
          dqf   = dqf + dqf_n
          IF ( switch_LEM_micro )  THEN
             dV_n   = dV(k) / ( n_end(i,j,k) - n_start(i,j,k) + 1.0 )
             dqf_n  = dqf_n / ( dV_n * rho_air(k) )
             particle(n)%var(iT) = particle(n)%var(iT) + dqf_n * fac_fus
          ENDIF
          nml_prof(k) = nml_prof(k) + particle(n)%ivar(iA)
          mel_prof(k) = mel_prof(k) + particle(n)%var(iM)
          micro_field(i,j,k,imel) = micro_field(i,j,k,imel) + particle(n)%var(iM)

       ENDIF

!--    Complete sublimation of an ice particle - convert it back to liquid
       IF (particle(n)%ivar(ITYPE) .LT. 0 .AND. (particle(n)%var(iM) / particle(n)%ivar(iA) .LE. (4.0 / 3.0 * pi * rho_ice * (0.01E-6)**3))) THEN
         particle(n)%ivar(iTYPE) = ABS(particle(n)%ivar(iTYPE)) !>0: liquid, <0 : ice, |2|: with IN, |1|: without IN
         particle(n)%ivar(iNmono) = 1
         particle(n)%var(iMrime) = 0.0
         IF (particle(n)%ivar(iTYPE) .EQ. 2) THEN
            particle(n)%ivar(iA) = INT(n_ice(MAXLOC(particle(n)%var(in_startind+1:in_startind+4), dim = 1)) * dV(k) / (np_xy * frac_ice), KIND=8)
         ELSE
            particle(n)%ivar(iA) = INT(n_aero(MAXLOC(particle(n)%var(iAero_portions), dim = 1)) * dV(k) / (np_xy * (1.0 - frac_ice)), KIND=8)
         ENDIF
         IF (TRIM(unit_aero) .EQ. "1/kg") THEN
            particle(n)%ivar(iA) = particle(n)%ivar(iA) * rho_air(k)
         ENDIF
         !To account for the fact that SIP doesn't increase the number of INP, revert ice multiplicity back to the number concentration of INP 
         !(based on the dominant INP mode in the ice particle) upon sublimation - B.Ascher 02.03.2026

         dqf_n = -particle(n)%var(iM)
         dqf = dqf + dqf_n
         IF ( switch_LEM_micro )  THEN
            dV_n   = dV(k) / ( n_end(i,j,k) - n_start(i,j,k) + 1.0 )
            dqf_n  = dqf_n / ( dV_n * rho_air(k) )
            particle(n)%var(iT) = particle(n)%var(iT) + dqf_n * fac_fus
         ENDIF

         nsub_prof(k) = nsub_prof(k) + particle(n)%ivar(iA)
      ENDIF
    ENDDO
    micro_field(i,j,k,imv)      = micro_field(i,j,k,imv) - ( dql + dqi ) / ( dV(k) * rho_air(k) )
    tabs_LCM(i,j,k)             = tabs_LCM(i,j,k) + ( fac_cond * dql + fac_sub * dqi + fac_fus * dqf ) / ( dV(k) * rho_air(k) )

    micro_field(i,j,k,iss)      = micro_field(i,j,k,imv) - qsatw( tabs_LCM(i,j,k), pres_LCM(k) ) !absolute supersaturation
    micro_field(i,j,k,isLES)    = micro_field(i,j,k,iss) / qsatw( tabs_LCM(i,j,k), pres_LCM(k) ) !relative supersaturation

    micro_field(i,j,k,ilastT)   = tabs_LCM(i,j,k)  ! temperature at last time step (advected quantity)
    micro_field(i,j,k,ilastp)   = pres_LCM(k)      ! pressure at last time step    (advected quantity)
!
!-- Save supersaturation after condensation, but before dynamics for estimating
!-- diabatic effects on supersaturation in micro_sgs_mixing.f90
    IF ( switch_LEM_micro )  THEN
       DO n = n_start(i,j,k), n_end(i,j,k)
          particle_tabs(n)      = particle(n)%var(iT) - particle(n)%var(iZ) * ggr / cp
          particle(n)%var(ieta) = REAL( particle(n)%var(iq) - qsatw( particle_tabs(n), pres_LCM(k) ), KIND=lcmp )
       ENDDO
    ENDIF

    CALL t_stopf ('micro_condensation')

 END SUBROUTINE micro_condensation

 REAL FUNCTION venteff(rad,ws,rho_air_int)

    USE micro_vars, ONLY: lcmp
    USE params,     ONLY: muelq

    IMPLICIT NONE

    REAL            :: rad, Re_p, rho_air_int
    REAL(KIND=lcmp) :: ws
!
!-- Compute ventilation effect for droplets > 40 microns (Rogers and Yau, 1989):
    IF ( rad > 40.0E-6 )  THEN  !.AND.  micro_field(i,j,k,iss) < 0.0 )  THEN
!
!--    Calculate droplet's Reynolds number
       Re_p = 2.0 * rad * ws / ( muelq / rho_air_int )

       IF ( Re_p > 2.5 )  THEN
          venteff = 0.78 + 0.28 * SQRT( Re_p )
       ELSE
          venteff = 1.0 + 0.09 * Re_p
       ENDIF
    ELSE
       venteff = 1.0
    ENDIF

    RETURN

 END FUNCTION venteff

 REAL FUNCTION eta_int( eta_lst, forcing, tau, dt )
!
!-- Get integral absolute supersaturation following Tzivion, Feingold, and Levin (1989, JAS)
    IMPLICIT NONE

    REAL, INTENT(IN) :: eta_lst, forcing, dt, tau

    eta_int = tau * ( 1.0 - EXP( -dt / tau ) ) * ( eta_lst - tau * forcing ) + tau * forcing * dt

 END FUNCTION eta_int

 REAL FUNCTION eta_new( eta_lst, forcing, tau, dt )
!
!-- Get new absolute supersaturation following Tzivion, Feingold, and Levin (1989, JAS)
    IMPLICIT NONE

    REAL, INTENT(IN) :: eta_lst, forcing, dt, tau

    eta_new = tau * forcing + ( eta_lst - tau * forcing ) * EXP( -dt / tau )

 END FUNCTION eta_new

 SUBROUTINE radius_liquid_euler(r_ini,dt_int, &
                                r0,C_pre,supersat,ventilation_effect, &
                                afactor,bfactor,r_aero,D_pre,radiation)
!
!-- Newton-Raphson scheme to solve droplet growth equation.
!-- Note that not dr/dt = ... is used, but dr^2/dt = ... based on Chen (1992)
    USE micro_vars, ONLY: rho_liq

    IMPLICIT NONE

    REAL, INTENT(INOUT)  :: r_ini
    REAL, INTENT(IN)     :: dt_int, r0, C_pre, supersat, ventilation_effect, &
                            afactor, bfactor, r_aero, D_pre, radiation

    INTEGER :: m
    REAL    :: r_eul, r_eul_old, rel_change, dr2dt, d2r2dtdr2, f, dfdr2, dt_eul, t_eul, G_pre

    G_pre = C_pre / rho_liq

    r_eul     = r_ini
    r_eul_old = r_ini

    dt_eul = dt_int
    t_eul  = 0.0

    DO WHILE ( t_eul .LT. dt_int - 1.0E-20 )

       DO m = 1, 500

          dr2dt      = 2.0 * G_pre * ventilation_effect &
                       * ( supersat - afactor / r_eul + bfactor / r_eul**3 - D_pre * radiation * r_eul ) &
                       * r_eul / ( r_eul + r0 )
          d2r2dtdr2  = G_pre * ventilation_effect * &
                       ( afactor * r_eul**3 - bfactor * ( 3.0 * r_eul + 2.0 * r0 ) &
                         - r_eul**3 * ( D_pre * radiation * r_eul * ( r_eul + 2.0 * r0 ) - r0 * supersat ) &
                       ) / ( r_eul**4 * ( r_eul + r0 )**2 )

          dt_eul = MIN( dt_int - t_eul, dt_int )!MIN( 0.5 * ABS( 1.0 / d2r2dtdr2 ), dt_int - t_eul, dt_int )
!
!--       To speed up the Newton-Raphson scheme, the square root is executed at every iteration
          f       = r_eul**2 - r_ini**2 - dt_eul * dr2dt
          dfdr2   = 1.0 - dt_eul * d2r2dtdr2 ! = 1.0 - ( 0.0 + dt_eul * d2r2dtdr2 )

          r_eul   = SQRT( MAX( r_eul**2 - f / dfdr2, r_aero**2 ) ) ! Newton-Raphson scheme (2nd order)

          rel_change = ABS( r_eul - r_eul_old ) / r_eul_old
          r_eul_old  = r_eul
          IF ( rel_change .LT. 1.0E-12 )  THEN
             EXIT
          ENDIF

       ENDDO

       t_eul = t_eul + dt_eul
       r_ini = r_eul

    ENDDO

    RETURN

 END SUBROUTINE radius_liquid_euler

 SUBROUTINE mass_ice_euler(m_ini,dt_int,C_pre,capacitance_pre,f_vent_mass_ice,supersat_ice)
!
!-- Simple Euler forward step
    USE params,     ONLY: pi
    USE micro_vars, ONLY: rho_ice

    IMPLICIT NONE

    REAL, INTENT(INOUT) :: m_ini
    REAL, INTENT(IN)    :: dt_int, C_pre, capacitance_pre, f_vent_mass_ice, supersat_ice

    REAL    :: dm23dt, massold
!
!-- In the following, we use a prognostic equation for m^(2/3) instead of m,
!-- where dmdt = 4.0 * pi * G_pre * CAP * f_vent * supersat_ice.
!-- Due to this choice, CAP / m_ini**(1.0/3.0) is approximately constant
!-- with time, and changes in the capacitance can be neglegted
    massold = m_ini
    dm23dt = 8.0 / 3.0 * pi * C_pre * ( capacitance_pre / m_ini**(1.0/3.0) ) * f_vent_mass_ice * supersat_ice
    m_ini  = MAX( ( m_ini**(2.0/3.0) + dt_int * dm23dt ), ( 4.0 / 3.0 * pi * rho_ice * 0.001E-6**3 )**(2.0/3.0) )**(3.0/2.0)
    RETURN

 END SUBROUTINE mass_ice_euler

 SUBROUTINE ice_growth_diskice(ssi,esati_pre,tkelv,ppa,therco_local,diffelq_local,&
   f_vent_mass_ice, f_vent_axes_ice, mfp, vapvel, scrit_a, scrit_c, gtp1, gtp2, gtp_stand, S_a, S_c, rhoa, lt, &
   m_ini,alen,clen,rho_apparent,cap_pre,m,slopelg,aresid,cresid,uresid,nup1,nup2, &
   ndwn1,ndwn2,iequiv,requiv,alpha_a,alpha_c)
   USE PARAMS, ONLY: cp, pi, lsub, rv, rgas
   USE MICRO_VARS, ONLY: rho_ice, dtn_LCM, switch_diskice_sizethresh, sizethresh_diskice
   
   IMPLICIT NONE
   REAL :: alen,clen,phi,alpha_a,alpha_c, requiv,capdel, ktprime, dvprime, geff, C_ice_pre
   REAL, EXTERNAL :: capacitance, esati, gamma_ice
   REAL, INTENT(IN) :: ssi, esati_pre, tkelv, ppa, therco_local, diffelq_local, &
                       f_vent_mass_ice, f_vent_axes_ice, &
                       mfp, vapvel, scrit_a, scrit_c, &
                       gtp1, gtp2, gtp_stand, S_a, S_c, rhoa, lt, &
                       m,slopelg, aresid, cresid, uresid, nup1, nup2, &
                       ndwn1, ndwn2, cap_pre
   INTEGER, INTENT(IN) :: iequiv
   REAL :: m_ini, dm23dt, dm, m_new, drho, rho_deposition, rho_apparent, V_ini, dV, V_new, da, daprelim, dcprelim, gamma_dep, gamma_temp
   REAL :: dlogV, dlna
   capdel = capacitance(alen+mfp, clen+mfp)
   !-- Calculates alpha_c and alpha_a
   call depcoeffparam(ssi,esati_pre,tkelv,ppa,diffelq_local,&
   mfp,vapvel,scrit_a,scrit_c,gtp1,gtp2,gtp_stand,S_a,S_c,rhoa,lt, &
   alen,clen,cap_pre,capdel,m,slopelg,aresid,cresid,uresid,nup1,nup2, &
   ndwn1,ndwn2,iequiv,requiv,alpha_a,alpha_c)
 
 
   dvprime = 2.0/3.0*(diffelq_local / (4.0*diffelq_local * cap_pre / (alpha_a * vapvel * alen*clen) + cap_pre/capdel) )&
            +1.0/3.0*(diffelq_local / (4.0*diffelq_local * cap_pre / (alpha_c * vapvel * alen**2) + cap_pre/capdel))
 
   ktprime = 2.0/3.0*(therco_local / (4.0 * therco_local * cap_pre / (0.96 * rhoa * cp*  vapvel * alen * clen) + cap_pre/capdel)) &
             +1.0/3.0*(therco_local / (4.0*therco_local*cap_pre / (0.96 * rhoa * cp * vapvel * alen **2) + cap_pre/capdel))
   geff  = 1.0 / ( rv * tkelv / ( 100.0 * esati_pre * dvprime ) +  &
             ( lsub / ( rv * tkelv ) - 1.0 ) * lsub /        &
             ( ktprime * tkelv )  )
   dm23dt = 8.0/3.0 * pi * geff * (cap_pre / m_ini ** (1.0/3.0)) * f_vent_mass_ice * ssi
   dm = m_ini
   m_ini  = MAX( ( m_ini**(2.0/3.0) + dtn_LCM * dm23dt ), ( 4.0 / 3.0 * pi * rho_ice * 0.01E-6**3 )**(2.0/3.0) )**(3.0/2.0)
   dm = m_ini - dm
   phi = clen / alen
   gamma_temp = gamma_ice(tkelv)
   gamma_dep = alpha_c / alpha_a * f_vent_axes_ice
   !
   !-- Deposition density following Chen and Lamb (1994)
   drho           = geff * ssi / dvprime ! vapor density difference, dmdt = 4 * pi * C_pre * diff_coeff * f_vent * drho = 4 * pi * C_pre * C_ice * f_vent * supersat_ice <==> drho = C_ice * supersat_ice / diff_coeff
   !-- Deposition density following Jensen and Harrington (2015)
   !For plates, rho_dep = rho_ice * gamma_dep if a-axis > 100 um, rho_dep = rho_ice otherwise
   !For columns, rho_dep = rho_ice / gamma_dep
   IF ( (gamma_dep .LT. 1.0  .AND.  alen .LT. 100.0E-6) .OR. (switch_diskice_sizethresh .and. min(alen, clen) .lt. sizethresh_diskice))  THEN !  Jensen and Harrington (2015) !Note that this threshold causes an abrubt "jump" in growth rates above 100 microns.
      rho_deposition = rho_ice
   ELSE
       rho_deposition = rho_ice * EXP( -3.0 * MAX( drho * 1000.0 - 0.05, 0.0 ) / gamma_temp )  ! (42) in Chen and Lamb (1994, JAS)
   ENDIF
   V_ini        = (m_ini - dm) / rho_apparent
 
   IF ( dm .LT. 0.0 )  THEN ! sublimation: density does not change, Eqn. (23) in Shima et al. (2020)
       dV = dm / rho_apparent
   ELSE                     ! deposition: considering current deposition density; Eqn. (41) in Chen and Lamb (1994, JAS)
       dV = dm / rho_deposition
   ENDIF
   
   V_new = MAX(V_ini + dV, 1.0e-24)
 
   IF (dm .ge. 0.0) THEN
      IF (switch_diskice_sizethresh .and. (MIN(alen, clen) .lt. sizethresh_diskice)) THEN
         alen        = ( ( m_ini) / ( 4.0 / 3.0 * pi * rho_ice ) )**0.333333333
         clen        = ( ( m_ini) / ( 4.0 / 3.0 * pi * rho_ice ) )**0.333333333
         rho_apparent = rho_ice
      ELSE
 
       daprelim = dV / V_ini * (1.0 / (2.0 / alen + gamma_dep / clen)) !Use a corrector scheme approach for now. Expensive, but accurate
       dcprelim = gamma_dep * daprelim
       da = 1.0 / 2.0 * daprelim + 1.0 / 2.0 * (dV / V_new) * (1.0 / (2.0 / (alen + daprelim) + gamma_dep / (clen + dcprelim)))
       alen = alen+da
       clen = clen+gamma_dep*da
       phi = clen/alen
       alen = (V_new / (4.0/3.0 * pi * phi))**(1.0/3.0) !Get rid of any small inconsistencies between mass volume and axes volume
       clen = phi * alen
       rho_apparent = m_ini / V_new
      ENDIF
   ELSE
     alen = (V_new / (4.0/3.0 * pi * phi))**(1.0/3.0)
     clen = phi * alen
   ENDIF
 END SUBROUTINE ice_growth_diskice
 REAL FUNCTION capacitance(aaxis,caxis)
!
!-- Capacities of ice. Eqns. (39) and (40) from Chen and Lamb (1994)
    IMPLICIT NONE

    REAL, INTENT(IN)  :: aaxis, caxis
    REAL              :: phi, epsilon

    phi = caxis / aaxis
!
!-- Oblate spheroids
    IF ( phi .LT. 1.0 )  THEN
       epsilon     = SQRT( 1.0 - phi**2 )
       capacitance = aaxis * epsilon / ASIN( epsilon )
!
!-- Prolate spheroids
    ELSEIF ( phi .GT. 1.0 )  THEN
        epsilon     = SQRT( 1.0 - phi**(-2) )
        capacitance = caxis * epsilon / LOG( ( 1.0 + epsilon ) * phi )
!
!-  Sphere
    ELSEIF ( phi .EQ. 1.0 )  THEN
       capacitance = aaxis
    ENDIF

    !capacitance = ( 2 * aaxis + caxis ) / 3.0 ! approximation from Shima et al. (2020)

    RETURN

 END FUNCTION capacitance

 REAL FUNCTION gamma_ice(T)
!
!-- Lookup table of inherent growth ratio gamma, using values of Fig. 3 in Chen and Lamb (1994)
    IMPLICIT NONE

    REAL, INTENT(IN) :: T

    REAL, DIMENSION(61), PARAMETER :: &
       T_values = (/ -30.0, -29.5, -29.0, -28.5, -28.0, -27.5, -27.0, -26.5, -26.0, -25.5, &
                     -25.0, -24.5, -24.0, -23.5, -23.0, -22.5, -22.0, -21.5, -21.0, -20.5, &
                     -20.0, -19.5, -19.0, -18.5, -18.0, -17.5, -17.0, -16.5, -16.0, -15.5, &
                     -15.0, -14.5, -14.0, -13.5, -13.0, -12.5, -12.0, -11.5, -11.0, -10.5, &
                     -10.0,  -9.5,  -9.0,  -8.5,  -8.0,  -7.5,  -7.0,  -6.5,  -6.0,  -5.5, &
                     -5.0,   -4.5,  -4.0,  -3.5,  -3.0,  -2.5,  -2.0,  -1.5,  -1.0,  -0.5, &
                      0.0 /)

    REAL, DIMENSION(61), PARAMETER :: &
       G_values = (/ 1.2845082,  1.3063233,  1.334114,   1.3740169,  1.4210837,  1.4697628,  1.5329632,  &
                     1.5988813,  1.6817352,  1.7688825,  1.8605458,  1.9405501,  2.0070235,  2.0325339,  &
                     1.9405501,  1.7540505,  1.5459258,  1.3229275,  1.1758163,  1.0276118,  0.90187675, &
                     0.79486565, 0.69760859, 0.60967931, 0.54188221, 0.47358129, 0.42269433, 0.37254004, &
                     0.33391297, 0.29929097, 0.2739657,  0.25504251, 0.24555938, 0.2465954,  0.25828424, &
                     0.27628232, 0.30309514, 0.34245537, 0.39849983, 0.50022382, 0.6604531,  1.0189954,  &
                     1.9901947,  2.9564446,  3.2161679,  3.2161679,  2.9814439,  2.4772788,  1.8683956,  &
                     1.2058968,  0.65491523, 0.3694163,  0.34535113, 0.46371622, 0.70945624, 0.78488925, &
                     0.82904522, 0.87937983, 0.91719354, 0.96066932, 0.9977689 /)

    INTEGER :: G_loc

    G_loc = FLOOR( ( ( T-273.15 ) + 30.5 ) * 2.0 )
    IF ( G_loc .LT. 1 ) THEN
       gamma_ice = G_values(1)
    ELSEIF ( G_loc .GE. 61 ) THEN
       gamma_ice = G_values(61)
    ELSE
       gamma_ice = G_values(G_loc) + ( G_values(G_loc+1) - G_values(G_loc) ) / 0.5 * ( ( T-273.15 ) - T_values(G_loc) )
    ENDIF

    RETURN

 END FUNCTION gamma_ice

 subroutine getscrit(temp,sa,sc)
  
   implicit none
   
   real temp,sa,sc
   real dt,x
   
   dt = temp-273.15
   x = max(-70.,min(dt,-1.0))
   if(dt.lt.-30.0)then
      ! original magee (high scrit) 
     !  sc = 1.8115 + 0.15585*x + 0.011569*x**2
     !  sa = 1.8115 + 0.15585*x + 0.011569*x**2
      
      ! with bailey and hallett
      sc = 3.7955 + 0.10614*x + 0.0075309*x**2
      sa = sc
      ! modified magee
      !sc = 3.3649 + 0.13266*x + 0.008959*x**2
      !sa = sc
   else if(dt.ge.-30.0.and.dt.le.-22.0)then
     !  sc = 753.63 + 105.97 * x + 5.5532 * x**2 + 0.12809 * x**3  &
     !       + 0.001103 * x**4
      sc = 1.1217 + 0.038098 * x - 0.083749 * x**2  &
           - 0.015734 * x**3 - 0.0010108 * x**4  &
           - 2.9148e-05 * x**5 - 3.1823e-07 * x**6
      sa = -0.71057 - 0.14775*x + 0.0041200*x**2
   ! else if (dt.gt.-22.0.and.dt.le.-1.0) then
   !    sc = 1.1217 + 0.038098 * x - 0.083749 * x**2  &
   !         - 0.015734 * x**3 - 0.0010108 * x**4  &
   !         - 2.9148e-05 * x**5 - 3.1823e-07 * x**6
   elseif (dt .gt. -22.0 .and. dt .le. -15.0) then
     sa = -5.2367 - 1.3184*x - 0.11066*x**2 - 0.0032303*x**3
     sc = 1.1217 + 0.038098 * x - 0.083749 * x**2  &
           - 0.015734 * x**3 - 0.0010108 * x**4  &
           - 2.9148e-05 * x**5 - 3.1823e-07 * x**6
   elseif (dt.gt.-15.0.and.dt.le. 0.0)then
     !sa = 0.755 + 0.03325*x + 0.00125*x**2   ! fit to Woods
     sc = 1.1217 + 0.038098 * x - 0.083749 * x**2  &
           - 0.015734 * x**3 - 0.0010108 * x**4  &
           - 2.9148e-05 * x**5 - 3.1823e-07 * x**6
     sa = 0.34572 - 0.0093029*x + 0.00030832*x**2 ! fit to Nelson & Night
   ! if(dt.ge.-30.0.and.dt.le.-22.0)then
   !   !  sa =  -0.71057 - 0.14775*x + 0.0042304*x**2
   !   sa = -0.71057 - 0.14775*x + 0.0041200*x**2
   ! elseif (dt.gt.-22.0.and.dt.le.-15.0)then
   !    sa = -5.2367 - 1.3184*x - 0.11066*x**2 - 0.0032303*x**3
   ! elseif (dt.gt.-15.0.and.dt.le.-10.0)then
      !sa = 0.755 + 0.03325*x + 0.00125*x**2   ! fit to Woods
     !  sa = 0.34572 - 0.0093029*x + 0.00030832*x**2 ! fit to Nelson & Night
   ! elseif (dt.gt.-10.0.and.dt.le.-1.0)then
      !          sa = 0.37445 + 0.029636*x + 0.011575*x**2 &      ! Wood
      !               + 0.00069444*x**3
     !  sa = 0.34572 - 0.0093029*x + 0.00030832*x**2 ! Nelson & Knight
   endif
      
   sa = sa/100.
   sc = sc/100.
   
   return
 end subroutine getscrit

subroutine diffparams(tkelv, ppa, therco_local, diffelq_local, esati_pre, ssi, mfp, vapvel, scrit_a, scrit_c, gtp1, gtp2, gtp_stand, S_a, S_c, rhoa, lt)
   use params, only: lsub, rv, rgas, rstar, pi, cp
   implicit none
   real, intent(in) :: tkelv, ppa, therco_local, diffelq_local, esati_pre, ssi, scrit_a, scrit_c
   real, external   :: diffcoeff
   real             :: mfp, vapvel, gtp1, gtp2, gtp_stand, S_a, S_c, rhoa, lt
   real, parameter :: mw = 0.01802, alphat = 0.96 !Molecular weight of water, thermal accommodation coefficient
 
   !-- Calculate coefficients used in Diskice growth scheme
   !ssi is similar, but not identical to supersat_ice. It is based on vapor pressure, rather than mixing ratio
   mfp = 6.6e-8*(1013.25/ppa)*(tkelv/273.15)        ! mean free path in m
   vapvel = sqrt(8.0*rstar*tkelv/(pi*mw))              ! molecular speed(m/s)
   ! scrit_a = supersat_a(tkelv) !This is an approximation based on interpolation which should be much faster than the getscrit routine- Ben A. Dec. 2024
   ! scrit_c = supersat_c(tkelv)
   !  Calculation of enthalpies -- commented this section out b/c Lsub is always extremely close to 2.83*10**6 - Ben A. Dec. 2024
   ! Lv = 2.5008e6*(273.15/tkelv)**(0.167+3.67e-4*tkelv)
   ! Lf = 0.3337e6+(2031.0-10.467*(tkelv-273.15))*(tkelv-273.15)
   ! Lsub = Lv+Lf
 
   gtp1 = rv*tkelv/(diffelq_local*100.*esati_pre)
   gtp2 = (lsub**2.0/(therco_local*rv*tkelv**2.0)) - &
      (lsub/(therco_local*tkelv))
   gtp_stand = 1.0/(gtp1+gtp2)                      ! effective diffusivity
   S_a = abs(ssi/scrit_a)                           ! critical sat ratio for a axis
   S_c = abs(ssi/scrit_c)                           ! critical sat ratio for c axis
 
   !pratio = 1013.15/ppa
   !rlambda = 6.6e-6 * pratio * (tkelv / 293.15)/100.
   !delv = rlambda
   !Commented the above line out because it seems the same as mfp above, though confusingly with 293.15 K instead of 273.15 K as in mfp.
 
   rhoa = 100.*ppa/(rgas*tkelv)
   lt = therco_local/(alphat*rhoa*cp)*sqrt(2.*pi/(rgas*tkelv))
 end subroutine diffparams

!----------------------------------------------------------------------
! SUBROUTINE GETFITTINGPARAMS
!
! This subroutine gets the various parameters used in the fits for the
! surface supersaturation based on the chosen value of M
!     JYH, Penn State University, May 15 2019
!
!  Input: m (parameter in Nelson and Baker (1996) approximation for
!            the deposition coefficient)
!
!  Output: slopelg,aresid,cresid,uresid,nup1,nup2,ndwn1,ndwn2
!          (Coefficients in surface supersaturation polynomial fit)
!
!----------------------------------------------------------------------
subroutine fittingparams(m,slopelg,aresid,cresid,uresid,nup1,nup2,ndwn1,ndwn2)

implicit none
real m   ! input, m in equation for deposition coefficient
real slopelg,aresid,cresid,uresid,nup1,nup2, &    ! output, fitting parameters
     ndwn1,ndwn2

slopelg = 0.832
slopelg = 0.1532 + 0.49078 * m - 0.18002 * m**2 + 0.04165 * m**3 - &
     0.0061541 * m**4 + 0.00058009 * m**5 - 3.3791e-05 * m**6 + &
     1.1092e-06 * m**7 - 1.5691e-08 * m**8
aresid = 0.46192 - 0.093575 * m + 0.012798 * m**2 - &
     0.0008756 * m**3 + 2.2466e-05 * m**4
cresid = aresid
uresid = 0.25872 - 0.057957 * m + 0.0079602 * m**2 - &
     0.00052905 * m**3 + 1.3145e-05 * m**4

nup1 = -0.69061 + 10.765 * m - 9.1666 * m**2 + 3.8888 * m**3 - &
     0.91779 * m**4 + 0.12544 * m**5 - 0.0098291 * m**6 +  &
     0.00040868 * m**7 - 6.9744e-06 * m**8

nup2 = 5.1239 - 19.204 * m + 15.236 * m**2 - 6.4802 * m**3 +  &
     1.5151 * m**4 - 0.20317 * m**5 + 0.015601 * m**6 - &
     0.00063752 * m**7 + 1.0736e-05 * m**8

ndwn1 = -0.86527 - 0.03202 * m + 0.027342 * m**2 - 0.04378 * m**3 &
     + 0.014734 * m**4 - 0.0021036 * m**5 + 0.0001431 * m**6 &
     - 4.317e-06 * m**7 + 3.9242e-08 * m**8

ndwn2 = -4.692 + 12.47 * m - 9.8936 * m**2 + 4.1881 * m**3 &
     - 1.0003 * m**4 + 0.14011 * m**5 - 0.011359 * m**6 &
     + 0.00049122 * m**7 - 8.7303e-06 * m**8      

return
end subroutine fittingparams

!-------------------------------------------------------------------------
! SUBROUTINE DEPCOEFFPARAM
!
! This subroutine finds the deposition coefficients for spheroids. 
! The theory follows that of Zhang and Harrington (2014)
! however the iteration to find the alpha values is replaced
! by a parameterization. This version allows for a choice of m.
!
!  Article: Harrington et al. (2021, J. Atmos. Sci.)
!
!   JYH, Penn State University, May 15, 2019
!
!  Input:
!  ------
!    ssi:  ambient ice supersaturation (units: fraction)
!    tkelv: Temperature  (units: Kelvin)  
!    ppa:   Atmospheric pressure (units: Pa)
!    alen:  spheroid a-semi-dimension (half diameter basal plane) (units: m)
!    clen:  spheroid c-semi-dimension (half prism plane height) (units: m)
!    iequiv:  = 1  (calculation done for equivalent volume sphere)
!             = any other value for spheroidal calculation
!    requiv:  equivalent volume sphere radius (units: m)
!    m:       parameter in Nelson and Baker (1996) approximation for
!            the deposition coefficient
!    Coefficients: slopelg,aresid,cresid,uresid,nup1,nup2,ndwn1,ndwn2
!           (Coefficients in surface supersaturation polynomial fit. 
!            Caculated in subroutine fittingparams)
!
!  Output:
!  -------
!      alpha_ap:  Deposition coefficient for a-axis
!      alpha_cp:  Deposition coefficient for c-axis
!                 For equivalent volume sphere, alpha_ap = alpha_cp
!      alpha_t:   Thermal accomodation coefficient, simply set to 0.96
!
!-------------------------------------------------------------------------

subroutine depcoeffparam(ssi,esati_pre,tkelv,ppa,diffelq_local,&
   mfp, vapvel, scrit_a, scrit_c, gtp1, gtp2, gtp_stand, S_a, S_c, rhoa, lt, &
   alen,clen,cap,capnew,m,slopelg,aresid,cresid,uresid,nup1,nup2, &
   ndwn1,ndwn2,iequiv,requiv,alpha_a,alpha_c)
 
   use params, only: lsub, pi, rstar, rv, rgas, cp
 
 implicit none
 integer iequiv
 ! real ssi,tkelv,pmb,C3a,Dv, esati_pre
 ! real mfp,vel,Pi,Mw,Runiv,m,ppa
 real :: scrit_a !The reason this isn't intent(in) is because it gets modified if iequiv is 1
 real, intent(in) :: ssi, esati_pre, tkelv, ppa, diffelq_local, mfp, vapvel, &
                    scrit_c, gtp1, gtp2, gtp_stand, S_a, S_c, rhoa, lt, &
                    alen,clen,m,slopelg,aresid,cresid,uresid,nup1,nup2,ndwn1,ndwn2,cap, capnew
 
 real C3a, &
     caplengtha,caplengthc,C3c,requiv, &
     alpha_a,alpha_c, &
     alenfix, clenfix, &
     slocal_a,slocal_c,suamb, &
     amn,cmn,capacitance,sisloc_c,sisloc_a,    &
     slocal_diff_a,slrata_interp,slocal_ap, &
     slocal_diff_c,slratc_interp,slocal_cp, &
     scrat_low,scrat_hi, &
     scratio_a,scratio_c,slrata_interp_lim,sldiffa_upratio, &
     slratc_interp_lim,sldiffc_upratio,To, &
     sa_temp,sc_temp
 
 alenfix = alen
 clenfix = clen
 if(iequiv.eq.1)then   ! assumed equiv vol sphere
 requiv = (alen**2 * clen)**(1.0/3.0)
 alenfix = requiv
 clenfix = requiv
 endif
 
 if(iequiv.eq.1) scrit_a = 2./3.*scrit_a + 1./3.*scrit_c 
 
 amn = alenfix + mfp
 cmn = clenfix + mfp
 caplengtha = alenfix*clen/cap
 caplengthc = alenfix**2/cap
 
 C3a = (vapvel*caplengtha*cap/capnew)/(4.*diffelq_local)       
 C3c = (vapvel*caplengthc*cap/capnew)/(4.*diffelq_local)
 
 ! Get alpha for a-axis using parameteriztion
 !_____________________________________________________
 
 sisloc_a = diffelq_local/(gtp_stand*rv*tkelv/(100.*esati_pre)) /(1./(1.+C3a))
 
 slocal_diff_a = ssi/sisloc_a
 scratio_a = min(max(slocal_diff_a/scrit_a,scrat_low),scrat_hi)
 slrata_interp = 10.**(slopelg*alog10(slocal_diff_a/scrit_a))
 slrata_interp_lim = max(min(slrata_interp,1.0) , &
   1.0/(sisloc_a))
 sldiffa_upratio = max(0.0,slrata_interp-1.) + 1.
 slrata_interp = max(slrata_interp , &
   1.0/(sisloc_a)) &
   + aresid*slrata_interp * &
   (slrata_interp/(1./sisloc_a))**(ndwn1) * &
   min((slrata_interp/(1./sisloc_a))**(ndwn2),1.0) &
   - uresid*slrata_interp * slrata_interp_lim**nup1 * &
   sldiffa_upratio**(nup2)
 slocal_ap = slocal_diff_a/slrata_interp
 alpha_a=(slocal_ap/scrit_a)**m * tanh((scrit_a/slocal_ap)**m)
 slocal_a = slocal_ap
 alpha_a = max(min(alpha_a,1.0),1.e-6)
 
 ! Get alpha for c-axis next from parameterization
 !_____________________________________________________
 
 if(iequiv.ne.1)then
 
 sisloc_c = diffelq_local/(gtp_stand*rv*tkelv/(100.*esati_pre)) /(1./(1.+C3c))
 
 slocal_diff_c = ssi/sisloc_c
 scratio_c = min(max(slocal_diff_c/scrit_c,scrat_low),scrat_hi)
 slratc_interp = 10.**(slopelg*alog10(slocal_diff_c/scrit_c))
 slratc_interp_lim = max(min(slratc_interp,1.0) , &
      1.0/(sisloc_c))
 sldiffc_upratio = max(0.0,slratc_interp-1.) + 1.
 slratc_interp = max(slratc_interp, &
      1./sisloc_c)     &
      + cresid*slratc_interp * &
      (slratc_interp/(1./sisloc_c))**(ndwn1) * &
      min((slratc_interp/(1./sisloc_c))**(ndwn2),1.0) &
 !           - 0.11*slratc_interp * slratc_interp_lim**nup1 * &
      - uresid*slratc_interp * slratc_interp_lim**nup1 * &
      sldiffc_upratio**(nup2)
 slocal_cp = slocal_diff_c/slratc_interp
 alpha_c=(slocal_cp/scrit_c)**m * tanh((scrit_c/slocal_cp)**m)
 alpha_c = max(min(alpha_c,1.0),1.e-6)
 slocal_c = slocal_cp
 else
 alpha_c = alpha_a
 slocal_c = slocal_a
 endif
 
 return
 end subroutine depcoeffparam

 SUBROUTINE acd_update(aaxis,caxis,rho_apparent,mass_int,dm_int,C_ice_int,supersat_ice_int,T_int,diffelq_local,f_vent_axes_ice_int)
   !
   !-- update crystal axes and density
       USE params,     ONLY: muelq, pi
       USE micro_vars, ONLY: rho_ice, sizethresh_chen
   
       IMPLICIT NONE
   
       REAL, INTENT(INOUT) ::  aaxis, caxis, rho_apparent
       REAL, INTENT(IN)    ::  mass_int, dm_int, C_ice_int, supersat_ice_int,     &
                               T_int, f_vent_axes_ice_int, diffelq_local
       REAL                ::  V_ini, dV_ice, dlogV_ice, mass_ini, gamma_ast,     &
                               gamma_ice, drho, rho_deposition, gamma_ice_int,    &
                               phi_pre

       IF ( MAX(aaxis,caxis) .LT. sizethresh_chen)  THEN
           gamma_ice_int = 1.0
       ELSE
           gamma_ice_int = gamma_ice(T_int)
       ENDIF
       IF ( T_int .LT. 243.15 )  THEN
          gamma_ice_int = 1.28
       ENDIF

       phi_pre = caxis/aaxis !Determine initial aspect ratio. Use during sublimation, when aspect ratio is preserved
   !
   !-- Deposition density following Chen and Lamb (1994)
       drho           = C_ice_int * supersat_ice_int / diffelq_local ! vapor density difference, dmdt = 4 * pi * C_pre * diff_coeff * f_vent * drho = 4 * pi * C_pre * C_ice * f_vent * supersat_ice <==> drho = C_ice * supersat_ice / diff_coeff
       ! drho           = C_ice_int * MIN( supersat_ice_int, esatw(T_int) / esati(T_int) - 1.0 ) / diffelq ! reduce to water saturation based on Miller and Young (1979) -> is done in the call of this routine
       IF ( gamma_ice_int .LT. 1.0  .AND.  aaxis .LT. 100.0E-6 )  THEN !  Jensen and Harrington (2015)
          rho_deposition = rho_ice
       ELSE
          rho_deposition = rho_ice * EXP( -3.0 * MAX( drho * 1000.0 - 0.05, 0.0 ) / gamma_ice_int )  ! (42) in Chen and Lamb (1994, JAS)
       ENDIF
   
       mass_ini     = mass_int - dm_int        ! mass before timestep
       V_ini        = mass_ini / rho_apparent  ! volume before timestep
   
       IF ( dm_int .LT. 0.0 )  THEN ! sublimation: density does not change, Eqn. (23) in Shima et al. (2020)
          dV_ice = dm_int / rho_apparent
       ELSE                     ! deposition: considering current deposition density; Eqn. (41) in Chen and Lamb (1994, JAS)
          dV_ice = dm_int / rho_deposition
       ENDIF
       dlogV_ice = dV_ice / V_ini
   
       gamma_ast = gamma_ice_int * f_vent_axes_ice_int
   
       IF ( MIN( aaxis, caxis ) .LT. sizethresh_chen .AND. (dm_int .gt. 0.0))  THEN   ! Eqns. (28) - (29) in Shima et al. (2020)
          aaxis        = ( ( mass_ini + dm_int ) / ( 4.0 / 3.0 * pi * rho_ice ) )**0.333333333
          caxis        = ( ( mass_ini + dm_int ) / ( 4.0 / 3.0 * pi * rho_ice ) )**0.333333333
          rho_apparent = rho_ice
       ELSEIF (MIN(aaxis, caxis) .GE. sizethresh_chen .AND. (dm_int .gt. 0.0)) THEN
          aaxis        = aaxis * EXP( dlogV_ice / ( 2.0 + gamma_ast ) )                 ! Eqn. (26) in Shima et al. (2020)
          caxis        = caxis * EXP( dlogV_ice * ( gamma_ast / ( 2.0 + gamma_ast ) ) ) ! Eqn. (27) in Shima et al. (2020)
          rho_apparent = ( mass_ini + dm_int ) / ( V_ini + dV_ice )                     ! Eqn. (25) in Shima et al. (2020)
       ELSE !dm_int .le. 0.0
          aaxis = (mass_int / rho_apparent / (4.0/3.0 * pi * phi_pre))**0.333333333
          caxis = phi_pre * aaxis
       ENDIF
   
       RETURN
   
    END SUBROUTINE acd_update

 SUBROUTINE ice_ventilation(f_vent_mass_ice,f_vent_axes_ice,aaxis,caxis,w_sedi,capacitance_pre,rho_air_int,diffelq_local)
!
!-- Get ventilation coefficient (Eqns. (29) and (30) of Chen and Lamb (1994, JAS))
    USE params, ONLY: muelq

    IMPLICIT NONE

    REAL, INTENT(OUT) ::  f_vent_mass_ice,f_vent_axes_ice
    REAL, INTENT(IN)  ::  aaxis, caxis, w_sedi, capacitance_pre, rho_air_int, diffelq_local
    REAL              ::  N_Re, N_Sc, XX

    N_Re   = MIN( ( 2.0 * MAX(aaxis,caxis) ) * w_sedi / ( muelq / rho_air_int ), 120.0 ) ! Reynolds number
    N_Sc   = ( muelq / rho_air_int ) / diffelq_local                                       ! Schmidt number
    XX     = N_Sc**(1.0/3.0) * N_Re**(1.0/2.0)

    IF ( XX .LT. 1.0 )  THEN
       f_vent_mass_ice = 1.0 + 0.14 * XX**2
       f_vent_axes_ice = ( 1.0 + 0.14 * XX**2 * ( caxis / capacitance_pre )**0.5 ) / &
                         ( 1.0 + 0.14 * XX**2 * ( aaxis / capacitance_pre )**0.5 )
    ELSE
       f_vent_mass_ice = 0.86 + 0.28 * XX
       f_vent_axes_ice = ( 0.86 + 0.28 * XX * ( caxis / capacitance_pre )**0.5 ) / &
                         ( 0.86 + 0.28 * XX * ( aaxis / capacitance_pre )**0.5 )
    ENDIF

    RETURN

 END SUBROUTINE ice_ventilation

!
!-- Surface tension between liquid water and air (in J/m2)
REAL FUNCTION sigma_air_liq(tabs)

IMPLICIT NONE

REAL, INTENT(IN) :: tabs
REAL             :: tabs_c

tabs_c = tabs - 273.15
!
!--    Pruppacher and Klett (1997), Eq. 5-12
sigma_air_liq = 75.93 + 0.115 * tabs_c    + 6.818e-2 * tabs_c**2 +      &
                        6.511e-3 * tabs_c**3 + 2.933e-4 * tabs_c**4 +      &
                        6.283e-6 * tabs_c**5 + 5.285e-8 * tabs_c**6
sigma_air_liq = sigma_air_liq * 1.0E-3
!
! !--    Straka (2009)
!        sigma_air_liq = 0.0761 - 0.000155 * tabs_c

 END FUNCTION sigma_air_liq
!
!-- Surface tension between ice and liquid water (in J/m2)
 REAL FUNCTION sigma_ice_liq(tabs)

    IMPLICIT NONE

    REAL, INTENT(IN) :: tabs
! !
! !--    Pruppacher and Klett (1997), Eq. 5-47
!        sigma_ice_liq = ( 28.0 + 0.25 * ( tabs - 273.15 ) ) * 1.0E-3
! !
! !--    DeMott and Rogers (1990)
!        sigma_ice_liq = ( 28.0 + 0.167 * ( tabs - 273.15 ) ) * 1.0E-3
!
!--    Chukin et al. (2010) -> Yields most realistic results for homogeneous freezing
    sigma_ice_liq = ( 28.5 + 0.25 * ( tabs - 273.15 ) ) * 1.0E-3

 END FUNCTION sigma_ice_liq
