 SUBROUTINE micro_sgs_mixing(i,j,k)

    USE micro_vars,   ONLY: n_start, n_end, particle, rho_liq, dV, iWsgs, iT,   &
                            diss_LCM, iA, iM, iLEM, lcmp, dz_LCM, iilst, ijlst, &
                            iklst, nstep_LCM_init, tabs_LCM, pres_LCM, rho_air, &
                            dtn_LCM, gpzs, tk_LCM, switch_sedimentation, iq, iWsedi,    &
                            smix_LCM, start_sedimentation, ieta, iX, iY, iZ,    &
                            i_global_LCM, j_global_LCM
    USE microphysics, ONLY: micro_field, imv, ilastp, ilastT
    USE micro_random, ONLY: randomu

    USE params,       ONLY: pi, muelq, rgas, lcond, cp, rv, ggr, fac_cond
    USE vars,         ONLY: qsatw, esatw, dtqsatw, ttend, qtend, w, fluxbq,     &
                            fluxbt, rhow, gamaz
    USE grid,         ONLY: rank, nstep
    USE rad,          ONLY: qrad
    USE domain,       ONLY: nx_gl, ny_gl


    IMPLICIT NONE

    INTEGER, INTENT(IN) :: i, j, k

    INTEGER, SAVE :: numbers(1:27) = (/  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
                                        11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
                                        21, 22, 23, 24, 25, 26, 27 /)

    INTEGER :: ii, jj, kk, m, n, n_time, n_tsteps, n_last, dn_LEM, &
               numbers_cube(-1:1,-1:1,-1:1), n_rearrangements, n_rsteps, nn, mm

    REAL    :: T_now, q_now, a0, dt_diff, dz_sgs, D_mol, D_turb, L_turb, eta,  &
               lambda, T_adi, q_adi, a1, dt_turb, tau_sgs, tau_ndg, tau_sed,   &
               dql, dql_n, dV_n, T_ndg, q_ndg, t_box, q_box, t_box_int,        &
               q_box_int, N_sed, A_box, eta_ndg, d2phi

    REAL, PARAMETER :: Pr = 1.0 ! as in SAM, see tke_full.f90

    REAL,    DIMENSION(:), ALLOCATABLE :: T_new, T_old, q_new, q_old, eta_new, eta_old
    INTEGER, DIMENSION(:), ALLOCATABLE :: id, id_old

    CALL t_startf ('micro_sgs_mixing')

    T_now = tabs_LCM(i,j,k)
    q_now = micro_field(i,j,k,imv)

    n_last = n_end(i,j,k) - n_start(i,j,k) + 1

    IF ( n_last .GT. 10 )  THEN !there must be at least a couple of SDs inside the grid box
!
!--    Determine how much temperature and humidity change (restrict qrad to particles with liquid only?)
       DO n = n_start(i,j,k), n_end(i,j,k)

          particle(n)%var(iT) = particle(n)%var(iT) + ( ttend(k) + qrad(i,j,k) ) * dtn_LCM !+ ( -ggr / cp ) * ( particle(n)%var(iW) + particle(n)%var(iWsgs) ) * dtn_LCM &
          particle(n)%var(iq) = particle(n)%var(iq) + qtend(k) * dtn_LCM
          IF ( k .eq. 1 )  THEN ! surface fluxes -> replaced by strong nudging at surface
             particle(n)%var(iT) = particle(n)%var(iT) + fluxbt(i,j) * rhow(1) / ( rho_air(1) * dz_LCM(1) ) * dtn_LCM
             particle(n)%var(iq) = particle(n)%var(iq) + fluxbq(i,j) * rhow(1) / ( rho_air(1) * dz_LCM(1) ) * dtn_LCM
          ENDIF

       ENDDO
!
!--    Set turbulence parameters
       dz_sgs = dz_LCM(1) / n_last
       eta    = 6.0 * dz_sgs

       L_turb = MAX( MIN( smix_LCM(i,j,k), dz_LCM(1) ), 0.001 ) ! L is the largest eddy size used in the LEM

       D_turb = MAX( Pr * tk_LCM(i,j,k), 0.0 )                                      ! turbulent diffusivity from LES SGS model (SAM uses Pr = 1.0)
       D_mol  = MAX( D_turb * ( eta / L_turb )**1.33333333, muelq / rho_air(k) )    ! scale molecular diffusion to take over unresolved eddies
       lambda = 54.0 / 5.0 * D_turb / L_turb**3 * ( L_turb / eta )**1.66666666      ! frequency of rearrangement events (Eq. 2.4, Krueger, 1993, JAS)

       dt_diff  = 0.2 * dz_sgs**2 / D_mol
       n_tsteps = CEILING( dtn_LCM / dt_diff ) ! number of sub-time steps
       dt_diff  = dtn_LCM / n_tsteps

       dt_turb  = ( lambda * dz_LCM(1) )**(-1)
       n_rsteps = CEILING( dt_diff / dt_turb ) ! number of rearrangements per sub-time step
       dt_turb  = dt_diff / n_rsteps
!
!--    Check if the time step is too long. If this is the case, turbulent rearrangements
!--    might be underestimated.
!       IF ( dt_diff / ( lambda * dz_LCM )**(-1) .GT. 1.0 )  THEN
!          PRINT*, '+++ time step >>dt_diff<< too long (turb. rearrangement) +++'
!          PRINT*, dt_diff, 1.0 / ( lambda * dz_LCM )
!          STOP
!       ENDIF

!
!--    Create arrays id and S_old based on the origins of all SDs inside the current
!--    grid box. First, rearrange numbers randomly and store them on a cube with
!--    sidelength 3. This is necessary to place entrained SDs at different locations
!--    at any timestep.
       CALL shuffle(numbers,27)
       numbers_cube = RESHAPE(numbers,(/3,3,3/))

       ALLOCATE( T_new(n_start(i,j,k):n_end(i,j,k)),   &
                 T_old(n_start(i,j,k):n_end(i,j,k)),   &
                 q_new(n_start(i,j,k):n_end(i,j,k)),   &
                 q_old(n_start(i,j,k):n_end(i,j,k)),   &
                 eta_new(n_start(i,j,k):n_end(i,j,k)), &
                 eta_old(n_start(i,j,k):n_end(i,j,k)), &
                 id(n_start(i,j,k):n_end(i,j,k)),      &
                 id_old(n_start(i,j,k):n_end(i,j,k)) )

       DO n = n_start(i,j,k), n_end(i,j,k)
          ii = particle(n)%ivar(iilst) - ( i + i_global_LCM )
          IF ( ii .GT. nx_gl / 2 )  THEN
             ii = ii - nx_gl
          ELSEIF ( ii .LE. -nx_gl / 2 )  THEN
             ii = ii + nx_gl
          ENDIF
          jj = particle(n)%ivar(ijlst) - ( j + j_global_LCM )
          IF ( jj .GT. ny_gl / 2 )  THEN
             jj = jj - ny_gl
          ELSEIF ( jj .LE. -ny_gl / 2 )  THEN
             jj = jj + ny_gl
          ENDIF
          kk = particle(n)%ivar(iklst) - k
          IF ( ABS(ii) .GT. 1  .OR.  ABS(jj) .GT. 1  .OR.  ABS(kk) .GT. 1 )  THEN
             PRINT*, "LCM ERROR: Particle moves too fast (micro_sgs_mixing.f90)."
             PRINT*, "           This is probably due to too strong collisional "
             PRINT*, "           growth. Consider running the model without the "
             PRINT*, "           SGS model (switch_LEM_micro = .FALSE.)."
             ! PRINT*, rank
             ! PRINT*, i+i_global_LCM, j+j_global_LCM, k
             ! PRINT*, particle(n)%ivar(iilst), particle(n)%ivar(ijlst), particle(n)%ivar(iklst)
             ! PRINT*, ii, jj, kk
             ! PRINT*, nx_gl, ny_gl
             ! PRINT*, i_global_LCM, j_global_LCM
             STOP
          ENDIF
!
!--       This works if the number of particles per grid box is smaller than
!--       100000 (which will probably be the case)
          id(n)      = numbers_cube(ii,jj,kk) * 100000 + particle(n)%ivar(iLEM)
          id_old(n)  = n
          T_old(n)   = particle(n)%var(iT)
          q_old(n)   = particle(n)%var(iq)
          eta_old(n) = particle(n)%var(ieta)
       ENDDO
!
!--    Order the SDs regarding their original grid box and their location in the
!--    original LEM array (i.e., if two SDs have been entrained which were previously
!--    neighbors in the LEM array, they will still be neighbors)
       CALL insertionsort( n_last, id, id_old, T_old, q_old, eta_old )
!
!--    Assign new LEM ids
       DO n = n_start(i,j,k), n_end(i,j,k)
          id(n) = id_old(n)
       ENDDO

       DO n_time = 1, n_tsteps
!
!--       Compute molecular (or, if dz_sgs is too large, eddy) diffusion
          d2phi                 = T_old(n_start(i,j,k)+1) - 2.0 * T_old(n_start(i,j,k)) + T_old(n_end(i,j,k))
          T_new(n_start(i,j,k)) = T_old(n_end(i,j,k)) + dt_diff * D_mol * d2phi / dz_sgs**2
          DO n = n_start(i,j,k)+1, n_end(i,j,k)-1
             d2phi    = T_old(n+1) - 2.0 * T_old(n) + T_old(n-1)
             T_new(n) = T_old(n) + dt_diff * D_mol * d2phi / dz_sgs**2
          ENDDO
          d2phi                 = T_old(n_start(i,j,k)) - 2.0 * T_old(n_end(i,j,k)) + T_old(n_end(i,j,k)-1)
          T_new(n_end(i,j,k))   = T_old(n_end(i,j,k)) + dt_diff * D_mol * d2phi / dz_sgs**2

          d2phi                 = q_old(n_start(i,j,k)+1) - 2.0 * q_old(n_start(i,j,k)) + q_old(n_end(i,j,k))
          q_new(n_start(i,j,k)) = q_old(n_start(i,j,k)) + dt_diff * D_mol * d2phi / dz_sgs**2
          DO n = n_start(i,j,k)+1, n_end(i,j,k)-1
             d2phi    = q_old(n+1) - 2.0 * q_old(n) + q_old(n-1)
             q_new(n) = q_old(n) + dt_diff * D_mol * d2phi / dz_sgs**2
          ENDDO
          d2phi               = q_old(n_start(i,j,k)) - 2.0 * q_old(n_end(i,j,k)) + q_old(n_end(i,j,k)-1)
          q_new(n_end(i,j,k)) = q_old(n_end(i,j,k)) + dt_diff * D_mol * d2phi / dz_sgs**2

          d2phi                   = eta_old(n_start(i,j,k)+1) - 2.0 * eta_old(n_start(i,j,k)) + eta_old(n_end(i,j,k))
          eta_new(n_start(i,j,k)) = eta_old(n_start(i,j,k)) + dt_diff * D_mol * d2phi / dz_sgs**2
          DO n = n_start(i,j,k)+1, n_end(i,j,k)-1
             d2phi      = eta_old(n+1) - 2.0 * eta_old(n) + eta_old(n-1)
             eta_new(n) = eta_old(n) + dt_diff * D_mol * d2phi / dz_sgs**2
          ENDDO
          d2phi                 = eta_old(n_start(i,j,k)) - 2.0 * eta_old(n_end(i,j,k)) + eta_old(n_end(i,j,k)-1)
          eta_new(n_end(i,j,k)) = eta_old(n_end(i,j,k)) + dt_diff * D_mol * d2phi / dz_sgs**2
!
!--       Compute probability for turbulent rearrangement (keep in mind that
!--       ( lambda * dx_total )**(-1) is the average time between two rearrangement
!--       events)
          DO n_rearrangements = 1, n_rsteps
             IF ( dt_diff / ( lambda * dz_LCM(1) )**(-1) .GT. randomu(rank+1) )  THEN
                CALL micro_sgs_rearrangement( id, T_new, q_new, eta_new, n_last, dz_sgs, eta, L_turb )
             ENDIF
          ENDDO

          T_old   = T_new
          q_old   = q_new
          eta_old = eta_new

       ENDDO
!
!--    After molecular diffusion and turbulent rearrangement is finished, the
!--    new values of T, q, and eta are stored to the particles
       DO m = n_start(i,j,k), n_end(i,j,k)
          DO n = n_start(i,j,k), n_end(i,j,k)
             IF ( m .EQ. id(n) )  THEN
                particle(m)%var(iT)    = T_new(n)
                particle(m)%var(iq)    = q_new(n)
                particle(m)%var(ieta)  = eta_new(n)
                particle(m)%ivar(iLEM) = n - n_start(i,j,k) + 1!this ensures that the SDs will maintain their order in the LEM array over several timesteps
                EXIT
             ENDIF
             IF ( n .EQ. n_end(i,j,k) )  THEN
                PRINT*, '>>>>SEVERE ERROR<<<<<'
                STOP
             ENDIF
          ENDDO
          particle(m)%ivar(iilst) = i + i_global_LCM
          particle(m)%ivar(ijlst) = j + j_global_LCM
          particle(m)%ivar(iklst) = k
       ENDDO

       DEALLOCATE( T_new, T_old, q_new, q_old, eta_new, eta_old, id, id_old )
!
!--    Nudging: ensure that SGS terms do not deviate too much from grid box mean
       dql = SUM( particle(n_start(i,j,k):n_end(i,j,k))%var(iM) ) / ( dV(k) * rho_air(k) )
       DO n = n_start(i,j,k), n_end(i,j,k)

          dV_n  = dV(k) / n_last
          dql_n = particle(n)%var(iM) / ( dV_n * rho_air(k) )

          tau_sgs = 900.0
          tau_sed = ( dz_LCM(1) / n_last ) / ( particle(n)%var(iWsedi) + 1.0E-20 )

          tau_ndg = 1.0 / ( 1.0 / tau_sed + 1.0 / tau_sgs )
!
!--       Nudging is only possible to a conserved variable. Therefore, iT and iq are nudged to
!--       the liquid water static energy and total water content.
          T_ndg   = particle(n)%var(iT)   - ( ( particle(n)%var(iT) - fac_cond * dql_n ) - ( tabs_LCM(i,j,k) + gamaz(k) - fac_cond * dql ) ) * dtn_LCM / tau_ndg
          q_ndg   = particle(n)%var(iq)   - ( ( particle(n)%var(iq) + dql_n ) - ( micro_field(i,j,k,imv) + dql ) ) * dtn_LCM / tau_ndg
          eta_ndg = particle(n)%var(ieta) - ( particle(n)%var(ieta) - ( micro_field(i,j,k,imv) - qsatw(micro_field(i,j,k,ilastT), micro_field(i,j,k,ilastp)) ) ) * dtn_LCM / tau_ndg

          particle(n)%var(iT)   = REAL( T_ndg,   KIND=lcmp )
          particle(n)%var(iq)   = REAL( q_ndg,   KIND=lcmp )
          particle(n)%var(ieta) = REAL( eta_ndg, KIND=lcmp )

       ENDDO

    ELSE
 !
 !--   If there are too less particles inside a grid box, the LEM is not applicable.
 !--   Therefore, the values are set to the grid box values immediately (similar to nudging with a timescale of 0 secs).
       dql = SUM( particle(n_start(i,j,k):n_end(i,j,k))%var(iM) ) / ( dV(k) * rho_air(k) )
       DO n = n_start(i,j,k), n_end(i,j,k)

          dV_n  = dV(k) / n_last
          dql_n = particle(n)%var(iM) / ( dV_n * rho_air(k) )

          particle(n)%var(iT)   = REAL( tabs_LCM(i,j,k) + gamaz(k), KIND=lcmp )!REAL( tabs_LCM(i,j,k), KIND=lcmp )
          particle(n)%var(iq)   = REAL( micro_field(i,j,k,imv) + dql - dql_n, KIND=lcmp )
          particle(n)%var(ieta) = REAL( micro_field(i,j,k,imv) - qsatw( micro_field(i,j,k,ilastT), micro_field(i,j,k,ilastp)), KIND=lcmp )

          particle(n)%ivar(iLEM)  = n - n_start(i,j,k) + 1
          particle(n)%ivar(iilst) = i + i_global_LCM
          particle(n)%ivar(ijlst) = j + j_global_LCM
          particle(n)%ivar(iklst) = k

       ENDDO

    ENDIF

    CALL t_stopf ('micro_sgs_mixing')

 END SUBROUTINE micro_sgs_mixing

 SUBROUTINE micro_sgs_rearrangement( id, T_old, q_old, eta_old, n_last, dz_sgs, eta, L_turb )
!
!-- To be called from micro_sgs_mixing. Calculates turbulent rearrangements in
!-- the LEM.
    USE micro_vars,   ONLY: dz_LCM
    USE micro_random, ONLY: randomu

    USE grid,         ONLY: rank
    USE params,       ONLY: ggr, cp

    IMPLICIT NONE

    INTEGER, INTENT(IN)    ::  n_last
    INTEGER, INTENT(INOUT) ::  id(n_last)
    REAL,    INTENT(IN)    ::  dz_sgs, eta, L_turb
    REAL,    INTENT(INOUT) ::  T_old(n_last), q_old(n_last), eta_old(n_last)

    INTEGER                ::  n_start, n_end, n_select, n, m, n_length
    REAL                   ::  length
!
!-- The rearrangement only works if eta is smaller than L(=dz_LCM), which is usually
!-- the case unless there are 6 or less SDs inside a grid box.
    IF ( n_last .LE. 6 )  RETURN
!
!-- If L_turb is smaller than eta, SGS turbulence cannot be resolved by LEM
    IF ( L_turb .LT. eta )  RETURN
!
!-- Randomly select starting position of rearrangement
    n_start = INT( randomu(rank+1) * REAL( n_last ) ) + 1
!
!-- Select randomly the segment's length from Eq. 2.3 (Krueger, 1993, JAS)
    IF ( L_turb / eta .LT. 1.01 )  THEN
       length = eta
    ELSE
       length = -99999999.99
       DO WHILE ( length .LT. eta  .OR.  length .GT. MIN( L_turb, dz_LCM(1) ) )
          length = ( eta**(-1.6666666666) -                                               &
                     randomu(rank+1) * ( eta**(-1.6666666666) - L_turb**(-1.6666666666) ) &
                   )**(-0.6)
       ENDDO
    ENDIF
    n_length = FLOOR( length / dz_sgs )
    n_length = MIN( MAX( NINT( n_length / 3.0 ) * 3, 6 ), FLOOR( n_last / 3.0 ) * 3 ) ! make sure that n_length is a multiple of 3
!
!-- Determine end position of rearrangement
    n_end = MODULO( n_start + n_length - 1, n_last )
    IF ( n_end .EQ. 0 ) n_end = n_last
!
!-- Start the remapping using the triplet map (see p. 3079 in Krueger (1993, JAS)):
    triplet_map: BLOCK

       TYPE rearrangement_type
          INTEGER :: id
          REAL    :: T, q, eta
       END TYPE rearrangement_type

       TYPE(rearrangement_type), DIMENSION(1:n_length) :: mapping, mapping_final
       LOGICAL, SAVE                                   :: switch_supersat_fluct_LEM = .FALSE.

       INTEGER :: mm, nn
!
!--    Store the selected super-droplets on a temporary array
       IF ( n_start .LE. n_end )  THEN
          mapping(1:n_length)%id  = id(n_start:n_end)
          mapping(1:n_length)%T   = T_old(n_start:n_end)
          mapping(1:n_length)%q   = q_old(n_start:n_end)
          mapping(1:n_length)%eta = eta_old(n_start:n_end)
       ELSE
          mapping(1:n_last-n_start+1)%id         = id(n_start:n_last)
          mapping(1:n_last-n_start+1)%T          = T_old(n_start:n_last)
          mapping(1:n_last-n_start+1)%q          = q_old(n_start:n_last)
          mapping(1:n_last-n_start+1)%eta        = eta_old(n_start:n_last)
          mapping(n_last-n_start+2:n_length)%id  = id(1:n_end)
          mapping(n_last-n_start+2:n_length)%T   = T_old(1:n_end)
          mapping(n_last-n_start+2:n_length)%q   = q_old(1:n_end)
          mapping(n_last-n_start+2:n_length)%eta = eta_old(1:n_end)
       ENDIF
!
!--    Relocation of particles (previous random version below): see parcel model with n_start = 1
       DO n = 1, n_length

          IF ( n .GE. 1  .AND.  n .LE. NINT( n_length / 3.0 ) )  THEN

             m = 3 * ( n - 1 ) + 1

          ELSEIF ( n .GT. NINT( n_length / 3.0 )  .AND.  n .LE. NINT( 2.0 * n_length / 3.0 ) )  THEN

             m = 2 * ( n_length - 1 ) - 3 * ( n - 1 ) + 1

          ELSEIF ( n .GT. NINT( 2.0 * n_length / 3.0 )  .AND.  n .LE. n_length )  THEN

             m = 3 * ( n - 1 ) - 2 * ( n_length - 1 ) + 1

          ENDIF

          mapping_final(n) = mapping(m)
!
!--       Add supersaturation fluctuations
          IF ( switch_supersat_fluct_LEM )  THEN
             mapping_final(n)%T = mapping_final(n)%T - ( n - m ) * dz_sgs * ggr / cp
          ENDIF

       ENDDO
!
!--    Finally, place the rearranged superdroplets into the original array
       IF ( n_start .LE. n_end )  THEN
          id(n_start:n_end)       = mapping_final(1:n_length)%id
          T_old(n_start:n_end)    = mapping_final(1:n_length)%T
          q_old(n_start:n_end)    = mapping_final(1:n_length)%q
          eta_old(n_start:n_end)  = mapping_final(1:n_length)%eta
       ELSE
          id(n_start:n_last)      = mapping_final(1:n_last-n_start+1)%id
          T_old(n_start:n_last)   = mapping_final(1:n_last-n_start+1)%T
          q_old(n_start:n_last)   = mapping_final(1:n_last-n_start+1)%q
          eta_old(n_start:n_last) = mapping_final(1:n_last-n_start+1)%eta
          id(1:n_end)             = mapping_final(n_last-n_start+2:n_length)%id
          T_old(1:n_end)          = mapping_final(n_last-n_start+2:n_length)%T
          q_old(1:n_end)          = mapping_final(n_last-n_start+2:n_length)%q
          eta_old(1:n_end)        = mapping_final(n_last-n_start+2:n_length)%eta
       ENDIF

    END BLOCK triplet_map

 END SUBROUTINE micro_sgs_rearrangement

 SUBROUTINE insertionsort( n_last, id, id_old, T_old, q_old, eta_old )

    IMPLICIT NONE

    INTEGER :: n, m, temp_id, temp_id_old
    REAL    :: temp_T, temp_q, temp_eta

    INTEGER, INTENT(IN)    :: n_last
    INTEGER, INTENT(INOUT) :: id(1:n_last), id_old(1:n_last)
    REAL,    INTENT(INOUT) :: T_old(1:n_last), q_old(1:n_last), eta_old(1:n_last)
!
!-- Perform an insertion sort
    DO n = 2, n_last
       m = n - 1
       temp_id     = id(n)
       temp_id_old = id_old(n)
       temp_T      = T_old(n)
       temp_q      = q_old(n)
       temp_eta    = eta_old(n)
       DO WHILE ( id(m) > temp_id )!.AND. m >= 1 )
           id(m+1)      = id(m)
           id_old(m+1)  = id_old(m)
           T_old(m+1)   = T_old(m)
           q_old(m+1)   = q_old(m)
           eta_old(m+1) = eta_old(m)
           m = m - 1
           IF ( m == 0 )  EXIT
       END DO
       id(m+1)      = temp_id
       id_old(m+1)  = temp_id_old
       T_old(m+1)   = temp_T
       q_old(m+1)   = temp_q
       eta_old(m+1) = temp_eta
    ENDDO

 END SUBROUTINE insertionsort

!  SUBROUTINE shellsort( n_last, id, S_old )
!
!     IMPLICIT NONE
!
!     INTEGER :: inc, n, m, temp_id
!     REAL    :: temp_S
!
!     INTEGER, INTENT(IN)    :: n_last
!     INTEGER, INTENT(INOUT) :: id(1:n_last)
!     REAL,    INTENT(INOUT) :: S_old(1:n_last)
!
!     inc = 1
! 1   inc = 3 * inc + 1
!     IF ( inc .LE. n_last )  GOTO 1
! 2   CONTINUE
!     inc = inc / 3
!     DO n = inc + 1, n_last
!        temp_id = id(n)
!        temp_S = S_old(n)
!        m = n
! 3      IF ( id(m-inc) .gt. temp_id )  THEN
!           id(m)    = id(m-inc)
!           S_old(m) = S_old(m-inc)
!           m = m - inc
!           IF ( m .LE. inc )  GOTO 4
!           GOTO 3
!        ENDIF
! 4      id(m)    = temp_id
!        S_old(m) = temp_S
!     ENDDO
!     IF ( inc .GT. 1 )  GOTO 2
!
!  END SUBROUTINE shellsort

 SUBROUTINE shuffle( a, n )

    USE micro_random, ONLY: randomu

    USE grid,         ONLY: rank

    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: a(1:n)
    INTEGER, INTENT(IN)    :: n
    INTEGER                :: i, randpos, temp
!
!-- Performs a Fisher-Yates shuffle
    DO i = n, 2, -1
       randpos = MIN( MAX( INT(randomu(rank+1) * i) + 1, 1 ), n )
       temp = a(randpos)
       a(randpos) = a(i)
       a(i) = temp
    ENDDO

 END SUBROUTINE shuffle
