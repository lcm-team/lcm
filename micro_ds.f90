MODULE micro_ds
    CONTAINS    
    !Drop shattering from heterogeneous freezing parameterized by Phillips et al. (2018)
    SUBROUTINE ds_phillips(i,j,k, n, temp, nfrag_tot, nfrag_b)
    USE params, ONLY: pi, lfus
    USE micro_vars, ONLY: rho_liq, rho_ice, n_start, n_end, particle, ds_prof, &
                            iA, iTYPE, iM, iAXS, iCXS, iDNS
    USE microphysics, ONLY: ids, micro_field
    USE micro_random, ONLY: randomu
    USE grid, ONLY: rank
    
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: i,j,k,n
    REAL, INTENT(IN)     :: temp !Diameter of the drop and temperature of the gridbox
    INTEGER, INTENT(OUT) :: nfrag_tot, nfrag_b
    INTEGER              :: nfrag_s !Number of 'large' and 'small' fragments which are emitted from the drop during freezing
    REAL                 :: chi, logdiam, beta, logzeta, zeta, logeta, eta_ds, t0_ds, &
                            zeta_big, eta_big, t0_big, xi, omega, tempc, diam !Variables used in droplet fragment number parameterization
    !REAL, PARAMETER      :: cpl = 4200 !J/kg/K, specific heat of water as defined in Phillips et al. (2018), Appendix B
    REAL                 :: realntot, realnbig, x_rand !Probabilistically assign fractions of fragments
    REAL                :: mindiff, tinymass, phi !Find difference in size between small ice fragments and other ice particles in the grid box
    INTEGER             :: MinN, p, oldmultiplicity

    nfrag_tot = 0
    nfrag_b = 0
    nfrag_s = 0
    tempc = temp-273.15
    diam = (6.*particle(n)%var(iM)/(pi*particle(n)%ivar(iA)*rho_liq))**(1./3.)*1000. !Diameter of the freezing drop in mm, calculated only from liquid mass (drops are large enough that aerosol size should be irrelevant)
    !First, find dependency of several variables on droplet diameter (in mm)
    xi = (MIN(diam, 0.06) - 0.05)/0.01 !Linearly interpolate droplet fragmentation number between 50um and 60um droplet diameters as in Eq. 1
    omega = (MAX(tempc, -6.) + 3.)/(-3.) !Linearly interpolate droplet fragmentation number based on temperature between -6 C and -3 C as in Eq. 1
    oldmultiplicity = particle(n)%ivar(iA)
    IF (diam .lt. 0.4) THEN
        chi = log10(MAX(diam, 0.06))
        beta = 0
        logzeta = 2.4268*chi**3 + 3.3274*chi**2 + 2.0783*chi + 1.2927
        logeta = 0.1242*chi**3 - 0.2316*chi**2 - 0.9874*chi - 0.0827
        t0_ds = -1.3999*chi**3 - 5.3285*chi**2 - 3.9847*chi - 15.0332
        zeta_big = -0.4651*chi**3 - 1.1072*chi**2 - 0.4539*chi + 0.5137
        eta_big = 28.5888*chi**3 + 49.8504*chi**2 + 22.4873*chi + 8.0481
        t0_big = 13.3588*chi**3 + 15.7432*chi**2 - 2.6543*chi - 18.4875
    ELSE
        chi = log10(MIN(diam, 1.6))
        beta = -0.1839*chi**2 - 0.2017*chi - 0.0512
        logzeta = 2.4268*chi**3 + 3.3274*chi**2 + 2.0783*chi + 1.2927
        logeta = 0.1242*chi**3 - 0.2316*chi**2 - 0.9874*chi - 0.0827
        t0_ds = -1.3999*chi**3 - 5.3285*chi**2 - 3.9847*chi - 15.0332
        zeta_big = -0.4651*chi**3 - 1.1072*chi**2 - 0.4539*chi + 0.5137
        eta_big = 28.5888*chi**3 + 49.8504*chi**2 + 22.4873*chi + 8.0481
        t0_big = 13.3588*chi**3 + 15.7432*chi**2 - 2.6543*chi - 18.4875

    ENDiF

    zeta = 10.**logzeta
    eta_ds = 10.**logeta

    realntot = (xi*omega)*((zeta*eta_ds**2)/((tempc-t0_ds)**2 + eta_ds**2) + beta*tempc)
    nfrag_tot = FLOOR(realntot)
    IF (randomu(rank+1) .lt. (realntot - nfrag_tot)) nfrag_tot = nfrag_tot + 1
    realnbig = MIN((xi*omega)*((zeta_big*eta_big**2)/((tempc-t0_big)**2 + eta_big**2)), realntot)
    nfrag_b = FLOOR(realnbig)
    IF (randomu(rank+1) .lt. (realnbig - nfrag_b)) nfrag_b = nfrag_b + 1

    END SUBROUTINE ds_phillips

    SUBROUTINE ds_lachapelle(i,j,k,n,temp,nfrag_tot) !Parameterization based on LaChapelle et al. (2024)
    USE params, ONLY: pi
    USE micro_vars, ONLY: rho_liq, rho_ice, particle, ds_prof, &
                            iA, iM
    USE microphysics, ONLY: ids, micro_field
    USE micro_random, ONLY: randomu
    USE grid, ONLY: rank

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,j,k,n
    REAL, INTENT(IN) :: temp
    INTEGER, INTENT(OUT) :: nfrag_tot

    REAL :: diam, nfrags, tempc
    !Diam is the spherical equivalent diameter of droplets, in micrometers
    !nfrags is the number of fragments produced by a shattering event. In this scheme, they're assumed to be monodispersive in size
    nfrag_tot = 0
    diam = 1.0E6*((6*particle(n)%var(iM)/(pi*rho_liq*particle(n)%ivar(iA)))**(1./3.))
    tempc = temp-273.16 !Temperature of the gridbox in celsius
    IF (diam .ge. 100.0) THEN!Only do fragmentation for droplets larger than 100 micron diameter
        nfrags = 2.5E-11 * (MIN(diam, 3500.))**4 !Eq. 1 of Lachapelle et al. (2024), limiting maximum size for fragment calculation to 3500 microns
        nfrags = nfrags * (MAX(tempc, -6.) + 3.)/(-3.)
        nfrag_tot = FLOOR(nfrags)
    ENDIF
    
    END SUBROUTINE ds_lachapelle

    SUBROUTINE ds_lauber(i,j,k,n,temp,nfrag_tot) !Parameterization based on Table 1 in Lauber et al. (2021)
    USE params, ONLY: pi
    USE micro_vars, ONLY: rho_liq, rho_ice, particle, ds_prof, &
                            iA, iM
    USE microphysics, ONLY: ids, micro_field
    USE micro_random, ONLY: randomu
    USE grid, ONLY: rank

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,j,k,n
    REAL, INTENT(IN) :: temp
    INTEGER, INTENT(OUT) :: nfrag_tot

    REAL :: diam, nfrags, psh, normfac, tempc
    tempc = temp-273.16 !Temperature of gridbox in degrees celsius
    !Diam is the spherical equivalent diameter of droplets, in meters
    !nfrags is the number of fragments produced by a shattering event. In this scheme, they're assumed to be monodispersive in size
    nfrag_tot = 0
    diam = ((6*particle(n)%var(iM)/(pi*rho_liq*particle(n)%ivar(iA)))**(1./3.))
    IF (diam .ge. 40.0E-6) THEN!Only do fragmentation for droplets larger than 40 micron diameter
        psh = 4.4E6*diam**2 !Eq. 4 of Lauber et al. (2021)
        IF (tempc .lt. -20.0) THEN
            psh = psh*(1 + (tempc+20)/5) !Scale fragmentation probability linearly between -20C and -25 based on Keinert et al. (2020)
        ENDIF
        nfrags = 9.0E4*diam*psh !Eq. 5 of Lauber et al. (2021)
        nfrag_tot = FLOOR(nfrags)
        !Note: Though Lauber et al. (2021) does also explore a parameterization with psh = 1 and nfrags = 2.2*10**5*diam, I've elected to stick to the first
        !estimates, because Keinert et al. (2020) did not show such high breakup probabilities even in moving airflow, and such high fragment numbers are at odds with previous studies.
    ENDIF

    END SUBROUTINE ds_lauber
END MODULE micro_ds

