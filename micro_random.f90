 MODULE micro_random

    USE micro_vars
    USE grid, ONLY: rank

 CONTAINS

 FUNCTION num_gauss( sigma ) RESULT( ng )
!
!-- Returns random number of the gauss distribution with zero mean and variance of sigma^2.
!-- Number within two standard deviation (2*sigma) is returned.

    IMPLICIT NONE

    REAL, INTENT(IN) :: sigma
    REAL             :: ng
    REAL, PARAMETER  :: nsigma = 2

    DO
       ng = gasdev( rank+1, sigma )
       IF ( abs(ng) .LE. nsigma * sigma ) EXIT
    ENDDO

 END FUNCTION num_gauss

 REAL FUNCTION gasdev(idum, sigma)
!
!-- gasdev_s in "Numerical Recipes in Fortran 90: the Art of Parallel Scientific Computing"
!-- gasdev in "Numerical Recipes in Fortran 77: the Art of Scientific Computing"
!-- Modified from subroutine to function. Returns a normally distributed deviate with zero
!-- mean and variance of sigma, using ran(idum) as the source of unit deviates.
!-- Input idum has to be positive. For LCM, idum = rank+1

    IMPLICIT NONE

    INTEGER, PARAMETER       :: K4B=selected_int_kind(9)
    INTEGER(K4B), INTENT(IN) :: idum
    REAL, INTENT(IN)         :: sigma   ! mean and standard deviation
    REAL                     :: rsq, v1, v2
    REAL, SAVE               :: gnext
    LOGICAL, SAVE            :: gaus_stored = .FALSE.

    IF (gaus_stored) THEN   ! We have an extra deviate handy, so return it, and unset the flag.
       gasdev = gnext
       gaus_stored = .FALSE.
    ELSE   ! We don't have an extra deviate handy, so pick two uniform numbers
       DO
          v1 = 2.0 * randomu(idum) - 1.0
          v2 = 2.0 * randomu(idum) - 1.0
          rsq = v1**2 + v2**2
          IF (rsq > 0.0 .AND. rsq < 1.0) EXIT   ! See if they are in the unit circle,
       ENDDO                                    ! otherwise try again
       rsq = SQRT(-2.0 * log(rsq) / rsq)        ! Now make the Box-Muller transformation to get
       gasdev = v1 * rsq * sigma                ! two normal deviates. Return one and save the
       gnext  = v2 * rsq * sigma                ! other for next time.
       gaus_stored = .TRUE.                     ! Set flag
    ENDIF

 END FUNCTION gasdev

 REAL FUNCTION randomu(idum)
!
!-- "ram" in "Numerical Recipes in Fortran 90: the Art of Parallel Scientific Computing"
!-- "Minimal" random number generator of Park and Miller combined with a Marsaglia shift sequence.
!-- Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values).
!-- This fully portable, scalar generator has the "traditional" (not F90) calling sequence with a
!-- random deviate as the returned function value: call with idum a nagative integer to initialize;
!-- thereafter, do not alter idum except to reinitialize. The period of this generator is about
!-- 3.1x10^18
    IMPLICIT NONE

    INTEGER, PARAMETER          :: K4B = SELECTED_INT_KIND(9)
    ! INTEGER(K4B), INTENT(INOUT) :: idum
    INTEGER(K4B), INTENT(IN)    :: idum
    INTEGER(K4B), PARAMETER     :: IA=16870, IM=2147483647, IQ=127773, IR=2836
    REAL, SAVE                  :: am
    INTEGER(K4B), SAVE          :: ix=-1, iy=-1, k

    IF (idum <= 0 .OR. iy < 0) THEN   ! initialize
       am = nearest(1.0,-1.0) / IM
       iy = ior(ieor(888889999,abs(idum)),1)
       ix = ieor(777755555,abs(idum))
       !!idum = abs(idum) + 1   ! set idum positive
    ENDIF
    ix = ieor(ix,ishft(ix,13))   ! Marsaglia shift sequence with period 2^23 - 1
    ix = ieor(ix,ishft(ix,-17))
    ix = ieor(ix,ishft(ix,5))
    k = iy / IQ   ! Park-Miller sequence by Schrage's method, period 2^31 - 2
    iy = IA * (iy - k * IQ) - IR * k
    IF (iy < 0) iy = iy + IM
    randomu = am * ior(iand(IM,ieor(ix,iy)),1)   ! Combine the two generators with masking to
                                                 ! ensure nonzero value
 END FUNCTION randomu

 END MODULE micro_random
