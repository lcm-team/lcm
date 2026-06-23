 SUBROUTINE micro_user_init

    USE grid, ONLY: rank, nrestart, masterproc, case, caseid

    IMPLICIT NONE

    INTEGER                :: stat, k
    CHARACTER(LEN=7), SAVE :: rank_char

    WRITE(rank_char,'(''_'',I6.6)') rank
    ! IF ( nrestart .EQ. 0 )  THEN
    !    OPEN(UNIT=100,FILE='./OUT_USER/particle_ana'//TRIM(rank_char)//'.dat',STATUS='REPLACE',FORM='UNFORMATTED')
    ! ELSE
    !    OPEN(UNIT=100,FILE='./OUT_USER/particle_ana'//TRIM(rank_char)//'.dat',STATUS='OLD',FORM='UNFORMATTED',POSITION='APPEND',IOSTAT=stat)
    !    IF ( stat .NE. 0 )  THEN
    !       OPEN(UNIT=100,FILE='./OUT_USER/particle_ana'//TRIM(rank_char)//'.dat',STATUS='NEW',FORM='UNFORMATTED')
    !    ENDIF
    ! ENDIF

    ! IF ( nrestart .EQ. 0 )  THEN
    !    OPEN(UNIT=101,FILE='./OUT_USER/particle_track'//TRIM(rank_char)//'.dat',STATUS='REPLACE',FORM='UNFORMATTED')
    ! ELSE
    !    OPEN(UNIT=101,FILE='./OUT_USER/particle_track'//TRIM(rank_char)//'.dat',STATUS='OLD',FORM='UNFORMATTED',POSITION='APPEND',IOSTAT=stat)
    !    IF ( stat .NE. 0 )  THEN
    !       OPEN(UNIT=101,FILE='./OUT_USER/particle_track'//TRIM(rank_char)//'.dat',STATUS='NEW',FORM='UNFORMATTED')
    !    ENDIF
    ! ENDIF

 END SUBROUTINE micro_user_init

 SUBROUTINE micro_user_init_particles(i,j,k)

    USE micro_vars,      ONLY: particle, n_start, n_end, lcmp
    USE micro_user_vars

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: i, j, k
    INTEGER             :: n

    ! DO n = n_start(i,j,k), n_end(i,j,k)
    !    particle(n)%var(...) = ...
    ! ENDDO

 END SUBROUTINE micro_user_init_particles

 SUBROUTINE micro_user_diagnose

    USE micro_vars,      ONLY: n_start, n_end, iM, iA, rho_liq, particle, get_r_particle
    USE micro_user_vars

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    INTEGER :: i, j, k, n
    REAL    :: radius

!     DO k = 1, nzm
!        DO j = 1, ny
!           DO i = 1, nx
! !
! !--          Compute diagnostic fields
!              IF ( n_end(i,j,k) - n_start(i,j,k) .LT. 0 )  CYCLE !if there are no particles inside that gridbox, do nothing
!
!              DO  n = n_start(i,j,k), n_end(i,j,k)
!                !radius = ( particle(n)%var(iM) / ( 4.0 / 3.0 * pi * particle(n)%ivar(iA) * rho_liq ) )**0.33333333333
!                radius = get_r_particle(n)
!
!              ENDDO
!
!           ENDDO
!        ENDDO
!     ENDDO

 END SUBROUTINE micro_user_diagnose
