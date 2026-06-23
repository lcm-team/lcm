 MODULE micro_task_util_MPI

 CONTAINS

 SUBROUTINE micro_particle_array2mpi_datatype()
!
!-- Define particle MPI data type
    USE micro_vars,      ONLY: newtype, nv0, nv1, niv0, lcmp
    USE micro_user_vars, ONLY: nvU, nivU

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    INTEGER :: buff_type(2), block_length(2), displacements(2), ierr
!
!-- particle_array: (int, ..., int, real, ..., real)
!-- buff_type:      (int, real)
!-- block_length:   (niv0+nivU for int, nv0+nv1 for real)
!-- displacement:   (0 bite for int, 4 bytes for first real)
    buff_type(1) = MPI_INTEGER8   ! 8 byte integer = MPI_INTEGER8; 4 byte integer = MPI_INTEGER4
    IF ( lcmp .EQ. 4 )  THEN
       buff_type(2) = MPI_REAL4   ! single precision, 4 byte real
    ELSEIF ( lcmp .EQ. 8 )  THEN
       buff_type(2) = MPI_REAL8   ! double precision, 8 byte real
    ELSE
       PRINT*,'LCM ERROR: Value of lcmp not allowed, it needs to be either 4 or 8.'
       CALL micro_task_abort
    ENDIF

    block_length(1)  = niv0 + nivU
    block_length(2)  = nv0 + nv1 + nvU
    displacements(1) = 0
    displacements(2) = 8 * ( niv0 + nivU )

    CALL MPI_TYPE_STRUCT( 2, block_length, displacements, buff_type, newtype, ierr )
    CALL MPI_TYPE_COMMIT(newtype, ierr)

 END SUBROUTINE micro_particle_array2mpi_datatype

 SUBROUTINE micro_task_recv_particle( buffer, length, tag, request )

    USE micro_vars, ONLY: newtype, particle_array

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    INTEGER, INTENT(IN)  :: length   ! buffers' length
    TYPE(particle_array), DIMENSION(length), INTENT(OUT) :: buffer  ! buffer of data
    INTEGER, INTENT(IN)  :: tag
    INTEGER, INTENT(OUT) :: request
    INTEGER              :: ierr

    CALL MPI_IRECV( buffer, length, newtype, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, request, ierr )

 END SUBROUTINE micro_task_recv_particle

 SUBROUTINE micro_task_bsend_particle( rank_to, buffer, length, tag )

    USE micro_vars, ONLY: newtype, particle_array

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    INTEGER, INTENT(IN) :: rank_to   ! receiving task's rank
    INTEGER, INTENT(IN) :: length    ! buffers' length
    TYPE(particle_array), DIMENSION(length), INTENT(IN) :: buffer   ! buffer of data
    INTEGER, INTENT(IN) :: tag       ! tag of the message
    INTEGER             :: ierr

    CALL MPI_SEND( buffer, length, newtype, rank_to, tag, MPI_COMM_WORLD, ierr )

 END SUBROUTINE micro_task_bsend_particle

 SUBROUTINE micro_task_abort()
!
!-- Abort process instead of using task_abort in task_util_MPI.f90

    USE grid, ONLY: dompi

    INCLUDE 'mpif.h'

    INTEGER:: ierr, rc

    IF (dompi) THEN
       CALL MPI_ABORT( MPI_COMM_WORLD, rc, ierr )
       CALL MPI_FINALIZE(ierr)
    ENDIF
    CALL EXIT(999) ! to avolid resubmission when finished

 END SUBROUTINE micro_task_abort

 SUBROUTINE micro_task_recv_integer( buffer, request )

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    INTEGER, INTENT(OUT) :: buffer   ! buffer of data
    INTEGER, INTENT(OUT) :: request
    INTEGER              :: ierr

    CALL MPI_IRECV( buffer,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,request,ierr )

 END SUBROUTINE micro_task_recv_integer

 SUBROUTINE micro_task_bsend_integer( rank_to, buffer, tag )

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    INTEGER, INTENT(IN) :: rank_to   ! receiving task's rank
    INTEGER, INTENT(IN) :: buffer    ! buffer of data
    INTEGER, INTENT(IN) :: tag       ! tag of the message
    INTEGER             :: ierr

    CALL MPI_SEND( buffer, 1, MPI_INTEGER, rank_to, tag, MPI_COMM_WORLD, ierr)

 END SUBROUTINE micro_task_bsend_integer

 SUBROUTINE micro_task_receive_LCM( buffer, length, request )

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    INTEGER, INTENT(IN)                  :: length   ! buffers' length
    REAL, DIMENSION(length), INTENT(OUT) :: buffer   ! buffer of data
    INTEGER, INTENT(OUT)                 :: request
    INTEGER                              :: ierr

    CALL MPI_IRECV(buffer,length,MPI_REAL8,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,request,ierr)

 END SUBROUTINE micro_task_receive_LCM

 SUBROUTINE micro_task_bsend_LCM( rank_to, buffer, length, tag )

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    INTEGER, INTENT(IN)                    :: rank_to ! receiving task's rank
    INTEGER, INTENT(IN)                    :: length  ! buffers' length
    REAL,    DIMENSION(length), INTENT(IN) :: buffer  ! buffer of data
    INTEGER, INTENT(IN)                    :: tag     ! tag of the message
    INTEGER                                :: ierr

    CALL MPI_SEND( buffer, length, MPI_REAL8, rank_to, tag, MPI_COMM_WORLD, ierr )

 END SUBROUTINE micro_task_bsend_LCM

 END MODULE micro_task_util_MPI
