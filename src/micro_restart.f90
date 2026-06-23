 MODULE micro_restart
!
!-- Writing and reading restart data
    USE micro_vars
    USE micro_task_util_MPI, ONLY: micro_task_abort

    USE grid,                ONLY: restart_sep, rank, nsubdomains, nstep, masterproc, nrestart, &
                                   dt, caseid_restart, case_restart, caseid, case, nx, ny, nzm
    USE vars,                ONLY: tabs, rho, pres, dtn

 CONTAINS

 SUBROUTINE micro_restart_read
!
!-- Read restart data
    IMPLICIT NONE

    TYPE(particle_array), ALLOCATABLE, DIMENSION(:) :: particle_temp
    CHARACTER (LEN=4) :: rankchar
    INTEGER           :: irank, ii, n, ierr
    INTEGER           :: lenstr
    EXTERNAL          :: lenstr

    CALL t_startf ('micro_read')

    IF(masterproc) PRINT*,'LCM: Reading restart file...'
!
!-- Allocate particle_temp (array) and particle (pointer-array)
    IF (restart_sep) THEN
       WRITE(rankchar,'(i4)') rank
       IF (nrestart .ne. 2) THEN
          OPEN(67, FILE='./RESTART/'//case(1:LEN_TRIM(case))//'_'//caseid(1:LEN_TRIM(caseid))// &
               '_'//rankchar(5-lenstr(rankchar):4)//'_restart_MICRO_LCM.bin', STATUS='old', &
               FORM='unformatted', iostat=ierr)
       ELSE
          OPEN(67, FILE='./RESTART/'//case_restart(1:LEN_TRIM(case_restart))//'_'//caseid_restart(1:LEN_TRIM(caseid_restart))//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart_MICRO_LCM.bin', STATUS='old', &
               FORM='unformatted', iostat=ierr)
       END IF
       IF (ierr/=0) THEN
         WRITE(*,*) 'Error', ierr, 'trying to open restart file. Cancelling.'
         call task_abort
       END IF

       READ(67) n_particles, history_length, filename_LCM, n_particles_alloc, n_ID 
       ALLOCATE( particle_temp(0:n_particles), particle(0:n_particles_alloc), STAT=ierr )
       IF ( ierr /= 0 ) THEN
          PRINT*,'LCM ERROR: Failed to allocate arrays in micro_restart_read on proc:', rank
          CALL micro_task_abort
       ENDIF
       READ(67) particle_temp
       CLOSE(67)
    ELSE
       WRITE(rankchar,'(i4)') nsubdomains
       IF (nrestart .ne. 2) THEN
          OPEN(67, FILE='./RESTART/'//case(1:LEN_TRIM(case))//'_'//caseid(1:LEN_TRIM(caseid))// &
               '_'//rankchar(5-lenstr(rankchar):4)//'_restart_MICRO_LCM.bin', STATUS='old', &
               FORM='unformatted', iostat=ierr)
       ELSE
          OPEN(67, FILE='./RESTART/'//case_restart(1:LEN_TRIM(case_restart))//'_'//caseid_restart(1:LEN_TRIM(caseid_restart))//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart_MICRO_LCM.bin', STATUS='old', &
               FORM='unformatted', iostat=ierr)
       END IF

       IF (ierr/=0) THEN
         WRITE(*,*) 'Error', ierr, 'trying to open restart file. Cancelling.'
         call task_abort
       END IF

       DO irank = 0, nsubdomains-1
          CALL task_barrier()
          IF (irank == rank) THEN
             DO ii = 0, irank-1   ! skip records
                READ(67)
                READ(67)
             ENDDO
             READ(67) n_particles, history_length, filename_LCM, n_particles_alloc, n_ID
             ALLOCATE( particle_temp(0:n_particles), particle(0:n_particles_alloc), STAT=ierr )
             IF ( ierr /= 0 ) THEN
                PRINT*,'LCM ERROR: Failed to allocate arrays in micro_restart_read on proc:', rank
                CALL micro_task_abort
             ENDIF
             READ(67) particle_temp
             CLOSE(67)
          ENDIF
       ENDDO
    ENDIF
!
!-- Store data from particle_temp (array) to particle (pointer-array)
!-- Deallocate particle_temp
    DO n = 1, n_particles
       particle(n)%ivar(:) = particle_temp(n)%ivar(:)
       particle(n)%var(:)  = particle_temp(n)%var(:)
    ENDDO
    DEALLOCATE( particle_temp, STAT=ierr )
    IF ( ierr /= 0 ) THEN
       PRINT*,'LCM ERROR: Failed to deallocate arrays in micro_restart_read on proc:', rank
       CALL micro_task_abort
    ENDIF
!
!-- Sort all particles by their grid box (this makes them easily accessible
!-- via n_start and n_end in a 3D loop structure)
    CALL micro_sort_by_grid_box
!
!-- Exchange LES fields
    CALL micro_coupling

    IF (rank == nsubdomains-1) THEN
       PRINT*,'Case:', caseid
       PRINT*,'Restart LCM at step:', nstep
       PRINT*,'Time:', nstep*dt
    ENDIF

    CALL task_barrier()

    CALL t_stopf ('micro_read')

 END SUBROUTINE micro_restart_read

 SUBROUTINE micro_restart_write
!
!-- Write restart data
    IMPLICIT NONE

    TYPE(particle_array), ALLOCATABLE, DIMENSION(:) :: particle_temp
    CHARACTER (LEN=4) :: rankchar
    INTEGER :: irank, n, ierr
    INTEGER :: lenstr
    EXTERNAL :: lenstr

    CALL t_startf ('micro_write')

    IF (masterproc) PRINT*,'LCM: Writting restart file...'
!
!-- Store data from particle (pointer-array) to particle_temp
    ALLOCATE( particle_temp(0:n_particles), STAT=ierr )
    IF ( ierr /= 0 ) THEN
       PRINT*,'LCM ERROR: Failed to allocate arrays in micro_restart_write on proc:', rank
       CALL micro_task_abort
    ENDIF
    DO n = 1, n_particles
       particle_temp(n)%ivar(:) = particle(n)%ivar(:)
       particle_temp(n)%var(:)  = particle(n)%var(:)
    ENDDO
!
!-- Write restart data
    IF (restart_sep) THEN
       WRITE(rankchar,'(i4)') rank
       OPEN(67, FILE='./RESTART/'//case(1:LEN_TRIM(case))//'_'//caseid(1:LEN_TRIM(caseid))// &
        '_'//rankchar(5-lenstr(rankchar):4)//'_restart_MICRO_LCM.bin', STATUS='unknown', &
        FORM='unformatted')
       WRITE(67) n_particles, history_length, filename_LCM, MAX(n_particles_alloc,n_particles), n_ID
       WRITE(67) particle_temp
       CLOSE(67)
    ELSE
       WRITE(rankchar,'(i4)') nsubdomains
       DO irank = 0, nsubdomains-1
          CALL task_barrier()
          IF (irank == rank) THEN
             IF (masterproc) THEN
              OPEN(67, FILE='./RESTART/'//case(1:LEN_TRIM(case))//'_'// &
               caseid(1:LEN_TRIM(caseid))//'_'//rankchar(5-lenstr(rankchar):4)// &
               '_restart_MICRO_LCM.bin', STATUS='unknown', FORM='unformatted')
           ELSE
              OPEN(67, FILE='./RESTART/'//case(1:LEN_TRIM(case))//'_'// &
               caseid(1:LEN_TRIM(caseid)) //'_'//rankchar(5-lenstr(rankchar):4)// &
               '_restart_MICRO_LCM.bin', STATUS='unknown', FORM='unformatted', POSITION='append')
           ENDIF
           WRITE(67) n_particles, history_length, filename_LCM, MAX(n_particles_alloc,n_particles), n_ID
           WRITE(67) particle_temp
           CLOSE(67)
          ENDIF
       ENDDO
    ENDIF
!
!-- Deallocate particle_temp
    DEALLOCATE( particle_temp, STAT=ierr )
    IF ( ierr /= 0 ) THEN
       PRINT*,'LCM ERROR: Failed to allocate arrays in micro_restart_write on proc:', rank
       CALL micro_task_abort
    ENDIF

    IF (masterproc) PRINT*,'LCM: Saved restart file. nstep=', nstep

    CALL task_barrier()

    CALL t_stopf ('micro_write')

 END SUBROUTINE micro_restart_write

 END MODULE micro_restart
