!TODO: exhange nx, ny by nx_LCM, ny_LCM

 SUBROUTINE micro_exchange_particles
!
!-- Exchange particles between subdomains
!
!-- This subroutine updates the position of each particle.
!-- 1. Check if each particle remains in the same subdomain or moves next
!--    subdomain. And make linked list with pointer for all directions
!--    including staying in the current rank.
!-- 2. MPI exchange (see task_exchange_particle)
!--    2.1 Store outgoing particle into array and send.
!--    2.2 Receive incoming particles.
!-- 3. Update particle pointer-array:
!--    3.1 Add incomping particles to the linked list.
!--    3.2 Store all particles into new particle.
    USE micro_task_util_MPI, ONLY: micro_task_abort
    USE micro_vars,          ONLY: particle, lcmp, particle,                     &
                                   particle_array, n_out, particle_out_nn_head,  &
                                   particle_out_nn_tail, particle_out_ne_head,   &
                                   particle_out_ne_tail, particle_out_ee_head,   &
                                   particle_out_ee_tail, particle_out_se_head,   &
                                   particle_out_se_tail, particle_out_ss_head,   &
                                   particle_out_ss_tail, particle_out_sw_head,   &
                                   particle_out_sw_tail, particle_out_ww_head,   &
                                   particle_out_ww_tail, particle_out_nw_head,   &
                                   particle_out_nw_tail, particle_tail, np_in,   &
                                   particle_in, particle_head, n_particles,      &
                                   np_out,  xwidth, ywidth, gpys, gpxs,          &
                                   n_particles_alloc

    USE grid,                ONLY: rank, nx, ny, RUN3D
    USE params,              ONLY: pi
    USE domain,              ONLY: YES3D

    IMPLICIT NONE

    INTEGER :: n, ierr
    TYPE(particle_array), ALLOCATABLE, DIMENSION(:) :: particle_temp

    CALL t_startf ('micro_exchange_particles')
!
!-- Initialize the empty list
    NULLIFY( particle_head, particle_tail, &
             particle_out_nn_head, particle_out_nn_tail, particle_out_ne_head, particle_out_ne_tail, &
             particle_out_ee_head, particle_out_ee_tail, particle_out_se_head, particle_out_se_tail, &
             particle_out_ss_head, particle_out_ss_tail, particle_out_sw_head, particle_out_sw_tail, &
             particle_out_ww_head, particle_out_ww_tail, particle_out_nw_head, particle_out_nw_tail )
!
!-- Prepare for loop
    n_out(:) = 0

    DO n = 1, n_particles
!
!--    STEP 1: Add particle to the list, particle_head or particle_out_##_head.
!--    Count n_out. Determine if each particle is going out or not, then make
!--    linked lists.
!--
!--    The south and west boundaries belong the current subdomain. The following
!--    diagram shows the indices for gpxs and gpys for the current rank to judge
!--    if a particle moves to a neighbour rank.
!--    _              _              _
!--    | north-west   | north        | north-east
!--    | ips <        | ips <=       |
!--    |              | ipe+1 >      | ipe+1 >=
!--    |              |              | jpe+1 >=
!--    | jpe+1 >=     | jpe+1 >=     |
!--    *------------| *------------| *------------| <=== jpe+1
!--    _              _              _
!--    | west         |              | east
!--    | ips <        |              |
!--    |              |    rank      | ipe+1 >=
!--    | jps >=       |              | jps >=
!--    | jpe+1 <      |              | jpe+1 <
!--    *------------| *------------| *------------| <=== jps
!--    _              _              _
!--    | south-west   | south        | south-east
!--    | ips <        | ips <=       |
!--    |              | ipe+1 >      | ipe+1 >=
!--    |              | jps <        | jps <
!--    | jps <        |              |
!--    *------------| *------------| *------------|
!--                   ^              ^
!--                   ^              ^
!--                  ips           ipe+1
       IF (RUN3D) THEN
!
!--       North
          IF ( (particle(n)%var(1) >= gpxs(1)).AND.(particle(n)%var(1) < gpxs(nx+1)) &
               .AND.(particle(n)%var(2) >= gpys(ny+1)) ) THEN
!
!--          Horizontal boundary position adjustment
             IF ( particle(n)%var(2) >= ywidth ) particle(n)%var(2) = particle(n)%var(2) - ywidth
!
!--          Add a new particle to the end of the list (append method)
             IF ( ASSOCIATED(particle_out_nn_head) ) THEN  ! List is not empty:
                particle_out_nn_tail%next => particle(n)   !  -> attach new particle at end of list (particle_out_nn_tail still points to the last particle, i.e., the pointer of the last element is changed from null to the current particle)
                NULLIFY(particle(n)%next)                  !  -> no successor of new last particle in list
                particle_out_nn_tail => particle(n)        !  -> reset particle_tail pointer
             ELSE                                          ! List is empty:
                particle_out_nn_head => particle(n)        !  -> start list with particle (head)
                particle_out_nn_tail => particle(n)        !  -> end list with the same particle (tail)
                NULLIFY(particle_out_nn_tail%next)         !  -> no successor of last particle in list
             ENDIF
!
!--          Count number of particle going to Northern subdomain
             n_out(6) = n_out(6) + 1
!
!--       North-east
          ELSE IF ( (particle(n)%var(1) >= gpxs(nx+1)).AND.(particle(n)%var(2) >= gpys(ny+1)) ) THEN
!
!--          Horizontal boundary position adjustment
             IF ( particle(n)%var(1) >= xwidth ) particle(n)%var(1) = particle(n)%var(1) - xwidth
             IF ( particle(n)%var(2) >= ywidth ) particle(n)%var(2) = particle(n)%var(2) - ywidth
!
!--          Add a new particle to the end of the list
             IF ( ASSOCIATED(particle_out_ne_head) ) THEN
                particle_out_ne_tail%next => particle(n)
                NULLIFY(particle(n)%next)
                particle_out_ne_tail => particle(n)
             ELSE
                 particle_out_ne_head => particle(n)
                 particle_out_ne_tail => particle(n)
                 NULLIFY(particle_out_ne_tail%next)
             ENDIF
             n_out(8) = n_out(8) + 1
!
!--       East
          ELSE IF ( (particle(n)%var(1) >= gpxs(nx+1)) &
                    .AND.(particle(n)%var(2) >= gpys(1)).AND.(particle(n)%var(2) < gpys(ny+1)) ) THEN
!
!--          Horizontal boundary position adjustment
             IF ( particle(n)%var(1) >= xwidth ) particle(n)%var(1) = particle(n)%var(1) - xwidth
!
!--          Add a new particle to the end of the list
             IF ( ASSOCIATED(particle_out_ee_head) ) THEN
                particle_out_ee_tail%next => particle(n)
                NULLIFY(particle(n)%next)
                particle_out_ee_tail => particle(n)
             ELSE
                particle_out_ee_head => particle(n)
                particle_out_ee_tail => particle(n)
                NULLIFY(particle_out_ee_tail%next)
             ENDIF
             n_out(1) = n_out(1) + 1
!
!--       South-east
          ELSE IF ( (particle(n)%var(1) >= gpxs(nx+1)).AND.(particle(n)%var(2) < gpys(1)) ) THEN
!
!--          Horizontal boundary position adjustment
             IF ( particle(n)%var(1) >= xwidth ) particle(n)%var(1) = particle(n)%var(1) - xwidth
             IF ( particle(n)%var(2) < 0.0 ) particle(n)%var(2) = ywidth + particle(n)%var(2)
!
!--          Add a new particle to the end of the list
             IF ( ASSOCIATED(particle_out_se_head) ) THEN
                particle_out_se_tail%next => particle(n)
                NULLIFY(particle(n)%next)
                particle_out_se_tail => particle(n)
             ELSE
                particle_out_se_head => particle(n)
                particle_out_se_tail => particle(n)
                NULLIFY(particle_out_se_tail%next)
             ENDIF
             n_out(3) = n_out(3) + 1
!
!--       South
          ELSE IF ( (particle(n)%var(1) >= gpxs(1)).AND.(particle(n)%var(1) < gpxs(nx+1)) &
                    .AND.(particle(n)%var(2) < gpys(1)) ) THEN
!
!--          Horizontal boundary position adjustment
             IF ( particle(n)%var(2) < 0.0 ) particle(n)%var(2) = ywidth + particle(n)%var(2)
!
!--          Add a new particle to the end of the list
             IF ( ASSOCIATED(particle_out_ss_head) ) THEN
                particle_out_ss_tail%next => particle(n)
                NULLIFY(particle(n)%next)
                particle_out_ss_tail => particle(n)
             ELSE
                particle_out_ss_head => particle(n)
                particle_out_ss_tail => particle(n)
                NULLIFY(particle_out_ss_tail%next)
             ENDIF
             n_out(5) = n_out(5) + 1
!
!--       South-west
          ELSE IF ( (particle(n)%var(1) < gpxs(1)).AND.(particle(n)%var(2) < gpys(1)) ) THEN
!
!--          Horizontal boundary position adjustment
             IF ( particle(n)%var(1) < 0.0 ) particle(n)%var(1) = xwidth + particle(n)%var(1)
             IF ( particle(n)%var(2) < 0.0 ) particle(n)%var(2) = ywidth + particle(n)%var(2)
!
!--          Add a new particle to the end of the list
             IF ( ASSOCIATED(particle_out_sw_head) ) THEN
                particle_out_sw_tail%next => particle(n)
                NULLIFY(particle(n)%next)
                particle_out_sw_tail => particle(n)
             ELSE
                particle_out_sw_head => particle(n)
                particle_out_sw_tail => particle(n)
                NULLIFY(particle_out_sw_tail%next)
             ENDIF
             n_out(7) = n_out(7) + 1
!
!--       West
          ELSE IF ( (particle(n)%var(1) < gpxs(1) ) &
                    .AND.(particle(n)%var(2) >= gpys(1)).AND.(particle(n)%var(2) < gpys(ny+1)) ) THEN
!
!--          Horizontal boundary position adjustment
             IF ( particle(n)%var(1) < 0.0 ) particle(n)%var(1) = xwidth + particle(n)%var(1)
!
!--          Add a new particle to the end of the list
             IF ( ASSOCIATED(particle_out_ww_head) ) THEN
                particle_out_ww_tail%next => particle(n)
                NULLIFY(particle(n)%next)
                particle_out_ww_tail => particle(n)
             ELSE
                particle_out_ww_head => particle(n)
                particle_out_ww_tail => particle(n)
                NULLIFY(particle_out_ww_tail%next)
             ENDIF
             n_out(2) = n_out(2) + 1
!
!--       North-west
          ELSE IF ( (particle(n)%var(1) < gpxs(1)).AND.(particle(n)%var(2) >= gpys(ny+1)) ) THEN
!
!--          Horizontal boundary position adjustment
             IF ( particle(n)%var(1) < 0.0 ) particle(n)%var(1) = xwidth + particle(n)%var(1)
             IF ( particle(n)%var(2) >= ywidth ) particle(n)%var(2) = particle(n)%var(2) - ywidth
!
!--          Add a new particle to the end of the list
             IF ( ASSOCIATED(particle_out_nw_head) ) THEN
                particle_out_nw_tail%next => particle(n)
                NULLIFY(particle(n)%next)
                particle_out_nw_tail => particle(n)
             ELSE
                particle_out_nw_head => particle(n)
                particle_out_nw_tail => particle(n)
                NULLIFY(particle_out_nw_tail%next)
             ENDIF
             n_out(4) = n_out(4) + 1
!
!--       Stay in this subdomain
          ELSE
!
!--          Add a new particle to the end of the list
             IF ( ASSOCIATED(particle_head) ) THEN
                particle_tail%next => particle(n)
                NULLIFY(particle(n)%next)
                particle_tail => particle(n)
             ELSE
                particle_head => particle(n)
                particle_tail => particle(n)
                NULLIFY(particle_tail%next)
             ENDIF

          ENDIF
!
!--    2D run:
       ELSE
!
!--       East
          IF ( particle(n)%var(1) >= gpxs(nx+1) ) THEN
!
!--          Horizontal boundary position adjustment
             IF ( particle(n)%var(1) >= xwidth ) particle(n)%var(1) = particle(n)%var(1) - xwidth
!
!--          Add a new particle to the end of the list
             IF ( ASSOCIATED(particle_out_ee_head) ) THEN
                particle_out_ee_tail%next => particle(n)
                NULLIFY(particle(n)%next)
                particle_out_ee_tail => particle(n)
             ELSE
                particle_out_ee_head => particle(n)
                particle_out_ee_tail => particle(n)
                NULLIFY(particle_out_ee_tail%next)
             ENDIF
             n_out(1) = n_out(1) + 1
!
!--       West
          ELSE IF ( particle(n)%var(1) < gpxs(1) ) THEN
!
!--          Horizontal boundary position adjustment
             IF ( particle(n)%var(1) < 0.0 ) particle(n)%var(1) = xwidth + particle(n)%var(1)
!
!--          Add a new particle to the end of the list
             IF ( ASSOCIATED(particle_out_ww_head) ) THEN
                particle_out_ww_tail%next => particle(n)
                NULLIFY(particle(n)%next)
                particle_out_ww_tail => particle(n)
             ELSE
                particle_out_ww_head => particle(n)
                particle_out_ww_tail => particle(n)
                NULLIFY(particle_out_ww_tail%next)
             ENDIF
             n_out(2) = n_out(2) + 1
!
!--       Stay in this subdomain
          ELSE
!
!--          Add a new particle to the end of the list
             IF ( ASSOCIATED(particle_head) ) THEN
                particle_tail%next => particle(n)
                NULLIFY(particle(n)%next)
                particle_tail => particle(n)
             ELSE
                particle_head => particle(n)
                particle_tail => particle(n)
                NULLIFY(particle_tail%next)
             ENDIF

          ENDIF

       ENDIF

    ENDDO
!
!-- STEP 2: MPI exchange
    np_out = n_out(1)+n_out(2)+n_out(3)+n_out(4)+n_out(5)+n_out(6)+n_out(7)+n_out(8)

    CALL task_exchange_particle
!
!-- STEP 3:   Update particle pointer-array
!-- STEP 3.1: Add incoming particles to the linked list
    IF (np_in > 0) THEN
       DO n = 1, np_in
          IF ( ASSOCIATED(particle_head) ) THEN
             particle_tail%next => particle_in(n)
             NULLIFY(particle_in(n)%next)
             particle_tail => particle_in(n)
          ELSE
             particle_head => particle_in(n)
             particle_tail => particle_in(n)
             NULLIFY(particle_tail%next)
          ENDIF
       ENDDO
    ENDIF
!
!-- STEP 3.2: Store particle_head (linked list data) into new particle
!--           (pointer-array) through particle_temp (array) for next time step.
!--           Do clean deallocation.
    n_particles = n_particles - np_out + np_in
    ! print*, 'motion', rank, n_particles_alloc, n_particles
    n_particles_alloc = MAX( n_particles_alloc, n_particles )

    ALLOCATE( particle_temp(0:n_particles_alloc), STAT=ierr )
    IF ( ierr /= 0 ) THEN
       PRINT*,'LCM ERROR: Failed to allocate arrays in micro_xyz on proc:', rank
       CALL micro_task_abort
    ENDIF

    IF ( n_particles > 0 .AND. .NOT.ASSOCIATED(particle_head) ) THEN
       PRINT*,'LCM: particle_head failed. Abort. rank:', rank
       CALL micro_task_abort
    ELSE
       DO n = 1, n_particles
          particle_temp(n)%ivar(:) = particle_head%ivar(:)
          particle_temp(n)%var(:)  = particle_head%var(:)
          particle_head => particle_head%next
       ENDDO
    ENDIF

    NULLIFY( particle_head, particle_tail, &
             particle_out_nn_head, particle_out_nn_tail, particle_out_ne_head, particle_out_ne_tail, &
             particle_out_ee_head, particle_out_ee_tail, particle_out_se_head, particle_out_se_tail, &
             particle_out_ss_head, particle_out_ss_tail, particle_out_sw_head, particle_out_sw_tail, &
             particle_out_ww_head, particle_out_ww_tail, particle_out_nw_head, particle_out_nw_tail )

    DEALLOCATE( particle, particle_in, STAT=ierr )

    IF ( ierr /= 0 ) THEN
       PRINT*,'LCM ERROR: Failed to deallocate arrays in micro_xyz on proc:', rank
       CALL micro_task_abort
    ENDIF

    ALLOCATE( particle(0:n_particles_alloc), STAT=ierr )

    IF ( ierr /= 0 ) THEN
       PRINT*,'LCM ERROR: Failed to allocate arrays in micro_xyz on proc:', rank
       CALL micro_task_abort
    ENDIF

    DO n = 1, n_particles
       particle(n)%ivar(:) = particle_temp(n)%ivar(:)
       particle(n)%var(:)  = particle_temp(n)%var(:)
    ENDDO

    DEALLOCATE( particle_temp, STAT=ierr )

    IF ( ierr /= 0 ) THEN
       PRINT*,'LCM ERROR: Failed to deallocate arrays in micro_xyz on proc:', rank
       CALL micro_task_abort
    ENDIF

    CALL t_stopf ('micro_exchange_particles')

 END SUBROUTINE micro_exchange_particles

 SUBROUTINE task_exchange_particle
!
!-- Non-blocking receives before blocking sends for efficiency on IBM SP.
!--
!-- Procedure:
!-- 1. 3D run or 2D run?
!-- 2. Set ranks for loop.
!-- 3. Obtain n_in, np_in, n_max. <= MPI exchange
!-- 4. Receive and send particles. <= MPI exchange
!--
!-- ----------------- ----------------- -----------------
!-- | i=4             | i=6             | i=8             |
!-- | ranks(4)=ranknw | ranks(6)=ranknn | ranks(8)=rankne |
!-- | tags(4)=3       | tags(6)=5       | tags(8)=7       |
!--  ----------------- ----------------- -----------------
!-- | i=2             |                 | i=1             |
!-- | ranks(2)=rankww |  current rank   | ranks(1)=rankee |
!-- | tags(2)=1       |                 | tags(1)=2       |
!--  ----------------- ----------------- -----------------
!-- | i=7             | i=5             | i=3             |
!-- | ranks(7)=ranksw | ranks(5)=rankss | ranks(3)=rankse |
!-- | tags(7)=8       | tags(5)=6       | tags(3)=4       |
!--  ----------------- ----------------- -----------------
!-- * i is the indices for neighbour ranks.
!-- * ranks are the actual neighbour rank indices.
!-- * tags converts to incoming index from outgoing index. If i=1 for the
!--   current rank, then i=2 for the right next rank.

    USE micro_vars

    USE grid, ONLY: nx, ny, nz, nxp1, nyp1, rank, rankee, rankne, ranknn, ranknw, &
                    rankse, rankss, ranksw, rankww, RUN3D, masterproc

    USE micro_task_util_MPI

    IMPLICIT NONE

    TYPE(particle_array),   ALLOCATABLE ::  buff_send(:,:)   ! buffer to send data
    TYPE(particle_array),   ALLOCATABLE ::  buff_recv(:,:)   ! buffer to receive data
    TYPE(particle_pointer), POINTER     ::  temp

    LOGICAL, DIMENSION(8) :: flag
    INTEGER, DIMENSION(8) :: ranks, n_in, reqs, tags, buff
    INTEGER               :: n_max, tag_in, count, status, i, j, iend, n, ierr
!
!-- STEP 1: 3D run or 2D run?
    IF ( RUN3D ) THEN
       iend = 8
    ELSE
       iend = 2
    ENDIF
!
!-- STEP 2: Prepare for MPI exchange
    ranks(1) = rankee   ! East
    ranks(2) = rankww   ! West
    ranks(3) = rankse   ! South east
    ranks(4) = ranknw   ! North west
    ranks(5) = rankss   ! South
    ranks(6) = ranknn   ! North
    ranks(7) = ranksw   ! South west
    ranks(8) = rankne   ! North East

    tags(1) = 2         ! rank tag for incoming particles
    tags(2) = 1
    tags(3) = 4
    tags(4) = 3
    tags(5) = 6
    tags(6) = 5
    tags(7) = 8
    tags(8) = 7
!
!-- STEP 3: Obtain n_in by MPI. Non-blocking receives n_in from each neighbors:
    flag(:) = .FALSE.
    DO i = 1, iend
       IF ( rank /= ranks(i) ) THEN
!
!--       i is the incoming rank index for the current rank. However, since number of
!--       communication is known (either 2 or 8) and each communication passes/receives only
!--       one number, MPI_ANY_TAG is used in task_rec_integer.
          CALL micro_task_recv_integer(buff(i), reqs(i))
       ELSE
!
!--       This flag indicates either no MPI communiation or complete MPI communication.
!--       flag has to be true to avoid calling task_test for MPI communication.
          flag(i) = .TRUE.
       ENDIF
    ENDDO
!
!-- Blocking sends n_out to each neighbors:
!-- Fill the data from the buffers that have the same rank (were not sent):
    count = 0
    DO i = 1, iend
       IF ( rank /= ranks(i) ) THEN
          CALL micro_task_bsend_integer(ranks(i), n_out(i), i)
       ELSE
          n_in(tags(i)) = n_out(i)   ! i is outgoing rank
          count = count + 1          ! tags(i) convers to incoming rank, e.g. i=1, tags(i)=2
       ENDIF
    ENDDO
!
!-- Monitor the progress of receiving
    DO WHILE ( count < iend )
       DO i = 1, iend
          IF ( .NOT.flag(i) ) THEN
             CALL task_test(reqs(i), flag(i), status, tag_in)
!
!--          tag_in is the same as i in task_bsend_integer.
!--          tags(tag_in) converts to the correct incoming rank index for the
!--          current rank.
             IF ( flag(i) ) THEN
                n_in(tags(tag_in)) = buff(i)
                count = count + 1
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    CALL task_barrier

    ! IF (masterproc) THEN
    !    PRINT*,'LCM: from E ', n_in(1)
    !    PRINT*,'LCM: from W ', n_in(2)
    !    PRINT*,'LCM: from SE', n_in(3)
    !    PRINT*,'LCM: from NW', n_in(4)
    !    PRINT*,'LCM: from S ', n_in(5)
    !    PRINT*,'LCM: from N ', n_in(6)
    !    PRINT*,'LCM: from SW', n_in(7)
    !    PRINT*,'LCM: from NE', n_in(8)
    ! ENDIF
!
!-- STEP 4: Obtain np_in, n_max and allocate arrays
    np_in = 0
    n_max = 0
    DO i = 1, iend
       np_in = np_in + n_in(i)
       IF ( n_in(i) > n_max ) n_max = n_in(i)
       IF ( n_out(i) > n_max ) n_max = n_out(i)
    ENDDO
    IF ( np_in == 0 ) np_in = 1   ! number of element has to be more than / equal to 1
    IF ( n_max == 0 ) n_max = 1   ! number of element has to be more than / equal to 1
    ALLOCATE( buff_recv(n_max,iend), buff_send(n_max,iend), particle_in(np_in), STAT=ierr )
    IF ( ierr /= 0 ) THEN
       PRINT*,'LCM ERROR: Failed to allocate arrays in task_exchange_particle on proc:', rank
       CALL micro_task_abort
    ENDIF
!
!-- STEP 5: Receive and send particles. Non-blocking receives first:
    flag(:) = .FALSE.
    DO i = 1, iend
       IF (( rank /= ranks(i) ).AND.( n_in(i) > 0 )) THEN
!
!--       i is the incoming rank index for the current rank.
!--       tags(i) is the incoming rank, which should be same as i for the destination rank in
!--       task_bsend_particle.
          CALL micro_task_recv_particle(buff_recv(1,i), n_in(i), tags(i), reqs(i))
       ELSE
!
!--       This flag indicates either no MPI communiation or complete MPI communication.
!--       flag has to be true to avoid calling task_test for MPI communication.
          flag(i) = .TRUE.
       ENDIF
    ENDDO
!
!-- Blocking sends:
    NULLIFY(temp)   ! empty

    IF ( n_out(1) > 0 ) THEN
       IF ( .NOT.ASSOCIATED(particle_out_ee_head) ) THEN
          IF (masterproc) PRINT*,'LCM: particle_out_ee_head failed. Abort. rank:', rank
          CALL micro_task_abort
       ELSE
          temp => particle_out_ee_head   ! in order not to mess head, use temporal pointer
          DO n = 1, n_out(1)
             buff_send(n,1)%ivar(:) = temp%ivar(:)
             buff_send(n,1)%var(:)  = temp%var(:)
             temp => temp%next
          ENDDO
       ENDIF
       NULLIFY(temp)   ! for clean deallocation
    ENDIF

    IF ( n_out(2) > 0 ) THEN
       IF ( .NOT.ASSOCIATED(particle_out_ww_head) ) THEN
          IF (masterproc) PRINT*,'LCM: particle_out_ww_head failed. Abort. rank:', rank
          CALL micro_task_abort
       ELSE
          temp => particle_out_ww_head
          DO n = 1, n_out(2)
             buff_send(n,2)%ivar(:) = temp%ivar(:)
             buff_send(n,2)%var(:)  = temp%var(:)
             temp => temp%next
          ENDDO
       ENDIF
       NULLIFY(temp)
    ENDIF

    IF ( RUN3D ) THEN

       IF ( n_out(3) > 0 ) THEN
          IF ( .NOT.ASSOCIATED(particle_out_se_head) ) THEN
             IF (masterproc) PRINT*,'LCM: particle_out_se_head failed. Abort. rank:', rank
             CALL micro_task_abort
          ELSE
             temp => particle_out_se_head
             DO n = 1, n_out(3)
                buff_send(n,3)%ivar(:) = temp%ivar(:)
                buff_send(n,3)%var(:)  = temp%var(:)
                temp => temp%next
             ENDDO
          ENDIF
          NULLIFY(temp)
       ENDIF

       IF ( n_out(4) > 0 ) THEN
          IF ( .NOT.ASSOCIATED(particle_out_nw_head) ) THEN
             IF (masterproc) PRINT*,'LCM: particle_out_nw_head failed. Abort. rank:', rank
             CALL micro_task_abort
          ELSE
             temp => particle_out_nw_head
             DO n = 1, n_out(4)
                buff_send(n,4)%ivar(:) = temp%ivar(:)
                buff_send(n,4)%var(:)  = temp%var(:)
                temp => temp%next
             ENDDO
          ENDIF
          NULLIFY(temp)
       ENDIF

       IF ( n_out(5) > 0 ) THEN
          IF ( .NOT.ASSOCIATED(particle_out_ss_head) ) THEN
             IF (masterproc) PRINT*,'LCM: particle_out_ss_head failed. Abort. rank:', rank
             CALL micro_task_abort
          ELSE
             temp => particle_out_ss_head
             DO n = 1, n_out(5)
                buff_send(n,5)%ivar(:) = temp%ivar(:)
                buff_send(n,5)%var(:)  = temp%var(:)
                temp => temp%next
             ENDDO
          ENDIF
          NULLIFY(temp)
       ENDIF

       IF ( n_out(6) > 0 ) THEN
          IF ( .NOT.ASSOCIATED(particle_out_nn_head) ) THEN
             IF (masterproc) PRINT*,'LCM: particle_out_nn_head failed. Abort. rank:', rank
             CALL micro_task_abort
          ELSE
             temp => particle_out_nn_head
             DO n = 1, n_out(6)
                buff_send(n,6)%ivar(:) = temp%ivar(:)
                buff_send(n,6)%var(:)  = temp%var(:)
                temp => temp%next
             ENDDO
          ENDIF
          NULLIFY(temp)
       ENDIF

       IF ( n_out(7) > 0 ) THEN
          IF ( .NOT.ASSOCIATED(particle_out_sw_head) ) THEN
             IF (masterproc) PRINT*,'LCM: particle_out_sw_head failed. Abort. rank:', rank
             CALL micro_task_abort
          ELSE
             temp => particle_out_sw_head
             DO n = 1, n_out(7)
                buff_send(n,7)%ivar(:) = temp%ivar(:)
                buff_send(n,7)%var(:)  = temp%var(:)
                temp => temp%next
             ENDDO
          ENDIF
          NULLIFY(temp)
       ENDIF

       IF ( n_out(8) > 0 ) THEN
          IF ( .NOT.ASSOCIATED(particle_out_ne_head) ) THEN
             IF (masterproc) PRINT*,'LCM: particle_out_ne_head failed. Abort. rank:', rank
             CALL micro_task_abort
          ELSE
             temp => particle_out_ne_head
             DO n = 1, n_out(8)
                buff_send(n,8)%ivar(:) = temp%ivar(:)
                buff_send(n,8)%var(:)  = temp%var(:)
                temp => temp%next
             ENDDO
          ENDIF
          NULLIFY(temp)
       ENDIF

    ENDIF

    DO i = 1, iend
       IF ( n_out(i) > 0 ) THEN
          IF ( rank /= ranks(i) ) THEN
!
!--          i is the destination rank index from the current rank
             CALL micro_task_bsend_particle(ranks(i), buff_send(1:n_out(i),i), n_out(i), i)
          ELSE
!
!--          Fill the buff_recv from buff_send that has the same rank
             DO j = 1, n_out(i)
                buff_recv(j,i) = buff_send(j,i)
             ENDDO
          ENDIF
       ENDIF
    ENDDO
!
!-- Fill particle_in from buff_recv that has the same rank:
    np_in = 0
    count = 0
    DO i = 1, iend
       IF ( flag(i) ) THEN   ! no MPI communication
          count = count + 1
          IF (( rank == ranks(tags(i)) ).AND.( n_in(tags(i)) > 0 )) THEN
!
!--          There are outgoing particles which will be periodically back to the same rank
             DO j = 1, n_in(tags(i))
                particle_in(np_in+j)%ivar(:) = buff_recv(j,i)%ivar(:)
                particle_in(np_in+j)%var(:)  = buff_recv(j,i)%var(:)
             ENDDO
             np_in = np_in + n_in(tags(i))
          ENDIF
       ENDIF
    ENDDO
!
!-- Monitor the progress of receiving; fill the data:
    DO WHILE ( count < iend )
       DO i = 1, iend
          IF ( .NOT.flag(i) ) THEN
!
!--          reqs(i) is the same reqs(i) in task_recv_particle.
!--          tag_in is the same as tags(i) in task_recv_particle.
             CALL task_test(reqs(i), flag(i), status, tag_in)
             IF ( flag(i) ) then
                count = count + 1
                DO j = 1, n_in(i)
                   particle_in(np_in+j)%ivar(:) = buff_recv(j,i)%ivar(:)
                   particle_in(np_in+j)%var(:)  = buff_recv(j,i)%var(:)
                ENDDO
                np_in = np_in + n_in(i)
               ENDIF
          ENDIF
       ENDDO
    ENDDO

    CALL task_barrier

    DEALLOCATE( buff_recv, buff_send, STAT=ierr )
    IF ( ierr /= 0 ) THEN
       PRINT*,'LCM ERROR: Failed to allocate arrays in task_exchange_particle on proc:', rank
       CALL micro_task_abort
    ENDIF

 END SUBROUTINE task_exchange_particle

 SUBROUTINE task_exchange_LCM(f, dimx1, dimx2, dimy1, dimy2, dimz, i_1, i_2, j_1, j_2, id)
!
!-- MPI CODE FOR GRID VALUE EXCHANGE (ghost points of u, v, w fields)
!--
!-- Non-blocking receives before blocking sends for efficiency on IBM SP.
!-- Sends and receives the boundary messages
!-- This subroutine is the same subroutine, "task_exchange.f90," but modified for LCM purpose.

    USE micro_vars

    USE grid, ONLY: nx, ny, nz, nxp1, nyp1, rank, rankee, rankne, ranknn, ranknw, &
                    rankse, rankss, ranksw, rankww, RUN3D, masterproc

    USE micro_task_util_MPI

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: dimx1, dimx2, dimy1, dimy2, dimz
    INTEGER, INTENT(IN) :: i_1, i_2, j_1, j_2
    REAL,    DIMENSION(dimx1:dimx2, dimy1:dimy2, dimz), INTENT(IN) :: f
    INTEGER, INTENT(IN) :: id   ! id of the sent field


    INTEGER, PARAMETER :: nghost  = MAX(ux1,ux2,vy1,vy2,sx1,sx2,sy1,sy2)
    INTEGER, PARAMETER :: bufflen = MAX(nx,ny) * nghost * nz
    REAL,    DIMENSION(bufflen)   :: buff_send   ! buff to send data
    REAL,    DIMENSION(bufflen,8) :: buff_recv ! buff to receive data
    INTEGER, DIMENSION(8) :: reqs_in, ids
    INTEGER, DIMENSION(8) :: i_start, i_end, j_start, j_end, ranks
    LOGICAL, DIMENSION(8) :: flag

    INTEGER i, j, k, n, rf, tag, m, mend, i1, i2, j1, j2, mask, count

    i1 = i_1 - 1
    i2 = i_2 - 1
    j1 = j_1 - 1
    j2 = j_2 - 1
!
!-- East
    i_start(1) = nx-i1
    i_end  (1) = nx
    j_start(1) = 1
    j_end  (1) = ny
!
!-- West
    i_start(2) = 1
    i_end  (2) = 1+i2
    j_start(2) = 1
    j_end  (2) = ny
!
!-- North
    i_start(7) = 1
    i_end  (7) = nx
    j_start(7) = ny-j1
    j_end  (7) = ny
!
!-- North East
    i_start(8) = nx-i1
    i_end  (8) = nx
    j_start(8) = ny-j1
    j_end  (8) = ny
!
!-- South East
    i_start(3) = nx-i1
    i_end  (3) = nx
    j_start(3) = 1
    j_end  (3) = 1+j2
!
!-- South
    i_start(4) = 1
    i_end  (4) = nx
    j_start(4) = 1
    j_end  (4) = 1+j2
!
!-- South West
    i_start(5) = 1
    i_end  (5) = 1+i2
    j_start(5) = 1
    j_end  (5) = 1+j2
!
!-- North West
    i_start(6) = 1
    i_end  (6) = 1+i2
    j_start(6) = ny-j1
    j_end  (6) = ny

    ranks(1) = rankee
    ranks(2) = rankww
    ranks(7) = ranknn
    ranks(8) = rankne
    ranks(3) = rankse
    ranks(4) = rankss
    ranks(5) = ranksw
    ranks(6) = ranknw

    mask=id*100000+i1*10000+i2*1000+j1*100+j2*10

    ids(1) = mask+7
    ids(2) = mask+3
    ids(7) = mask+5
    ids(8) = mask+6
    ids(3) = mask+8
    ids(4) = mask+1
    ids(5) = mask+2
    ids(6) = mask+4
!
!-- Send/receive buffs to/from neighbors
    IF ( RUN3D ) THEN
       mend = 8
    ELSE
       mend = 2
    ENDIF
!
!-  Non-blocking receives first:
    DO m = 1, mend
       IF ( rank /= ranks(m) ) THEN
          CALL micro_task_receive_LCM(buff_recv(:,m), bufflen, reqs_in(m))
          flag(m) = .FALSE.
       ELSE
          flag(m) = .TRUE.
       ENDIF
    ENDDO
!
!-- Blocking sends:
    DO m = 1, mend
       n = 0
       DO k = 1, dimz
          DO j = j_start(m), j_end(m)
             DO i = i_start(m), i_end(m)
                n = n+1
                buff_send(n) = f(i,j,k)
             ENDDO
          ENDDO
       ENDDO
       IF ( rank /= ranks(m) ) THEN
          CALL micro_task_bsend_LCM(ranks(m), buff_send, n, ids(m))
       ELSE
          DO i = 1, n
             buff_recv(i,m) = buff_send(i)
          ENDDO
       ENDIF
    ENDDO
!
!-- Fill the data from the buffers that have the same rank (were not sent):
    count = 0
    DO m = 1, mend
       IF ( flag(m) ) THEN
          count = count+1
          tag = ids(m)
          CALL task_dispatch_LCM(bufflen, buff_recv(1:bufflen,m), tag)
       ENDIF
    ENDDO
!
!-- Monitor the progress of receiving; fill the data:
    DO WHILE (count < mend)
       DO m = 1, mend
          IF ( .NOT.flag(m) ) THEN
             CALL task_test(reqs_in(m), flag(m), rf, tag)
             IF ( flag(m) ) THEN
                count = count + 1
                CALL task_dispatch_LCM(bufflen, buff_recv(1:bufflen,m), tag)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    CALL task_barrier()

 END SUBROUTINE task_exchange_LCM

 SUBROUTINE task_dispatch_LCM(blen, buff, tag)
!
!-- Dispatches the messages according to the field sent
    USE micro_vars

    USE grid, ONLY: nx, ny, nz, nxp1, nyp1, rank, rankee, rankne, ranknn, ranknw, &
                    rankse, rankss, ranksw, rankww, RUN3D, masterproc

    USE micro_task_util_MPI

    USE grid, ONLY: nzm

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: blen
    REAL, DIMENSION(blen), INTENT(IN) :: buff
    INTEGER, INTENT(IN) :: tag
    INTEGER :: field, n

    field = tag / 100000
!
!-- resolved velocity components, field index: 1-3
    SELECT CASE(field)
       CASE(1)
          CALL task_assign_LCM(u_LCM, x1u, x2u, y1u, y2u, nzm, buff, blen, tag)
       CASE(2)
          CALL task_assign_LCM(v_LCM, x1v, x2v, y1v, y2v, nzm, buff, blen, tag)
       CASE(3)
          CALL task_assign_LCM(w_LCM, x1w, x2w, y1w, y2w, nzm, buff, blen, tag)
    END SELECT

 END SUBROUTINE task_dispatch_LCM

 SUBROUTINE task_assign_LCM(f, dimx1, dimx2, dimy1, dimy2, dimz, buff, blen, tag)
!
!-- This routine assigns the ghost points after MPI exchange
    USE micro_vars

    USE grid, ONLY: nx, ny, nz, nxp1, nyp1, rank, rankee, rankne, ranknn, ranknw, &
                    rankse, rankss, ranksw, rankww, RUN3D, masterproc

    USE micro_task_util_MPI

    IMPLICIT NONE

    INTEGER, INTENT(IN)                  :: dimx1, dimx2, dimy1, dimy2, dimz
    REAL,    DIMENSION(dimx1:dimx2,dimy1:dimy2,dimz), INTENT(INOUT) :: f
    INTEGER, INTENT(IN)                  :: blen
    REAL,    DIMENSION(blen), INTENT(IN) :: buff
    INTEGER, INTENT(IN)                  :: tag
    INTEGER                              :: i, j, k, n, proc, i1, i2, j1, j2
!
!-- Decode the tag:
    i    = tag / 100000
    i1   = (tag - i * 100000) / 10000
    i2   = (tag - i * 100000 - i1 * 10000)/1000
    j1   = (tag - i * 100000 - i1 * 10000 - i2 * 1000) / 100
    j2   = (tag - i * 100000 - i1 * 10000 - i2 * 1000 - j1 * 100) / 10
    proc = tag - i * 100000 - i1 * 10000 - i2 * 1000 - j1 * 100 - j2 * 10

    SELECT CASE(proc)
       CASE(1)   ! From "North":
          n = 0
          DO k = 1, dimz
             DO j = nyp1, nyp1+j2
                DO i = 1, nx
                   n = n + 1
                   f(i,j,k) = buff(n)
                ENDDO
             ENDDO
          ENDDO
       CASE(2)   ! From "North-East":
          n = 0
          DO k = 1, dimz
             DO j = nyp1, nyp1+j2
                DO i = nxp1, nxp1+i2
                   n = n + 1
                   f(i,j,k) = buff(n)
                ENDDO
             ENDDO
          ENDDO
       CASE(3)   ! From "East":
          n = 0
          DO k = 1, dimz
             DO j = 1, ny
                DO i = nxp1, nxp1+i2
                   n = n + 1
                   f(i,j,k) = buff(n)
                ENDDO
             ENDDO
          ENDDO
       CASE(4)   ! From "South-East":
          n = 0
          DO k = 1, dimz
             DO j = -j1, 0
                DO i = nxp1, nxp1+i2
                   n = n + 1
                   f(i,j,k) = buff(n)
                ENDDO
             ENDDO
          ENDDO
       CASE(5)   ! From "South":
          n = 0
          DO k = 1, dimz
             DO j = -j1, 0
                DO i = 1, nx
                   n = n + 1
                   f(i,j,k) = buff(n)
                ENDDO
             ENDDO
          ENDDO
       CASE(6)   ! From "South-West":
          n = 0
          DO k = 1, dimz
             DO j = -j1, 0
                DO i = -i1, 0
                   n = n + 1
                   f(i,j,k) = buff(n)
                ENDDO
             ENDDO
          ENDDO
       CASE(7)   ! From "West":
          n = 0
          DO k = 1, dimz
             DO j = 1, ny
                DO i = -i1, 0
                   n = n + 1
                   f(i,j,k) = buff(n)
                ENDDO
             ENDDO
          ENDDO
       CASE(8)   ! From "North-West":
          n = 0
          DO k = 1, dimz
             DO j = nyp1, nyp1+j2
                DO i = -i1, 0
                   n = n + 1
                   f(i,j,k) = buff(n)
                ENDDO
             ENDDO
          ENDDO
    END SELECT

END SUBROUTINE task_assign_LCM
