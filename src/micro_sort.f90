 SUBROUTINE micro_sort_by_grid_box
!
!-- Droplets by their physical location, i.e., by their grid box
    USE micro_vars,   ONLY: n_particles, particle, n_start, n_end, particle_array, &
                            n_particles_alloc, nx_LCM, ny_LCM, nz_LCM, idx, idy,   &
                            nx_sb, ny_sb, gpxs, gpys, gpzs, gpzc, iX, iY, iZ,      &
                            dz_LCM, iID, nzl, search_sorted_array

    USE grid,         ONLY: RUN3D, rank, z

    IMPLICIT NONE

    INTEGER(KIND=8),      DIMENSION(:), ALLOCATABLE  :: n_boxid
    INTEGER,              DIMENSION(:), ALLOCATABLE  :: n_index
    TYPE(particle_array), DIMENSION(:), ALLOCATABLE  :: particle_temp

    INTEGER(KIND=8)  :: n_boxid_temp
    INTEGER          :: n_index_temp

    INTEGER          :: n, m, i, j, k, i_last, j_last, k_last, inc
    REAL             :: alpha, beta, gamma

    CALL t_startf ('micro_sort_by_grid_box')

    IF ( n_particles .GE. 1 )  THEN ! Do this only if particles are available

       n_particles_alloc = MAX( n_particles_alloc, n_particles )
       ALLOCATE(n_boxid(0:n_particles_alloc))
       ALLOCATE(n_index(0:n_particles_alloc))
!
!--    Assign n_boxid an individual id based on the grid box in which each particle
!--    is located
!
!--    i index
       n_index(1:n_particles) = MAX( MIN( FLOOR( particle(1:n_particles)%var(iX) * idx ) + 1 - nx_sb, nx_LCM ), 1 )
       n_boxid(1:n_particles) = INT(n_index(1:n_particles),KIND=8)
!
!--    j index
       IF (RUN3D) THEN
          n_index(1:n_particles) = MAX( MIN( FLOOR( particle(1:n_particles)%var(iY) * idy ) + 1 - ny_sb, ny_LCM ), 1 )
          n_boxid(1:n_particles) = 10000_8 * INT(n_index(1:n_particles),KIND=8) + n_boxid(1:n_particles)
       ELSE
          n_boxid(1:n_particles) = 10000_8 + n_boxid(1:n_particles)
       ENDIF
!
!--    k index
       DO n = 1, n_particles
          !n_index(n) = minloc(abs(z-particle(n)%var(iZ)), dim = 1)
          ! use a binary search to get z index (works for non-uniform grids, faster than minloc)
          IF (particle(n)%var(iZ) < gpzs(nzl+1)) THEN
            ! further reduce efforts by only searching levels where particles were placed
            n_index(n) = search_sorted_array(gpzs, particle(n)%var(iZ), nzl+1)
          ELSE
            n_index(n) = search_sorted_array(gpzs, particle(n)%var(iZ), nz_LCM+1)
          END IF
       ENDDO
       n_boxid(1:n_particles) = 100000000_8 * INT(n_index(1:n_particles),KIND=8) + n_boxid(1:n_particles)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!--    Now, sort particles by their grid box using n_boxid
       ALLOCATE(particle_temp(0:n_particles_alloc))

       DO n = 1, n_particles
          particle_temp(n)%ivar(:) = particle(n)%ivar(:)
          particle_temp(n)%var(:)  = particle(n)%var(:)
       ENDDO

       DO n = 1, n_particles
          n_index(n) = n
       ENDDO
   ! !
   ! !-- Do insertion sort
       ! DO n = 2, n_particles
       !    m = n - 1
       !    n_boxid_temp = n_boxid(n)
       !    n_index_temp = n_index(n)
       !    DO WHILE ( n_boxid(m) > n_boxid_temp )!.AND. m >= 1 )
       !       n_boxid(m+1) = n_boxid(m)
       !       n_index(m+1) = n_index(m)
       !       m = m - 1
       !       IF ( m == 0 )  EXIT
       !    END DO
       !    n_boxid(m+1) = n_boxid_temp
       !    n_index(m+1) = n_index_temp
       ! END DO

!
!--    Do Shell sort (from numrical recipes), more than 40 % faster than insertion sort
       inc = 1
1      inc = 100 * inc + 1
       IF ( inc .LE. n_particles )  GOTO 1
2      CONTINUE
       inc = inc / 100
       DO n = inc + 1, n_particles
          n_boxid_temp = n_boxid(n)
          n_index_temp = n_index(n)
          m = n
3         IF ( n_boxid(m-inc) .gt. n_boxid_temp ) THEN
             n_boxid(m) = n_boxid(m-inc)
             n_index(m) = n_index(m-inc)
             m = m - inc
             IF ( m .LE. inc )  GOTO 4
             GOTO 3
          ELSE IF (n_boxid(m-inc) .eq. n_boxid_temp) THEN
             IF ( particle(n_index(m-inc))%ivar(iID) .gt. particle(n_index_temp)%ivar(iID)) THEN
               n_boxid(m) = n_boxid(m-inc)
               n_index(m) = n_index(m-inc)
               m = m - inc
               IF ( m .LE. inc )  GOTO 4
               GOTO 3
             END IF
          ENDIF
4         n_boxid(m) = n_boxid_temp
          n_index(m) = n_index_temp
       ENDDO
       IF ( inc .GT. 1 )  GOTO 2

!
!--    Rearrange particles using the sorted tags
       DO n = 1, n_particles
          particle(n)%ivar(:) = particle_temp(n_index(n))%ivar(:)
          particle(n)%var(:)  = particle_temp(n_index(n))%var(:)
       ENDDO

       DEALLOCATE(particle_temp)
       DEALLOCATE(n_index)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!--    Finally, determine the arrays n_start, n_end. The defaults of 0 and -1 avoid
!--    the execution of n-loops if there are no super-droplets inside a grid box
       n_start = 0
       n_end   = -1

       ! CALL get_ijk( particle(1)%var(1), particle(1)%var(2), particle(1)%var(3),     &
       !               i, j, k, alpha, beta, gamma )

       k = INT( n_boxid(1) / 100000000_8, KIND=4)
       j = INT( ( n_boxid(1) - 100000000_8 * INT(k,KIND=8) ) / 10000_8, KIND=4)
       i = INT( n_boxid(1) - 100000000_8 * INT(k,KIND=8) - 10000_8 * INT(j,KIND=8), KIND=4)

       n_start(i,j,k) = 1

       k_last = k
       j_last = j
       i_last = i
       DO n = 2, n_particles
          ! CALL get_ijk( particle(n)%var(1), particle(n)%var(2), particle(n)%var(3),  &
          !               i, j, k, alpha, beta, gamma )

          k = INT( n_boxid(n) / 100000000_8, KIND=4)
          j = INT( ( n_boxid(n) - 100000000_8 * INT(k,KIND=8) ) / 10000_8, KIND=4)
          i = INT( n_boxid(n) - 100000000_8 * INT(k,KIND=8) - 10000_8 * INT(j,KIND=8), KIND=4)

          IF ( n_boxid(n) .NE. n_boxid(n-1) )  THEN
             n_end(i_last,j_last,k_last) = n-1
             n_start(i,j,k)              = n
          ENDIF
          k_last = k
          j_last = j
          i_last = i
       ENDDO
       n_end(i_last,j_last,k_last) = n_particles

       DEALLOCATE(n_boxid)

    ELSE

       n_start = 0
       n_end   = -1

    ENDIF

    CALL t_stopf ('micro_sort_by_grid_box')

 END SUBROUTINE micro_sort_by_grid_box

 SUBROUTINE get_ijk( xp, yp, zp, i, j, k, alpha, beta, gamma )

    USE micro_vars, ONLY: idx, idy, nx_sb, ny_sb, gpxs, gpys, gpzs, gpzc, lcmp, &
                          dz_LCM, nx_LCM, ny_LCM, nz_LCM, search_sorted_array, nzl

    USE grid,       ONLY: RUN3D, z

    IMPLICIT NONE

!
!-- Calculates the indices of the current grid box (i,j,k) in which the considered
!-- particle is located. Additionally, the nomralized distances the the grid box
!-- bounderies in x, y, and z direction (alpha, beta, gamma) are calculated.
    REAL,    INTENT(IN)  :: xp, yp, zp
    REAL,    INTENT(OUT) :: alpha, beta, gamma
    INTEGER, INTENT(OUT) :: i, j, k
!
!-- x direction
    i     = MAX( MIN( FLOOR( xp * idx ) + 1 - nx_sb, nx_LCM ), 1 )
    alpha = MAX( MIN( ( xp - gpxs(i) ) * idx, 1.0 ), 0.0 )
!
!-- y direction
    IF (RUN3D) THEN
       j    = MAX( MIN( FLOOR( yp * idy ) + 1 - ny_sb, ny_LCM ), 1 )
       beta = MAX( MIN( ( yp - gpys(j) ) * idy, 1.0 ), 0.0 )
    ELSE
       j    = 1
       beta = 0.0
    ENDIF
!
   !  k = MAX( MIN( FLOOR( zp / dz_LCM ) + 1, nz_LCM ), 1 )
    !k = minloc(abs(z-zp), dim = 1)
    ! use a binary search to get z index (works for non-uniform grids, faster than minloc)
    IF (zp<gpzs(nzl+1)) THEN
      ! further reduce efforts by only searching levels where particles were placed
      k = search_sorted_array(gpzs, zp, nzl+1)
    ELSE
      k = search_sorted_array(gpzs, zp, nz_LCM+1)
    END IF
    gamma = MAX( MIN( ( zp - gpzs(k) ) / ( gpzs(k+1) - gpzs(k) ), 1.0 ), 0.0 )

 END SUBROUTINE get_ijk
