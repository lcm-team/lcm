program test_minimal_les_host_initialization
   use iso_fortran_env, only: int64, real64
   use lcm_initialization_state, only: lcm_initial_particle_count_global, &
                                       lcm_initial_particle_count_local
   use lcm_particle_state, only: cell_particles
   use minimal_les_host, only: minimal_host_config, minimal_host_dx, &
                               minimal_host_dy, minimal_host_np_xy, &
                               minimal_host_nx_local, minimal_host_ny_local, &
                               minimal_host_nz, minimal_host_nzl, &
                               minimal_host_u, minimal_host_v, minimal_host_w, &
                               minimal_host_z_face, setup_minimal_les_host, &
                               advance_minimal_les_transport
   use minimal_les_host_references, only: pi, reference_aerosol_radius, &
                                          reference_bin_aerosol,        &
                                          reference_in_particle_weight, &
                                          reference_random_aerosol_weight, &
                                          reference_usm_aerosol

   implicit none

   integer, parameter :: liquid_particle_type = 1
   integer, parameter :: ice_nucleating_particle_type = 2
   integer(int64), parameter :: expected_count = &
      int(minimal_host_nx_local, int64) * &
      int(minimal_host_ny_local, int64) * &
      int(minimal_host_nzl, int64) * &
      int(minimal_host_np_xy, int64)

   character(len=32) :: case_name
   real(real64), allocatable :: initial_x(:,:,:,:)
   real(real64), allocatable :: initial_y(:,:,:,:)
   real(real64), allocatable :: initial_z(:,:,:,:)

   call get_command_argument(1, case_name)
   if (len_trim(case_name) == 0) case_name = "bin"

   call setup_minimal_les_host(trim(case_name))

   call assert_true(lcm_initial_particle_count_local == expected_count, &
                    "minimal host local particle count")
   call assert_true(lcm_initial_particle_count_global == expected_count, &
                    "minimal host global particle count")
   call assert_particle_cells(trim(case_name))
   call store_initial_particle_positions()
   call advance_minimal_les_transport(1.0_real64)
   call assert_particle_velocities_and_positions(1.0_real64)

contains

   subroutine assert_particle_cells(case_name)
      character(len=*), intent(in) :: case_name
      integer :: i, j, k, n
      integer :: total_count

      total_count = 0
      do k = 1, minimal_host_nz
         do j = 1, minimal_host_ny_local
            do i = 1, minimal_host_nx_local
               if (k <= minimal_host_nzl) then
                  call assert_true(cell_particles(i,j,k)%count == &
                                   minimal_host_np_xy, &
                                   "minimal host initialized cell count")
                  do n = 1, cell_particles(i,j,k)%count
                     total_count = total_count + 1
                     call assert_particle(i, j, k, n, case_name)
                  enddo
               else
                  call assert_true(cell_particles(i,j,k)%count == 0, &
                                   "minimal host cell above nzl is empty")
               endif
            enddo
         enddo
      enddo

      call assert_true(total_count == expected_count, &
                       "minimal host summed particle count")
      call assert_unique_particle_ids()
   end subroutine assert_particle_cells

   subroutine assert_particle(i, j, k, n, case_name)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n
      character(len=*), intent(in) :: case_name

      call assert_particle_metadata(i, j, k, n)
      call assert_particle_position(i, j, k, n)

      select case (trim(case_name))
      case ("bin")
         call assert_exact_aerosol_state(i, j, k, n, reference_mode="bin")
         call assert_liquid_state(i, j, k, n)
      case ("USM")
         call assert_exact_aerosol_state(i, j, k, n, reference_mode="USM")
         call assert_liquid_state(i, j, k, n)
      case ("random")
         call assert_random_aerosol_state(i, j, k, n)
         call assert_liquid_state(i, j, k, n)
      case ("bin_ice")
         call assert_ice_state(i, j, k, n)
      case default
         error stop "unknown minimal LES host initialization assertion case"
      end select
   end subroutine assert_particle

   subroutine assert_particle_metadata(i, j, k, n)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n

      call assert_true(cell_particles(i,j,k)%particles(n)%cell_particle_index == n, &
                       "minimal host cell-local particle index")
      call assert_true(cell_particles(i,j,k)%particles(n)%last_i == i, &
                       "minimal host last i index")
      call assert_true(cell_particles(i,j,k)%particles(n)%last_j == j, &
                       "minimal host last j index")
      call assert_true(cell_particles(i,j,k)%particles(n)%last_k == k, &
                       "minimal host last k index")
      call assert_true(cell_particles(i,j,k)%particles(n)%w_sedi == 0.0_real64, &
                       "minimal host sedimentation velocity")
      call assert_true(cell_particles(i,j,k)%particles(n)%spray_marker == -999, &
                       "minimal host default spray marker")
      call assert_true(cell_particles(i,j,k)%particles(n)%track_marker == -999, &
                       "minimal host default tracking marker")
   end subroutine assert_particle_metadata

   subroutine assert_particle_position(i, j, k, n)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n

      call assert_true(cell_particles(i,j,k)%particles(n)%x >= &
                       real(i - 1, real64) * minimal_host_dx, &
                       "minimal host particle x lower bound")
      call assert_true(cell_particles(i,j,k)%particles(n)%x < &
                       real(i, real64) * minimal_host_dx, &
                       "minimal host particle x upper bound")
      call assert_true(cell_particles(i,j,k)%particles(n)%y >= &
                       real(j - 1, real64) * minimal_host_dy, &
                       "minimal host particle y lower bound")
      call assert_true(cell_particles(i,j,k)%particles(n)%y < &
                       real(j, real64) * minimal_host_dy, &
                       "minimal host particle y upper bound")
      call assert_true(cell_particles(i,j,k)%particles(n)%z >= &
                       minimal_host_z_face(k), &
                       "minimal host particle z lower bound")
      call assert_true(cell_particles(i,j,k)%particles(n)%z < &
                       minimal_host_z_face(k + 1), &
                       "minimal host particle z upper bound")
   end subroutine assert_particle_position

   subroutine store_initial_particle_positions()
      integer :: i
      integer :: j
      integer :: k
      integer :: n

      allocate(initial_x(minimal_host_nx_local, minimal_host_ny_local, &
                         minimal_host_nzl, minimal_host_np_xy))
      allocate(initial_y(minimal_host_nx_local, minimal_host_ny_local, &
                         minimal_host_nzl, minimal_host_np_xy))
      allocate(initial_z(minimal_host_nx_local, minimal_host_ny_local, &
                         minimal_host_nzl, minimal_host_np_xy))

      do k = 1, minimal_host_nzl
         do j = 1, minimal_host_ny_local
            do i = 1, minimal_host_nx_local
               do n = 1, cell_particles(i,j,k)%count
                  initial_x(i,j,k,n) = cell_particles(i,j,k)%particles(n)%x
                  initial_y(i,j,k,n) = cell_particles(i,j,k)%particles(n)%y
                  initial_z(i,j,k,n) = cell_particles(i,j,k)%particles(n)%z
               enddo
            enddo
         enddo
      enddo
   end subroutine store_initial_particle_positions

   subroutine assert_particle_velocities_and_positions(time_step)
      real(real64), intent(in) :: time_step

      integer :: i
      integer :: j
      integer :: k
      integer :: n
      real(real64) :: expected_u
      real(real64) :: expected_v
      real(real64) :: expected_w

      do k = 1, minimal_host_nzl
         do j = 1, minimal_host_ny_local
            do i = 1, minimal_host_nx_local
               do n = 1, cell_particles(i,j,k)%count
                  call expected_particle_velocity( &
                     initial_x(i,j,k,n), initial_y(i,j,k,n), &
                     initial_z(i,j,k,n), &
                     expected_u, expected_v, expected_w)
                  call assert_near(cell_particles(i,j,k)%particles(n)%u, &
                                   expected_u, &
                                   "minimal host prepared particle u")
                  call assert_near(cell_particles(i,j,k)%particles(n)%v, &
                                   expected_v, &
                                   "minimal host prepared particle v")
                  call assert_near(cell_particles(i,j,k)%particles(n)%w, &
                                   expected_w, &
                                   "minimal host prepared particle w")
                  call assert_near(cell_particles(i,j,k)%particles(n)%x, &
                                   initial_x(i,j,k,n) + expected_u * time_step, &
                                   "minimal host particle x after micro_xyz")
                  call assert_near(cell_particles(i,j,k)%particles(n)%y, &
                                   initial_y(i,j,k,n) + expected_v * time_step, &
                                   "minimal host particle y after micro_xyz")
                  call assert_near(cell_particles(i,j,k)%particles(n)%z, &
                                   initial_z(i,j,k,n) + expected_w * time_step, &
                                   "minimal host particle z after micro_xyz")
               enddo
            enddo
         enddo
      enddo
   end subroutine assert_particle_velocities_and_positions

   subroutine expected_particle_velocity(x, y, z, expected_u, expected_v, &
                                         expected_w)
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      real(real64), intent(in) :: z
      real(real64), intent(out) :: expected_u
      real(real64), intent(out) :: expected_v
      real(real64), intent(out) :: expected_w

      integer :: iteration
      real(real64) :: predicted_u
      real(real64) :: predicted_v
      real(real64) :: predicted_w
      real(real64) :: unow
      real(real64) :: vnow
      real(real64) :: wnow
      real(real64) :: xp
      real(real64) :: yp
      real(real64) :: zp

      call interpolate_expected_velocity(x, y, z, expected_u, expected_v, &
                                         expected_w)
      unow = expected_u
      vnow = expected_v
      wnow = expected_w

      do iteration = 1, minimal_host_config%n_uvw_iteration
         xp = x + expected_u
         yp = y + expected_v
         zp = z + expected_w
         call interpolate_expected_velocity(xp, yp, zp, predicted_u, predicted_v, &
                                            predicted_w)
         expected_u = 0.5_real64 * (unow + predicted_u)
         expected_v = 0.5_real64 * (vnow + predicted_v)
         expected_w = 0.5_real64 * (wnow + predicted_w)
      enddo
   end subroutine expected_particle_velocity

   subroutine interpolate_expected_velocity(x, y, z, expected_u, expected_v, &
                                            expected_w)
      real(real64), intent(in) :: x
      real(real64), intent(in) :: y
      real(real64), intent(in) :: z
      real(real64), intent(out) :: expected_u
      real(real64), intent(out) :: expected_v
      real(real64), intent(out) :: expected_w

      integer :: i
      integer :: j
      integer :: k
      real(real64) :: alpha
      real(real64) :: beta
      real(real64) :: gamma

      i = max(min(floor(x / minimal_host_dx) + 1, minimal_host_nx_local), 1)
      j = max(min(floor(y / minimal_host_dy) + 1, minimal_host_ny_local), 1)
      k = expected_vertical_index(z)

      alpha = max(min((x - real(i - 1, real64) * minimal_host_dx) / &
                      minimal_host_dx, 1.0_real64), 0.0_real64)
      beta = max(min((y - real(j - 1, real64) * minimal_host_dy) / &
                     minimal_host_dy, 1.0_real64), 0.0_real64)
      gamma = max(min((z - minimal_host_z_face(k)) / &
                      (minimal_host_z_face(k + 1) - minimal_host_z_face(k)), &
                      1.0_real64), 0.0_real64)

      expected_u = minimal_host_u(i,j,k) * (1.0_real64 - alpha) + &
                   minimal_host_u(i + 1,j,k) * alpha
      expected_v = minimal_host_v(i,j,k) * (1.0_real64 - beta) + &
                   minimal_host_v(i,j + 1,k) * beta
      expected_w = minimal_host_w(i,j,k) * (1.0_real64 - gamma) + &
                   minimal_host_w(i,j,k + 1) * gamma
   end subroutine interpolate_expected_velocity

   integer function expected_vertical_index(z) result(index)
      real(real64), intent(in) :: z

      integer :: kk

      index = minimal_host_nz
      do kk = 1, minimal_host_nz
         if (z < minimal_host_z_face(kk + 1)) then
            index = kk
            exit
         endif
      enddo
      index = max(min(index, minimal_host_nz), 1)
   end function expected_vertical_index

   subroutine assert_exact_aerosol_state(i, j, k, n, reference_mode)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n
      character(len=*), intent(in) :: reference_mode

      integer(int64) :: expected_weight(minimal_host_np_xy)
      real(real64) :: expected_mass_aero(minimal_host_np_xy)
      real(real64) :: expected_portions(4, minimal_host_np_xy)

      select case (trim(reference_mode))
      case ("bin")
         call reference_bin_aerosol(k, expected_weight, expected_mass_aero, &
                                    expected_portions)
      case ("USM")
         call reference_usm_aerosol(k, expected_weight, expected_mass_aero, &
                                    expected_portions)
      case default
         error stop "unknown exact aerosol reference mode"
      end select

      call assert_true(cell_particles(i,j,k)%particles(n)%weight == &
                       expected_weight(n), &
                       "minimal host exact aerosol weight")
      call assert_near_relative(cell_particles(i,j,k)%particles(n)%mass_aero, &
                                expected_mass_aero(n), &
                                "minimal host exact aerosol mass")
      call assert_near_array(cell_particles(i,j,k)%particles(n)%aerosol_portions, &
                             expected_portions(:, n), &
                             "minimal host exact aerosol portions")
      call assert_liquid_particle_type(i, j, k, n)
   end subroutine assert_exact_aerosol_state

   subroutine assert_random_aerosol_state(i, j, k, n)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n

      call assert_true(cell_particles(i,j,k)%particles(n)%weight == &
                       reference_random_aerosol_weight(k), &
                       "minimal host exact random aerosol weight")
      call assert_true(cell_particles(i,j,k)%particles(n)%mass_aero > 0.0_real64, &
                       "minimal host random aerosol mass")
      call assert_near_array(cell_particles(i,j,k)%particles(n)%aerosol_portions, &
                             [1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64], &
                             "minimal host random aerosol portions")
      call assert_liquid_particle_type(i, j, k, n)
   end subroutine assert_random_aerosol_state

   subroutine assert_liquid_particle_type(i, j, k, n)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n

      call assert_true(cell_particles(i,j,k)%particles(n)%particle_type == &
                       liquid_particle_type, &
                       "minimal host liquid particle type")
      call assert_true(maxval(abs( &
                       cell_particles(i,j,k)%particles(n)%in_surface_area)) == &
                       0.0_real64, &
                       "minimal host no IN surface area")
      call assert_true(cell_particles(i,j,k)%particles(n)%freezing_temperature == &
                       233.15_real64, &
                       "minimal host default freezing temperature")
   end subroutine assert_liquid_particle_type

   subroutine assert_liquid_state(i, j, k, n)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n
      real(real64) :: dry_radius
      real(real64) :: dry_equivalent_liquid_mass

      dry_radius = reference_aerosol_radius( &
         cell_particles(i,j,k)%particles(n)%weight, &
         cell_particles(i,j,k)%particles(n)%mass_aero)
      dry_equivalent_liquid_mass = &
         real(cell_particles(i,j,k)%particles(n)%weight, real64) * &
         4.0_real64 / 3.0_real64 * pi * 1000.0_real64 * dry_radius**3

      call assert_true(cell_particles(i,j,k)%particles(n)%mass_water >= &
                       dry_equivalent_liquid_mass, &
                       "minimal host liquid mass covers dry radius")
      call assert_true(cell_particles(i,j,k)%particles(n)%n_monomer == 1_int64, &
                       "minimal host monomer count")
      call assert_true(cell_particles(i,j,k)%particles(n)%temperature > &
                       250.0_real64 .and. &
                       cell_particles(i,j,k)%particles(n)%temperature < &
                       330.0_real64, &
                       "minimal host plausible particle temperature")
      call assert_true(cell_particles(i,j,k)%particles(n)%qv > 0.0_real64, &
                       "minimal host positive particle qv")
      call assert_true(abs(cell_particles(i,j,k)%particles(n)%eta) < &
                       1.0_real64, &
                       "minimal host finite particle eta")
   end subroutine assert_liquid_state

   subroutine assert_ice_state(i, j, k, n)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n

      call assert_true(cell_particles(i,j,k)%particles(n)%particle_type == &
                       ice_nucleating_particle_type, &
                       "minimal host IN particle type")
      call assert_true(cell_particles(i,j,k)%particles(n)%weight == &
                       reference_in_particle_weight(k), &
                       "minimal host exact IN particle weight")
      call assert_true(cell_particles(i,j,k)%particles(n)%mass_aero > 0.0_real64, &
                       "minimal host positive IN aerosol mass")
      call assert_true(cell_particles(i,j,k)%particles(n)%in_surface_area(3) > &
                       0.0_real64, &
                       "minimal host active IN surface area")
      call assert_true(maxval(abs( &
                       cell_particles(i,j,k)%particles(n)%in_surface_area([1,2,4]))) == &
                       0.0_real64, &
                       "minimal host inactive IN surface areas")
      call assert_true(cell_particles(i,j,k)%particles(n)%freezing_temperature >= &
                       233.15_real64 .and. &
                       cell_particles(i,j,k)%particles(n)%freezing_temperature <= &
                       273.15_real64, &
                       "minimal host IN freezing temperature range")
   end subroutine assert_ice_state

   subroutine assert_unique_particle_ids()
      integer :: i, j, k, n
      integer :: ii, jj, kk, nn
      integer(int64) :: particle_id

      do k = 1, minimal_host_nzl
         do j = 1, minimal_host_ny_local
            do i = 1, minimal_host_nx_local
               do n = 1, cell_particles(i,j,k)%count
                  particle_id = cell_particles(i,j,k)%particles(n)%id
                  call assert_true(particle_id >= 0_int64, &
                                   "minimal host nonnegative particle id")
                  do kk = k, minimal_host_nzl
                     do jj = 1, minimal_host_ny_local
                        do ii = 1, minimal_host_nx_local
                           do nn = 1, cell_particles(ii,jj,kk)%count
                              if (same_particle_position(i, j, k, n, &
                                                         ii, jj, kk, nn)) cycle
                              call assert_true( &
                                 particle_id /= &
                                 cell_particles(ii,jj,kk)%particles(nn)%id, &
                                 "minimal host unique particle id")
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine assert_unique_particle_ids

   logical function same_particle_position(i, j, k, n, ii, jj, kk, nn) &
      result(is_same)
      integer, intent(in) :: i, j, k, n
      integer, intent(in) :: ii, jj, kk, nn

      is_same = i == ii .and. j == jj .and. k == kk .and. n == nn
   end function same_particle_position

   subroutine assert_near_array(actual, expected, message)
      real(real64), intent(in) :: actual(:)
      real(real64), intent(in) :: expected(:)
      character(len=*), intent(in) :: message

      call assert_true(size(actual) == size(expected), message // " shape")
      call assert_true(maxval(abs(actual - expected)) <= 1.0e-12_real64, message)
   end subroutine assert_near_array

   subroutine assert_near(actual, expected, message)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: message

      call assert_true(abs(actual - expected) <= 1.0e-12_real64, message)
   end subroutine assert_near

   subroutine assert_near_relative(actual, expected, message)
      real(real64), intent(in) :: actual
      real(real64), intent(in) :: expected
      character(len=*), intent(in) :: message
      real(real64) :: scale

      scale = max(1.0e-30_real64, abs(expected))
      call assert_true(abs(actual - expected) / scale <= 1.0e-10_real64, message)
   end subroutine assert_near_relative

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         error stop "assertion failed: " // message
      endif
   end subroutine assert_true

end program test_minimal_les_host_initialization
