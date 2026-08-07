module minimal_les_host_references
   use iso_fortran_env, only: int64, real64
   use minimal_les_host, only: minimal_host_config, minimal_host_density, &
                               minimal_host_dx, minimal_host_dy, &
                               minimal_host_dz, minimal_host_np_xy

   implicit none

   private

   public :: reference_bin_aerosol
   public :: reference_usm_aerosol
   public :: reference_random_aerosol_weight
   public :: reference_in_particle_weight
   public :: reference_cell_volume
   public :: reference_aerosol_radius
   public :: pi

   real(real64), parameter :: pi = 3.141592653589793238462643383279502884_real64

contains

   subroutine reference_bin_aerosol(k, expected_weight, expected_mass_aero, &
                                    expected_portions)
      integer, intent(in) :: k
      integer(int64), intent(out) :: expected_weight(minimal_host_np_xy)
      real(real64), intent(out) :: expected_mass_aero(minimal_host_np_xy)
      real(real64), intent(out) :: expected_portions(4, minimal_host_np_xy)

      integer :: m
      integer :: m_max
      integer :: n
      real(real64) :: dlogr
      real(real64) :: n_dummy
      real(real64) :: r_aero
      real(real64) :: r_aero_max
      real(real64) :: r_aero_min
      real(real64) :: r_l
      real(real64) :: r_r

      expected_weight = 0_int64
      expected_mass_aero = 0.0_real64
      expected_portions = 0.0_real64

      call reference_aerosol_radius_bounds(1, r_aero_min, r_aero_max)
      dlogr = (log(r_aero_max) - log(r_aero_min)) / &
              real(minimal_host_np_xy, real64)

      do n = 1, minimal_host_np_xy
         r_l = exp(log(r_aero_min) + real(n - 1, real64) * dlogr)
         r_r = exp(log(r_aero_min) + real(n, real64) * dlogr)
         r_aero = sqrt(r_l * r_r)

         n_dummy = reference_log_normal_pdf( &
            real(minimal_host_config%n_aero(1), real64), &
            real(minimal_host_config%sigma_aero(1), real64), &
            real(minimal_host_config%rm_aero(1), real64), r_aero)
         expected_weight(n) = int(n_dummy * r_aero * dlogr * &
                                  reference_cell_volume(k) * &
                                  minimal_host_density(k), int64)
         expected_mass_aero(n) = reference_aerosol_mass(expected_weight(n), r_aero)
         expected_portions(1, n) = 1.0_real64
      enddo

      do n = 1, minimal_host_np_xy
         if (expected_weight(n) == 0_int64) then
            m_max = n
            do m = 1, minimal_host_np_xy
               if (expected_weight(m) > expected_weight(m_max) .and. &
                   expected_weight(m) >= 2_int64) then
                  m_max = m
               endif
            enddo

            if (expected_weight(m_max) == 0_int64) then
               error stop "reference initialization of aerosols not possible"
            endif

            expected_weight(n) = expected_weight(m_max) / 2_int64
            expected_weight(m_max) = expected_weight(m_max) - expected_weight(n)
            expected_mass_aero(n) = expected_mass_aero(m_max) / &
                                    real(expected_weight(m_max) + &
                                         expected_weight(n), real64) * &
                                    real(expected_weight(n), real64)
            expected_mass_aero(m_max) = expected_mass_aero(m_max) - &
                                        expected_mass_aero(n)

            if (expected_weight(n) == 0_int64 .or. &
                expected_weight(m_max) == 0_int64) then
               error stop "reference initialization of aerosols not possible"
            endif
         endif
      enddo
   end subroutine reference_bin_aerosol

   subroutine reference_usm_aerosol(k, expected_weight, expected_mass_aero, &
                                    expected_portions)
      integer, intent(in) :: k
      integer(int64), intent(out) :: expected_weight(minimal_host_np_xy)
      real(real64), intent(out) :: expected_mass_aero(minimal_host_np_xy)
      real(real64), intent(out) :: expected_portions(4, minimal_host_np_xy)

      integer :: n
      real(real64) :: p_usm
      real(real64) :: r_aero
      real(real64) :: r_aero_max
      real(real64) :: r_aero_min
      real(real64) :: r_l
      real(real64) :: r_r

      expected_weight = 0_int64
      expected_mass_aero = 0.0_real64
      expected_portions = 0.0_real64

      call reference_aerosol_radius_bounds(1, r_aero_min, r_aero_max)
      p_usm = 1.0_real64 / real(minimal_host_np_xy, real64)
      r_l = r_aero_min

      do n = 1, minimal_host_np_xy - 1
         r_r = reference_percentile_radius(1, p_usm * real(n, real64), &
                                           r_l, r_aero_max)
         r_aero = sqrt(r_l * r_r)
         r_l = r_r
         expected_weight(n) = reference_equal_mode_weight(k)
         expected_mass_aero(n) = reference_aerosol_mass(expected_weight(n), r_aero)
         expected_portions(1, n) = 1.0_real64
      enddo

      n = minimal_host_np_xy
      r_aero = sqrt(r_l * r_aero_max)
      expected_weight(n) = int(reference_cell_volume(k) * &
                               sum(real(minimal_host_config%n_aero, real64)) / &
                               real(minimal_host_np_xy, real64), int64)
      expected_mass_aero(n) = reference_aerosol_mass(expected_weight(n), r_aero)
   end subroutine reference_usm_aerosol

   integer(int64) function reference_random_aerosol_weight(k) result(weight)
      integer, intent(in) :: k

      weight = reference_equal_mode_weight(k)
   end function reference_random_aerosol_weight

   integer(int64) function reference_in_particle_weight(k) result(weight)
      integer, intent(in) :: k

      weight = int(real(minimal_host_config%n_ice(3), real64) * &
                   reference_cell_volume(k) * minimal_host_density(k) / &
                   (real(minimal_host_np_xy, real64) * &
                    real(minimal_host_config%frac_ice, real64)), int64)
   end function reference_in_particle_weight

   real(real64) function reference_cell_volume(k) result(volume)
      integer, intent(in) :: k

      volume = minimal_host_dx * minimal_host_dy * minimal_host_dz(k)
   end function reference_cell_volume

   real(real64) function reference_aerosol_radius(weight, mass_aero) result(radius)
      integer(int64), intent(in) :: weight
      real(real64), intent(in) :: mass_aero
      real(real64) :: volume

      if (weight <= 0_int64) then
         error stop "reference aerosol radius requires positive weight"
      endif
      volume = mass_aero / (real(minimal_host_config%rho_aero(1), real64) * &
                            real(weight, real64))
      radius = (volume * 3.0_real64 / (4.0_real64 * pi))**(1.0_real64 / 3.0_real64)
   end function reference_aerosol_radius

   integer(int64) function reference_equal_mode_weight(k) result(weight)
      integer, intent(in) :: k

      weight = int(reference_cell_volume(k) * minimal_host_density(k) * &
                   real(minimal_host_config%n_aero(1), real64) / &
                   real(minimal_host_np_xy, real64), int64)
   end function reference_equal_mode_weight

   real(real64) function reference_aerosol_mass(weight, radius) result(mass)
      integer(int64), intent(in) :: weight
      real(real64), intent(in) :: radius

      mass = real(weight, real64) * 4.0_real64 / 3.0_real64 * pi * &
             radius**3 * real(minimal_host_config%rho_aero(1), real64)
   end function reference_aerosol_mass

   subroutine reference_aerosol_radius_bounds(i_mode, r_aero_min, r_aero_max)
      integer, intent(in) :: i_mode
      real(real64), intent(out) :: r_aero_min
      real(real64), intent(out) :: r_aero_max

      r_aero_min = reference_percentile_radius(i_mode, 5.0e-5_real64, &
                                               1.0e-9_real64, 1.0e-5_real64)
      r_aero_max = reference_percentile_radius(i_mode, 1.0_real64 - 5.0e-5_real64, &
                                               r_aero_min, 1.0e-5_real64)
   end subroutine reference_aerosol_radius_bounds

   real(real64) function reference_percentile_radius(i_mode, percentile, &
                                                     r_left, r_right) result(radius)
      integer, intent(in) :: i_mode
      real(real64), intent(in) :: percentile
      real(real64), intent(in) :: r_left
      real(real64), intent(in) :: r_right

      real(real64) :: ak
      real(real64) :: bk
      real(real64) :: ck
      real(real64) :: f_ak
      real(real64) :: f_bk
      real(real64) :: f_ck

      ak = r_left
      f_ak = reference_log_normal_cdf( &
         real(minimal_host_config%n_aero(i_mode), real64), &
         real(minimal_host_config%sigma_aero(i_mode), real64), &
         real(minimal_host_config%rm_aero(i_mode), real64), ak) / &
         real(minimal_host_config%n_aero(i_mode), real64) - percentile

      bk = r_right
      f_bk = reference_log_normal_cdf( &
         real(minimal_host_config%n_aero(i_mode), real64), &
         real(minimal_host_config%sigma_aero(i_mode), real64), &
         real(minimal_host_config%rm_aero(i_mode), real64), bk) / &
         real(minimal_host_config%n_aero(i_mode), real64) - percentile

      do
         ck = (ak * f_bk - bk * f_ak) / (f_bk - f_ak)
         f_ck = reference_log_normal_cdf( &
            real(minimal_host_config%n_aero(i_mode), real64), &
            real(minimal_host_config%sigma_aero(i_mode), real64), &
            real(minimal_host_config%rm_aero(i_mode), real64), ck) / &
            real(minimal_host_config%n_aero(i_mode), real64) - percentile

         if (f_ak * f_ck > 0.0_real64) then
            ak = ck
            f_ak = f_ck
         else
            bk = ck
            f_bk = f_ck
         endif

         if (abs(f_ck) < 1.0e-12_real64) exit
      enddo
      radius = ck
   end function reference_percentile_radius

   real(real64) function reference_log_normal_pdf(n, s, m, x) result(pdf)
      real(real64), intent(in) :: n
      real(real64), intent(in) :: s
      real(real64), intent(in) :: m
      real(real64), intent(in) :: x

      pdf = 0.0_real64
      if (n /= 0.0_real64 .and. s /= 0.0_real64 .and. &
          m /= 0.0_real64 .and. x > 0.0_real64) then
         pdf = n / (sqrt(2.0_real64 * pi) * log(s) * x) * &
               exp(-1.0_real64 * log(x / m)**2 / &
                   (2.0_real64 * log(s)**2))
      endif
   end function reference_log_normal_pdf

   real(real64) function reference_log_normal_cdf(n, s, m, x) result(cdf)
      real(real64), intent(in) :: n
      real(real64), intent(in) :: s
      real(real64), intent(in) :: m
      real(real64), intent(in) :: x

      cdf = 0.0_real64
      if (n /= 0.0_real64 .and. s /= 0.0_real64 .and. &
          m /= 0.0_real64 .and. x > 0.0_real64) then
         cdf = n * 0.5_real64 * &
               (1.0_real64 + erf(log(x / m) / &
                (sqrt(2.0_real64) * log(s))))
      endif
   end function reference_log_normal_cdf

end module minimal_les_host_references
