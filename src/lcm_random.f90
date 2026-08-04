!> Standalone LCM random-number utilities.
!!
!! This module keeps random-number generation independent from SAM-specific
!! process state such as grid%rank. Callers pass the intended stream identifier
!! explicitly.
module lcm_random
   implicit none

   private

   public :: lcm_random_uniform
   public :: lcm_random_gaussian
   public :: lcm_random_truncated_gaussian
   public :: randomu

contains

   real function lcm_random_truncated_gaussian(seed_id, sigma) result(value)
      integer, intent(in) :: seed_id
      real, intent(in) :: sigma

      real, parameter :: nsigma = 2.0

      do
         value = lcm_random_gaussian(seed_id, sigma)
         if (abs(value) <= nsigma * sigma) exit
      enddo
   end function lcm_random_truncated_gaussian

   real function lcm_random_gaussian(seed_id, sigma) result(value)
      integer, intent(in) :: seed_id
      real, intent(in) :: sigma

      real :: radius_squared
      real :: v1
      real :: v2
      real, save :: cached_value
      logical, save :: has_cached_value = .false.

      if (has_cached_value) then
         value = cached_value
         has_cached_value = .false.
      else
         do
            v1 = 2.0 * lcm_random_uniform(seed_id) - 1.0
            v2 = 2.0 * lcm_random_uniform(seed_id) - 1.0
            radius_squared = v1**2 + v2**2
            if (radius_squared > 0.0 .and. radius_squared < 1.0) exit
         enddo

         radius_squared = sqrt(-2.0 * log(radius_squared) / radius_squared)
         value = v1 * radius_squared * sigma
         cached_value = v2 * radius_squared * sigma
         has_cached_value = .true.
      endif
   end function lcm_random_gaussian

   real function lcm_random_uniform(seed_id) result(value)
      integer, intent(in) :: seed_id

      integer, parameter :: k4b = selected_int_kind(9)
      integer(k4b), parameter :: ia = 16870
      integer(k4b), parameter :: im = 2147483647
      integer(k4b), parameter :: iq = 127773
      integer(k4b), parameter :: ir = 2836
      integer(k4b) :: k
      integer(k4b), save :: ix = -1
      integer(k4b), save :: iy = -1
      real, save :: am

      if (seed_id <= 0 .or. iy < 0) then
         am = nearest(1.0, -1.0) / im
         iy = ior(ieor(888889999, abs(seed_id)), 1)
         ix = ieor(777755555, abs(seed_id))
      endif

      ix = ieor(ix, ishft(ix, 13))
      ix = ieor(ix, ishft(ix, -17))
      ix = ieor(ix, ishft(ix, 5))
      k = iy / iq
      iy = ia * (iy - k * iq) - ir * k
      if (iy < 0) iy = iy + im

      value = am * ior(iand(im, ieor(ix, iy)), 1)
   end function lcm_random_uniform

   real function randomu(seed_id) result(value)
      integer, intent(in) :: seed_id

      value = lcm_random_uniform(seed_id)
   end function randomu

end module lcm_random
