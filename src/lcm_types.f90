module lcm_types
   use iso_fortran_env, only: real64

   implicit none

   private

   public :: lcm_grid_t
   public :: lcm_init_grid

   type :: lcm_grid_t
      integer :: nx_local = 0
      integer :: ny_local = 0
      integer :: nz = 0

      integer :: nx_global = 0
      integer :: ny_global = 0

      integer :: i_start = 0
      integer :: i_end = 0
      integer :: j_start = 0
      integer :: j_end = 0

      real(real64) :: dx = 0.0_real64
      real(real64) :: dy = 0.0_real64

      real(real64), allocatable :: z_center(:)
      real(real64), allocatable :: z_face(:)
      real(real64), allocatable :: dz(:)

      integer :: halo_x = 0
      integer :: halo_y = 0
      integer :: halo_z = 0

      integer :: subdomain_x_index = 0
      integer :: subdomain_y_index = 0
      integer :: n_subdomains_x = 1
      integer :: n_subdomains_y = 1

      integer :: global_i_start = 0
      integer :: global_j_start = 0

      integer :: west_rank = -1
      integer :: east_rank = -1
      integer :: south_rank = -1
      integer :: north_rank = -1

      logical :: periodic_x = .false.
      logical :: periodic_y = .false.
   end type lcm_grid_t

contains

   subroutine lcm_init_grid(grid, nx_local, ny_local, nz, nx_global, ny_global, &
                            i_start, i_end, j_start, j_end, dx, dy, z_center, &
                            z_face, dz, halo_x, halo_y, halo_z,                &
                            subdomain_x_index, subdomain_y_index,             &
                            n_subdomains_x, n_subdomains_y, global_i_start,   &
                            global_j_start, west_rank, east_rank, south_rank, &
                            north_rank, periodic_x, periodic_y)
      type(lcm_grid_t), intent(out) :: grid
      integer, intent(in) :: nx_local
      integer, intent(in) :: ny_local
      integer, intent(in) :: nz
      integer, intent(in) :: nx_global
      integer, intent(in) :: ny_global
      integer, intent(in) :: i_start
      integer, intent(in) :: i_end
      integer, intent(in) :: j_start
      integer, intent(in) :: j_end
      real(real64), intent(in) :: dx
      real(real64), intent(in) :: dy
      real(real64), intent(in) :: z_center(:)
      real(real64), intent(in) :: z_face(:)
      real(real64), intent(in) :: dz(:)
      integer, intent(in) :: halo_x
      integer, intent(in) :: halo_y
      integer, intent(in) :: halo_z
      integer, intent(in) :: subdomain_x_index
      integer, intent(in) :: subdomain_y_index
      integer, intent(in) :: n_subdomains_x
      integer, intent(in) :: n_subdomains_y
      integer, intent(in) :: global_i_start
      integer, intent(in) :: global_j_start
      integer, intent(in) :: west_rank
      integer, intent(in) :: east_rank
      integer, intent(in) :: south_rank
      integer, intent(in) :: north_rank
      logical, intent(in) :: periodic_x
      logical, intent(in) :: periodic_y

      grid%nx_local = nx_local
      grid%ny_local = ny_local
      grid%nz = nz

      grid%nx_global = nx_global
      grid%ny_global = ny_global

      grid%i_start = i_start
      grid%i_end = i_end
      grid%j_start = j_start
      grid%j_end = j_end

      grid%dx = dx
      grid%dy = dy

      allocate(grid%z_center(size(z_center)))
      allocate(grid%z_face(size(z_face)))
      allocate(grid%dz(size(dz)))

      grid%z_center = z_center
      grid%z_face = z_face
      grid%dz = dz

      grid%halo_x = halo_x
      grid%halo_y = halo_y
      grid%halo_z = halo_z

      grid%subdomain_x_index = subdomain_x_index
      grid%subdomain_y_index = subdomain_y_index
      grid%n_subdomains_x = n_subdomains_x
      grid%n_subdomains_y = n_subdomains_y

      grid%global_i_start = global_i_start
      grid%global_j_start = global_j_start

      grid%west_rank = west_rank
      grid%east_rank = east_rank
      grid%south_rank = south_rank
      grid%north_rank = north_rank

      grid%periodic_x = periodic_x
      grid%periodic_y = periodic_y
   end subroutine lcm_init_grid

end module lcm_types
