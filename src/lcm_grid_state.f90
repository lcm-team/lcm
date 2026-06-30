!> LCM-owned grid state for the current structured-grid implementation.
module lcm_grid_state
   use lcm_types, only: lcm_grid_t

   implicit none

   private

   public :: t_lcm_grid_state
   public :: grid_state
   public :: lcm_set_grid

   type :: t_lcm_grid_state
      type(lcm_grid_t) :: grid
      logical :: is_set = .false.

      integer :: ux1 = 3
      integer :: ux2 = 4
      integer :: uy1 = 0
      integer :: uy2 = 0
      integer :: vx1 = 0
      integer :: vx2 = 4
      integer :: vy1 = 0
      integer :: vy2 = 0
      integer :: wx1 = 4
      integer :: wx2 = 4
      integer :: wy1 = 0
      integer :: wy2 = 0
      integer :: sx1 = 4
      integer :: sx2 = 4
      integer :: sy1 = 0
      integer :: sy2 = 0

      integer :: x1u = 0
      integer :: x2u = 0
      integer :: y1u = 0
      integer :: y2u = 0
      integer :: x1v = 0
      integer :: x2v = 0
      integer :: y1v = 0
      integer :: y2v = 0
      integer :: x1w = 0
      integer :: x2w = 0
      integer :: y1w = 0
      integer :: y2w = 0
      integer :: x1s = 0
      integer :: x2s = 0
      integer :: y1s = 0
      integer :: y2s = 0

      integer :: nx_sb = 0
      integer :: ny_sb = 0
      integer :: i_global_LCM = 0
      integer :: j_global_LCM = 0

      real, allocatable :: gpxc(:)
      real, allocatable :: gpyc(:)
      real, allocatable :: gpzc(:)
      real, allocatable :: gpxs(:)
      real, allocatable :: gpys(:)
      real, allocatable :: gpzs(:)
      real, allocatable :: dV(:)

      real :: xwidth = 0.0
      real :: ywidth = 0.0
      real :: idx = 0.0
      real :: idy = 0.0
      real :: ztop = 0.0
      real :: zbottom = 0.0
   end type t_lcm_grid_state

   type(t_lcm_grid_state), protected, save :: grid_state

contains

   subroutine lcm_set_grid(host_grid)
      type(lcm_grid_t), intent(in) :: host_grid

      if (grid_state%is_set) then
         error stop "LCM grid has already been initialized"
      endif

      grid_state%grid = host_grid
      grid_state%is_set = .true.

      call initialize_geometry()
   end subroutine lcm_set_grid

   subroutine initialize_geometry()
      integer :: i, j
      integer :: y_dim

      y_dim = grid_state%grid%y_dimension

      grid_state%uy1 = 3 * y_dim + 1
      grid_state%uy2 = 3 * y_dim + 1
      grid_state%vx1 = 3 + y_dim
      grid_state%vy1 = 2 * y_dim + 1
      grid_state%vy2 = 3 * y_dim + 1
      grid_state%wy1 = 3 * y_dim + 1
      grid_state%wy2 = 3 * y_dim + 1
      grid_state%sy1 = 3 * y_dim + 1
      grid_state%sy2 = 3 * y_dim + 1

      grid_state%x1u = 1 - grid_state%ux1
      grid_state%y1u = 1 - grid_state%uy1 * y_dim
      grid_state%x1v = 1 - grid_state%vx1
      grid_state%y1v = 1 - grid_state%vy1 * y_dim
      grid_state%x1w = 1 - grid_state%wx1
      grid_state%y1w = 1 - grid_state%wy1 * y_dim
      grid_state%x1s = 1 - grid_state%sx1
      grid_state%y1s = 1 - grid_state%sy1 * y_dim

      grid_state%x2u = grid_state%grid%nx_local + grid_state%ux2
      grid_state%y2u = grid_state%grid%ny_local + grid_state%uy2 * y_dim
      grid_state%x2v = grid_state%grid%nx_local + grid_state%vx2
      grid_state%y2v = grid_state%grid%ny_local + grid_state%vy2 * y_dim
      grid_state%x2w = grid_state%grid%nx_local + grid_state%wx2
      grid_state%y2w = grid_state%grid%ny_local + grid_state%wy2 * y_dim
      grid_state%x2s = grid_state%grid%nx_local + grid_state%sx2
      grid_state%y2s = grid_state%grid%ny_local + grid_state%sy2 * y_dim

      grid_state%idx = 1.0 / real(grid_state%grid%dx)
      grid_state%idy = 1.0 / real(grid_state%grid%dy)

      grid_state%nx_sb = grid_state%grid%global_i_start - 1
      grid_state%ny_sb = grid_state%grid%global_j_start - 1
      grid_state%i_global_LCM = grid_state%grid%global_i_start
      grid_state%j_global_LCM = grid_state%grid%global_j_start

      allocate(grid_state%gpxc(-1:grid_state%grid%nx_local+1))
      allocate(grid_state%gpyc(-1:grid_state%grid%ny_local+1))
      allocate(grid_state%gpzc(grid_state%grid%nz))
      allocate(grid_state%gpxs(0:grid_state%grid%nx_local+1))
      allocate(grid_state%gpys(0:grid_state%grid%ny_local+1))
      allocate(grid_state%gpzs(grid_state%grid%nz+1))
      allocate(grid_state%dV(grid_state%grid%nz))

      grid_state%gpzc = real(grid_state%grid%z_center)
      grid_state%gpzs = real(grid_state%grid%z_face)
      grid_state%dV = real(grid_state%grid%dx * grid_state%grid%dy * grid_state%grid%dz)

      do i = 0, grid_state%grid%nx_local + 1
         grid_state%gpxs(i) = real(i - 1 + grid_state%nx_sb) * real(grid_state%grid%dx)
      enddo
      do j = 0, grid_state%grid%ny_local + 1
         grid_state%gpys(j) = real(j - 1 + grid_state%ny_sb) * real(grid_state%grid%dy)
      enddo
      do i = -1, grid_state%grid%nx_local + 1
         grid_state%gpxc(i) = (real(i + grid_state%nx_sb) - 0.5) * real(grid_state%grid%dx)
      enddo
      do j = -1, grid_state%grid%ny_local + 1
         grid_state%gpyc(j) = (real(j + grid_state%ny_sb) - 0.5) * real(grid_state%grid%dy)
      enddo

      grid_state%xwidth = real(grid_state%grid%nx_global) * real(grid_state%grid%dx)
      grid_state%ywidth = real(grid_state%grid%ny_global) * real(grid_state%grid%dy)
      grid_state%ztop = grid_state%gpzs(grid_state%grid%nz + 1) - 0.001
      grid_state%zbottom = grid_state%gpzs(1) + 0.001
   end subroutine initialize_geometry

end module lcm_grid_state
