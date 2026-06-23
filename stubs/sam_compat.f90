module domain
  implicit none
  integer, parameter :: YES3D = 1
  integer, parameter :: nx_gl = 4
  integer, parameter :: ny_gl = 4
  integer, parameter :: nsubdomains_x = 1
  integer, parameter :: nsubdomains_y = 1
  integer, parameter :: naqspc = 1
end module domain

module grid
  implicit none
  integer, parameter :: nx = 4
  integer, parameter :: ny = 4
  integer, parameter :: nz = 5
  integer, parameter :: nzm = 4
  integer, parameter :: nxp1 = nx + 1
  integer, parameter :: nyp1 = ny + 1
  integer, parameter :: dimx1_s = 1
  integer, parameter :: dimx2_s = nx
  integer, parameter :: dimy1_s = 1
  integer, parameter :: dimy2_s = ny
  integer :: rank = 0
  integer :: rankee = 0
  integer :: rankne = 0
  integer :: ranknn = 0
  integer :: ranknw = 0
  integer :: rankse = 0
  integer :: rankss = 0
  integer :: ranksw = 0
  integer :: rankww = 0
  integer, parameter :: nsubdomains = 1
  integer, parameter :: nsubdomains_x = 1
  integer, parameter :: nsubdomains_y = 1
  integer :: nstep = 0
  integer :: nstat = 1
  integer :: nrestart = 0
  integer :: nrestart_skip = 0
  integer :: nstop = 1
  integer :: nelapse = 1
  integer :: icycle = 1
  integer :: ncycle = 1
  integer :: nprint = 1
  integer :: nsave3D_mod = 1
  integer :: nsave3Dstart_mod = 1
  integer :: nsave3Dend_mod = 0
  integer :: nsave2D_mod = 1
  integer :: nsave2Dstart_mod = 1
  integer :: nsave2Dend_mod = 0
  integer :: nsave2D = 1
  integer :: nsave2Dstart = 1
  integer :: nsave2Dend = 0
  integer :: nsave3D = 1
  integer :: nsave3Dstart = 1
  integer :: nsave3Dend = 0
  integer :: nstatis = 1
  integer, parameter :: NADVS = 0
  logical :: restart_sep = .false.
  real :: dx = 1.0
  real :: dy = 1.0
  real :: dz = 1.0
  real :: dt = 1.0
  real :: dtn = 1.0
  real :: dtfactor = 1.0
  real :: time = 0.0
  real :: z(nzm) = [0.5, 1.5, 2.5, 3.5]
  real :: zi(nz) = [0.0, 1.0, 2.0, 3.0, 4.0]
  real :: adz(nzm) = 1.0
  logical :: RUN3D = .true.
  logical :: masterproc = .true.
  logical :: dompi = .false.
  logical :: compute_reffc = .true.
  logical :: compute_reffi = .true.
  character(len=128) :: case = "standalone"
  character(len=128) :: caseid = "standalone"
  character(len=128) :: case_restart = "standalone"
  character(len=128) :: caseid_restart = "standalone"
end module grid

module params
  implicit none
  real, parameter :: pi = 3.14159265358979323846
  real, parameter :: cp = 1004.0
  real, parameter :: rv = 461.5
  real, parameter :: rgas = 287.0
  real, parameter :: rstar = 8.314462618
  real, parameter :: ggr = 9.81
  real, parameter :: lcond = 2.5e6
  real, parameter :: lsub = 2.834e6
  real, parameter :: lfus = 3.34e5
  real, parameter :: fac_cond = lcond / cp
  real, parameter :: fac_sub = lsub / cp
  real, parameter :: fac_fus = lfus / cp
  real, parameter :: epsv = 0.622
  real, parameter :: muelq = 1.8e-5
  real, parameter :: therco = 2.4e-2
  real, parameter :: diffelq = 2.4e-5
  logical :: doprecip = .true.
  logical :: docloud = .true.
  logical :: dochemistry = .false.
  logical :: dorestart_last_only = .false.
end module params

module vars
  use grid, only: nx, ny, nzm
  implicit none
  real :: u(nx, ny, nzm) = 0.0
  real :: v(nx, ny, nzm) = 0.0
  real :: w(nx, ny, nzm) = 0.0
  real :: t(nx, ny, nzm) = 300.0
  real :: tabs(nx, ny, nzm) = 300.0
  real :: ttend(nx, ny, nzm) = 0.0
  real :: qtend(nx, ny, nzm) = 0.0
  real :: dtn = 1.0
  real :: rho(nzm) = 1.0
  real :: pres(nzm) = 100000.0
  real :: tabs0(nzm) = 300.0
  real :: q0(nzm) = 0.0
  real :: gamaz(nzm) = 0.0
  real :: total_water_prec = 0.0
  real :: precsfc = 0.0
  real :: prec_xy(nx, ny) = 0.0
  real :: fluxbq(nx, ny) = 0.0
  real :: fluxtq(nx, ny) = 0.0
  real :: fluxlq(ny, nzm) = 0.0
  real :: fluxrq(ny, nzm) = 0.0
  real :: fluxqq(nx, nzm) = 0.0
  real :: fluxhq(nx, nzm) = 0.0
  real :: qv(nx, ny, nzm) = 0.0
  real :: qcl(nx, ny, nzm) = 0.0
  real :: qci(nx, ny, nzm) = 0.0
  real :: qpl(nx, ny, nzm) = 0.0
  real :: qpi(nx, ny, nzm) = 0.0
contains
  real function esatw(tkelv)
    real, intent(in) :: tkelv
    esatw = 610.94 * exp(17.625 * (tkelv - 273.15) / max(tkelv - 30.11, 1.0))
  end function esatw

  real function qsatw(tkelv, ppa)
    use params, only: epsv
    real, intent(in) :: tkelv, ppa
    real :: e
    e = esatw(tkelv)
    qsatw = epsv * e / max(ppa - e, 1.0)
  end function qsatw

  real function dtqsatw(tkelv, ppa)
    real, intent(in) :: tkelv, ppa
    dtqsatw = (qsatw(tkelv + 0.01, ppa) - qsatw(tkelv - 0.01, ppa)) / 0.02
  end function dtqsatw

  real function esati(tkelv)
    real, intent(in) :: tkelv
    esati = 610.94 * exp(22.587 * (tkelv - 273.15) / max(tkelv + 0.7, 1.0))
  end function esati

  real function qsati(tkelv, ppa)
    use params, only: epsv
    real, intent(in) :: tkelv, ppa
    real :: e
    e = esati(tkelv)
    qsati = epsv * e / max(ppa - e, 1.0)
  end function qsati

  real function dtqsati(tkelv, ppa)
    real, intent(in) :: tkelv, ppa
    dtqsati = (qsati(tkelv + 0.01, ppa) - qsati(tkelv - 0.01, ppa)) / 0.02
  end function dtqsati
end module vars

module rad
  use grid, only: nx, ny, nzm
  implicit none
  real :: qrad(nx, ny, nzm) = 0.0
end module rad

module hbuffer
  implicit none
contains
  subroutine hbuf_put(name, values, scale)
    character(len=*), intent(in) :: name
    real, intent(in) :: values(:)
    real, intent(in) :: scale
  end subroutine hbuf_put
end module hbuffer

module hbuf_1Dtimeseries
  implicit none
  integer :: ncondavg_1Dtimeseries = 1
  real :: condavg_factor_1Dtimeseries = 1.0
end module hbuf_1Dtimeseries

module chemistry
  implicit none
  integer :: flag_chem_init = 0
  integer :: nstep_chem_init = 0
  integer :: iHp = 1
  integer :: iOHm = 1
  real :: ionic_strength_threshold = 0.0
  real :: aqspc_molmasses(1) = 1.0
  character(len=16) :: aqspc_names(1) = ["stub"]
contains
  subroutine chem_supply_caretaking()
  end subroutine chem_supply_caretaking

  subroutine chem_set_initial_concentrations()
  end subroutine chem_set_initial_concentrations

  subroutine chem_monster_proc()
  end subroutine chem_monster_proc

  subroutine chem_proc()
  end subroutine chem_proc

  real function calc_ionic_strength()
    calc_ionic_strength = 0.0
  end function calc_ionic_strength
end module chemistry

subroutine t_startf(name)
  character(len=*), intent(in) :: name
end subroutine t_startf

subroutine t_stopf(name)
  character(len=*), intent(in) :: name
end subroutine t_stopf

subroutine task_rank_to_index(rank, i_index, j_index)
  integer, intent(in) :: rank
  integer, intent(out) :: i_index, j_index
  i_index = 0
  j_index = 0
end subroutine task_rank_to_index

subroutine task_exchange(field, x1, x2, y1, y2, nz, wx1, wx2, wy1, wy2, tag)
  real, intent(inout) :: field(:,:,:)
  integer, intent(in) :: x1, x2, y1, y2, nz, wx1, wx2, wy1, wy2, tag
end subroutine task_exchange
