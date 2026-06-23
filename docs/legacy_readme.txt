================================================================================
-------------------- Lagrangian Cloud Microphysics for SAM ---------------------
================================================================================

In this readme, all necessary steps to run Lagrangian microphysics in SAM are
outlined. For further information, do not hesitate to contact Fabian Hoffmann
(mail.fabian.hoffmann@gmail.com).

Please refer to the following publications for more information on this model:

 * Hoffmann et al. (2015), https://doi.org/10.1016/j.atmosres.2014.12.008
 * Hoffmann (2017), https://doi.org/10.5194/acp-17-8343-2017
 * Hoffmann et al. (2019), https://doi.org/10.1175/JAS-D-18-0087.1
 * Hoffmann and Feingold (2019), https://doi.org/10.1175/JAS-D-18-0318.1

Please note that some parts of this code originate from the Lagrangian particle
model previously implemented in PALM (https://palm.muk.uni-hannover.de). Further
note that the underlying particle model has been been graciously provided by Tak
Yamaguchi (CIRES/NOAA).

The Lagrangian cloud microphysics operate largely autonomous. Only the following
changes are required:

 * SGS_TKE/tke_full.f90: Provide three-dimensional fields for the kinetic energy
                         dissipation rate (diss_LCM), subgrid-scale turbulence
                         kinetic energy (tke_LCM), subgrid-scale turbulent
                         diffusion coefficient (tk_LCM), and the subgrid-scale
                         mixing length (smix_LCM).

 * write_fields3D.f90:   The loop for the output of microphysical fields needs
                         to run from 1 to nmicro_fields + 17 (and not from 1 to
                         nmicro_fields) to output additional quantities
                         calculated in the model.

Further output quantities are only available in combination with Tak Yamaguchi's
version of SAM (preprocessing option SAMTAK).

CHANGES: symbolic link in RUNDATA removed: lst->lst_SAMTAK
CHANGES: symbolic link in RUNDATA/RADAEROSIMPLE removed. aer_opts_lw->RADAEROSIMPLE_rh28bin_nca80bin/aer_opts_lw
CHANGES: symbolic link in RUNDATA/RADAEROSIMPLE removed. aer_opts_sw->RADAEROSIMPLE_rh28bin_nca80bin/aer_opts_sw 

================================================================================
-------------------- Namelist Parameters and Descriptions ----------------------
================================================================================

nstep_LCM_init  = 1           ! start LCM at a prescribed time step
n_uvw_iteration = 1           ! number of iteration to find velocity at t=t+1/2.

--------------------------------------------------------------------------------

np_xy           = 20          ! particles per grid box
nzl             = 100         ! highest grid level for placing particles

n_aero          = 30.0E6      ! prescribe number concentration of cloud droplets/aerosols (1/m3)
rm_aero         = 0.08E-6     ! geometric mean radius of aerosol size distribution (m)
sigma_aero      = 1.4         ! geometric stddev of aerosol size distribution (unitless)

init_aero       = 'random'    ! mode to initialize aerosols (bin, USM, or random)

--------------------------------------------------------------------------------

switch_sgs_velocities = .FALSE.     ! SGS velocity fluctuations
switch_LEM_micro      = .TRUE.      ! SGS mixing model

--------------------------------------------------------------------------------

switch_steady_pichamber           = .FALSE. ! turn on re-injection of fallen-out droplets as aerosols, creates steady state pi chamber
switch_walls_pichamber            = .FALSE. ! turn on lateral walls of pi chamber
switch_reflective_walls_pichamber = .FALSE. ! turn on reflective walls of pi chamber (if false, droplets are absorbed)

--------------------------------------------------------------------------------

switch_sedimentation = .FALSE.       ! sedimentation
start_sedimentation  = 1             ! number of timesteps to initiate sedimentation after model starts

switch_collection = .FALSE.          ! collection growth
start_collection  = 1                ! number of timesteps to initiate collection after model starts
switch_linear_sampling_collection = .TRUE.  ! LSM method to speed up collision process (Shima et al. 2009)

switch_analytical_micro = .FALSE.    ! default is .FALSE., .TRUE. enables analytical
                                     ! solution of diffusional growth equation
                                     ! (no radiation, solute, curvature effects)
switch_radiation_micro  = .FALSE.    ! radiation effects of diffusional growth
switch_FFC_micro        = .FALSE.    ! consider film-forming compounds in diffusional growth equation

r_separate_micro        = 25.0E-6    ! separation radius between cloud and rain droplets

--------------------------------------------------------------------------------

switch_particle_source = .FALSE.  ! add a particle source
switch_phase_sources = .FALSE.    ! If true, keep track of separate CCN and IN sources
n_start_source_aero  = 10              ! prescribe timestep when aerosol particle source becomes active
switch_moving_sources = .FALSE.        ! if true, allow particle sources to move in time.
                                       ! can be useful for galilean domains

nz_b_source_aero     = 30              ! bottom boundary of aerosol source
nz_t_source_aero     = 31              ! top boundary of aerosol source
ny_s_source_aero     = 1               ! south boundary of aerosol source
ny_n_source_aero     = 1               ! north boundary of aerosol source
nx_w_source_aero     = 90              ! west boundary of aerosol source
nx_e_source_aero     = 90              ! east boundary of aerosol source

vx_source_aero       = 0               ! eastward velocity of aerosol source
vy_source_aero       = 0               ! northward velocity of aerosol source
vz_source_aero       = 0               ! upward velocity of aerosol source

n_source_aero        = 1.0E6           ! aerosol source strength (1/m^3/s)
p_source_aero        = 0               ! Probability of emitting an aerosol particle per gridbox per timestep
r_source_aero   = 0.1E-6          ! geometric mean of aerosol source (m)
sigma_source_aero = 2.0           ! geometric standard deviation of aerosol source

frac_in_source_particles          ! fraction of emitted superparticles which should contain IN
frac_in_source_weight             ! multiply IN-containing particles by a multiplication factor
                                  ! to account for IN being much rarer than CCN (so should be <<1)

type_source_in = 3                ! which species of IN should be added by the source volume?
                                  ! Default: Snowmax. Otherwise, specify this even when using switch_phase_sources = .FALSE.

! IN source flags only have an effect if switch_phase_sources is .TRUE. If not, 
n_start_source_in  = 10              ! prescribe timestep when IN particle source becomes active

nz_b_source_in     = 30              ! bottom boundary of IN source
nz_t_source_in     = 31              ! top boundary of IN source
ny_s_source_in     = 1               ! south boundary of IN source
ny_n_source_in     = 1               ! north boundary of IN source
nx_w_source_in     = 90              ! west boundary of IN source
nx_e_source_in     = 90              ! east boundary of IN source

vx_source_in       = 0               ! eastward velocity of in source
vy_source_in       = 0               ! northward velocity of in source
vz_source_in       = 0               ! upward velocity of in source

n_source_in        = 1.0E6           ! IN source strength (1/m^3/s)
p_source_in        = 0               ! Probability of emitting an IN particle per gridbox per timestep
r_source_in   = 0.1E-6          ! geometric mean of IN source (m)
sigma_source_in = 2.0           ! geometric standard deviation of IN source

--------------------------------------------------------------------------------

spec_nsave      = 30          ! interval of spectra output
spec_nsavestart = 1           ! start output of spectra after spec_nsavestart timesteps

--------------------------------------------------------------------------------

switch_ice       = .FALSE.    ! ice switch
switch_inmass    = .FALSE.    ! whether to use behavior which varies heterogeneous nucleation rate based on IN type or mass
                              ! FALSE means that all ice particles use a SNOWMAX parameterization which nucleates at high temperatures
switch_habit_ice = .TRUE.     ! predict i habit switch
habitparam       = "Chen"     ! Choose which ice habit paramterization to use. Currently available options are "Chen" and "DISKICE"
                              ! "Chen" is a parameterization from Chen & Lamb (1994) assuming diffusion-limited growth
                              ! and habit development based only on temperature
                              ! "DISKICE" accounts for kinetic and surface effects and prescribes habit based on
                              ! temperature, vapor pressure, and background atmospheric pressure
sizethresh_chen  = 10.0E-6    ! If using habitparam = "Chen", the size threshold which ice particles must reach to develop habits
sizethresh_diskice = 10.0E-6  ! Same as above, but for using habitparam = "DISKICE"
switch_diskice_sizethresh = .TRUE. 
                              ! DISKICE is more numerically stable than Chen & Lamb, 
                              ! so a size threshold for habit development is optional
                              ! If .FALSE., sizethresh_diskice will be disregarded and ice crytal can immediately develop habits
n_ice            = 0.001E6    ! number concentration of ice nuclei (1/m3)
frac_ice         = 0.0        ! fraction of INITIAL (potential) ice particles of all Lagrangian particles per grid box (-)
frac_ice_source  = 0.0        ! fraction of (potential) ice particles in source (-)
start_ice        = 0          ! start ice physics at timestep

================================================================================
-------------------- Building and Running SAM on a Laptop ----------------------
================================================================================

use homebrew to install gcc, openmpi, netcdf, ...

REMEMBER: set NOTIMERS=-DDISABLE_TIMERS in the Makefile since the c routines do
          not work anymore

./Clean     <- delete compiled files
./RM_OUTPUT <- delete output
./Build     <- build SAM

cd UTIL
make -f Makefile <- create utilities

mpirun -np 4 --oversubscribe ./SAM  <- run SAM on 4 virtual cores
mpirun -np 2 ./SAM                  <- run SAM on two cores 2
