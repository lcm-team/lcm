!> Regression values from SAM micro_sgs_mixing with its first-temperature-point
!! diffusion error corrected. The reference was evaluated in real64 with SAM
!! constants, LCM's saturation function and the same seeded uniform stream.
!! No SAM modules are needed to run these tests.
program test_lem_sgs_sam_reference
   use iso_fortran_env, only: real64, int64
   use minimal_les_host, only: setup_minimal_les_host, minimal_host_qv
   use lcm_sgs_mixing, only: lcm_attach_lem_environment, lcm_apply_sgs_mixing
   use lcm_particle_state, only: lcm_reset_particle_cells, lcm_add_particle, cell_particles
   use lcm_parameters, only: fac_cond
   use lcm_random, only: lcm_random_uniform
   use lcm_types, only: lcm_particle_t
   implicit none

   integer, parameter :: np = 24
   integer, parameter :: retained_indices(np,3) = reshape([ &
      24,6,18,12,23,5,17,11,22,4,16,10,21,3,15,9,20,2,14,8,19,1,13,7, &
      18,6,24,12,17,5,23,11,16,4,22,10,15,3,21,9,14,2,20,8,13,1,19,7, &
      18,6,24,12,17,5,23,11,16,4,22,10,15,3,21,9,14,2,20,8,13,1,19,7], [np,3])
   integer, parameter :: upper_indices(np) = [ &
      12,24,11,23,10,22,9,21,8,20,7,19,6,18,5,17,4,16,3,15,2,14,1,13]
   real(real64), parameter :: retained_moments(3,3) = reshape([ &
      9.0300592313992733e4_real64, 2.7012536879394697_real64, 1.7383338980031934e-3_real64, &
      9.0301195119872937e4_real64, 2.7009863735876953_real64, 5.5049478232465218e-3_real64, &
      9.0301805273770136e4_real64, 2.7010183181148308_real64, 9.1523961687710667e-3_real64], [3,3])
   real(real64), parameter :: upper_moments(3) = [ &
      9.0293807654826218e4_real64, 2.6987361206290683_real64, 8.4596783171007645e-4_real64]
   ! Allow roundoff from SAM's decimal exponents versus exact rational exponents.
   real(real64), parameter :: tolerance(3) = [1.e-6_real64, 3.e-10_real64, 1.e-10_real64]
   real(real64), target :: diffusivity(2,2,2), length(2,2,2), total_water(2,2,2)
   real(real64), target :: liquid_static(2,2,2), heat_source(2,2,2), vapor_source(2,2,2)
   real(real64) :: carrier, cell_air_mass, liquid_mean, initial_mass(np), moments(3), discarded
   integer :: n, k, step, steps, origin, expected_indices(np)
   logical :: upper_layer
   character(len=32) :: mode
   type(lcm_particle_t) :: p

   call get_command_argument(1, mode)
   upper_layer = trim(mode) == 'first_layer_spacing'
   if (.not. upper_layer .and. trim(mode) /= 'retained_shuffle') error stop 'Unknown SGS reference case'
   call setup_minimal_les_host()
   call lcm_reset_particle_cells()
   k = 1
   steps = 3
   cell_air_mass = 1.2_real64 * 250000.0_real64
   if (upper_layer) then
      k = 2
      steps = 1
      ! The fixture has dz=[100,120] m. Carrier mass uses the actual second cell;
      ! SAM's SGS line length and sedimentation relaxation still use 100 m.
      cell_air_mass = 1.05_real64 * 300000.0_real64
   endif
   carrier = cell_air_mass / np
   do n = 1, np
      p = lcm_particle_t()
      p%id = int(n,int64)
      p%weight = 1_int64
      p%mass_aero = 1.e-9_real64
      p%aerosol_portions(1) = 1.0_real64
      p%mass_water = .001_real64 * carrier * (1.0_real64 + .3_real64 * sin(real(n,real64)))
      initial_mass(n) = p%mass_water
      p%temperature = 301.0_real64
      if (upper_layer) p%temperature = p%temperature + .5_real64 * sin(real(n,real64))
      p%qv = .010_real64 - p%mass_water / carrier + .0003_real64 * cos(real(n,real64))
      p%eta = .0001_real64 * sin(real(n,real64))
      p%w_sedi = .02_real64
      p%cell_particle_index = np - n + 1
      origin = 1
      if (mod(n,2) == 0) origin = 2
      p%last_i = origin
      p%last_j = 1
      p%last_k = k
      call lcm_add_particle(1,1,k,p)
   enddo
   liquid_mean = sum(initial_mass) / cell_air_mass
   diffusivity = 20.0_real64
   length = 100.0_real64
   total_water = minimal_host_qv + liquid_mean
   liquid_static = 300.0_real64 - fac_cond * liquid_mean
   heat_source = .01_real64
   vapor_source = 1.e-7_real64
   heat_source(:,:,1) = heat_source(:,:,1) + .02_real64 * 1.3_real64 / (1.2_real64 * 100.0_real64)
   vapor_source(:,:,1) = vapor_source(:,:,1) + 2.e-7_real64 * 1.3_real64 / (1.2_real64 * 100.0_real64)
   call lcm_attach_lem_environment(diffusivity, length, total_water, liquid_static, heat_source, vapor_source)
   discarded = real(lcm_random_uniform(-123),real64)

   do step = 1, steps
      if (.not. upper_layer) then
         ! Model repeated entrainment from four origin cells. This distinguishes
         ! SAM's retained permutation from a permutation reset at every call.
         do n = 1, np
            cell_particles(1,1,k)%particles(n)%last_i = 1 + mod(n,2)
            cell_particles(1,1,k)%particles(n)%last_j = 1 + mod(n/2,2)
         enddo
      endif
      call lcm_apply_sgs_mixing(.5_real64)
      expected_indices = retained_indices(:,step)
      if (upper_layer) expected_indices = upper_indices
      moments = 0.0_real64
      do n = 1, np
         p = cell_particles(1,1,k)%particles(n)
         if (p%cell_particle_index /= expected_indices(n)) error stop 'SAM SGS parcel ordering mismatch'
         if (p%mass_water /= initial_mass(n)) error stop 'SGS changed particle liquid mass'
         moments = moments + real(n,real64) * [p%temperature, p%qv, p%eta]
      enddo
      if (upper_layer) then
         if (any(abs(moments - upper_moments) > tolerance)) error stop 'SAM first-layer SGS spacing mismatch'
      else
         if (any(abs(moments - retained_moments(:,step)) > tolerance)) error stop 'SAM repeated SGS mixing mismatch'
      endif
   enddo
   print *, 'PASS: SAM SGS reference ', trim(mode)
end program test_lem_sgs_sam_reference
