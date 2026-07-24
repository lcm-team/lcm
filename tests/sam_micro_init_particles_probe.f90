program sam_micro_init_particles_probe
   use iso_fortran_env, only: int64
   use grid, only: dimx1_s, dimx2_s, dimy1_s, dimy2_s, nzm
   use micro_vars
   use microphysics, only: imv, micro_field, nmicro_fields, nmicro_fields_diag

   implicit none

   integer :: n
   integer :: output_unit
   character(len=1024) :: output_path
   character(len=32) :: aerosol_mode

   interface
      subroutine micro_init_particles(i, j, k)
         integer, intent(in) :: i
         integer, intent(in) :: j
         integer, intent(in) :: k
      end subroutine micro_init_particles
   end interface

   call get_command_argument(1, output_path)
   if (len_trim(output_path) == 0) then
      error stop "SAM-LCM probe output path is required"
   endif
   call get_command_argument(2, aerosol_mode)
   if (len_trim(aerosol_mode) == 0) aerosol_mode = "bin"

   call setup_sam_lcm_state(trim(aerosol_mode))
   call micro_init_particles(1, 1, 1)

   open(newunit=output_unit, file=trim(output_path), status="replace", &
        action="write")
   do n = 1, np_xy
      write(output_unit,'(i0,1x,i0,1x,es24.16,1x,4(es24.16,1x),i0)') &
         n, particle(n)%ivar(iA), particle(n)%var(iN), &
         particle(n)%var(iAero_portions(1)), &
         particle(n)%var(iAero_portions(2)), &
         particle(n)%var(iAero_portions(3)), &
         particle(n)%var(iAero_portions(4)), &
         particle(n)%ivar(iTYPE)
   enddo
   close(output_unit)

contains

   subroutine setup_sam_lcm_state(aerosol_mode)
      character(len=*), intent(in) :: aerosol_mode

      np_xy = 8
      nzl = 1
      nx_LCM = 1
      ny_LCM = 1
      nz_LCM = 1
      dx_LCM = 50.0
      dy_LCM = 50.0
      init_aero = aerosol_mode
      unit_aero = "1/kg"
      switch_ice = .false.
      frac_ice = real(0.0, kind=4)
      n_aero = real([100.0e6, 0.0, 0.0, 0.0], kind=4)
      rm_aero = real([0.04e-6, 0.0, 0.0, 0.0], kind=4)
      sigma_aero = real([1.4, 0.0, 0.0, 0.0], kind=4)
      shares_aero = real([1.0, 0.0, 0.0, 0.0], kind=4)
      rho_aero = real([2170.0, 0.0, 0.0, 0.0], kind=4)
      n_init = np_xy
      n_particles = np_xy
      n_particles_alloc = np_xy

      allocate(dV(1))
      allocate(rho_air(1))
      allocate(pres_LCM(1))
      allocate(dz_LCM(1))
      allocate(tabs_LCM(1,1,1))
      allocate(n_start(1,1,1))
      allocate(n_end(1,1,1))
      allocate(particle(1:np_xy))
      allocate(inv_cdf(4, 0:n_inv_cdf))
      allocate(micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nz_LCM, &
                           nmicro_fields + nmicro_fields_diag))

      dV = dx_LCM * dy_LCM * 100.0
      rho_air = 1.2
      pres_LCM = 100000.0
      dz_LCM = 100.0
      tabs_LCM = 280.0
      micro_field = 0.0
      micro_field(1,1,1,imv) = 0.010
      n_start = 1
      n_end = np_xy
      do n = 1, np_xy
         particle(n)%ivar = 0_int64
         particle(n)%var = 0.0_lcmp
      enddo
   end subroutine setup_sam_lcm_state

end program sam_micro_init_particles_probe
