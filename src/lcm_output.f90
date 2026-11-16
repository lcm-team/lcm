!> Contains netcdf-output routines
!!
module lcm_output
   use iso_fortran_env, only: real32, real64
   use mpi_f08, only: MPI_Comm_rank, MPI_Comm_size, MPI_Barrier, MPI_Wtime,   &
                      MPI_Bcast, MPI_Scatter, MPI_INT, MPI_REAL, MPI_LOGICAL

   use netcdf, only: nf90_noerr, nf90_strerror, nf90_inq_libvers,             &
                     NF90_NETCDF4, NF90_GLOBAL, NF90_FLOAT, NF90_FILL,        &
                     nf90_create, nf90_close, nf90_def_dim, nf90_put_att,     &
                     nf90_def_var, nf90_def_var_fill, nf90_enddef, nf90_put_var

   use lcm_parallel_state, only: parallel_state
   use lcm_grid_state,     only: grid_state
   use lcm_field_state,    only: field_state

   implicit none

   private

   public :: lcm_write_global_3d_fields

contains

   subroutine lcm_write_global_3d_fields(filename, current_lcm_step)
      character(*), intent(in) :: filename
      integer, intent(in)      :: current_lcm_step
      !real(real64), intent(in) :: field_3d(:,:,:)
      integer :: ncid, dimid_x, dimid_y, dimid_z, varid_temp, varid_qv
      integer :: k, len_x, len_y, len_z
      real(real32) :: fill_value

      if (parallel_state%parallel%rank /= 0) return

      len_x = grid_state%grid%nx_local
      len_y = grid_state%grid%ny_local
      len_z = grid_state%grid%nz
      fill_value = -999.0

      print *, 'writing to netcdf4 file now from rank 0...'
      print *, 'filename: ', filename

      ! create file and dimensions

      call nc_wrap(nf90_create(path=TRIM(filename), cmode=NF90_NETCDF4, ncid=ncid))

      call nc_wrap(nf90_def_dim(ncid=ncid, name='nx_c', len=len_x, dimid=dimid_x))
      call nc_wrap(nf90_def_dim(ncid=ncid, name='ny_c', len=len_y, dimid=dimid_y))
      call nc_wrap(nf90_def_dim(ncid=ncid, name='nz_c', len=len_z, dimid=dimid_z))


      ! write global attributes

      call nc_wrap(nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, &
                                name='Conventions',           &
                                values='CF-1.13'))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, &
                                name='title',                 &
                                values='Test output file created from LCM'))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, &
                                name='institution',           &
                                values='Institut für Meteorologie / Freie Universität Berlin'))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, &
                                name='source',                &
                                values='DALES-LCM'))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, &
                                name='lcm_step',              &
                                values=current_lcm_step))


      ! create variables and their attributes

      call nc_wrap(nf90_def_var(ncid=ncid, varid=varid_temp,   &
                                name='temperature',            &
                                xtype=NF90_FLOAT,              &
                                ! dimension order has to be given flipped around (Fortran order):
                                dimids=[dimid_z, dimid_y, dimid_x]))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=varid_temp,   &
                                name='standard_name',          &
                                values='air_temperature'))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=varid_temp,   &
                                name='long_name',              &
                                values='grid cell air temperature'))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=varid_temp,   &
                                name='units',                  &
                                values='K'))
      call nc_wrap(nf90_def_var_fill(ncid, varid_temp, NF90_FILL, fill_value))

      call nc_wrap(nf90_def_var(ncid=ncid, varid=varid_qv,     &
                                name='qv',                     &
                                xtype=NF90_FLOAT,              &
                                ! dimension order has to be given flipped around (Fortran order):
                                dimids=[dimid_z, dimid_y, dimid_x]))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=varid_qv,     &
                                name='standard_name',          &
                                values='specific_humidity'))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=varid_qv,     &
                                name='long_name',              &
                                values='grid cell specific humidity'))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=varid_qv,     &
                                name='units',                  &
                                values='kg/kg'))
      call nc_wrap(nf90_def_var_fill(ncid, varid_qv, NF90_FILL, fill_value))


      call nc_wrap(nf90_enddef(ncid=ncid))

      ! write variable values to created variables

      ! write 3d field level by level in order to have transposable 2d arrays
      do k = 1, len_z
         print *, 'k=', k, &
                  'min=', minval(field_state%temperature(:,:,k)),                      &
                  'max=', maxval(field_state%temperature(:,:,k)),                      &
                  'mean=', sum(field_state%temperature(:,:,k))/                        &
                            (real(len_x, kind=real64)*real(len_y, kind=real64))

         call nc_wrap(nf90_put_var(ncid=ncid, varid=varid_temp,                        &
                                   ! make sure to give flipped around 2d array (Fortran order):
                                   values=transpose(field_state%temperature(:,:,k)),   &
                                   start=[k, 1, 1],                                    &
                                   count=[1, len_y, len_x]))
      enddo
      do k = 1, len_z
         print *, 'k=', k, &
                  'min=', minval(field_state%qv(:,:,k)),                      &
                  'max=', maxval(field_state%qv(:,:,k)),                      &
                  'mean=', sum(field_state%qv(:,:,k))/                        &
                            (real(len_x, kind=real64)*real(len_y, kind=real64))

         call nc_wrap(nf90_put_var(ncid=ncid, varid=varid_qv,                 &
                                   ! make sure to give flipped around 2d array (Fortran order):
                                   values=transpose(field_state%qv(:,:,k)),   &
                                   start=[k, 1, 1],                           &
                                   count=[1, len_y, len_x]))
      enddo


      call nc_wrap(nf90_close(ncid=ncid))
      print *, 'finished writing'

   end subroutine lcm_write_global_3d_fields


   subroutine nc_handle_error(status)
      integer, intent(in) :: status

      if(status /= nf90_noerr) then
         error stop trim(nf90_strerror(status))
      endif
   end subroutine nc_handle_error

   ! if error return code, prints netcdf error message and stops
   subroutine nc_wrap(status)
      integer, intent(in) :: status

      if(status /= nf90_noerr) then
         error stop trim(nf90_strerror(status))
      endif
   end subroutine nc_wrap

end module lcm_output
