!> Contains netCDF4-output routines
!>
module lcm_output
   use iso_fortran_env, only: real32, real64
   use mpi_f08, only: MPI_Wtime, MPI_Gather, MPI_Cart_coords, MPI_Barrier,          &
                      MPI_SUCCESS, MPI_REAL, MPI_INFO_NULL
   use netcdf, only: nf90_strerror, nf90_inq_varid, nf90_var_par_access,            &
                     nf90_create, nf90_open, nf90_close, nf90_def_dim,              &
                     nf90_put_att, nf90_def_var, nf90_def_var_fill, nf90_put_var,   &
                     NF90_NOERR, NF90_NETCDF4, NF90_WRITE, NF90_COLLECTIVE,         &
                     NF90_GLOBAL, NF90_FLOAT, NF90_FILL

   use lcm_types,          only: lcm_grid_t
   use lcm_grid_state,     only: grid_state
   use lcm_field_state,    only: field_state
   use lcm_parallel_state, only: parallel_state, lcm_parallel_t

   implicit none

   private

   public :: lcm_write_global_3d_fields

contains

   !> This is the main routine for global 3d output of eulerian fields.
   !> First communicates the data from the subdomains intranode-wise to the node heads.
   !> Then from root process creates the netCDF4 file and defines the dims, vars, attrs.
   !> At last I/O writes fields data to disk (single precision) in parallel from all node heads.
   !> The result is a single file for the global fields for one time step.
   subroutine lcm_write_global_3d_fields(current_lcm_step, time_step)
      integer, intent(in)      :: current_lcm_step
      real(real64), intent(in) :: time_step

      character(99) :: filename = ''
      integer :: ncid, dimid_x_c, dimid_x_e, dimid_y_c, dimid_y_e, dimid_z_fl, dimid_z_hl
      integer :: varid_temp, varid_qv, ierr
      integer :: k, rank_i, len_node_x, len_global_x, len_global_y, len_global_z
      real(real32) :: fill_value
      real(real64) :: start_time, end_time
      real(real32), allocatable :: temperature_3d_node(:,:,:)
      real(real32), allocatable :: qv_3d_node(:,:,:)

      if (parallel_state%intranode%dimensions(2) * grid_state%grid%ny_local      &
            /= grid_state%grid%ny_global) then
         error stop 'intranode dimensions are not set up correctly'
      endif

      !> do intranode MPI communication of subdomain data

      start_time = MPI_Wtime()
      call get_intranode_3d_field(parallel_state%intranode, grid_state%grid,     &
                                  field_state%temperature, temperature_3d_node)
      call get_intranode_3d_field(parallel_state%intranode, grid_state%grid,     &
                                  field_state%qv, qv_3d_node)
      end_time = MPI_Wtime()
      if (parallel_state%intranode%rank == 0) then
         print *, 'finished gathering intranode data in', end_time - start_time, 's'
      endif

      if (parallel_state%intranode%rank /= 0) return
      !> from here on only internode processes (=node heads) execute

      !write (filename,'(a,i0,a)') 'lcm_output_3d_fields_step', current_lcm_step, '.nc'
      filename = 'lcm_output_3d_fields_step5.nc'

      print *, 'writing to netcdf4 file now from global rank', parallel_state%global%rank,   &
               'filename: ', trim(filename)

      len_node_x   = parallel_state%intranode%dimensions(1) * grid_state%grid%nx_local
      len_global_x = grid_state%grid%nx_global
      len_global_y = grid_state%grid%ny_global
      len_global_z = grid_state%grid%nz

      fill_value = -999.0

      !> root process creates the netCDF4 file and defines the dims, vars, attrs...

      if (parallel_state%internode%rank == 0) then
         call nc_wrap(nf90_create(path=trim(filename), cmode=NF90_NETCDF4, ncid=ncid))

         !> write dimensions

         call nc_wrap(nf90_def_dim(ncid=ncid, name='nx_c', len=len_global_x,   dimid=dimid_x_c))
         call nc_wrap(nf90_def_dim(ncid=ncid, name='nx_e', len=len_global_x+1, dimid=dimid_x_e))
         call nc_wrap(nf90_def_dim(ncid=ncid, name='ny_c', len=len_global_y,   dimid=dimid_y_c))
         call nc_wrap(nf90_def_dim(ncid=ncid, name='ny_e', len=len_global_y+1, dimid=dimid_y_e))
         call nc_wrap(nf90_def_dim(ncid=ncid, name='nz_fl', len=len_global_z,   dimid=dimid_z_fl))
         call nc_wrap(nf90_def_dim(ncid=ncid, name='nz_hl', len=len_global_z+1, dimid=dimid_z_hl))


         !> write global attributes

         call nc_def_global_attr(ncid, 'Conventions', 'CF-1.13')
         call nc_def_global_attr(ncid, 'title', 'Test output file created from LCM')
         call nc_def_global_attr(ncid, 'institution',                   &
           'Institut für Meteorologie / Freie Universität Berlin')
         call nc_def_global_attr(ncid, 'source', 'DALES-LCM')
         call nc_def_global_attr(ncid, 'lcm_step', content_int=current_lcm_step)
         call nc_def_global_attr(ncid, 'time_step_sec', content_real64=time_step)


         !> create variables and their attributes

         call nc_def_var(ncid=ncid, varid=varid_temp,                &
                         xtype=NF90_FLOAT,                           &
                         fill_value=fill_value,                      &
                         ndims=3,                                    &
                         dimids=[dimid_x_c, dimid_y_c, dimid_z_fl],  &
                         var_name='temperature',                     &
                         standard_name='air_temperature',            &
                         long_name='grid cell air temperature',      &
                         units='K')

         call nc_def_var(ncid=ncid, varid=varid_qv,                  &
                         xtype=NF90_FLOAT,                           &
                         fill_value=fill_value,                      &
                         ndims=3,                                    &
                         dimids=[dimid_x_c, dimid_y_c, dimid_z_fl],  &
                         var_name='qv',                              &
                         standard_name='specific_humidity',          &
                         long_name='grid cell specific humidity',    &
                         units='kg/kg')

         call nc_wrap(nf90_close(ncid))
      endif

      !> synchronize all internode processes to make sure one process does not attempt below
      !> to write before the above creation of the file by the root process
      call MPI_Barrier(parallel_state%internode%communicator, ierr)
      if (ierr /= MPI_SUCCESS) then
         error stop 'LCM failed to synchronize internode processes after netcdf file creation'
      endif


      !> write variable values to created variables in parallel with all node heads

      start_time = MPI_Wtime()
      call nc_wrap(nf90_open(path=trim(filename), mode=NF90_WRITE, ncid=ncid,    &
                             comm=parallel_state%internode%communicator%MPI_Val, &
                             info=MPI_INFO_NULL%MPI_Val))

      call nc_wrap(nf90_inq_varid(ncid, 'temperature', varid_temp))
      call nc_wrap(nf90_var_par_access(ncid, varid_temp, NF90_COLLECTIVE))

      call nc_write_var3d(ncid, varid_temp, temperature_3d_node,        &
                          parallel_state%internode%rank, len_node_x,    &
                          len_global_x, len_global_y, len_global_z)

      call MPI_Barrier(parallel_state%internode%communicator, ierr)

      call nc_wrap(nf90_inq_varid(ncid, 'qv', varid_qv))
      call nc_wrap(nf90_var_par_access(ncid, varid_qv, NF90_COLLECTIVE))

      call nc_write_var3d(ncid, varid_qv, qv_3d_node,                   &
                          parallel_state%internode%rank, len_node_x,    &
                          len_global_x, len_global_y, len_global_z)

      call nc_wrap(nf90_close(ncid))
      end_time = MPI_Wtime()
      print *, 'finished writing intranode data to disk in', end_time - start_time, 's'

   end subroutine lcm_write_global_3d_fields

   !> wrapper routine for defining global netcdf attributes
   !> when calling this routine, if not specified with a keyword,
   !> the default is string content
   subroutine nc_def_global_attr(ncid, name, content_str, content_int,  &
                                 content_real32, content_real64)
      integer, intent(in)                :: ncid
      character(*), intent(in)           :: name
      character(*), intent(in), optional :: content_str
      integer, intent(in), optional      :: content_int
      real(real32), intent(in), optional :: content_real32
      real(real64), intent(in), optional :: content_real64

      if (present(content_str)) then
         call nc_wrap(nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, &
                                   name=trim(name),              &
                                   values=content_str))
      elseif (present(content_int)) then
         call nc_wrap(nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, &
                                   name=trim(name),              &
                                   values=content_int))
      elseif (present(content_real32)) then
         call nc_wrap(nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, &
                                   name=trim(name),              &
                                   values=content_real32))
      elseif (present(content_real64)) then
         call nc_wrap(nf90_put_att(ncid=ncid, varid=NF90_GLOBAL, &
                                   name=trim(name),              &
                                   values=content_real64))
      endif
   end subroutine nc_def_global_attr

   !> wrapper routine for defining a netcdf variable and filling its attributes
   subroutine nc_def_var(ncid, varid, xtype, fill_value, ndims, dimids,  &
                         var_name, standard_name, long_name, units)
      integer, intent(in)                :: ncid
      integer, intent(out)               :: varid
      integer, intent(in)                :: xtype
      real(real32), intent(in)           :: fill_value
      integer, intent(in)                :: ndims
      integer, intent(in)                :: dimids(ndims)
      character(*), intent(in)           :: var_name
      character(*), intent(in)           :: standard_name
      character(*), intent(in)           :: long_name
      character(*), intent(in)           :: units

      call nc_wrap(nf90_def_var(ncid=ncid, varid=varid,   &
                                name=var_name,            &
                                xtype=xtype,              &
                                !> dimensions have to be flipped (Fortran <--> C order):
                                dimids=dimids(ndims:1:-1)))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=varid,   &
                                name='standard_name',     &
                                values=standard_name))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=varid,   &
                                name='long_name',         &
                                values=long_name))
      call nc_wrap(nf90_put_att(ncid=ncid, varid=varid,   &
                                name='units',             &
                                values=units))
      call nc_wrap(nf90_def_var_fill(ncid, varid, NF90_FILL, fill_value))
   end subroutine nc_def_var

   !> wrapper routine for writing 3d values to a netcdf variable
   !> optionally prints out minimal level-wise stats
   subroutine nc_write_var3d(ncid, varid, values_3d, position_node_x, len_node_x,               &
                             len_global_x, len_global_y, len_global_z, lprint_level_stats)
      integer, intent(in)           :: ncid
      integer, intent(in)           :: varid
      integer, intent(in)           :: position_node_x, len_node_x
      integer, intent(in)           :: len_global_x, len_global_y, len_global_z
      real(real32), intent(in)      :: values_3d(:,:,:)
      logical, intent(in), optional :: lprint_level_stats
      integer :: k

      call nc_wrap(nf90_put_var(ncid=ncid, varid=varid,                                         &
                                !> dimensions have to be flipped (Fortran <--> C order):
                                values=reshape(values_3d,                                       &
                                               shape=[len_global_z, len_global_y, len_node_x],  &
                                               order=[3, 2, 1]),                                &
                                start=[1, 1, position_node_x*len_node_x+1],                     &
                                count=[len_global_z, len_global_y, len_node_x]))

      if (present(lprint_level_stats)) then
         if (lprint_level_stats) then
            do k = 1, len_global_z
               print *, 'k=', k,                               &
                        'min=', minval(values_3d(:,:,k)),      &
                        'max=', maxval(values_3d(:,:,k)),      &
                        'mean=', sum(values_3d(:,:,k))/        &
                           (real(len_node_x, kind=real32)      &
                            *real(len_global_y, kind=real32))
            enddo
         endif
      endif
   end subroutine nc_write_var3d

   !> wrapper routine for handling of netcdf errors
   !> if error detected, prints out error message and stops LCM
   subroutine nc_wrap(status)
      integer, intent(in) :: status

      if(status /= NF90_NOERR) then
         error stop trim(nf90_strerror(status))
      endif
   end subroutine nc_wrap

   !> This routine communicates the data from all subdomains inside a node to the node head.
   !> Works with MPI_Gather() with a 1d buffer and a following unpacking algorithm.
   !> Optimized for low memory usage
   subroutine get_intranode_3d_field(intranode, localgrid, field3d_local, field3d_node)
      type(lcm_parallel_t), intent(in) :: intranode
      type(lcm_grid_t), intent(in)     :: localgrid
      real(real64), intent(in)               :: field3d_local(:,:,:)
      real(real32), allocatable, intent(out) :: field3d_node(:,:,:)
      real(real32), allocatable              :: node_buffer1d(:)
      real(real32), allocatable              :: dummy_recv(:)
      integer      :: n_local, ierr, rank_i, start_x, end_x, start_y, end_y, offset_flat
      integer      :: coords_i(2)

      n_local = size(field3d_local)

      if (intranode%rank == 0) then
         allocate(node_buffer1d(n_local*intranode%size))
         allocate(field3d_node(intranode%dimensions(1)*localgrid%nx_local,       &
                               intranode%dimensions(2)*localgrid%ny_local,       &
                               localgrid%nz))

         call MPI_Gather(real(field3d_local, kind=real32), n_local, MPI_REAL,    &
                         node_buffer1d, n_local, MPI_REAL,                       &
                         0, intranode%communicator, ierr)
         if (ierr /= MPI_SUCCESS) then
            error stop 'LCM failed to gather the intranode field3d on intranode rank 0'
         endif

         !> unpack the 1d buffer data into the large intranode 3d array
         do rank_i = 0, intranode%size - 1
            call MPI_Cart_coords(intranode%communicator, rank_i, 2, coords_i, ierr)
            if (ierr /= MPI_SUCCESS) then
               error stop 'LCM failed to get the intranode cartesian ' //        &
                          'coordinates on intranode rank 0'
            endif

            start_x = coords_i(1) * localgrid%nx_local + 1
            end_x = (coords_i(1) + 1) * localgrid%nx_local

            start_y = coords_i(2) * localgrid%ny_local + 1
            end_y = (coords_i(2) + 1) * localgrid%ny_local
            offset_flat = rank_i * n_local

            field3d_node(start_x:end_x, start_y:end_y, :) = reshape(             &
               node_buffer1d(offset_flat + 1 : offset_flat + n_local),           &
               shape=[localgrid%nx_local, localgrid%ny_local, localgrid%nz],     &
               order=[1, 2, 3] )
         enddo

      else
         allocate(dummy_recv(n_local))
         call MPI_Gather(real(field3d_local, kind=real32), n_local, MPI_REAL,    &
                         dummy_recv, n_local, MPI_REAL,                          &
                         0, intranode%communicator, ierr)
         if (ierr /= MPI_SUCCESS) then
            error stop 'LCM failed to gather the intranode field3d on intranode rank /= 0'
         endif
      endif
   end subroutine get_intranode_3d_field

end module lcm_output
