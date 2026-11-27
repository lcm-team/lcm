!> SAM-derived particle-air SGS mixing, independent of host-model modules.
module lcm_sgs_mixing
 use iso_fortran_env,only:real64
 use ieee_arithmetic,only:ieee_is_finite
 use lcm_grid_state,only:grid_state
 use lcm_field_state,only:field_state
 use lcm_parallel_state,only:parallel_state
 use lcm_particle_state,only:cell_particles,qsatw_liquid
 use lcm_microphysics_state,only:lcm_get_cell_mean_history
 use lcm_parameters,only:fac_cond,molecular_viscosity_air,gravity,specific_heat_dry_air
 use lcm_random,only:lcm_random_uniform
 implicit none
 private
 public::lcm_set_lem_environment,lcm_validate_lem_environment,lcm_apply_sgs_mixing
 public::lem_diffuse_periodic,lem_triplet_map
 real(real64),allocatable:: diffusivity(:,:,:),mixing_length(:,:,:)
 real(real64),allocatable:: total_water(:,:,:),liquid_static_energy(:,:,:)
 real(real64),allocatable:: temperature_tendency(:,:,:),vapor_tendency(:,:,:)
contains
 !> Copy the SGS environment; the host supplies updated values before each advance.
 subroutine lcm_set_lem_environment(eddy_diffusivity,eddy_length,total_water_dry, &
   liquid_static_temperature,static_temperature_tendency,vapor_mixing_ratio_tendency)
  real(real64),intent(in)::eddy_diffusivity(:,:,:),eddy_length(:,:,:),total_water_dry(:,:,:)
  real(real64),intent(in)::liquid_static_temperature(:,:,:)
  real(real64),intent(in)::static_temperature_tendency(:,:,:),vapor_mixing_ratio_tendency(:,:,:)
  integer::shape_expected(3)
  if(.not.grid_state%is_set)error stop 'LEM environment requires grid setup'
  shape_expected=[grid_state%grid%nx_local,grid_state%grid%ny_local,grid_state%grid%nz]
  if(any(shape(eddy_diffusivity)/=shape_expected).or.any(shape(eddy_length)/=shape_expected).or. &
    any(shape(total_water_dry)/=shape_expected).or.any(shape(liquid_static_temperature)/=shape_expected).or. &
    any(shape(static_temperature_tendency)/=shape_expected).or. &
    any(shape(vapor_mixing_ratio_tendency)/=shape_expected))error stop 'LEM environment shape mismatch'
  if(.not.allocated(diffusivity))then
   allocate(diffusivity,mold=eddy_diffusivity)
   allocate(mixing_length,mold=eddy_length)
   allocate(total_water,mold=total_water_dry)
   allocate(liquid_static_energy,mold=liquid_static_temperature)
   allocate(temperature_tendency,mold=static_temperature_tendency)
   allocate(vapor_tendency,mold=vapor_mixing_ratio_tendency)
  endif
  diffusivity=eddy_diffusivity;mixing_length=eddy_length
  total_water=total_water_dry;liquid_static_energy=liquid_static_temperature
  temperature_tendency=static_temperature_tendency;vapor_tendency=vapor_mixing_ratio_tendency
 end subroutine
 subroutine lcm_validate_lem_environment()
  if(.not.allocated(diffusivity))error stop 'LEM advance requires supplied SGS environment'
  if(.not.all(ieee_is_finite(diffusivity)).or..not.all(ieee_is_finite(mixing_length)).or. &
     .not.all(ieee_is_finite(total_water)).or..not.all(ieee_is_finite(liquid_static_energy)).or. &
     .not.all(ieee_is_finite(temperature_tendency)).or..not.all(ieee_is_finite(vapor_tendency))) &
    error stop 'Nonfinite LEM environment'
  if(any(diffusivity<0).or.any(mixing_length<0).or.any(total_water<0))error stop 'Invalid LEM environment'
 end subroutine
 subroutine lcm_apply_sgs_mixing(dt)
  real(real64),intent(in)::dt
  integer::i,j,k
  call lcm_validate_lem_environment()
  do k=1,grid_state%grid%nz
   do j=1,grid_state%grid%ny_local
    do i=1,grid_state%grid%nx_local
     if(cell_particles(i,j,k)%count>0)call mix_cell(i,j,k,dt)
    enddo
   enddo
  enddo
 end subroutine
 subroutine mix_cell(i,j,k,dt)
  integer,intent(in)::i,j,k
  real(real64),intent(in)::dt
  integer::n,m,l,Np,ii,jj,kk,gi,gj,cube(-1:1,-1:1,-1:1)
  ! SAM retains this permutation between cells and timesteps before reshuffling.
  integer,save::origin_permutation(27)=[1,2,3,4,5,6,7,8,9,10,11,12,13,14, &
    15,16,17,18,19,20,21,22,23,24,25,26,27]
  integer::nsub,nrearr,s,t,tmp,keytmp,idtmp,stream
  integer,allocatable::idx(:),keys(:)
  real(real64),allocatable::phi(:,:),nextphi(:,:)
  real(real64)::rho,V,carrier,rl_mean,rl,tprev,pprev,qprev,eta_target
  real(real64)::dz,ds,eta,Lt,eddy_D,molecular_D,lambda,h,phi_tmp(3),tau,tau_sed,eddy_size
  ! Cells are not TARGETs; use direct cell access to keep storage ownership unchanged.
  Np=cell_particles(i,j,k)%count
  if(Np>=100000)error stop 'SAM LEM ordering supports fewer than 100000 particles per cell'
  rho=field_state%density(k,1,1)
  if(field_state%dry_air_density_is_set)rho=field_state%dry_air_density(i,j,k)
  V=real(grid_state%dV(k),real64)
  if(.not.ieee_is_finite(rho).or..not.ieee_is_finite(V).or..not.ieee_is_finite(dt).or. &
     rho<=0.or.V<=0.or.dt<=0)error stop 'Invalid LEM mixing geometry or timestep'
  carrier=rho*V/real(Np,real64)
  rl_mean=sum(cell_particles(i,j,k)%particles(1:Np)%mass_water)/(rho*V)
  call lcm_get_cell_mean_history(i,j,k,tprev,pprev,qprev)
  eta_target=field_state%qv(i,j,k)-qsatw_liquid(tprev,pprev)
  gi=grid_state%grid%global_i_start+i-1;gj=grid_state%grid%global_j_start+j-1
  if(Np>10)then
   allocate(phi(3,Np),nextphi(3,Np),idx(Np),keys(Np))
   do n=1,Np
    cell_particles(i,j,k)%particles(n)%temperature=cell_particles(i,j,k)%particles(n)%temperature &
      +temperature_tendency(i,j,k)*dt
    cell_particles(i,j,k)%particles(n)%qv=cell_particles(i,j,k)%particles(n)%qv+vapor_tendency(i,j,k)*dt
   enddo
   stream=parallel_state%parallel%rank+1
   do n=27,2,-1
    m=min(n,max(1,int(lcm_random_uniform(stream)*n)+1))
    tmp=origin_permutation(n)
    origin_permutation(n)=origin_permutation(m)
    origin_permutation(m)=tmp
   enddo
   cube=reshape(origin_permutation,shape(cube))
   do n=1,Np
    ii=cell_particles(i,j,k)%particles(n)%last_i-gi
    jj=cell_particles(i,j,k)%particles(n)%last_j-gj
    if(ii>grid_state%grid%nx_global/2)ii=ii-grid_state%grid%nx_global
    if(2*ii<=-grid_state%grid%nx_global)ii=ii+grid_state%grid%nx_global
    if(jj>grid_state%grid%ny_global/2)jj=jj-grid_state%grid%ny_global
    if(2*jj<=-grid_state%grid%ny_global)jj=jj+grid_state%grid%ny_global
    kk=cell_particles(i,j,k)%particles(n)%last_k-k
    if(max(abs(ii),abs(jj),abs(kk))>1)error stop 'Particle crossed more than one cell before LEM mixing'
    keys(n)=cube(ii,jj,kk)*100000+cell_particles(i,j,k)%particles(n)%cell_particle_index
    idx(n)=n
    phi(:,n)=[cell_particles(i,j,k)%particles(n)%temperature, &
      cell_particles(i,j,k)%particles(n)%qv,cell_particles(i,j,k)%particles(n)%eta]
   enddo
   ! Stable insertion sort, with the bound checked before array access.
   do n=2,Np
    keytmp=keys(n);idtmp=idx(n);phi_tmp=phi(:,n);m=n-1
    do while(m>=1)
     if(keys(m)<=keytmp)exit
     keys(m+1)=keys(m);idx(m+1)=idx(m);phi(:,m+1)=phi(:,m);m=m-1
    enddo
    keys(m+1)=keytmp;idx(m+1)=idtmp;phi(:,m+1)=phi_tmp
   enddo
   ! Preserve SAM's first-layer LEM line length, including on stretched grids.
   ! Carrier-air mass above still uses this cell's actual volume, as in SAM.
   dz=grid_state%grid%dz(1);ds=dz/real(Np,real64);eta=6*ds
   Lt=max(min(mixing_length(i,j,k),dz),.001_real64)
   eddy_D=max(diffusivity(i,j,k),0.0_real64)
   molecular_D=max(eddy_D*(eta/Lt)**(4.0_real64/3.0_real64),molecular_viscosity_air/rho)
   lambda=54.0_real64/5*eddy_D/Lt**3*(Lt/eta)**(5.0_real64/3.0_real64)
   nsub=max(1,ceiling(dt/(.2_real64*ds**2/molecular_D)));h=dt/real(nsub,real64)
   nrearr=max(1,ceiling(h*lambda*dz))
   do s=1,nsub
    call lem_diffuse_periodic(phi,nextphi,h*molecular_D/ds**2)
    phi=nextphi
    if(lambda>0)then
     do t=1,nrearr
      if(h*lambda*dz>real(lcm_random_uniform(stream),real64))then
       if(Np>6.and.Lt>=eta)then
        m=min(Np,int(lcm_random_uniform(stream)*Np)+1)
        if(Lt/eta<1.01_real64)then
         eddy_size=eta
        else
         do
          eddy_size=(eta**(-5.0_real64/3.0_real64)-real(lcm_random_uniform(stream),real64)* &
             (eta**(-5.0_real64/3.0_real64)-Lt**(-5.0_real64/3.0_real64)))**(-.6_real64)
          if(eddy_size>=eta.and.eddy_size<=min(Lt,dz))exit
         enddo
        endif
        l=min(max(nint(real(floor(eddy_size/ds),real64)/3)*3,6),(Np/3)*3)
        call lem_triplet_map(phi,idx,m,l)
       endif
      endif
     enddo
    endif
   enddo
   do m=1,Np
    n=idx(m)
    cell_particles(i,j,k)%particles(n)%temperature=phi(1,m)
    cell_particles(i,j,k)%particles(n)%qv=phi(2,m)
    cell_particles(i,j,k)%particles(n)%eta=phi(3,m)
    cell_particles(i,j,k)%particles(n)%cell_particle_index=m
    rl=cell_particles(i,j,k)%particles(n)%mass_water/carrier
    tau_sed=ds/(cell_particles(i,j,k)%particles(n)%w_sedi+1.e-20_real64)
    tau=1/(1/tau_sed+1/900.0_real64)
    if(dt>tau)error stop 'Fixed timestep exceeds explicit SAM LEM relaxation timescale'
    cell_particles(i,j,k)%particles(n)%temperature=phi(1,m) &
      -((phi(1,m)-fac_cond*rl)-liquid_static_energy(i,j,k))*dt/tau
    cell_particles(i,j,k)%particles(n)%qv=phi(2,m)-((phi(2,m)+rl)-total_water(i,j,k))*dt/tau
    cell_particles(i,j,k)%particles(n)%eta=phi(3,m)-(phi(3,m)-eta_target)*dt/tau
   enddo
  else
   do n=1,Np
    rl=cell_particles(i,j,k)%particles(n)%mass_water/carrier
    cell_particles(i,j,k)%particles(n)%temperature=liquid_static_energy(i,j,k)+fac_cond*rl_mean
    cell_particles(i,j,k)%particles(n)%qv=total_water(i,j,k)-rl
    cell_particles(i,j,k)%particles(n)%eta=eta_target
    cell_particles(i,j,k)%particles(n)%cell_particle_index=n
   enddo
  endif
  do n=1,Np
   cell_particles(i,j,k)%particles(n)%last_i=gi;cell_particles(i,j,k)%particles(n)%last_j=gj
   cell_particles(i,j,k)%particles(n)%last_k=k
   if(.not.ieee_is_finite(cell_particles(i,j,k)%particles(n)%temperature).or. &
      .not.ieee_is_finite(cell_particles(i,j,k)%particles(n)%qv).or. &
      .not.ieee_is_finite(cell_particles(i,j,k)%particles(n)%eta).or. &
      cell_particles(i,j,k)%particles(n)%temperature-cell_particles(i,j,k)%particles(n)%z* &
      gravity/specific_heat_dry_air<=0.or.cell_particles(i,j,k)%particles(n)%qv<0)error stop 'Invalid particle thermodynamics after LEM mixing'
  enddo
 end subroutine
 subroutine lem_diffuse_periodic(old,new,coefficient)
  real(real64),intent(in)::old(:,:),coefficient
  real(real64),intent(out)::new(:,:)
  integer::n,Np,left,right
  Np=size(old,2)
  do n=1,Np
   left=modulo(n-2,Np)+1;right=modulo(n,Np)+1
   new(:,n)=old(:,n)+coefficient*(old(:,left)-2*old(:,n)+old(:,right))
  enddo
 end subroutine
 subroutine lem_triplet_map(phi,idx,start,length)
  real(real64),intent(inout)::phi(:,:)
  integer,intent(inout)::idx(:)
  integer,intent(in)::start,length
  real(real64)::tmp(size(phi,1),length)
  integer::ids(length),n,m,dest,src,Np
  Np=size(idx)
  do n=1,length
   if(n<=length/3)then
    m=3*(n-1)+1
   else if(n<=2*length/3)then
    m=2*(length-1)-3*(n-1)+1
   else
    m=3*(n-1)-2*(length-1)+1
   endif
   src=modulo(start+m-2,Np)+1
   tmp(:,n)=phi(:,src);ids(n)=idx(src)
  enddo
  do n=1,length
   dest=modulo(start+n-2,Np)+1;phi(:,dest)=tmp(:,n);idx(dest)=ids(n)
  enddo
 end subroutine
end module
