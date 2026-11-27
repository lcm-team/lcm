program test_lem_sgs
 use iso_fortran_env,only:real64,int64
 use minimal_les_host
 use lcm_sgs_mixing
 use lcm_particle_state,only:lcm_reset_particle_cells,lcm_add_particle,cell_particles
 use lcm_parameters,only:fac_cond,pi,rho_liq
 use lcm_types,only:lcm_particle_t
 implicit none
 real(real64)::old(3,12),new(3,12),before(3),m,carrier,rl,expected,Tsum,qsum,varold,varnew
 real(real64)::D(2,2,2),L(2,2,2),rt(2,2,2),hl(2,2,2),Ts(2,2,2),qs(2,2,2)
 integer::idx(12),i,j,n,np,cas
 type(lcm_particle_t)::p
 ! Independent periodic stencil and triplet-map invariants.
 do i=1,12
  old(:,i)=[real(i*i,real64),sin(real(i,real64)),1.0_real64];idx(i)=i
 enddo
 call lem_diffuse_periodic(old,new,.2_real64)
 call require(maxval(abs(sum(new,dim=2)-sum(old,dim=2)))<1.e-12_real64,'periodic diffusion conserves sums')
 call require(all(new(3,:)==1.0_real64),'constant field unchanged')
 call require(abs(new(1,1)-(old(1,1)+.2_real64*(old(1,12)-2*old(1,1)+old(1,2))))<1.e-12_real64, &
   'first periodic element uses its own old value')
 before=sum(old,dim=2)
 call lem_triplet_map(old,idx,10,9)
 call require(maxval(abs(sum(old,dim=2)-before))<1.e-12_real64,'wrapped triplet conserves sums')
 do i=1,12
  call require(count(idx==i)==1,'triplet is a permutation')
  call require(old(1,i)==real(idx(i)*idx(i),real64),'triplet retains parcel identity')
 enddo
 call setup_minimal_les_host()
 D=0;L=100;rt=.010_real64;hl=300.0_real64;Ts=0;qs=0
 call lcm_set_lem_environment(D,L,rt,hl,Ts,qs)
 do cas=1,5
  call lcm_reset_particle_cells()
  select case(cas)
  case(1);np=4
  case(2);np=10
  case(3);np=11
  case default;np=24
  end select
  D=0;Ts=0;qs=0
  if(cas==4)D=1000
  if(cas==5)then
   Ts=.01_real64;qs=1.e-6_real64
  endif
  carrier=1.2_real64*250000.0_real64/real(np,real64)
  m=4*pi/3*rho_liq*(2.e-6_real64)**3*1.e10_real64;rl=m/carrier
  varold=0
  do n=1,np
   p=lcm_particle_t();p%id=int(n,int64);p%weight=10000000000_int64;p%mass_water=m
   p%mass_aero=1.e-9_real64;p%aerosol_portions(1)=1
   p%temperature=hl(1,1,1)+fac_cond*rl+1
   p%qv=rt(1,1,1)-rl+.001_real64
   p%eta=.0001_real64;p%cell_particle_index=n;p%last_i=1;p%last_j=1;p%last_k=1
   if(cas==4)then
    p%temperature=hl(1,1,1)+fac_cond*rl+sin(2*pi*real(n,real64)/real(np,real64))
    p%qv=rt(1,1,1)-rl+.001_real64*sin(2*pi*real(n,real64)/real(np,real64))
    varold=varold+(p%temperature-hl(1,1,1)-fac_cond*rl)**2
   endif
   call lcm_add_particle(1,1,1,p)
  enddo
  call lcm_set_lem_environment(D,L,rt,hl,Ts,qs)
  ! Host changes after copying must not alter the SGS environment in LCM.
  D=-D;Ts=-Ts;qs=-qs;rt=-rt;hl=-hl;L=-L
  call lcm_apply_sgs_mixing(.5_real64)
  D=-D;Ts=-Ts;qs=-qs;rt=-rt;hl=-hl;L=-L
  Tsum=0;qsum=0;varnew=0
  do n=1,np
   p=cell_particles(1,1,1)%particles(n)
   call require(p%mass_water==m,'mixing does not change liquid mass')
   Tsum=Tsum+p%temperature-fac_cond*rl;qsum=qsum+p%qv+rl
   if(cas<=2)then
    call require(abs(p%qv-(rt(1,1,1)-rl))<1.e-14_real64,'sparse total-water reconstruction')
    call require(abs(p%temperature-(hl(1,1,1)+fac_cond*rl))<1.e-12_real64,'sparse thermal reconstruction')
   else if(cas/=4)then
    expected=rt(1,1,1)-rl+.001_real64+qs(1,1,1)*.5_real64 &
      -(.001_real64+qs(1,1,1)*.5_real64)*.5_real64/900
    call require(abs(p%qv-expected)<1.e-14_real64,'SAM explicit conserved-water relaxation')
    expected=hl(1,1,1)+fac_cond*rl+1+Ts(1,1,1)*.5_real64 &
      -(1+Ts(1,1,1)*.5_real64)*.5_real64/900
    call require(abs(p%temperature-expected)<1.e-11_real64,'SAM explicit liquid-static relaxation')
   endif
   varnew=varnew+(p%temperature-hl(1,1,1)-fac_cond*rl)**2
  enddo
  if(cas==4)then
   call require(abs(Tsum/np-hl(1,1,1))<1.e-11_real64,'mixing preserves mean liquid-static energy')
   call require(abs(qsum/np-rt(1,1,1))<1.e-14_real64,'mixing preserves mean total water')
   call require(varnew<varold,'diffusion decreases variance')
  endif
 enddo
 print *,'PASS: periodic diffusion, wrapped triplet, sparse fallback, conserved relaxation, forcing, turbulent mean conservation.'
contains
 subroutine require(ok,message)
 logical,intent(in)::ok
 character(len=*),intent(in)::message
 if(.not.ok)then
  print *,'FAIL: ',message
  error stop 1
 endif
 end subroutine
end program
