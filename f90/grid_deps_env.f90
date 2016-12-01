subroutine dep_curr_env(coord,momenta,wghts,curr,leftX,Rgrid,dx_inv,dr_inv,&
                    kx0,np,nx,nr,nkO)
implicit none
integer, intent(in)        :: np,nx,nr,nkO
real (kind=8), intent(in)  :: coord(3,np),momenta(3,np),wghts(np),leftX,Rgrid(0:nr),&
                              dx_inv,dr_inv,kx0
complex(kind=8),intent(inout):: curr(0:nx,0:nr,-nkO:nkO,3)
integer         :: ix,ir,ip,k,l,iO
real(kind=8)    :: xp,yp,zp,rp,gp,veloc(3),S0(0:1,2), inv9=1.0d0/9.0d0
complex(kind=8) :: ii=(0.0d0,1.0d0),curr_p(0:1,0:1),phaseO(0:nkO),&
                   phase_m,wp

!f2py intent(in) :: coord,momenta,wghts,leftX,Rgrid,dx_inv,dr_inv,kx0
!f2py intent(in,out) :: curr
!f2py intent(hide) :: np,nx,nr,nkO

do ip=1,np
  wp = wghts(ip)
  if(wp == 0.0) CYCLE
  xp = coord(1,ip)
  yp = coord(2,ip)
  zp = coord(3,ip)
  veloc = momenta(:,ip)
  rp = DSQRT(yp*yp+zp*zp)
  if ((rp>=Rgrid(nr)) .or. (SUM(ABS(veloc)).eq. 0.0d0)) CYCLE
  gp = DSQRT(1.0+veloc(1)*veloc(1)+veloc(2)*veloc(2)+veloc(3)*veloc(3))
  veloc = veloc/gp

  wp = wp*(dcos(xp*kx0) - ii*dsin(xp*kx0))

  ix = FLOOR((xp-leftX)*dx_inv)
  ir = FLOOR((rp-Rgrid(0))*dr_inv)

  S0 = 0.0d0
  S0(1,1) = (xp-leftX)*dx_inv - ix
  S0(0,1) = 1.0d0 -S0(1,1)
  S0(1,2) = (rp-Rgrid(ir))*dr_inv
  S0(0,2) = 1.0d0 -S0(1,2)

  if (rp>0.0) then
    phase_m = (yp-ii*zp)/rp
  else
    phase_m = 0.0d0
  endif

  phaseO(0) = 1.0
  if (nkO>0) then
    do iO = 1,nkO
      phaseO(iO) = phaseO(iO-1)*phase_m
    enddo
  endif

  do k = 0,1
    curr_p(:,k) = S0(:,1)*S0(k,2)
  enddo
  curr_p = curr_p*wp

  do l=1,3
    do iO = 0,nkO
      curr(ix:ix+1,ir:ir+1,iO,l) = curr(ix:ix+1,ir:ir+1,iO,l)+ phaseO(iO)*veloc(l)*curr_p
      if (iO>0) curr(ix:ix+1,ir:ir+1,-iO,l) = curr(ix:ix+1,ir:ir+1,-iO,l)+ &
        CONJG(phaseO(iO))*veloc(l)*curr_p
    enddo
  enddo
enddo
 
!$omp parallel do default(shared) private(l) schedule(static)
do l=1,3
  curr(:,1,0,l) = curr(:,1,0,l) - curr(:,0,0,l)
  if (nkO>0) then
    curr(:,1,1:nkO,l) = curr(:,2,1:nkO,l)*inv9
    curr(:,1,-nkO:-1,l) = curr(:,2,-nkO:-1,l)*inv9
  endif
  curr(:,0,:,l) = 0.0
enddo
!$omp end parallel do

end subroutine

subroutine dep_dens_env(coord,wghts,dens,leftX,Rgrid,dx_inv,dr_inv,&
                    kx0,np,nx,nr,nkO)
implicit none
integer, intent(in)        :: np,nx,nr,nkO
real (kind=8), intent(in)  :: coord(3,np),wghts(np),leftX,Rgrid(0:nr),&
                              dx_inv,dr_inv,kx0
complex(kind=8),intent(inout):: dens(0:nx,0:nr,-nkO:nkO)
integer         :: ix,ir,ip,k,iO
real(kind=8)    :: xp,yp,zp,rp,S0(0:1,2), inv9=1.0d0/9.0d0
complex(kind=8) :: ii=(0.0d0,1.0d0),dens_p(0:1,0:1),phaseO(0:nkO),phase_m,wp

!f2py intent(in) :: coord,wghts,leftX,Rgrid,dx_inv,dr_inv,kx0
!f2py intent(in,out) :: dens
!f2py intent(hide) :: np,nx,nr,nkO

do ip=1,np
  wp = wghts(ip)
  if(wp == 0.0) CYCLE
  xp = coord(1,ip)
  yp = coord(2,ip)
  zp = coord(3,ip)
  rp = DSQRT(yp*yp+zp*zp)
  if(rp>=Rgrid(nr)) CYCLE

  wp = wp*(dcos(xp*kx0) - ii*dsin(xp*kx0))
  ix = FLOOR((xp-leftX)*dx_inv)
  ir = FLOOR((rp-Rgrid(0))*dr_inv)

  S0 = 0.0d0
  S0(1,1) = (xp-leftX)*dx_inv - ix
  S0(0,1) = 1.0d0 -S0(1,1)
  S0(1,2) = (rp-Rgrid(ir))*dr_inv
  S0(0,2) = 1.0d0 -S0(1,2)

  if (rp>0.0) then
    phase_m = (yp-ii*zp)/rp
  else
    phase_m = 0.0d0
  endif

  phaseO(0) = 1.0
  if (nkO>0) then
    do iO = 1,nkO
      phaseO(iO) = phaseO(iO-1)*phase_m
    enddo
  endif

  do k = 0,1
    dens_p(:,k) = S0(:,1)*S0(k,2)*wp
  enddo
  dens_p = dens_p*wp

  do iO = 0,nkO
    dens(ix:ix+1,ir:ir+1,iO) = dens(ix:ix+1,ir:ir+1,iO)+ phaseO(iO)*dens_p
    if (iO>0) dens(ix:ix+1,ir:ir+1,-iO) = dens(ix:ix+1,ir:ir+1,-iO)+ CONJG(phaseO(iO))*dens_p
  enddo
enddo

dens(:,1,0) = dens(:,1,0) - dens(:,0,0)
if (nkO>0) then
  dens(:,1,  1:nkO) = dens(:,2,1:nkO)*inv9
  dens(:,1,-nkO:-1) = dens(:,2,-nkO:-1)*inv9
endif
dens(:,0,:) = 0.0

end subroutine

subroutine proj_fld_env(coord,wghts,Fld,Fld_tot,leftX,Rgrid,dx_inv,&
                             dr_inv,kx0,np,nx,nr,nkO)
implicit none
integer, intent(in)          :: np,nx,nr,nkO
real (kind=8),    intent(in) :: coord(3,np),wghts(np),leftX,Rgrid(0:nr),dx_inv,dr_inv,kx0
complex (kind=8), intent(in) :: Fld(0:nx,0:nr,-nkO:nkO,6)
real (kind=8), intent(inout) :: Fld_tot(6,np)
integer         :: ix,ir,ip,iO,k,l
real(kind=8)    :: xp,yp,zp,wp,rp,S0(0:1,2),Fld_p(6)
complex(kind=8) :: ii=(0.0d0,1.0d0),phase_p,phaseO(0:nkO),projcomp(0:1,0:1),projcompO(0:1,0:1)

!f2py intent(in)     :: coord,wghts,Fld,leftX,Rgrid,dx_inv,dr_inv,kx0
!f2py intent(in,out) :: Fld_tot
!f2py intent(hide)   :: np,nx,nr,nkO

!$omp parallel default(shared) private(ix,ir,ip,k,l,iO,xp,yp,zp,rp,wp,S0,Fld_p,&
!$omp phaseO,phase_p,projcomp,projcompO)
!$omp do schedule(static)
do ip=1,np
  wp = wghts(ip)
  if(wp == 0.0) CYCLE
  xp = coord(1,ip)
  yp = coord(2,ip)
  zp = coord(3,ip)
  rp = DSQRT(yp*yp+zp*zp)

  if(rp>=Rgrid(nr)) CYCLE

  ix = FLOOR((xp-leftX)*dx_inv)
  ir = FLOOR((rp-Rgrid(0))*dr_inv)

  S0 = 0.0d0
  S0(1,1) = (xp-leftX)*dx_inv - ix
  S0(0,1) = 1.0d0 -S0(1,1)
  S0(1,2) = (rp-Rgrid(ir))*dr_inv
  S0(0,2) = 1.0d0 -S0(1,2)

  if (rp>0.0) then
    phase_p = (yp+ii*zp)/rp
  else
    phase_p = 1.0d0
  endif

  phaseO(0) = 1.0d0
  if (nkO>0) then
    do iO = 1,nkO
      phaseO(iO) = phaseO(iO-1)*phase_p
    enddo
  endif

  do k=0,1
    projcomp(:,k) = S0(k,2)*S0(:,1)
  enddo
  projcomp = projcomp*(dcos(xp*kx0) + ii*dsin(xp*kx0))

  Fld_p = 0.0d0
  do iO=0,nkO
    projcompO = projcomp*phaseO(iO)
    do l=1,6
      Fld_p(l) = Fld_p(l) + SUM(DBLE(projcompO*Fld(ix:ix+1,ir:ir+1,iO,l)))
    enddo
    if (iO>0) then
      projcompO = projcomp*CONJG(phaseO(iO))
      do l=1,6
        Fld_p(l) = Fld_p(l) + SUM(DBLE(projcompO*Fld(ix:ix+1,ir:ir+1,-iO,l)))
      enddo
    endif
  enddo
  Fld_tot(:,ip) = Fld_tot(:,ip) + Fld_p
enddo
!$omp end do      
!$omp end parallel
end subroutine

subroutine eb_corr_axis_env(eb_spc,nx,nr,nkO)
implicit none 
integer, intent(in)        :: nx,nr,nkO
complex(kind=8),intent(inout):: eb_spc(0:nx,0:nr,-nkO:nkO,6) 
real (kind=8) :: pi=4.d0*DATAN(1.d0),pi_inv
integer :: l,ix,ir,iO

!f2py intent(in,out) :: eb_spc
!f2py intent(hide) :: nx,nr,nkO

pi_inv = 1./pi

!$omp parallel default(shared) private(ir,iO,l)
do l=1,6
  do iO = -nkO,nkO
    !$omp do schedule(static)
    do ir=0,nr
      eb_spc(:,ir,iO,l) = pi_inv*eb_spc(:,ir,iO,l)
    enddo
    !$omp end do
  enddo
enddo
!$omp end parallel

!$omp parallel default(shared) private(ix,l,iO)
do l=1,6
  do iO = -nkO,nkO
    if (nkO==0) then
      !$omp do schedule(static)
      do ix=0,nx
        eb_spc(ix,0,iO,l) = eb_spc(ix,1,iO,l)
      enddo
      !$omp end do
    else
      !$omp do schedule(static)
      do ix=0,nx
        eb_spc(ix,0,iO,l) = -eb_spc(ix,1,iO,l)
      enddo
      !$omp end do
    endif
  enddo      
enddo
!$omp end parallel
end subroutine
