subroutine dep_curr(coord,coord_cntr,curr,leftX,Rgrid,dx_inv,dr_inv,&
                    np,nx,nr,nkO)
implicit none
integer, intent(in)        :: np,nx,nr,nkO
real (kind=8), intent(in)  :: coord_cntr(3,np),coord(8,np),leftX,Rgrid(0:nr),&
                              dx_inv,dr_inv
complex(kind=8),intent(inout):: curr(0:nx,0:nr,0:nkO,3)
integer         :: ix,ir,ip,k,l,iO
real(kind=8)    :: xp,yp,zp,rp,wp,S0(0:1,2),veloc(3),curr_p(0:1,0:1)
complex(kind=8) :: ii=(0.0d0,1.0d0),phaseO(0:nkO),phase_m

!f2py intent(in) :: coord,coord_cntr,leftX,Rgrid,dx_inv,dr_inv
!f2py intent(in,out) :: curr
!f2py intent(hide) :: np,nx,nr,nkO

do ip=1,np
  wp = coord(8,ip)/coord(7,ip)
  if(wp == 0.0) CYCLE
  xp = coord_cntr(1,ip)
  yp = coord_cntr(2,ip)
  zp = coord_cntr(3,ip)
  veloc = coord(4:6,ip)
  rp = DSQRT(yp*yp+zp*zp)
  if ((rp>=Rgrid(nr)) .or. (SUM(ABS(veloc)).eq. 0.0d0)) CYCLE
  veloc = veloc*wp

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

  phaseO(0) = 1.0d0
  if (nkO>0) then
    do iO = 1,nkO
      phaseO(iO) = phaseO(iO-1)*phase_m
    enddo
  endif

  do k = 0,1
    curr_p(:,k) = S0(:,1)*S0(k,2)
  enddo

  do l=1,3
    do iO = 0,nkO
      curr(ix:ix+1,ir:ir+1,iO,l) = curr(ix:ix+1,ir:ir+1,iO,l)+ phaseO(iO)*veloc(l)*curr_p
    enddo
  enddo
enddo

if (Rgrid(0)<0) then
  curr(:,1,0,:) = curr(:,1,0,:) + curr(:,0,0,:)
  if (nkO>0) curr(:,1,1:nkO,:) = curr(:,2,1:nkO,:)/9.0d0
  curr(:,0,:,:) = 0.0
endif

end subroutine

subroutine dep_dens(coord,dens,leftX,Rgrid,dx_inv,dr_inv,&
                    np,nx,nr,nkO)
implicit none
integer, intent(in)        :: np,nx,nr,nkO
real (kind=8), intent(in)  :: coord(8,np),leftX,Rgrid(0:nr),&
                              dx_inv,dr_inv
complex(kind=8),intent(inout):: dens(0:nx,0:nr,0:nkO)
integer         :: ix,ir,ip,k,iO
real(kind=8)    :: xp,yp,zp,rp,wp,S0(0:1,2), dens_p(0:1,0:1)
complex(kind=8) :: ii=(0.0d0,1.0d0),phaseO(0:nkO),phase_m

!f2py intent(in) :: coord,leftX,Rgrid,dx_inv,dr_inv
!f2py intent(in,out) :: dens
!f2py intent(hide) :: np,nx,nr,nkO

do ip=1,np
  wp = coord(8,ip)
  if (wp == 0.0) CYCLE
  xp = coord(1,ip)
  yp = coord(2,ip)
  zp = coord(3,ip)
  rp = DSQRT(yp*yp+zp*zp)
  if (rp>=Rgrid(nr)) CYCLE

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

  phaseO(0) = 1.0d0
  if (nkO>0) then
    do iO = 1,nkO
      phaseO(iO) = phaseO(iO-1)*phase_m
    enddo
  endif

  do k = 0,1
    dens_p(:,k) = S0(:,1)*S0(k,2)*wp
  enddo

  do iO = 0,nkO
    dens(ix:ix+1,ir:ir+1,iO) = dens(ix:ix+1,ir:ir+1,iO)+ phaseO(iO)*dens_p
  enddo
enddo

if (Rgrid(0)<0) then
  dens(:,1,0) = dens(:,1,0) + dens(:,0,0)
  if (nkO>0) dens(:,1,1:nkO) = dens(:,2,1:nkO)/9.0d0
  dens(:,0,:) = 0.0
endif

end subroutine

subroutine proj_fld(coord,Fld,Fld_tot,leftX,Rgrid,dx_inv,&
                             dr_inv,np,nx,nr,nkO)
implicit none
integer, intent(in)          :: np,nx,nr,nkO
real (kind=8),    intent(in) :: coord(8,np),leftX,Rgrid(0:nr),dx_inv,dr_inv
complex (kind=8), intent(in) :: Fld(0:nx,0:nr,0:nkO,6)
real (kind=8), intent(inout) :: Fld_tot(6,np)
integer         :: ix,ir,ip,iO,k,l
real(kind=8)    :: xp,yp,zp,wp,rp,S0(0:1,2),Fld_p(6)
complex(kind=8) :: ii=(0.0d0,1.0d0),phase_p,phaseO(0:nkO),projcomp(0:1,0:1,0:nkO)

!f2py intent(in)     :: coord,Fld,leftX,Rgrid,dx_inv,dr_inv
!f2py intent(in,out) :: Fld_tot
!f2py intent(hide)   :: np,nx,nr,nkO

!$omp parallel shared(Fld_tot,coord,dx_inv,dr_inv,leftX,Rgrid,ii,np,nx,nr,nkO,Fld)
!$omp do schedule(static) private(ix,ir,ip,k,l,iO,xp,yp,zp,rp,wp,S0,Fld_p,phaseO,phase_p,projcomp)
do ip=1,np
  wp = coord(8,ip)
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
    phase_p = 0.0d0
  endif

  phaseO(0) = 1.0
  if (nkO>0) then
    do iO = 1,nkO
      phaseO(iO) = phaseO(iO-1)*phase_p
    enddo
  endif

  projcomp = 0.0
  do iO = 0,nkO
    do k=0,1
      projcomp(:,k,iO) = projcomp(:,k,iO)+ S0(k,2)*S0(:,1)*phaseO(iO)
    enddo
  enddo

  Fld_p = 0.0d0
  do l=1,6
    Fld_p(l) = Fld_p(l) + SUM(DBLE(projcomp*Fld(ix:ix+1,ir:ir+1,:,l)))
  enddo

  Fld_tot(:,ip) = Fld_tot(:,ip)+Fld_p
enddo
!$omp end do      
!$omp end parallel
end subroutine

subroutine dep_curr_env(coord,coord_cntr,curr,leftX,Rgrid,dx_inv,dr_inv,&
                    kx0,np,nx,nr,nkO)
implicit none
integer, intent(in)        :: np,nx,nr,nkO
real (kind=8), intent(in)  :: coord_cntr(3,np),coord(8,np),leftX,Rgrid(0:nr),&
                              dx_inv,dr_inv,kx0
complex(kind=8),intent(inout):: curr(0:nx,0:nr,-nkO:nkO,3)
integer         :: ix,ir,ip,k,l,iO
real(kind=8)    :: xp,yp,zp,rp,S0(0:1,2)
complex(kind=8) :: ii=(0.0d0,1.0d0),curr_p(0:1,0:1),phaseO(0:nkO),&
                   phase_m,wp

!f2py intent(in) :: coord,coord_cntr,leftX,Rgrid,dx_inv,dr_inv,kx0
!f2py intent(in,out) :: curr
!f2py intent(hide) :: np,nx,nr,nkO

do ip=1,np
  wp = coord(8,ip)/coord(7,ip)
  if(wp == 0.0) CYCLE
  xp = coord_cntr(1,ip)
  yp = coord_cntr(2,ip)
  zp = coord_cntr(3,ip)
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
    curr_p(:,k) = S0(:,1)*S0(k,2)
  enddo
  curr_p = curr_p*wp

  do l=1,3
    do iO = 0,nkO
      curr(ix:ix+1,ir:ir+1,iO,l) = curr(ix:ix+1,ir:ir+1,iO,l)+ phaseO(iO)*coord(3+l,ip)*curr_p
      if (iO>0) curr(ix:ix+1,ir:ir+1,-iO,l) = curr(ix:ix+1,ir:ir+1,-iO,l)+ &
        CONJG(phaseO(iO))*coord(3+l,ip)*curr_p
    enddo
  enddo
enddo
 
if (Rgrid(0)<0) then
  curr(:,1,nkO,:) = curr(:,1,nkO,:) + curr(:,0,nkO,:)
  if (nkO>0) then
    curr(:,1,1:nkO,:) = curr(:,2,1:nkO,:)/9.0d0
    curr(:,1,-nkO:-1,:) = curr(:,2,-nkO:-1,:)/9.0d0
  endif
  curr(:,0,:,:) = 0.0
endif

end subroutine

subroutine dep_dens_env(coord,dens,leftX,Rgrid,dx_inv,dr_inv,&
                    kx0,np,nx,nr,nkO)
implicit none
integer, intent(in)        :: np,nx,nr,nkO
real (kind=8), intent(in)  :: coord(8,np),leftX,Rgrid(0:nr),&
                              dx_inv,dr_inv,kx0
complex(kind=8),intent(inout):: dens(0:nx,0:nr,-nkO:nkO)
integer         :: ix,ir,ip,k,iO
real(kind=8)    :: xp,yp,zp,rp,S0(0:1,2)
complex(kind=8) :: ii=(0.0d0,1.0d0),dens_p(0:1,0:1),phaseO(0:nkO),phase_m,wp

!f2py intent(in) :: coord,leftX,Rgrid,dx_inv,dr_inv,kx0
!f2py intent(in,out) :: dens
!f2py intent(hide) :: np,nx,nr,nkO

do ip=1,np
  wp = coord(8,ip)
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

if (Rgrid(0)<0) then
  dens(:,1,0) = dens(:,1,0) + dens(:,0,0)
  if (nkO>0) then
    dens(:,1,1:nkO) = dens(:,2,1:nkO)/9.0d0
    dens(:,1,-nkO:-1:1) = dens(:,2,-nkO:-1:-1)/9.0d0
  endif
  dens(:,0,:) = 0.0
endif

end subroutine

subroutine proj_fld_env(coord,Fld,Fld_tot,leftX,Rgrid,dx_inv,&
                             dr_inv,kx0,np,nx,nr,nkO)
implicit none
integer, intent(in)          :: np,nx,nr,nkO
real (kind=8),    intent(in) :: coord(8,np),leftX,Rgrid(0:nr),dx_inv,dr_inv,kx0
complex (kind=8), intent(in) :: Fld(0:nx,0:nr,-nkO:nkO,6)
real (kind=8), intent(inout) :: Fld_tot(6,np)
integer         :: ix,ir,ip,iO,k,l
real(kind=8)    :: xp,yp,zp,wp,rp,S0(0:1,2),Fld_p(6)!,pi2m1=0.125d0/DATAN(1.d0)
complex(kind=8) :: ii=(0.0d0,1.0d0),phase_p,phaseO(0:nkO),projcomp(0:1,0:1),projcompO(0:1,0:1)

!f2py intent(in)     :: coord,Fld,leftX,Rgrid,dx_inv,dr_inv,kx0
!f2py intent(in,out) :: Fld_tot
!f2py intent(hide)   :: np,nx,nr,nkO

!$omp parallel default(shared) private(ix,ir,ip,k,l,iO,xp,yp,zp,rp,wp,S0,Fld_p,phaseO,phase_p,projcomp,projcompO)
!$omp do schedule(static)
do ip=1,np
  wp = coord(8,ip)
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
  Fld_tot(:,ip) = Fld_tot(:,ip)+Fld_p
enddo
!$omp end do      
!$omp end parallel
end subroutine


