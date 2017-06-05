! This file is a part of CHIMERA software
! CHIMERA is a simulation code for FEL and laser plasma simulations
! Copyright (C)  2016 Igor A. Andriyash <igor.andriyash@gmail.com>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

subroutine dep_curr(coord,momenta,wghts,curr,leftX,Rgrid,dx_inv,dr_inv,&
                    np,nx,nr,nkO)
implicit none
integer, intent(in)           :: np,nx,nr,nkO
real (kind=8), intent(in)     :: coord(3,np),momenta(3,np),wghts(np)
real (kind=8), intent(in)     :: leftX,Rgrid(0:nr),dx_inv,dr_inv
complex(kind=8),intent(inout) :: curr(0:nx,0:nr,0:nkO,3)

integer         :: ix,ir,ip,k,l,iO
real(kind=8)    :: xp,yp,zp,rp,wp,gp,S0(0:1,2),veloc(3),curr_p(0:1,0:1)
complex(kind=8) :: ii=(0.0d0,1.0d0),phaseO(0:nkO),phase_m

!f2py intent(in) :: coord,momenta,wghts,leftX,Rgrid,dx_inv,dr_inv
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
  gp = DSQRT(1+veloc(1)*veloc(1)+veloc(2)*veloc(2)+veloc(3)*veloc(3))
  veloc = veloc*wp/gp

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
      curr(ix:ix+1,ir:ir+1,iO,l) = curr(ix:ix+1,ir:ir+1,iO,l) &
                                 + phaseO(iO)*veloc(l)*curr_p
    enddo
  enddo
enddo

!$omp parallel do default(shared) private(l) schedule(static)
do l=1,3
  curr(:,1,:,l) = curr(:,1,:,l) - curr(:,0,:,l)
  curr(:,0,:,l) = 0.0
enddo
!$omp end parallel do

end subroutine

subroutine dep_dens(coord,wghts,dens,leftX,Rgrid,dx_inv,dr_inv,&
                    np,nx,nr,nkO)
implicit none
integer, intent(in)           :: np,nx,nr,nkO
real (kind=8), intent(in)     :: coord(3,np),wghts(np)
real (kind=8), intent(in)     :: leftX,Rgrid(0:nr),dx_inv,dr_inv
complex(kind=8),intent(inout) :: dens(0:nx,0:nr,0:nkO)
integer         :: ix,ir,ip,k,iO
real(kind=8)    :: xp,yp,zp,rp,wp,S0(0:1,2), dens_p(0:1,0:1)
complex(kind=8) :: ii=(0.0d0,1.0d0),phaseO(0:nkO),phase_m

!f2py intent(in) :: coord,wghts,leftX,Rgrid,dx_inv,dr_inv
!f2py intent(in,out) :: dens
!f2py intent(hide) :: np,nx,nr,nkO

do ip=1,np
  wp = wghts(ip)
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

dens(:,1,:) = dens(:,1,:) - dens(:,0,:)
dens(:,0,:) = 0.0

end subroutine

subroutine proj_fld(coord,wghts,Fld,Fld_tot,leftX,Rgrid,dx_inv,&
                             dr_inv,np,nx,nr,nkO)
implicit none
integer, intent(in)          :: np,nx,nr,nkO
real (kind=8),    intent(in) :: coord(3,np),wghts(np)
real (kind=8),    intent(in) :: leftX,Rgrid(0:nr),dx_inv,dr_inv
complex (kind=8), intent(in) :: Fld(0:nx,0:nr,0:nkO,6)
real (kind=8), intent(inout) :: Fld_tot(6,np)
integer         :: ix,ir,ip,iO,k,l
real(kind=8)    :: xp,yp,zp,wp,rp,S0(0:1,2),Fld_p(6)
complex(kind=8) :: ii=(0.0d0,1.0d0),phase_p,phaseO(0:nkO)
complex(kind=8) :: projcomp(0:1,0:1,0:nkO)

!f2py intent(in)     :: coord,wghts,Fld,leftX,Rgrid,dx_inv,dr_inv
!f2py intent(in,out) :: Fld_tot
!f2py intent(hide)   :: np,nx,nr,nkO

!$omp parallel default(shared) private(xp,yp,zp,rp,wp,S0,Fld_p, &
!$omp    phaseO,phase_p,projcomp,ix,ir,ip,k,l,iO)
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

  Fld_tot(:,ip) = Fld_tot(:,ip) + Fld_p
enddo
!$omp end do      
!$omp end parallel
end subroutine

subroutine eb_corr_axis(eb_spc,nx,nr,nkO)
implicit none 
integer, intent(in)        :: nx,nr,nkO
complex(kind=8),intent(inout):: eb_spc(0:nx,0:nr,nkO,6) 
real (kind=8) :: pi=4.d0*DATAN(1.d0), pi2_inv,pi_inv
integer :: l,ix,ir,iO

!f2py intent(in,out) :: eb_spc
!f2py intent(hide) :: nx,nr,nkO

pi_inv = 1./pi
pi2_inv = 0.5*pi_inv

!$omp parallel default(shared) private(ir,iO,l)
do l=1,6
  !$omp do schedule(static)
  do ir=0,nr
    eb_spc(:,ir,1,l) = pi2_inv*eb_spc(:,ir,1,l)

    if (nkO>1) then
      do iO = 2,nkO
        eb_spc(:,ir,iO,l) = pi_inv*eb_spc(:,ir,iO,l)
      enddo
    endif
  enddo
  !$omp end do
enddo
!$omp end parallel

!$omp parallel default(shared) private(ix,l,iO)
do l=1,6
  !$omp do schedule(static)
  do ix=0,nx
    eb_spc(ix,0,1,l) = eb_spc(ix,1,1,l)
  enddo
  !$omp end do
  if (nkO>1) then
    do iO = 2,nkO
      !$omp do schedule(static)
      do ix=0,nx
        eb_spc(ix,0,iO,l) = -eb_spc(ix,1,iO,l)
      enddo
      !$omp end do
    enddo
  endif
enddo
!$omp end parallel
end subroutine
