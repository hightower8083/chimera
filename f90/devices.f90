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

subroutine Undul_Mapped(coord,Fld,t,a0,params,np,nx)

implicit none
integer, intent(in)        :: np,nx
real (kind=8), intent(in)  :: t,a0(2,nx),params(3)
real (kind=8), intent(in)  :: coord(3,np)
real (kind=8), intent(inout)::Fld(6,np)
real (kind=8)              :: S0(-1:1),Xleft,Xright,lambda,dx,dx_inv,&
                              ddx,xp,yp,ku,pi=4.d0*DATAN(1.d0)
integer                    :: ip,ix

!f2py intent(in)  :: coord,t,a0,params
!f2py intent(in,out) :: Fld
!f2py intent(hide):: np,nx

lambda = params(1)
Xleft  = params(2)
dx     = params(3)

ku     = 2.d0*pi/lambda
dx_inv = 1.0d0/dx
Xright= Xleft+nx*dx

!$omp parallel default(private) shared(coord,Fld,t,a0,params,np,nx, &
!$omp                                  lambda,Xleft,dx,ku,dx_inv,Xright)
!$omp do schedule(static)
do ip=1,np
  xp = coord(1,ip)
  if ((xp<Xleft+dx) .or. (xp>Xright-dx)) CYCLE
  yp = coord(2,ip)

  S0 = 0.0d0
  ix  = FLOOR((xp-Xleft)*dx_inv +0.5d0)
  ddx = (xp-Xleft)*dx_inv - ix
  S0(-1) = 0.5d0*(0.5d0-ddx)**2
  S0( 0) = 0.75d0-ddx**2
  S0( 1) = 0.5d0*(0.5d0+ddx)**2

  Fld(5,ip) = Fld(5,ip) + SUM(S0*a0(1,ix-1:ix+1))*cosh(ku*yp)
  Fld(4,ip) = Fld(4,ip) + SUM(S0*a0(2,ix-1:ix+1))*sinh(ku*yp)
enddo
!$omp end do
!$omp end parallel

end subroutine

subroutine Undul_Mapped_Tap(coord,Fld,t,a0,params,np,nx)

implicit none
integer, intent(in)        :: np,nx
real (kind=8), intent(in)  :: t,a0(2,nx),params(5)
real (kind=8), intent(in)  :: coord(3,np)
real (kind=8), intent(inout)::Fld(6,np)
real (kind=8)              :: S0(-1:1),Xleft,Xright,lambda,dx,Lx,taper, &
                              x_shift, dx_inv,ddx,xp,yp,ku,ampl_tap, &
                              pi=4.d0*DATAN(1.d0)
integer                    :: ip,ix

!f2py intent(in)  :: coord,t,a0,params
!f2py intent(in,out) :: Fld
!f2py intent(hide):: np,nx

lambda = params(1)
Xleft  = params(2)
dx     = params(3)
Lx     = params(4)
taper  = params(5)

ku     = 2.d0*pi/lambda
dx_inv = 1.0d0/dx

Xright= Xleft+nx*dx
x_shift = 0.5*(nx*dx - Lx)

!$omp parallel default(private) shared(coord,Fld,t,a0,params,np,nx, &
!$omp   lambda,Xleft,dx,ku,dx_inv,Xright, Lx,taper,x_shift)
!$omp do schedule(static)
do ip=1,np
  xp = coord(1,ip)
  if ((xp<Xleft+dx) .or. (xp>Xright-dx)) CYCLE
  yp = coord(2,ip)

  S0 = 0.0d0
  ix  = FLOOR((xp-Xleft)*dx_inv +0.5d0)
  ddx = (xp-Xleft)*dx_inv - ix
  S0(-1) = 0.5d0*(0.5d0-ddx)**2
  S0( 0) = 0.75d0-ddx**2
  S0( 1) = 0.5d0*(0.5d0+ddx)**2

  ampl_tap = 1+taper*(coord(1,ip)-x_shift-Xleft-0.5*Lx)/(0.5*Lx)

  Fld(5,ip) = Fld(5,ip) + ampl_tap*SUM(S0*a0(1,ix-1:ix+1))*cosh(ku*yp)
  Fld(4,ip) = Fld(4,ip) + ampl_tap*SUM(S0*a0(2,ix-1:ix+1))*sinh(ku*yp)
enddo
!$omp end do
!$omp end parallel

end subroutine

subroutine Undul_Analytic_Taper(coord,Fld,t,params,np)
implicit none
integer, intent(in) :: np
real (kind=8), intent(in) :: t,params(5)
real (kind=8), intent(in) :: coord(3,np)
real (kind=8), intent(inout) :: Fld(6,np)
real (kind=8)   :: pi=4.d0*DATAN(1.d0), a0, X0,lambda,Lx,ku, ampl,taper
integer         :: ip

!f2py intent(in) :: coord,t,params
!f2py intent(in,out) :: Fld
!f2py intent(hide) :: np

a0     = params(1)
lambda = params(2)
X0     = params(3)
Lx     = params(4)
taper =  params(5)
ku     = 2.d0*pi/lambda

!$omp parallel default(private) shared(coord,Fld,t,params,np,a0, &
!$omp                                  lambda,X0,Lx,taper,ku)
!$omp do schedule(static)
do ip=1,np
  if(coord(1,ip)<=X0 .or. coord(1,ip)>=X0+Lx) then
    ampl = 0.0d0
  elseif(coord(1,ip)>X0 .AND. coord(1,ip)<X0+lambda) then
    ampl = (coord(1,ip)-X0)/lambda
  elseif(coord(1,ip)>X0+Lx-lambda .AND. coord(1,ip)<X0+Lx) then
    ampl = (X0+Lx-coord(1,ip))/lambda
  else
    ampl = 1.d0
  endif

  ampl = ampl*(1+taper*(coord(1,ip)-X0-0.5*Lx)/(0.5*Lx))

  ampl = ampl*a0
  Fld(5,ip) = Fld(5,ip) + ampl*sin(ku*(coord(1,ip)-X0))*cosh(ku*coord(2,ip))
  Fld(4,ip) = Fld(4,ip) + ampl*cos(ku*(coord(1,ip)-X0))*sinh(ku*coord(2,ip))
enddo
!$omp end do
!$omp end parallel

end subroutine

subroutine Undul_Analytic(coord,Fld,t,params,np)
implicit none
integer, intent(in) :: np
real (kind=8), intent(in) :: t,params(4)
real (kind=8), intent(in) :: coord(3,np)
real (kind=8), intent(inout) :: Fld(6,np)
real (kind=8)   :: pi=4.d0*DATAN(1.d0), a0, X0,lambda,Lx,ku, ampl
integer         :: ip

!f2py intent(in) :: coord,t,params
!f2py intent(in,out) :: Fld
!f2py intent(hide) :: np

a0     = params(1)
lambda = params(2)
X0     = params(3)
Lx     = params(4)

ku     = 2.d0*pi/lambda

!$omp parallel default(private) shared(coord,Fld,t,params,np,a0,&
!$omp                                  lambda,X0,Lx,ku)
!$omp do schedule(static)
do ip=1,np
  if(coord(1,ip)<=X0 .or. coord(1,ip)>=X0+Lx) then
    ampl = 0.0d0
  elseif(coord(1,ip)>X0 .AND. coord(1,ip)<X0+lambda) then
    ampl = (coord(1,ip)-X0)/lambda
  elseif(coord(1,ip)>X0+Lx-lambda .AND. coord(1,ip)<X0+Lx) then
    ampl = (X0+Lx-coord(1,ip))/lambda
  else
    ampl = 1.d0
  endif

  ampl = ampl*a0
  Fld(5,ip) = Fld(5,ip) + ampl*sin(ku*(coord(1,ip)-X0))*cosh(ku*coord(2,ip))
  Fld(4,ip) = Fld(4,ip) + ampl*cos(ku*(coord(1,ip)-X0))*sinh(ku*coord(2,ip))
enddo
!$omp end do
!$omp end parallel

end subroutine

subroutine PlaneWave(coord,Fld,t,params,np)
implicit none
integer, intent(in) :: np
real (kind=8), intent(in) :: t,params(7)
real (kind=8), intent(in) :: coord(3,np)
real (kind=8), intent(inout) :: Fld(6,np)
real (kind=8) :: pi=4.d0*DATAN(1.d0),sinth,costh,k0,&
                 theta,ampl,a0,X0,Lx,lambda,ramp,phi0
integer       :: ip

!f2py intent(in) :: coord, t,params
!f2py intent(in,out) :: Fld
!f2py intent(hide) :: np

a0     = params(1)
lambda = params(2)
X0     = params(3)
Lx     = params(4)
ramp   = params(5)
theta  = params(6)
phi0   = params(7)

sinth = sin(theta)
costh = cos(theta)
k0=2.d0*pi/lambda

!$omp parallel default(shared) private(ip,ampl) 
!$omp do schedule(static)
do ip=1,np
  if(coord(1,ip)<=X0 .or. coord(1,ip)>=X0+Lx) then
    ampl = 0.0d0
  elseif(coord(1,ip)>X0 .AND. coord(1,ip)<X0+ramp) then
    ampl = (coord(1,ip)-X0)/ramp
  elseif(coord(1,ip)>X0+Lx-ramp .AND. coord(1,ip)<X0+Lx) then
    ampl = (X0+Lx-coord(1,ip))/ramp
  else
    ampl = 1.d0
  endif

  ampl = a0*ampl*sin(k0*(coord(1,ip)*costh + coord(2,ip)*sinth-t)+phi0)

  Fld(3,ip) = Fld(3,ip) + ampl
  Fld(4,ip) = Fld(4,ip) + ampl*sinth
  Fld(5,ip) = Fld(5,ip) - ampl*costh
enddo
!$omp end do
!$omp end parallel
end subroutine

subroutine GaussBeam(coord, Fld, time, a0, params, np)

implicit none
integer, intent(in) :: np
real (kind=8), intent(in) :: time, params(8), a0
real (kind=8), dimension(3, np), intent(in) :: coord
real (kind=8), dimension(6, np), intent(inout) :: Fld

real (kind=8) :: axis, lambda, E_field, pi=4.d0*DATAN(1.d0), &
                 xp, yp, zp, x0, y0, z0, &
                 Lx2_inv, Ly2_inv, Lz2_inv
integer       :: ip

!f2py intent(in) :: coord, time, params, a0
!f2py intent(in,out) :: Fld
!f2py intent(hide) :: np

lambda = params(1)
axis = params(2)

x0 = params(3)
y0 = params(4)
z0 = params(5)

Lx2_inv = 1.0 / (params(6)*params(6))
Ly2_inv = 1.0 / (params(7)*params(7))
Lz2_inv = 1.0 / (params(8)*params(8))

!$omp parallel default(shared) private(ip, xp, yp, zp, E_field)
!$omp do schedule(static)
do ip=1,np
    xp = coord(1, ip) - x0 - axis*time
    yp = coord(2, ip) - y0
    zp = coord(3, ip) - z0

    E_field = a0 * EXP( -xp*xp*Lx2_inv - yp*yp*Ly2_inv - zp*zp*Lz2_inv ) &
              * sin(2.d0*pi/lambda*xp)

    Fld(3,ip) = Fld(3,ip) + E_field
    Fld(5,ip) = Fld(5,ip) - axis*E_field
enddo
!$omp end do
!$omp end parallel
end subroutine


!######################################################
!############# OLD STUFF ##############################
!######################################################


subroutine Pulse_circ(coord,Fld,time,a0,param_spc,params,num_part)

implicit none
integer, intent(in) :: num_part
real (kind=8), intent(in) :: time,params(7),a0,param_spc(4)
real (kind=8), dimension(3,num_part), intent(in) :: coord
real (kind=8), dimension(6,num_part), intent(inout) :: Fld
real (kind=8), dimension(num_part) :: env, xx, yy,zz
real (kind=8):: pi=3.141592653589793d0, angl,lambda

!f2py intent(in) :: coord, time, params,a0,param_spc
!f2py intent(in,out) :: Fld
!f2py intent(hide) :: num_part

lambda = param_spc(1)

angl = params(4)

xx = (coord(1,:)- params(5))*COS(angl) - (coord(2,:)- params(6))*SIN(angl)
yy = (coord(1,:)- params(5))*SIN(angl) + (coord(2,:)- params(6))*COS(angl)
zz =  coord(3,:)- params(7)

env = a0*EXP(-((xx-time)/params(1))**2-(yy/params(2))**2-(zz/params(3))**2)

!env = env*sin(2*pi/lambda*(xx-time))

Fld(2,:) = Fld(2,:) + env*cos(2.d0*pi/lambda*(xx-time))
Fld(3,:) = Fld(3,:) + env*sin(2.d0*pi/lambda*(xx-time))

Fld(5,:) = Fld(5,:) + Fld(3,:)
Fld(6,:) = Fld(6,:) - Fld(2,:)

end subroutine

subroutine Pulse(coord,Fld,time,a0,param_spc,params,num_part)

implicit none
integer, intent(in) :: num_part
real (kind=8), intent(in) :: time,params(7),a0,param_spc(4)
real (kind=8), dimension(3,num_part), intent(in) :: coord
real (kind=8), dimension(6,num_part), intent(inout) :: Fld
real (kind=8), dimension(num_part) :: env, xx, yy,zz
real (kind=8):: pi=3.141592653589793d0, angl,lambda

!f2py intent(in) :: coord, time, params,a0,param_spc
!f2py intent(in,out) :: Fld
!f2py intent(hide) :: num_part

lambda = param_spc(1)

angl = params(4)

xx = (coord(1,:)- params(5))*COS(angl) - (coord(2,:)- params(6))*SIN(angl)
yy = (coord(1,:)- params(5))*SIN(angl) + (coord(2,:)- params(6))*COS(angl)
zz =  coord(3,:)- params(7)

env = a0*EXP(-((xx-time)/params(1))**2-(yy/params(2))**2-(zz/params(3))**2)

env = env*sin(2.d0*pi/lambda*(xx-time))

Fld(3,:) = Fld(3,:) + env
Fld(4,:) = Fld(4,:) - env*SIN(angl)
Fld(5,:) = Fld(5,:) - env*COS(angl)
end subroutine

