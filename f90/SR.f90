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

subroutine sr_calc_far_tot(spect,coords,momenta_prv,momenta_nxt,wghts,dt,&
                         omega,SinTh,CosTh,SinPh,CosPh,nt,np,nom,nth,nph)
implicit none
integer, intent(in) :: nt,np,nth,nph,nom
real(kind=8), dimension(3,nt,np), intent(in) :: coords,momenta_prv,momenta_nxt
real(kind=8), intent(in) :: wghts(np), dt, omega(nom)
real(kind=8), intent(in) :: SinTh(nth),CosTh(nth),SinPh(nph),CosPh(nph)
real(kind=8), intent(inout) :: spect(nom,nth,nph)

real (kind=8), allocatable    :: spect_loc(:,:,:)
integer          :: ip,it,ith,iph,iom
real (kind=8)    :: coord(3),accel(3),veloc_prv(3),veloc_nxt(3),veloc(3)
real (kind=8)    :: wp,gp_inv
real (kind=8)    :: C1, C2, C41, C42, C43, C3,C2_inv, C2_inv2
real (kind=8)    :: dt_inv,dt2p,pi=4.d0*DATAN(1.d0)
real (kind=8)    :: sin_th,cos_th,sin_ph,cos_ph,omg
complex (kind=8) :: ii=(0.0d0,1.0d0), integral1(nom)
complex (kind=8) :: integral2(nom),integral3(nom)

!f2py intent(in) :: coords,momenta_prv,momenta_nxt,wghts,dt
!f2py intent(in) :: omega,SinTh,CosTh,SinPh,CosPh
!f2py intent(in,out) :: spect
!f2py intent(hide) :: nt,np,nom,nth,nph

dt_inv = 1.0d0/dt
dt2p = dt/(2.0d0*pi)

!$omp parallel default(shared) private(spect_loc,ip,it,ith,iph,iom,coord,  &
!$omp                                  accel,veloc_prv,veloc_nxt,veloc,wp, &
!$omp                                  gp_inv, C1,C2, C3, C41, C42, C43,   &
!$omp                                  C2_inv, C2_inv2, sin_th,cos_th,     &
!$omp                                  sin_ph,cos_ph,omg,integral1,        &
!$omp                                  integral2,integral3)

allocate(spect_loc(nom,nth,nph))
spect_loc = 0.0d0

!$omp do schedule(static)
do ip=1,np
  wp = ABS(wghts(ip))
  do iph=1,nph
    sin_ph = SinPh(iph)
    cos_ph = CosPh(iph)
    do ith=1,nth
      sin_th = SinTh(ith)
      cos_th = CosTh(ith)
      integral1 = 0.0d0
      integral2 = 0.0d0
      integral3 = 0.0d0
      do it=1,nt
        coord = coords(:,it,ip)
        veloc_prv = momenta_prv(:,it,ip)
        veloc_nxt = momenta_nxt(:,it,ip)

        gp_inv = 1.0d0/SQRT(1.0d0+ SUM(veloc_prv*veloc_prv) )
        veloc_prv = veloc_prv*gp_inv
        gp_inv = 1.0d0/SQRT(1.0d0+ SUM(veloc_nxt*veloc_nxt) )
        veloc_nxt = veloc_nxt*gp_inv

        accel = (veloc_nxt - veloc_prv)*dt_inv
        veloc = 0.5d0*(veloc_nxt + veloc_prv)

        C2 = 1.0d0 - (veloc(3)*sin_th*cos_ph + &
          veloc(2)*sin_th*sin_ph + veloc(1)*cos_th)
        C2_inv = 1.0d0/C2
        C2_inv2 = C2_inv*C2_inv

        C1 = accel(3)*sin_th*cos_ph + accel(2)*sin_th*sin_ph + accel(1)*cos_th
        C3 = 2.0d0 * pi * (it*dt-(coord(3)*sin_th*cos_ph &
             + coord(2)*sin_th*sin_ph + coord(1)*cos_th))

        C41 = ( C1*(       cos_th-veloc(1)) - C2*accel(1) )*C2_inv2
        C42 = ( C1*(sin_ph*sin_th-veloc(2)) - C2*accel(2) )*C2_inv2
        C43 = ( C1*(sin_th*cos_ph-veloc(3)) - C2*accel(3) )*C2_inv2

        do iom=1,nom
          omg = omega(iom)
          if (2.0d0 * omg * dt * C2 <= 1.0d0) then
            integral1(iom) = integral1(iom)+C41*dt*exp(ii*omg*C3)
            integral2(iom) = integral2(iom)+C42*dt*exp(ii*omg*C3)
            integral3(iom) = integral3(iom)+C43*dt*exp(ii*omg*C3)
          endif
        enddo
      enddo
      spect_loc(:,ith,iph) = spect_loc(:,ith,iph) &
         + wp*(ABS(integral1)**2+ABS(integral2)**2+ABS(integral3)**2)
    enddo
  enddo
enddo
!$omp end do 

do iph=1,nph
  do ith=1,nth
    do iom=1,nom
      !$omp atomic
      spect(iom,ith,iph) = spect(iom,ith,iph)+spect_loc(iom,ith,iph) 
    enddo
  enddo
enddo

deallocate(spect_loc)
!$omp barrier
!$omp end parallel
end subroutine

subroutine sr_calc_far_comp(spect,coords,momenta_prv,momenta_nxt,wghts, &
                              comp,dt,omega,SinTh,CosTh,SinPh,CosPh,nt,   &
                              np,nom,nth,nph)
implicit none
integer, intent(in) :: comp,nt,np,nth,nph,nom
real(kind=8), dimension(3,nt,np), intent(in) :: coords,momenta_prv,momenta_nxt
real(kind=8), intent(in) :: wghts(np), dt, omega(nom)
real(kind=8), intent(in) :: SinTh(nth),CosTh(nth),SinPh(nph),CosPh(nph)
real(kind=8), intent(inout) :: spect(nom,nth,nph)

real (kind=8), allocatable    :: spect_loc(:,:,:)
integer          :: ip,it,ith,iph,iom
real (kind=8)    :: coord(3),accel(3),veloc_prv(3),veloc_nxt(3),veloc(3)
real (kind=8)    :: wp,gp_inv
real (kind=8)    :: C1, C2, C4, C3,C2_inv, C2_inv2
real (kind=8)    :: dt_inv,dt2p,pi=4.d0*DATAN(1.d0)
real (kind=8)    :: sin_th,cos_th,sin_ph,cos_ph,omg
complex (kind=8) :: ii=(0.0d0,1.0d0), integral(nom)

!f2py intent(in) :: coords,momenta_prv,momenta_nxt,wghts,comp,dt
!f2py intent(in) :: omega,SinTh,CosTh,SinPh,CosPh
!f2py intent(in,out) :: spect
!f2py intent(hide) :: nt,np,nom,nth,nph

dt_inv = 1.0d0/dt
dt2p = dt/(2.0d0*pi)

!$omp parallel default(shared) private(spect_loc,ip,it,ith,iph,iom,coord, &
!$omp                                  accel,veloc_prv,veloc_nxt,veloc,   &
!$omp                                  wp,gp_inv,C1,C2, C4, C3,C2_inv,    &
!$omp                                  C2_inv2,sin_th,cos_th,sin_ph,cos_ph,&
!$omp                                  omg,integral)

allocate(spect_loc(nom,nth,nph))
spect_loc = 0.0d0

!$omp do schedule(static)
do ip=1,np
  wp = ABS(wghts(ip))
  do iph=1,nph
    sin_ph = SinPh(iph)
    cos_ph = CosPh(iph)
    do ith=1,nth
      sin_th = SinTh(ith)
      cos_th = CosTh(ith)
      integral = 0.0d0
      do it=1,nt
        coord = coords(:,it,ip)
        veloc_prv = momenta_prv(:,it,ip)
        veloc_nxt = momenta_nxt(:,it,ip)

        gp_inv = 1.0d0/SQRT(1.0d0+ SUM(veloc_prv*veloc_prv) )
        veloc_prv = veloc_prv*gp_inv
        gp_inv = 1.0d0/SQRT(1.0d0+ SUM(veloc_nxt*veloc_nxt) )
        veloc_nxt = veloc_nxt*gp_inv

        accel = (veloc_nxt - veloc_prv)*dt_inv
        veloc = 0.5d0*(veloc_nxt + veloc_prv)

        C2 = 1.0d0 - (veloc(3)*sin_th*cos_ph + &
          veloc(2)*sin_th*sin_ph + veloc(1)*cos_th)
        C2_inv = 1.0d0/C2
        C2_inv2 = C2_inv*C2_inv

        C1 = accel(3)*sin_th*cos_ph + accel(2)*sin_th*sin_ph + accel(1)*cos_th
        C3 = 2.0d0 * pi * (it*dt - (coord(3)*sin_th*cos_ph &
             + coord(2)*sin_th*sin_ph + coord(1)*cos_th))

        if    (comp==1) then
          C4 = ( C1*(       cos_th-veloc(1)) - C2*accel(1) )*C2_inv2
        elseif(comp==2) then
          C4 = ( C1*(sin_ph*sin_th-veloc(2)) - C2*accel(2) )*C2_inv2
        elseif(comp==3) then
          C4 = ( C1*(sin_th*cos_ph-veloc(3)) - C2*accel(3) )*C2_inv2
        else
          C4 = 0.0
        endif

        do iom=1,nom
          omg = omega(iom)
          if (2.0d0*omg*dt*C2 <= 1.0d0) then
            integral(iom) = integral(iom)+C4*dt*exp(ii*omg*C3)
          endif
        enddo
      enddo
      spect_loc(:,ith,iph) = spect_loc(:,ith,iph)+ wp*ABS(integral)**2
    enddo
  enddo
enddo
!$omp end do 

do iph=1,nph
  do ith=1,nth
    do iom=1,nom
      !$omp atomic
      spect(iom,ith,iph) = spect(iom,ith,iph)+spect_loc(iom,ith,iph) 
    enddo
  enddo
enddo

deallocate(spect_loc)
!$omp barrier
!$omp end parallel
end subroutine

subroutine sr_calc_near_comp(spect, coords, momenta, wghts, &
                              comp, dt, omega, Xgrid, Ygrid, z_scr, nt, &
                              np, nom, nx, ny)
implicit none
integer, intent(in) :: comp, nt, np, nx, ny, nom
real(kind=8), dimension(3,nt,np), intent(in) :: coords,momenta
real(kind=8), intent(in) :: wghts(np), dt, omega(nom)
real(kind=8), intent(in) :: Xgrid(nx), Ygrid(ny), z_scr
real(kind=8), intent(inout) :: spect(nom,nx,ny)

real (kind=8), allocatable    :: spect_loc(:,:,:)
integer          :: ip, it, ix, iy, iom
real (kind=8)    :: coord(3), veloc(3), pointing(3)
real (kind=8)    :: wp, gp_inv, x_scr, y_scr, omg, R_0, R_inv, kotelnikov
real (kind=8)    :: pi2, pi2_inv, dt_inv, dt_2p, pi=4.0d0*DATAN(1.0d0)
complex (kind=8) :: ii=(0.0d0,1.0d0), integral(nom) 
complex (kind=8) :: arg_phase, arg_amp1, arg_amp2,phase_prv

!f2py intent(in) :: coords, momenta, wghts, comp, dt
!f2py intent(in) :: z_scr, omega, Xgrid, Ygrid
!f2py intent(in,out) :: spect
!f2py intent(hide) :: nt, np, nom, nx, ny

pi2 = 2.0d0*pi
pi2_inv = 1.0d0/pi2

dt_inv = 1.0d0/dt
dt_2p = dt*pi2_inv

!$omp parallel default(shared) private(spect_loc,ip,it,ix,iy,iom,wp, &
!$omp                                  x_scr,y_scr,integral,coord, pointing, &
!$omp                                  R_0,R_inv,veloc,gp_inv,arg_phase, &
!$omp                                  phase_prv, arg_amp1,arg_amp2,omg, &
!$omp                                  kotelnikov)

allocate(spect_loc(nom,nx,ny))
spect_loc = 0.0d0

!$omp do schedule(static)
do ip=1,np
  wp = ABS(wghts(ip))
  do ix=1,nx
    x_scr = Xgrid(ix)
    do iy=1,ny
      y_scr =  Ygrid(iy)
      integral = 0.0d0
      phase_prv = 0.0d0
      do it=1,nt
        coord = coords(:,it,ip)
        pointing(1) = z_scr - coord(1)
        pointing(2) = y_scr - coord(2)
        pointing(3) = x_scr - coord(3)
        R_0 = SQRT(SUM(pointing*pointing))
        R_inv = 1.0d0 / R_0
        pointing = pointing * R_inv

        veloc = momenta(:,it,ip)
        gp_inv = 1.0d0 / SQRT(1.0d0+SUM(veloc*veloc))
        veloc = veloc * gp_inv

        arg_phase = ii * pi2 * (it*dt+R_0)
        arg_amp1 = ii * dt * R_inv * (veloc(comp)-pointing(comp))
        arg_amp2 = dt * R_inv * R_inv * pi2_inv * pointing(comp)

        kotelnikov =  0.0d0
!ABS(arg_phase - phase_prv)
!        phase_prv = arg_phase
        !dt * (1.0d0-SUM(veloc*pointing))

        do iom=1,nom
          omg = omega(iom)
          if (kotelnikov*omg < pi2) then
            integral(iom) = integral(iom) + (arg_amp1*omg+arg_amp2) &
                                           *exp(arg_phase*omg)
          endif
        enddo
      enddo
      spect_loc(:,ix,iy) = spect_loc(:,ix,iy) + wp*ABS(integral)**2
    enddo
  enddo
enddo
!$omp end do 

do iy=1,ny
  do ix=1,nx
    do iom=1,nom
      !$omp atomic
      spect(iom,ix,iy) = spect(iom,ix,iy) + spect_loc(iom,ix,iy) 
    enddo
  enddo
enddo

deallocate(spect_loc)
!$omp barrier
!$omp end parallel
end subroutine
