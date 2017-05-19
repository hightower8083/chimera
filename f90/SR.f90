subroutine sr_calc_incoh(spect,coords,momenta_prv,momenta_nxt,wghts,comp,Step,dt,&
                         omega,SinTh,CosTh,SinPh,CosPh,np,nom,nth,nph)

implicit none
integer, intent(in) :: comp,Step, np,nth,nph,nom
real(kind=8), dimension(3,np), intent(in) :: coords,momenta_prv,momenta_nxt
real(kind=8), intent(in) :: wghts(np), dt, omega(nom)
real(kind=8), intent(in) :: SinTh(nth),CosTh(nth),SinPh(nph),CosPh(nph)
complex(kind=8), intent(inout) :: spect(nom,nth,nph)

integer          :: ip, ith,iph,iom
real (kind=8)    :: wps(np),accels(3,np),velocs(3,np)
real (kind=8)    :: coord(3),accel(3),veloc_prv(3),veloc_nxt(3),veloc(3),wp,gp_inv
real (kind=8)    :: dt_inv,dt2p,C1, C2, C4, C3,C2_inv, C2_inv2,pi=4.d0*DATAN(1.d0)
real (kind=8)    :: sin_th,cos_th,sin_ph,cos_ph,omg
complex (kind=8) :: ii=(0.0d0,1.0d0), integral

!f2py intent(in) :: coords,momenta_prv,momenta_nxt,wghts,comp,Step,dt
!f2py intent(in) :: omega,SinTh,CosTh,SinPh,CosPh
!f2py intent(in,out) :: spect
!f2py intent(hide) :: np,nom,nth,nph

dt_inv = 1.0d0/dt
dt2p = dt/(2.0*pi)

wps = 0.0
accels = 0.0
velocs = 0.0

!$omp parallel do default(shared) private(veloc_prv,gp_inv,veloc_nxt,ip) schedule(static)
do ip=1,np
  wps(ip) = SQRT(ABS(wghts(ip)))
  veloc_prv = momenta_prv(:,ip)
  gp_inv = 1.0d0/SQRT(1.0d0+ SUM(veloc_prv*veloc_prv) )
  veloc_prv = veloc_prv*gp_inv

  veloc_nxt = momenta_nxt(:,ip)
  gp_inv = 1.0d0/SQRT(1.0d0+ SUM(veloc_nxt*veloc_nxt) )
  veloc_nxt = veloc_nxt*gp_inv

  accels(:,ip) = (veloc_nxt - veloc_prv)*dt_inv
  velocs(:,ip) = 0.5d0*(veloc_nxt + veloc_prv)
enddo
!$omp end parallel do

do iph=1,nph
  sin_ph = SinPh(iph)
  cos_ph = CosPh(iph)
  do ith=1,nth
    sin_th = SinTh(ith)
    cos_th = CosTh(ith)
    do iom=1,nom
      omg = omega(iom)

      integral = 0.0d0

      !$omp parallel do schedule(static) default(shared) private(ip,coord,wp,&
      !$omp accel,veloc,C1,C2,C4,C3,C2_inv,C2_inv2) REDUCTION(+:integral)
      do ip=1,np
        coord = coords(:,ip)
        wp = wps(ip)
        accel = accels(:,ip)
        veloc = velocs(:,ip)

        C2 = 1.0 - (veloc(3)*sin_th*cos_ph + &
          veloc(2)*sin_th*sin_ph + veloc(1)*cos_th)
        C2_inv = 1.0d0/C2
        C2_inv2 = C2_inv*C2_inv

        C1 = accel(3)*sin_th*cos_ph + accel(2)*sin_th*sin_ph + accel(1)*cos_th
        C3 = 2.0*pi*(Step*dt - (coord(3)*sin_th*cos_ph + coord(2)*sin_th*sin_ph + coord(1)*cos_th))

        if(comp==0) then
          C4 = ( C1*(       cos_th-veloc(1)) - C2*accel(1) )*C2_inv2
        elseif(comp==1) then
          C4 = ( C1*(sin_ph*sin_th-veloc(2)) - C2*accel(2) )*C2_inv2
        elseif(comp==2) then
          C4 = ( C1*(sin_th*cos_ph-veloc(3)) - C2*accel(3) )*C2_inv2
        else
          C4 = 0.0
        endif

!        C4 = ( C1*(sin_th*cos_ph-veloc(3)) - C2*accel(3) )*C2_inv2

!        if (omg <= 0.5d0*dt_inv*C2_inv) then
!          integral = integral + C4*dt2p*exp(ii*omg*C3)*wp
!        endif
        integral = integral + C4*dt2p*exp(ii*omg*C3)*wp
      enddo
      !$omp end parallel do 
      spect(iom,ith,iph) = spect(iom,ith,iph)+ integral      
    enddo
  enddo
enddo
end subroutine


subroutine sr_calc_incoh_track(spect,coords,momenta_prv,momenta_nxt,wghts,dt,&
                         omega,SinTh,CosTh,SinPh,CosPh,nt,np,nom,nth,nph)
implicit none
integer, intent(in) :: nt,np,nth,nph,nom
real(kind=8), dimension(3,nt,np), intent(in) :: coords,momenta_prv,momenta_nxt
real(kind=8), intent(in) :: wghts(np), dt, omega(nom)
real(kind=8), intent(in) :: SinTh(nth),CosTh(nth),SinPh(nph),CosPh(nph)
real(kind=8), intent(inout) :: spect(nom,nth,nph)

!real (kind=8), allocatable    :: spect_loc(nom,nth,nph)
real (kind=8), allocatable    :: spect_loc(:,:,:)
integer          :: ip,it,ith,iph,iom
real (kind=8)    :: coord(3),accel(3),veloc_prv(3),veloc_nxt(3),veloc(3),wp,gp_inv
real (kind=8)    :: C1, C2, C4, C3,C2_inv, C2_inv2
real (kind=8)    :: dt_inv,dt2p,pi=4.d0*DATAN(1.d0)
real (kind=8)    :: sin_th,cos_th,sin_ph,cos_ph,omg
complex (kind=8) :: ii=(0.0d0,1.0d0), integral(nom)

!f2py intent(in) :: coords,momenta_prv,momenta_nxt,wghts,dt
!f2py intent(in) :: omega,SinTh,CosTh,SinPh,CosPh
!f2py intent(in,out) :: spect
!f2py intent(hide) :: nt,np,nom,nth,nph

dt_inv = 1.0d0/dt
dt2p = dt/(2.0*pi)

!$omp parallel default(shared) private(spect_loc,ip,it,ith,iph,iom,coord,accel, &
!$omp                                  veloc_prv,veloc_nxt,veloc,wp,gp_inv, &
!$omp                                  C1,C2, C4, C3,C2_inv, C2_inv2, &
!$omp                                  sin_th,cos_th,sin_ph,cos_ph,omg,integral)

allocate(spect_loc(nom,nth,nph))
spect_loc = 0.0

!$omp do schedule(static)
do ip=1,np
  wp = SQRT(ABS(wghts(ip)))
  do iph=1,nph
    sin_ph = SinPh(iph)
    cos_ph = CosPh(iph)
    do ith=1,nth
      sin_th = SinTh(ith)
      cos_th = CosTh(ith)
      integral = 0.0
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

        C2 = 1.0 - (veloc(3)*sin_th*cos_ph + &
          veloc(2)*sin_th*sin_ph + veloc(1)*cos_th)
        C2_inv = 1.0d0/C2
        C2_inv2 = C2_inv*C2_inv

        C1 = accel(3)*sin_th*cos_ph + accel(2)*sin_th*sin_ph + accel(1)*cos_th
        C3 = 2.0*pi*(it*dt - (coord(3)*sin_th*cos_ph + coord(2)*sin_th*sin_ph + coord(1)*cos_th))

        C4 = ( C1*(sin_th*cos_ph-veloc(3)) - C2*accel(3) )*C2_inv2

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
