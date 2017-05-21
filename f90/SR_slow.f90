subroutine sr_calc_incoh(spect,coords,momenta_prv,momenta_nxt,wghts,comp,Step,dt,&
                         omega,SinTh,CosTh,SinPh,CosPh,np,nom,nth,nph)

implicit none
integer, intent(in) :: comp,Step, np,nth,nph,nom
real(kind=8), dimension(3,np), intent(in) :: coords,momenta_prv,momenta_nxt
real(kind=8), intent(in) :: wghts(np), dt, omega(nom)
real(kind=8), intent(in) :: SinTh(nth),CosTh(nth),SinPh(nph),CosPh(nph)
complex(kind=8), intent(inout) :: spect(nom,nth,nph)

integer          :: ip, ith,iph,iom
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
      !$omp veloc_prv,veloc_nxt,accel,veloc,C1,C2,C4,C3,C2_inv,C2_inv2) REDUCTION(+:integral)
      do ip=1,np
        coord = coords(:,ip)
        wp = ABS(wghts(ip))
!        wp = SQRT(ABS(wghts(ip)))

        veloc_prv = momenta_prv(:,ip)
        veloc_nxt = momenta_nxt(:,ip)

        accel = dt_inv*(veloc_nxt - veloc_prv)
        veloc = 0.5d0 *(veloc_nxt + veloc_prv)
        gp_inv = 1.0d0/SQRT(1.0d0+ SUM(veloc*veloc) )

        veloc = veloc*gp_inv
        accel = accel*gp_inv

        C2 = 1.0 - (veloc(3)*sin_th*cos_ph + &
          veloc(2)*sin_th*sin_ph + veloc(1)*cos_th)
        C2_inv = 1.0d0/C2
        C2_inv2 = C2_inv*C2_inv

        C1 = accel(3)*sin_th*cos_ph + accel(2)*sin_th*sin_ph + accel(1)*cos_th
        C3 = 2.0*pi*(Step*dt - (coord(3)*sin_th*cos_ph + coord(2)*sin_th*sin_ph + coord(1)*cos_th))
        C4 = ( C1*(sin_th*cos_ph-veloc(3)) - C2*accel(3) )*C2_inv2

        integral = integral + C4*dt2p*exp(ii*omg*C3)*wp
      enddo
      !$omp end parallel do 
      spect(iom,ith,iph) = spect(iom,ith,iph)+ integral      
    enddo
  enddo
enddo
end subroutine
