subroutine intens_profO(PWR_RO,Fld,NO,nx,nr,nkO)
implicit none
integer, intent(in)        :: nx,nr,nkO,NO
complex (kind=8),intent(in):: Fld(0:nx,0:nr,-nkO:nkO,3)
real (kind=8), intent(out) :: PWR_RO(NO,1:nr)
integer         :: ix,ikO,ir,iO,l
real(kind=8)    :: pi=4.d0*DATAN(1.d0)
complex(kind=8) :: ii=(0.0d0,1.0d0),phase_p,phase_m,phaseO(-nkO:nkO),Os(NO)
!f2py intent(in)  :: Fld,NO
!f2py intent(out) :: PWR_RO
!f2py intent(hide):: nx,nr,nkO

PWR_RO  = 0.0d0

Os(1) = 1.0
phase_p = DCOS(2.0d0*pi/DBLE(NO-1)) + ii*DSIN(2.0d0*pi/DBLE(NO-1))
do iO = 2,NO
  Os(iO) = Os(iO-1)*phase_p
enddo

do iO = 1,NO
  phase_p = Os(iO)
  phase_m = 1./Os(iO)
  phaseO = 1.0d0
  if (nkO>0) then
    do ikO = 1,nkO
      phaseO(ikO) = phaseO(ikO-1)*phase_p
      phaseO(-ikO) = phaseO(-ikO+1)*phase_m
    enddo
  endif

  do l=1,3
    do ir=1,nr
      do ix=0,nx
        PWR_RO(iO,ir) = PWR_RO(iO,ir) + ABS(SUM(phaseO*Fld(ix,ir,:,l)))**2
      enddo
    enddo
  enddo
enddo
end subroutine


subroutine myfftgramm(gramm,signal,wind,sigm,ns)

use, intrinsic :: iso_c_binding
include "fftw3.f03"

type(C_PTR) :: plan
integer, intent(in)        :: ns,wind
real (kind=8), intent(in)  :: signal(ns),sigm
real (kind=8), intent(out) :: gramm(wind,ns)
integer                    :: i,wind_2
real (kind=8)              :: signal_ext(ns+wind),filtr(wind),sigm2
complex(C_DOUBLE_COMPLEX) :: Afft(wind),Aifft(wind)

!f2py intent(in) :: signal,wind,sigm
!f2py intent(out) :: gramm
!f2py intent(hide) :: ns

signal_ext = 0.
Afft = 0.0
Aifft=0.0
gramm = 0.
wind_2 = wind/2
sigm2 = 1.d0/sigm**2
plan = fftw_plan_dft_1d(wind,Aifft,Afft,FFTW_FORWARD,FFTW_ESTIMATE)
signal_ext(wind_2:wind_2+ns-1) = signal

do i=1,wind
  filtr(i) = EXP(-0.5d0*DBLE((i-wind_2)**2)/DBLE(wind**2)*sigm2)
enddo

do i=1,ns
  Aifft = signal_ext(i:i+wind-1)
  call fftw_execute_dft(plan,Aifft,Afft)  
  gramm(:,i) = ABS(Afft)*filtr
enddo

end subroutine

subroutine get_amplitude1D(traject,period, amplitude, num_nodes)
implicit none
integer, intent(in)        :: num_nodes,period
real (kind=8), intent(in)  :: traject(num_nodes)
real (kind=8), intent(out) :: amplitude(num_nodes)
integer                    :: node,half_period
real (kind=8), dimension(num_nodes+period):: traj_extend

!f2py intent(in) :: traject,period
!f2py intent(out) :: amplitude
!f2py intent(hide) :: num_nodes

traj_extend = 0.
amplitude = 0.
half_period = period/2

traj_extend(half_period:half_period+num_nodes-1) = traject(:)
do node=1,num_nodes
  amplitude(node) = SQRT(2.*SUM(traj_extend(node:node+period)**2)/period)
enddo
end subroutine


subroutine get_amplitude1D2(traject,period, amplitude, num_nodes)
implicit none
integer, intent(in)        :: num_nodes,period
real (kind=8), intent(in)  :: traject(num_nodes)
real (kind=8), intent(out) :: amplitude(num_nodes)
integer                    :: node,half_period
real (kind=8), dimension(num_nodes+period):: traj_extend

!f2py intent(in) :: traject,period
!f2py intent(out) :: amplitude
!f2py intent(hide) :: num_nodes

traj_extend = 0.
amplitude = 0.
half_period = period/2

traj_extend(half_period:half_period+num_nodes-1) = traject(:)
do node=1,num_nodes
  amplitude(node) = MAXVAL(ABS(traj_extend(node:node+period)))
enddo
end subroutine


subroutine get_smooth1D(traject,period, avrgd, num_nodes)
implicit none
integer, intent(in)        :: num_nodes,period
real (kind=8), intent(in)  :: traject(num_nodes)
real (kind=8), intent(out) :: avrgd(num_nodes)
integer                    :: node,half_period
real (kind=8), dimension(num_nodes+period):: traj_extend

!f2py intent(in) :: traject,period
!f2py intent(out) :: avrgd
!f2py intent(hide) :: num_nodes

traj_extend = 0.
avrgd = 0.
half_period = INT(0.5*period)

traj_extend(half_period:half_period+num_nodes) = traject(:)

traj_extend(1:half_period)=traject(half_period:1:-1)
traj_extend(num_nodes+half_period:period+num_nodes)=traject(num_nodes:num_nodes-half_period:-1)

do node=1,num_nodes
  avrgd(node) = SUM(traj_extend(node:node+period))/(period+1.0)
enddo

end subroutine

subroutine get_strength(traject,period, strength, num_nodes,num_part)

implicit none
integer, intent(in)        :: num_nodes,num_part,period
real (kind=8), intent(in)  :: traject(0:num_nodes,num_part)
real (kind=8), intent(out) :: strength(num_part)
integer                    :: np, step
real (kind=8)              :: c1, c2, c3, c4
!real (kind=8), dimension(0:num_nodes+period):: time
real (kind=8), dimension(0:num_nodes+period):: traj_extend
real (kind=8), dimension(0:num_nodes):: y

!f2py intent(in) :: traject,period
!f2py intent(out) :: strength
!f2py intent(hide) :: num_nodes,num_part

strength = 0.

do np =1,num_part
  traj_extend(period/2:num_nodes+period/2) = traject(:,np)
  c1 = sum(traject(0:period,np))/period
  c2 = traject(0,np)
  c3 = sum(traject(num_nodes-period:num_nodes,np))/period
  c4 = traject(num_nodes,np)
  do step=0,period/2
    traj_extend(step) = c1 + (c2-c1)*EXP(-6.*(step-period/2)**2/period**2)
  enddo
  do step=num_nodes+period/2, num_nodes+period
    traj_extend(step) = c3 + (c4-c3)*EXP(-6.*(step-num_nodes-period/2)**2/period**2)
  enddo
  y = 0.
  do step = period/2,num_nodes+period/2
    y(step-period/2) = sum(traj_extend(step-period/2:step+period/2))/period
  enddo
  strength(np) =  sqrt( 2.*(sum((traject(:,np)-y)**2)/num_nodes - (sum(traject(:,np)-y)/num_nodes)**2))
enddo

end subroutine

SUBROUTINE DENSITY_2x(x,y,wght,grid,bins_x,bins_y,dens, N_part)

IMPLICIT NONE
INTEGER, INTENT(IN)        :: bins_x, bins_y, N_part
REAL(kind=8), INTENT(IN), dimension(N_part) :: x, y,wght
REAL(kind=8), INTENT(IN)   :: grid(4)
REAL(kind=8), INTENT(OUT)  :: dens(-2:bins_x+2,-2:bins_y+2)
REAL(kind=8)               :: dx, dy, dlt_xg_inv, dlt_yg_inv, xg(-2:bins_x+2), yg(-2:bins_y+2),x_max,y_max,&
                              origx, origy, dlt_xg, dlt_yg
INTEGER                    :: jp, kx, ky, j
REAL(kind=8)               :: S0(-3:3,1:2)

!f2py intent(in) :: x,y,wght, grid, bins_x, bins_y
!f2py intent(out) :: dens
!f2py intent(hide) :: N_part

dens = 0.d0
origx = grid(1)
origy = grid(3)
dlt_xg = (grid(2) - grid(1))/bins_x
dlt_yg = (grid(4) - grid(3))/bins_y
dlt_xg_inv = 1.0d0/dlt_xg
dlt_yg_inv= 1.0d0/dlt_yg
x_max=origx+dlt_xg*bins_x
y_max=origy+dlt_yg*bins_y

DO j  =  -2,bins_x+2
   xg(j)  =  dlt_xg*REAL(j) + origx
ENDDO
DO j  =  -2,bins_y+2
   yg(j)  =  dlt_yg*REAL(j) + origy
ENDDO

DO jp = 1, N_part
 if((x(jp)>=origx).AND.(x(jp)<=x_max).AND.(y(jp)>=origy).AND.(y(jp)<=y_max)) then
   S0(:,:) = 0.d0
   kx     = FLOOR((x(jp)-origx  )*dlt_xg_inv+0.5d0)
   ky     = FLOOR((y(jp)-origy  )*dlt_yg_inv+0.5d0)
   dx     =       (x(jp)-xg(kx ))*dlt_xg_inv
   dy     =       (y(jp)-yg(ky ))*dlt_yg_inv
   S0(-1,1) = (0.25d0-0.5d0*dx+(dx**3)/3.d0)
   S0( 0,1) = (0.5d0-dabs(dx**3)/3.d0)
   S0( 1,1) = (0.25d0+0.5d0*dx-(dx**3)/3.d0)
   if(dx >= 0.0) then
       S0(2,1)  = dabs(dx**3)/3.d0
   else
       S0(-2,1) = dabs(dx**3)/3.d0
   endif
   S0(-1,2) = (0.25d0-0.5d0*dy+(dy**3)/3.d0)
   S0( 0,2) = (0.5d0-dabs(dy**3)/3.d0)
   S0( 1,2) = (0.25d0+0.5d0*dy-(dy**3)/3.d0)
   if(dy >= 0.0) then
      S0(2,2)  = dabs(dy**3)/3.d0
   else
      S0(-2,2) = dabs(dy**3)/3.d0
   endif
   DO j = -2, 2
      dens(kx-2:kx+2,ky+j) = dens(kx-2:kx+2,ky+j) + wght(jp)*S0(-2:2,1)*S0(j,2)
   ENDDO
 endif
ENDDO

dens = dens * dlt_xg_inv * dlt_yg_inv

END SUBROUTINE
