subroutine fb_vec_in(vec_fb,vec,leftX,kx,In,nkx,nr,nkO,nkr)
use, intrinsic :: iso_c_binding
implicit none 
include "fftw3.f03"

integer, intent(in) :: nkx,nr,nkO,nkr
complex(kind=8),intent(in) :: vec(nkx,0:nr,nkO,3)
real (kind=8), intent(in)  :: leftX,kx(nkx),In(nr,nkr,nkO)
complex(kind=8),intent(inout):: vec_fb(nkx,nkr,nkO,3)

integer         :: l,ir,ik,iO
complex(kind=8) :: ii=(0.0d0,1.0d0),shiftX(nkx)
complex(C_DOUBLE_COMPLEX) :: Afft(nkx),Aifft(nkx)
type(C_PTR) :: plan_in

!f2py intent(in) :: vec,leftX,kx,In
!f2py intent(in,out) :: vec_fb
!f2py intent(hide) :: nkx,nr,nkO,nkr

plan_in = fftw_plan_dft_1d(nkx,Aifft,Afft,FFTW_FORWARD,&
                           FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
vec_fb = 0.0d0
shiftX = dcos(leftX*kx) - ii*dsin(leftX*kx)

!$omp parallel default(shared), private(iO,l,Aifft,Afft,ik,ir)
do l=1,3
  do iO = 1,nkO
    !$omp do schedule(static)
    do ik=1,nkr
      Aifft = 0.0
      do ir=1,nr
        Aifft = Aifft + In(ir,ik,iO)*vec(:,ir,iO,l)
      enddo
      call fftw_execute_dft(plan_in, Aifft,Afft)
      vec_fb(:,ik,iO,l) = vec_fb(:,ik,iO,l) + Afft*shiftX
    enddo
    !$omp end do 
  enddo
enddo
!$omp end parallel 
call fftw_destroy_plan(plan_in)
end subroutine

subroutine fb_scl_in(scl_fb,scl,leftX,kx,In,nkx,nr,nkO,nkr)
use, intrinsic :: iso_c_binding
implicit none 
include "fftw3.f03"
integer, intent(in) :: nkx,nr,nkO,nkr
real (kind=8), intent(in)  :: leftX,kx(nkx),In(nr,nkr,nkO)
complex(kind=8),intent(in) :: scl(nkx,0:nr,nkO)
complex(kind=8),intent(inout):: scl_fb(nkx,nkr,nkO)
integer         :: ir,ik,iO
complex(kind=8) :: ii=(0.0d0,1.0d0),shiftX(nkx)
complex(C_DOUBLE_COMPLEX) :: Afft(nkx),Aifft(nkx)
type(C_PTR) :: plan_in

!f2py intent(in) :: scl,leftX,kx,In
!f2py intent(in,out) :: scl_fb
!f2py intent(hide) :: nkx,nr,nkO,nkr

plan_in = fftw_plan_dft_1d(nkx,Aifft,Afft,FFTW_FORWARD,&
                           FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
scl_fb = 0.0d0
shiftX = dcos(leftX*kx) - ii*dsin(leftX*kx)

!$omp parallel default(shared), private(iO,Aifft,Afft,ik,ir)
do iO = 1,nkO
  !$omp do schedule(static)
  do ik=1,nkr
    Aifft = 0.0
    do ir=1,nr
      Aifft = Aifft + In(ir,ik,iO)*scl(:,ir,iO)
    enddo
    call fftw_execute_dft(plan_in, Aifft,Afft)
    scl_fb(:,ik,iO) = scl_fb(:,ik,iO) +Afft*shiftX
  enddo
  !$omp end do 
enddo
!$omp end parallel 
call fftw_destroy_plan(plan_in)
end subroutine

subroutine fb_vec_out(vec,vec_fb,leftX,kx,Out,nkx,nr,nkO,nkr)
use, intrinsic :: iso_c_binding
implicit none 
include "fftw3.f03"
integer, intent(in)        :: nkx,nr,nkO,nkr
real (kind=8), intent(in)  :: leftX,kx(nkx),Out(nkr,nr,nkO)
complex(kind=8),intent(in) :: vec_fb(nkx,nkr,nkO,3)
complex(kind=8),intent(out):: vec(nkx,0:nr,nkO,3) 
integer :: l,ir,ik,iO
complex(kind=8) :: ii=(0.0d0,1.0d0),shiftX(nkx)
complex(C_DOUBLE_COMPLEX) :: Afft(nkx),Aifft(nkx)
type(C_PTR) :: plan_out

!f2py intent(in) :: vec_fb,leftX,kx,Out
!f2py intent(out) :: vec
!f2py intent(hide) :: nkx,nr,nkO,nkr

vec = 0.0d0
plan_out = fftw_plan_dft_1d(nkx,Afft,Aifft,FFTW_BACKWARD,&
                            FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
shiftX = dcos(leftX*kx) + ii*dsin(leftX*kx)

!$omp parallel default(shared), private(iO,l,Aifft,Afft,ik,ir)
do l=1,3
  do iO = 1,nkO
    !$omp do schedule(static)
    do ir=1,nr
      Afft = 0.0
      do ik=1,nkr
        Afft = Afft + Out(ik,ir,iO)*vec_fb(:,ik,iO,l)
      enddo
      Afft = Afft*shiftX
      call fftw_execute_dft(plan_out, Afft,Aifft)
      vec(:,ir,iO,l) = vec(:,ir,iO,l) + Aifft
    enddo
    !$omp end do 
  enddo
enddo
!$omp end parallel
call fftw_destroy_plan(plan_out)
end subroutine

subroutine fb_scl_out(scl,scl_fb,leftX,kx,Out,nkx,nr,nkO,nkr)
use, intrinsic :: iso_c_binding
implicit none 
include "fftw3.f03"
integer, intent(in)        :: nkx,nr,nkO,nkr
real (kind=8), intent(in)  :: leftX,kx(nkx),Out(nkr,nr,nkO)
complex(kind=8),intent(in) :: scl_fb(nkx,nkr,nkO)
complex(kind=8),intent(out):: scl(nkx,0:nr,nkO)
integer :: ir,ik,iO
complex(kind=8) :: ii=(0.0d0,1.0d0),shiftX(nkx)
complex(C_DOUBLE_COMPLEX) :: Afft(nkx),Aifft(nkx)
type(C_PTR) :: plan_out

!f2py intent(in) :: scl_fb,leftX,kx,Out
!f2py intent(out) :: scl
!f2py intent(hide) :: nkx,nr,nkO,nkr

scl = 0.0d0
plan_out = fftw_plan_dft_1d(nkx,Afft,Aifft,FFTW_BACKWARD,&
                            FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
shiftX = dcos(leftX*kx) + ii*dsin(leftX*kx)

!$omp parallel default(shared), private(iO,Aifft,Afft,ik,ir)
do iO = 1,nkO
  !$omp do schedule(static)
  do ir=1,nr
    Afft = 0.0
    do ik=1,nkr
      Afft = Afft + Out(ik,ir,iO)*scl_fb(:,ik,iO)
    enddo
    Afft = Afft*shiftX
    call fftw_execute_dft(plan_out, Afft,Aifft)
    scl(:,ir,iO) = scl(:,ir,iO) + Aifft
  enddo
  !$omp end do 
enddo
!$omp end parallel
call fftw_destroy_plan(plan_out)
end subroutine

subroutine fb_eb_out(eb_spc,e_fb,b_fb,leftX,kx,Out,nkx,nr,nkO,nkr)
use, intrinsic :: iso_c_binding
implicit none 
include "fftw3.f03"
integer, intent(in)        :: nkx,nr,nkO,nkr
real (kind=8), intent(in)  :: leftX,kx(nkx),Out(nkr,nr,nkO)
complex(kind=8),intent(in) :: e_fb(nkx,nkr,nkO,6),b_fb(nkx,nkr,nkO,3)
complex(kind=8),intent(inout):: eb_spc(nkx,0:nr,nkO,6) 
integer :: l,ir,ik,iO
complex(kind=8) :: ii=(0.0d0,1.0d0),shiftX(nkx)
complex(C_DOUBLE_COMPLEX) :: Afft(nkx),Aifft(nkx)
type(C_PTR) :: plan_out

!f2py intent(in) :: e_fb,b_fb,leftX,kx,Out
!f2py intent(in,out) :: eb_spc
!f2py intent(hide) :: nkx,nr,nkO,nkr,nkO

eb_spc = 0.0d0
plan_out = fftw_plan_dft_1d(nkx,Afft,Aifft,FFTW_BACKWARD,&
                            FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
shiftX = dcos(leftX*kx) + ii*dsin(leftX*kx)
!$omp parallel default(shared), private(iO,l,Aifft,Afft,ik,ir)
do l=1,3
  do iO = 1,nkO
    !$omp do schedule(static)
    do ir=1,nr
      Afft = 0.0
      do ik=1,nkr
        Afft = Afft + Out(ik,ir,iO)*e_fb(:,ik,iO,l)
      enddo
      Afft = Afft*shiftX
      call fftw_execute_dft(plan_out, Afft,Aifft)
      eb_spc(:,ir,iO,l) = eb_spc(:,ir,iO,l) + Aifft
      Afft = 0.0
      do ik=1,nkr
        Afft = Afft + Out(ik,ir,iO)*b_fb(:,ik,iO,l)
      enddo
      Afft = Afft*shiftX
      call fftw_execute_dft(plan_out, Afft,Aifft)
      eb_spc(:,ir,iO,l+3) = eb_spc(:,ir,iO,l+3) + Aifft
    enddo
    !$omp end do 
  enddo
enddo
!$omp end parallel
call fftw_destroy_plan(plan_out)
end subroutine

subroutine fb_filtr(vec,leftX,kx,filtr,nkx,nkr,nkO)
use, intrinsic :: iso_c_binding
implicit none
include "fftw3.f03"
integer, intent(in) :: nkx,nkr,nkO
complex(kind=8),intent(inout) :: vec(nkx,nkr,nkO,3)
real (kind=8), intent(in)  :: leftX,kx(nkx),filtr(nkx)
integer         :: l,ik,iO
complex(kind=8) :: ii=(0.0d0,1.0d0),shiftX(nkx),shiftX_inv(nkx)
complex(C_DOUBLE_COMPLEX) :: Afft(nkx),Aifft(nkx)
type(C_PTR) :: plan_in, plan_out

!f2py intent(in) :: leftX,kx,filtr
!f2py intent(in,out) :: vec
!f2py intent(hide) :: nkx,nkr,nkO

shiftX = dcos(leftX*kx) + ii*dsin(leftX*kx)
shiftX_inv = 1.0d0/(DBLE(nkx)*shiftX)

plan_in = fftw_plan_dft_1d(nkx,Aifft,Afft,FFTW_FORWARD,&
                           FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
plan_out = fftw_plan_dft_1d(nkx,Afft,Aifft,FFTW_BACKWARD,&
                            FFTW_ESTIMATE+FFTW_DESTROY_INPUT)

!$omp parallel default(shared), private(iO,l,Aifft,Afft,ik)
do l=1,3
  do iO=1,nkO
    !$omp do schedule(static)
    do ik=1,nkr
      Afft = vec(:,ik,iO,l)*shiftX
      call dfftw_execute_dft(plan_out,Afft,Aifft)
      call dfftw_execute_dft(plan_in,Aifft*filtr,Afft)
      vec(:,ik,iO,l) = Afft*shiftX_inv
    enddo
    !$omp end do
  enddo
enddo
!$omp end parallel
call fftw_destroy_plan(plan_out)
call fftw_destroy_plan(plan_in)
end subroutine
