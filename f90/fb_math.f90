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

plan_in = fftw_plan_dft_1d(nkx,Aifft,Afft,FFTW_FORWARD,FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
vec_fb = 0.0d0
shiftX = dcos(leftX*kx) - ii*dsin(leftX*kx)

!$omp parallel shared(vec_fb,vec,In,plan_in,shiftX,nkr,nr,nkO), private(iO,l)
do l=1,3
  do iO = 1,nkO
    !$omp do schedule(static) private(Aifft,Afft,ik,ir)
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

plan_in = fftw_plan_dft_1d(nkx,Aifft,Afft,FFTW_FORWARD,FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
scl_fb = 0.0d0
shiftX = dcos(leftX*kx) - ii*dsin(leftX*kx)

!$omp parallel shared(scl_fb,scl,In,plan_in,shiftX,nkr,nr,nkO), private(iO)
do iO = 1,nkO
  !$omp do schedule(static) private(Aifft,Afft,ik,ir)
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
plan_out = fftw_plan_dft_1d(nkx,Afft,Aifft,FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
shiftX = dcos(leftX*kx) + ii*dsin(leftX*kx)

!$omp parallel shared(vec_fb,vec,Out,plan_out,shiftX,nkr,nr,nkO), private(iO,l)
do l=1,3
  do iO = 1,nkO
    !$omp do schedule(static) private(Aifft,Afft,ik,ir)
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

subroutine fb_eb_out(eb_spc,e_fb,b_fb,leftX,kx,Out,nkx,nr,nkO,nkr)
use, intrinsic :: iso_c_binding
implicit none 
include "fftw3.f03"
integer, intent(in)        :: nkx,nr,nkO,nkr
real (kind=8), intent(in)  :: leftX,kx(nkx),Out(nkr,nr,nkO)
complex(kind=8),intent(in) :: e_fb(nkx,nkr,nkO,3),b_fb(nkx,nkr,nkO,3)
complex(kind=8),intent(inout):: eb_spc(nkx,0:nr,nkO,6) 
integer :: l,ir,ik,iO
complex(kind=8) :: ii=(0.0d0,1.0d0),shiftX(nkx)
complex(C_DOUBLE_COMPLEX) :: Afft(nkx),Aifft(nkx)
type(C_PTR) :: plan_out

!f2py intent(in) :: e_fb,b_fb,leftX,kx,Out
!f2py intent(in,out) :: eb_spc
!f2py intent(hide) :: nkx,nr,nkO,nkr

eb_spc = 0.0d0
plan_out = fftw_plan_dft_1d(nkx,Afft,Aifft,FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
shiftX = dcos(leftX*kx) + ii*dsin(leftX*kx)

!$omp parallel shared(eb_spc,e_fb,b_fb,Out,plan_out,shiftX,nkr,nr,nkO), private(iO,l)
do l=1,3
  do iO = 1,nkO
    !$omp do schedule(static) private(Aifft,Afft,ik,ir)
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
plan_out = fftw_plan_dft_1d(nkx,Afft,Aifft,FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
shiftX = dcos(leftX*kx) + ii*dsin(leftX*kx)

!$omp parallel shared(scl_fb,scl,Out,plan_out,shiftX,nkr,nr,nkO), private(iO)
do iO = 1,nkO
  !$omp do schedule(static) private(Aifft,Afft,ik,ir)
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

plan_in = fftw_plan_dft_1d(nkx,Aifft,Afft,FFTW_FORWARD,FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
plan_out = fftw_plan_dft_1d(nkx,Afft,Aifft,FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_DESTROY_INPUT)

!$omp parallel shared(vec,plan_in,plan_out,shiftX,shiftX_inv,filtr,nkx,nkr,nkO),private(iO,l)
do l=1,3
  do iO=1,nkO
    !$omp do schedule(static) private(Aifft,Afft,ik)
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

subroutine fb_rot(vec_fb_loc,vec_fb,DpS2S,DmS2S,kx,nkx,nkr,nkO,nkr_loc)
implicit none
integer, intent(in) :: nkx,nkr,nkO,nkr_loc
complex(kind=8),intent(in) :: vec_fb(nkx,nkr,0:nkO,3)
real (kind=8), intent(in)  :: DpS2S(nkr,nkr_loc,0:nkO),&
                              DmS2S(nkr,nkr_loc,0:nkO),kx(nkx)
complex(kind=8),intent(inout):: vec_fb_loc(nkx,nkr_loc,0:nkO,3)
integer         :: ik_loc,ik,iO
complex(kind=8) :: sliceKR(nkx), ii=(0.0d0,1.0d0),vec_fb_ext(nkx,nkr,3)

!f2py intent(in) :: vec_fb,DpS2S,DmS2S,kx
!f2py intent(in,out) :: vec_fb_loc
!f2py intent(hide) :: nkx,nkr,nkO,nkr_loc

vec_fb_loc = 0.0d0
vec_fb_ext = 0.0d0

vec_fb_ext(1    ,:,:) = -CONJG(vec_fb(1       ,:,1,:))
vec_fb_ext(2:nkx,:,:) = -CONJG(vec_fb(nkx:2:-1,:,1,:))

!$omp parallel shared(vec_fb,vec_fb_loc,DpS2S,DmS2S,kx,ii,nkr,nkr_loc,nkO,vec_fb_ext), private(iO)
do iO=0,nkO
  !$omp do schedule(static) private(sliceKR,ik,ik_loc)
  do ik_loc=1,nkr_loc
    vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) - ii*kx*vec_fb(:,ik_loc,iO,3)
    vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) + ii*kx*vec_fb(:,ik_loc,iO,2)

    if (iO<nkO) then
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR - DpS2S(ik,ik_loc,iO)*(ii*vec_fb(:,ik,iO+1,2) - vec_fb(:,ik,iO+1,3))
      enddo
      vec_fb_loc(:,ik_loc,iO,1) = vec_fb_loc(:,ik_loc,iO,1) + sliceKR

      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR + DpS2S(ik,ik_loc,iO)*vec_fb(:,ik,iO+1,1)
      enddo
      vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) + ii*sliceKR
      vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) -    sliceKR
    endif

    if (iO>0) then
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR - DmS2S(ik,ik_loc,iO)*(ii*vec_fb(:,ik,iO-1,2)+vec_fb(:,ik,iO-1,3))
      enddo
      vec_fb_loc(:,ik_loc,iO,1) = vec_fb_loc(:,ik_loc,iO,1) + sliceKR
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*vec_fb(:,ik,iO-1,1)
      enddo
      vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) + ii*sliceKR
      vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) +    sliceKR
    else
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR - DmS2S(ik,ik_loc,iO)*(ii*vec_fb_ext(:,ik,2)+vec_fb_ext(:,ik,3))
      enddo
      vec_fb_loc(:,ik_loc,iO,1) = vec_fb_loc(:,ik_loc,iO,1) + sliceKR
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*vec_fb_ext(:,ik,1)
      enddo
      vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) + ii*sliceKR
      vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) +    sliceKR
    endif
  enddo
  !$omp end do 
enddo
!$omp end parallel
end subroutine

subroutine fb_grad(vec_fb_loc,scl_fb,DpS2S,DmS2S,kx,nkx,nkr,nkO,nkr_loc)
implicit none
integer, intent(in) :: nkx,nkr,nkO,nkr_loc
complex(kind=8),intent(in) :: scl_fb(nkx,nkr,0:nkO)
real (kind=8), intent(in)  :: DpS2S(nkr,nkr_loc,0:nkO),&
                              DmS2S(nkr,nkr_loc,0:nkO),kx(nkx)
complex(kind=8),intent(inout):: vec_fb_loc(nkx,nkr_loc,0:nkO,3)
integer         :: ik_loc,ik,iO
complex(kind=8) :: sliceKR(nkx), ii=(0.0d0,1.0d0),scl_fb_ext(nkx,nkr)

!f2py intent(in) :: scl_fb,DpS2S,DmS2S,kx
!f2py intent(in,out) :: vec_fb_loc
!f2py intent(hide) :: nkx,nkr,nkO,nkr_loc

vec_fb_loc = 0.0d0

scl_fb_ext(1    ,:) = -CONJG(scl_fb(1       ,:,1))
scl_fb_ext(2:nkx,:) = -CONJG(scl_fb(nkx:2:-1,:,1))

!$omp parallel shared(scl_fb,vec_fb_loc,DpS2S,DmS2S,kx,ii,nkr,nkr_loc,nkO,scl_fb_ext), private(iO)
do iO=0,nkO
  !$omp do schedule(static) private(sliceKR,ik,ik_loc)
  do ik_loc=1,nkr_loc
    vec_fb_loc(:,ik_loc,iO,1) = vec_fb_loc(:,ik_loc,iO,1) + ii*kx*scl_fb(:,ik_loc,iO)

    if (iO>0) then
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*scl_fb(:,ik,iO-1)
      enddo
      vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) - sliceKR
      vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) + ii*sliceKR
    else
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*scl_fb_ext(:,ik)
      enddo
      vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) - sliceKR
      vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) + ii*sliceKR
    endif
    if (iO<nkO) then
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR + DpS2S(ik,ik_loc,iO)*scl_fb(:,ik,iO+1)
      enddo
      vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) +    sliceKR
      vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) + ii*sliceKR
    endif
  enddo
  !$omp end do
enddo
!$omp end parallel
end subroutine

subroutine fb_div(scl_fb_loc,vec_fb,DpS2S,DmS2S,kx,nkx,nkr,nkO,nkr_loc)
implicit none
integer, intent(in) :: nkx,nkr,nkO,nkr_loc
complex(kind=8),intent(in) :: vec_fb(nkx,nkr,0:nkO,3)
real (kind=8), intent(in)  :: DpS2S(nkr,nkr_loc,0:nkO),&
                              DmS2S(nkr,nkr_loc,0:nkO),kx(nkx)
complex(kind=8),intent(inout):: scl_fb_loc(nkx,nkr_loc,0:nkO)
integer         :: ik_loc,ik,iO
complex(kind=8) :: sliceKR(nkx), ii=(0.0d0,1.0d0),vec_fb_extYZ(nkx,nkr,2)

!f2py intent(in) :: vec_fb,DpS2S,DmS2S,kx
!f2py intent(in,out) :: scl_fb_loc
!f2py intent(hide) :: nkx,nkr,nkO,nkr_loc

scl_fb_loc = 0.0d0
vec_fb_extYZ(1    ,:,:) = -CONJG(vec_fb(1       ,:,1,2:3))
vec_fb_extYZ(2:nkx,:,:) = -CONJG(vec_fb(nkx:2:-1,:,1,2:3))

!$omp parallel shared(vec_fb,scl_fb_loc,DpS2S,DmS2S,kx,ii,nkr,nkr_loc,nkO,vec_fb_extYZ), private(iO)
do iO=0,nkO
  !$omp do schedule(static) private(sliceKR,ik,ik_loc)
  do ik_loc=1,nkr_loc
    scl_fb_loc(:,ik_loc,iO) = scl_fb_loc(:,ik_loc,iO) + ii*kx*vec_fb(:,ik_loc,iO,1)

    sliceKR = 0.0d0
    if (iO>0) then
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*(ii*vec_fb(:,ik,iO-1,3)-vec_fb(:,ik,iO-1,2))
      enddo
    else
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*(ii*vec_fb_extYZ(:,ik,2)-vec_fb_extYZ(:,ik,1))
      enddo
    endif
    if (iO<nkO) then
      do ik=1,nkr
        sliceKR = sliceKR +  DpS2S(ik,ik_loc,iO)*(ii*vec_fb(:,ik,iO+1,3)+vec_fb(:,ik,iO+1,2))
      enddo
    endif
    scl_fb_loc(:,ik_loc,iO) = scl_fb_loc(:,ik_loc,iO) + sliceKR
  enddo
  !$omp end do
enddo
!$omp end parallel
end subroutine

subroutine fb_grad_env(vec_fb_loc,scl_fb,DpS2S,DmS2S,kx,nkx,nkr,nkO,nkr_loc)
implicit none
integer, intent(in) :: nkx,nkr,nkO,nkr_loc
complex(kind=8),intent(in) :: scl_fb(nkx,nkr,-nkO:nkO)
real (kind=8), intent(in)  :: DpS2S(nkr,nkr_loc,-nkO:nkO),&
                              DmS2S(nkr,nkr_loc,-nkO:nkO),kx(nkx)
complex(kind=8),intent(inout):: vec_fb_loc(nkx,nkr_loc,-nkO:nkO,3)
integer         :: ik_loc,ik,iO
complex(kind=8) :: sliceKR(nkx), ii=(0.0d0,1.0d0)

!f2py intent(in) :: scl_fb,DpS2S,DmS2S,kx
!f2py intent(in,out) :: vec_fb_loc
!f2py intent(hide) :: nkx,nkr,nkO,nkr_loc

vec_fb_loc = 0.0d0

!$omp parallel shared(scl_fb,vec_fb_loc,DpS2S,DmS2S,kx,ii,nkr,nkr_loc,nkO), private(iO)
do iO=-nkO,nkO
  !$omp do schedule(static) private(sliceKR,ik,ik_loc)
  do ik_loc=1,nkr_loc
    vec_fb_loc(:,ik_loc,iO,1) = vec_fb_loc(:,ik_loc,iO,1) + ii*kx*scl_fb(:,ik_loc,iO)

    if (iO>-nkO) then
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*scl_fb(:,ik,iO-1)
      enddo
      vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) -    sliceKR
      vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) + ii*sliceKR
    endif
    if (iO<nkO) then
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR + DpS2S(ik,ik_loc,iO)*scl_fb(:,ik,iO+1)
      enddo
      vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) +    sliceKR
      vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) + ii*sliceKR
    endif
  enddo
  !$omp end do
enddo
!$omp end parallel
end subroutine

subroutine fb_div_env(scl_fb_loc,vec_fb,DpS2S,DmS2S,kx,nkx,nkr,nkO,nkr_loc)
implicit none
integer, intent(in) :: nkx,nkr,nkO,nkr_loc
complex(kind=8),intent(in) :: vec_fb(nkx,nkr,-nkO:nkO,3)
real (kind=8), intent(in)  :: DpS2S(nkr,nkr_loc,-nkO:nkO),&
                              DmS2S(nkr,nkr_loc,-nkO:nkO),kx(nkx)
complex(kind=8),intent(inout):: scl_fb_loc(nkx,nkr_loc,-nkO:nkO)
integer         :: ik_loc,ik,iO
complex(kind=8) :: sliceKR(nkx), ii=(0.0d0,1.0d0)

!f2py intent(in) :: vec_fb,DpS2S,DmS2S,kx
!f2py intent(in,out) :: scl_fb_loc
!f2py intent(hide) :: nkx,nkr,nkO,nkr_loc

scl_fb_loc = 0.0d0

!$omp parallel shared(vec_fb,scl_fb_loc,DpS2S,DmS2S,kx,ii,nkr,nkr_loc,nkO), private(iO)
do iO=-nkO,nkO
  !$omp do schedule(static) private(sliceKR,ik,ik_loc)
  do ik_loc=1,nkr_loc
    scl_fb_loc(:,ik_loc,iO) = scl_fb_loc(:,ik_loc,iO) + ii*kx*vec_fb(:,ik_loc,iO,1)

    sliceKR = 0.0d0
    if (iO>-nkO) then      
      do ik=1,nkr
        sliceKR = sliceKR - DmS2S(ik,ik_loc,iO)*(vec_fb(:,ik,iO-1,2) - ii*vec_fb(:,ik,iO-1,3))
      enddo
    endif
    if (iO<nkO) then
      do ik=1,nkr
        sliceKR = sliceKR + DpS2S(ik,ik_loc,iO)*(vec_fb(:,ik,iO+1,2)+ii*vec_fb(:,ik,iO+1,3))
      enddo
    endif
    scl_fb_loc(:,ik_loc,iO) = scl_fb_loc(:,ik_loc,iO) + sliceKR
  enddo
  !$omp end do
enddo
!$omp end parallel
end subroutine

subroutine fb_rot_env(vec_fb_loc,vec_fb,DpS2S,DmS2S,kx,nkx,nkr,nkO,nkr_loc)
implicit none
integer, intent(in) :: nkx,nkr,nkO,nkr_loc
complex(kind=8),intent(in) :: vec_fb(nkx,nkr,-nkO:nkO,3)
real (kind=8), intent(in)  :: DpS2S(nkr,nkr_loc,-nkO:nkO),&
                              DmS2S(nkr,nkr_loc,-nkO:nkO),kx(nkx)
complex(kind=8),intent(inout):: vec_fb_loc(nkx,nkr_loc,-nkO:nkO,3)
integer         :: ik_loc,ik,iO
complex(kind=8) :: sliceKR(nkx), ii=(0.0d0,1.0d0)

!f2py intent(in) :: vec_fb,DpS2S,DmS2S,kx
!f2py intent(in,out) :: vec_fb_loc
!f2py intent(hide) :: nkx,nkr,nkO,nkr_loc

vec_fb_loc = 0.0d0

!$omp parallel shared(vec_fb,vec_fb_loc,DpS2S,DmS2S,kx,ii,nkr,nkr_loc,nkO), private(iO)
do iO=-nkO,nkO
  !$omp do schedule(static) private(sliceKR,ik,ik_loc)
  do ik_loc=1,nkr_loc
    vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) - ii*kx*vec_fb(:,ik_loc,iO,3)
    vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) + ii*kx*vec_fb(:,ik_loc,iO,2)

    if (iO<nkO) then
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR - DpS2S(ik,ik_loc,iO)*(ii*vec_fb(:,ik,iO+1,2)-vec_fb(:,ik,iO+1,3))
      enddo
      vec_fb_loc(:,ik_loc,iO,1) = vec_fb_loc(:,ik_loc,iO,1) + sliceKR
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR + DpS2S(ik,ik_loc,iO)*vec_fb(:,ik,iO+1,1)
      enddo
      vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) + ii*sliceKR
      vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) -    sliceKR
    endif
    if (iO>-nkO) then
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR - DmS2S(ik,ik_loc,iO)*(ii*vec_fb(:,ik,iO-1,2)+vec_fb(:,ik,iO-1,3))
      enddo
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*vec_fb(:,ik,iO-1,1)
      enddo
      vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) + ii*sliceKR
      vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) +    sliceKR
    endif
  enddo
  !$omp end do
enddo
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
