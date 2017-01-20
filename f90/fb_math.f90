subroutine fb_rot(vec_fb_loc,vec_fb,DpS2S,DmS2S,kx,nkx,nkr,nkO,nkr_loc)
implicit none
integer, intent(in) :: nkx,nkr,nkO,nkr_loc
complex(kind=8),intent(in) :: vec_fb(nkx,nkr,0:nkO,3)
real (kind=8), intent(in)  :: DpS2S(nkr,nkr_loc,0:nkO+1)
real (kind=8), intent(in)  :: DmS2S(nkr,nkr_loc,0:nkO+1),kx(nkx)
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

!$omp parallel default(shared), private(iO,sliceKR,ik,ik_loc)
do iO=0,nkO
  !$omp do schedule(static)
  do ik_loc=1,nkr_loc
    vec_fb_loc(:,ik_loc,iO,2) = vec_fb_loc(:,ik_loc,iO,2) &
                              - ii*kx*vec_fb(:,ik_loc,iO,3)
    vec_fb_loc(:,ik_loc,iO,3) = vec_fb_loc(:,ik_loc,iO,3) &
                              + ii*kx*vec_fb(:,ik_loc,iO,2)

    if (iO<nkO) then
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR - DpS2S(ik,ik_loc,iO)*(ii*vec_fb(:,ik,iO+1,2) &
                                                  - vec_fb(:,ik,iO+1,3))
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
        sliceKR = sliceKR - DmS2S(ik,ik_loc,iO)*(ii*vec_fb(:,ik,iO-1,2) & 
                                                 +vec_fb(:,ik,iO-1,3))
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
        sliceKR = sliceKR - DmS2S(ik,ik_loc,iO)*(ii*vec_fb_ext(:,ik,2) &
                                                 +vec_fb_ext(:,ik,3))
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
real (kind=8), intent(in)  :: DpS2S(nkr,nkr_loc,0:nkO+1),&
                              DmS2S(nkr,nkr_loc,0:nkO+1),kx(nkx)
complex(kind=8),intent(inout):: vec_fb_loc(nkx,nkr_loc,0:nkO,3)
integer         :: ik_loc,ik,iO
complex(kind=8) :: sliceKR(nkx), ii=(0.0d0,1.0d0),scl_fb_ext(nkx,nkr)

!f2py intent(in) :: scl_fb,DpS2S,DmS2S,kx
!f2py intent(in,out) :: vec_fb_loc
!f2py intent(hide) :: nkx,nkr,nkO,nkr_loc

vec_fb_loc = 0.0d0

scl_fb_ext(1    ,:) = -CONJG(scl_fb(1       ,:,1))
scl_fb_ext(2:nkx,:) = -CONJG(scl_fb(nkx:2:-1,:,1))

!$omp parallel default(shared), private(iO,sliceKR,ik,ik_loc)
do iO=0,nkO
  !$omp do schedule(static)
  do ik_loc=1,nkr_loc
    vec_fb_loc(:,ik_loc,iO,1) = vec_fb_loc(:,ik_loc,iO,1) &
                              + ii*kx*scl_fb(:,ik_loc,iO)

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
real (kind=8), intent(in)  :: DpS2S(nkr,nkr_loc,0:nkO+1),&
                              DmS2S(nkr,nkr_loc,0:nkO+1),kx(nkx)
complex(kind=8),intent(inout):: scl_fb_loc(nkx,nkr_loc,0:nkO)
integer         :: ik_loc,ik,iO
complex(kind=8) :: sliceKR(nkx), ii=(0.0d0,1.0d0),vec_fb_extYZ(nkx,nkr,2)

!f2py intent(in) :: vec_fb,DpS2S,DmS2S,kx
!f2py intent(in,out) :: scl_fb_loc
!f2py intent(hide) :: nkx,nkr,nkO,nkr_loc

scl_fb_loc = 0.0d0
vec_fb_extYZ(1    ,:,:) = -CONJG(vec_fb(1       ,:,1,2:3))
vec_fb_extYZ(2:nkx,:,:) = -CONJG(vec_fb(nkx:2:-1,:,1,2:3))

!$omp parallel default(shared), private(iO,sliceKR,ik,ik_loc)
do iO=0,nkO
  !$omp do schedule(static)
  do ik_loc=1,nkr_loc
    scl_fb_loc(:,ik_loc,iO) = scl_fb_loc(:,ik_loc,iO) & 
                            + ii*kx*vec_fb(:,ik_loc,iO,1)

    sliceKR = 0.0d0
    if (iO>0) then
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*(ii*vec_fb(:,ik,iO-1,3) & 
                                                 -vec_fb(:,ik,iO-1,2))
      enddo
    else
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*(ii*vec_fb_extYZ(:,ik,2) &
                                                 -vec_fb_extYZ(:,ik,1))
      enddo
    endif
    if (iO<nkO) then
      do ik=1,nkr
        sliceKR = sliceKR +  DpS2S(ik,ik_loc,iO)*(ii*vec_fb(:,ik,iO+1,3) &
                                                  +vec_fb(:,ik,iO+1,2))
      enddo
    endif
    scl_fb_loc(:,ik_loc,iO) = scl_fb_loc(:,ik_loc,iO) + sliceKR
  enddo
  !$omp end do
enddo
!$omp end parallel
end subroutine

subroutine fb_graddiv(vec_fb,DpS2S,DmS2S,kx,nkx,nkr,nkO,nkr_loc)
implicit none
integer, intent(in) :: nkx,nkr,nkO,nkr_loc
complex(kind=8),intent(inout) :: vec_fb(nkx,nkr,0:nkO,3)
real (kind=8), intent(in)  :: DpS2S(nkr,nkr_loc,0:nkO+1),&
                              DmS2S(nkr,nkr_loc,0:nkO+1),kx(nkx)
complex(kind=8) :: scl_fb_loc(nkx,nkr_loc,0:nkO+1)
integer         :: ik_loc,ik,iO
complex(kind=8) :: sliceKR(nkx), ii=(0.0d0,1.0d0),aux_comp(nkx,nkr,2)

!f2py intent(in) :: DpS2S,DmS2S,kx
!f2py intent(in,out) :: vec_fb
!f2py intent(hide) :: nkx,nkr,nkO,nkr_loc

scl_fb_loc = 0.0d0
aux_comp  = 0.0d0
if (nkO>0) then
  aux_comp(1    ,:,:) = -CONJG(vec_fb(1       ,:,1,2:3))
  aux_comp(2:nkx,:,:) = -CONJG(vec_fb(nkx:2:-1,:,1,2:3))
endif

!$omp parallel default(shared), private(iO,sliceKR,ik,ik_loc)
do iO=0,nkO+1
  !$omp do schedule(static)
  do ik_loc=1,nkr_loc
    if (iO<nkO+1) then
      scl_fb_loc(:,ik_loc,iO) = scl_fb_loc(:,ik_loc,iO)&
                              + ii*kx*vec_fb(:,ik_loc,iO,1)
    endif
    sliceKR = 0.0d0
    if (iO>0) then
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*(ii*vec_fb(:,ik,iO-1,3) &
                                                 -vec_fb(:,ik,iO-1,2))
      enddo
    else
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*(ii*aux_comp(:,ik,2) &
                                                 -aux_comp(:,ik,1))
      enddo
    endif
    if (iO<nkO) then
      do ik=1,nkr
        sliceKR = sliceKR +  DpS2S(ik,ik_loc,iO)*(ii*vec_fb(:,ik,iO+1,3) &
                                                  +vec_fb(:,ik,iO+1,2))
      enddo
    endif
    scl_fb_loc(:,ik_loc,iO) = scl_fb_loc(:,ik_loc,iO) + sliceKR
  enddo
  !$omp end do
enddo
!$omp end parallel

vec_fb = 0.0d0
aux_comp = 0.0d0

if (nkO>0) then
  aux_comp(1    ,:,1) = -CONJG(scl_fb_loc(1       ,:,1))
  aux_comp(2:nkx,:,1) = -CONJG(scl_fb_loc(nkx:2:-1,:,1))
endif

!$omp parallel default(shared), private(iO,sliceKR,ik,ik_loc)
do iO=0,nkO
  !$omp do schedule(static)
  do ik_loc=1,nkr_loc
    vec_fb(:,ik_loc,iO,1) = vec_fb(:,ik_loc,iO,1) &
                          + ii*kx*scl_fb_loc(:,ik_loc,iO)
    if (iO>0) then
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*scl_fb_loc(:,ik,iO-1)
      enddo
      vec_fb(:,ik_loc,iO,2) = vec_fb(:,ik_loc,iO,2) - sliceKR
      vec_fb(:,ik_loc,iO,3) = vec_fb(:,ik_loc,iO,3) + ii*sliceKR
    else
      sliceKR = 0.0d0
      do ik=1,nkr
        sliceKR = sliceKR + DmS2S(ik,ik_loc,iO)*aux_comp(:,ik,1)
      enddo
      vec_fb(:,ik_loc,iO,2) = vec_fb(:,ik_loc,iO,2) - sliceKR
      vec_fb(:,ik_loc,iO,3) = vec_fb(:,ik_loc,iO,3) + ii*sliceKR
    endif
    sliceKR = 0.0d0
    do ik=1,nkr
      sliceKR = sliceKR + DpS2S(ik,ik_loc,iO)*scl_fb_loc(:,ik,iO+1)
    enddo
    vec_fb(:,ik_loc,iO,2) = vec_fb(:,ik_loc,iO,2) +    sliceKR
    vec_fb(:,ik_loc,iO,3) = vec_fb(:,ik_loc,iO,3) + ii*sliceKR
  enddo
  !$omp end do
enddo
!$omp end parallel
end subroutine
