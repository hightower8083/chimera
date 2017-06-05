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

subroutine maxwell_push_with_spchrg(EG_fb,j_fb,grad_rho_n_fb,&
                                    grad_rho_np1_fb,&
                                    C1,C2,nkx,nkr,nkO)
implicit none
integer, intent(in) :: nkx,nkr,nkO
complex(kind=8), intent(inout) ::           EG_fb(nkx,nkr,nkO,6)
complex(kind=8), intent(in)    ::            j_fb(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    ::   grad_rho_n_fb(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    :: grad_rho_np1_fb(nkx,nkr,nkO,3)
real (kind=8), intent(in)      :: C1(nkx,nkr,nkO,5)
real (kind=8), intent(in)      :: C2(nkx,nkr,nkO,5)
integer         :: l,ikr,iO
complex(kind=8) :: slice_kx(nkx)

!f2py intent(in) :: j_fb,grad_rho_n_fb,grad_rho_np1_fb,C1,C2
!f2py intent(in,out) :: EG_fb
!f2py intent(hide) :: nkx,nkO,nkr

!$omp parallel default(shared) private(iO,l,ikr,slice_kx)
do iO=1,nkO
  !$omp do schedule(static)
  do ikr=1,nkr
    do l=1,3

      slice_kx = C1(:,ikr,iO,1)*EG_fb(:,ikr,iO,l) &
               + C1(:,ikr,iO,2)*EG_fb(:,ikr,iO,3+l) &
               + C1(:,ikr,iO,3)*j_fb(:,ikr,iO,l) & 
               + C1(:,ikr,iO,4)*grad_rho_n_fb(:,ikr,iO,l) & 
               + C1(:,ikr,iO,5)*grad_rho_np1_fb(:,ikr,iO,l)

      EG_fb(:,ikr,iO,l+3) = C2(:,ikr,iO,1)*EG_fb(:,ikr,iO,l) &
                          + C2(:,ikr,iO,2)*EG_fb(:,ikr,iO,3+l) &
                          + C2(:,ikr,iO,3)*j_fb(:,ikr,iO,l) &
                          + C2(:,ikr,iO,4)*grad_rho_n_fb(:,ikr,iO,l) &
                          + C2(:,ikr,iO,5)*grad_rho_np1_fb(:,ikr,iO,l)

      EG_fb(:,ikr,iO,l) = slice_kx
    enddo
  enddo
  !$omp end do
enddo
!$omp end parallel
end subroutine

subroutine maxwell_push_wo_spchrg(EG_fb,j_fb,C1,C2,nkx,nkr,nkO)
implicit none
integer, intent(in)            :: nkx,nkr,nkO
complex(kind=8), intent(inout) :: EG_fb(nkx,nkr,nkO,6)
complex(kind=8), intent(in)    :: j_fb(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    :: C1(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    :: C2(nkx,nkr,nkO,3)

integer         :: l,ikr,iO
complex(kind=8) :: slice_kx(nkx)

!f2py intent(in) :: j_fb,C1,C2
!f2py intent(in,out) :: EG_fb
!f2py intent(hide) :: nkx,nkO,nkr

!$omp parallel default(shared) private(l,ikr,iO,slice_kx)
do iO=1,nkO
  !$omp do schedule(static)
  do ikr=1,nkr
    do l=1,3
      slice_kx(:) = C1(:,ikr,iO,1)*EG_fb(:,ikr,iO,l) &
                  + C1(:,ikr,iO,2)*EG_fb(:,ikr,iO,3+l) &
                  + C1(:,ikr,iO,3)*j_fb(:,ikr,iO,l)

      EG_fb(:,ikr,iO,l+3) = C2(:,ikr,iO,1)*EG_fb(:,ikr,iO,l) &
                          + C2(:,ikr,iO,2)*EG_fb(:,ikr,iO,3+l) &
                          + C2(:,ikr,iO,3)*j_fb(:,ikr,iO,l)

      EG_fb(:,ikr,iO,l) = slice_kx
    enddo
  enddo
  !$omp end do
enddo
!$omp end parallel
end subroutine

subroutine maxwell_init_push(EG_fb,j_fb,grad_rho_n_fb,C1,C2,nkx,nkr,nkO)
implicit none
integer, intent(in) :: nkx,nkr,nkO
complex(kind=8), intent(inout) ::         EG_fb(nkx,nkr,nkO,6)
complex(kind=8), intent(in)    ::          j_fb(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    :: grad_rho_n_fb(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    :: C1(nkx,nkr,nkO,2)
complex(kind=8), intent(in)    :: C2(nkx,nkr,nkO,2)
integer                        :: l,ikr,iO

!f2py intent(in) :: j_fb,grad_rho_n_fb,C1,C2
!f2py intent(in,out) :: EG_fb
!f2py intent(hide) :: nkx,nkO,nkr

!$omp parallel default(shared) private(iO,l,ikr)
do iO=1,nkO
  !$omp do schedule(static)
  do ikr=1,nkr
    do l=1,3
      EG_fb(:,ikr,iO,l)  = EG_fb(:,ikr,iO,l) &
                         + C1(:,ikr,iO,1)*j_fb(:,ikr,iO,l) &
                         + C1(:,ikr,iO,2)*grad_rho_n_fb(:,ikr,iO,l)

      EG_fb(:,ikr,iO,l+3)= EG_fb(:,ikr,iO,l+3) & 
                         + C2(:,ikr,iO,1)*j_fb(:,ikr,iO,l) &
                         + C2(:,ikr,iO,2)*grad_rho_n_fb(:,ikr,iO,l)
    enddo
  enddo
  !$omp end do
enddo
!$omp end parallel
end subroutine

subroutine poiss_corr(j_fb,grad_div_j_fb,&
                      grad_rho_n_fb,grad_rho_np1_fb,&
                      dt_inv,w2_inv,nkx,nkr,nkO)
implicit none
integer, intent(in)            :: nkx,nkr,nkO
complex(kind=8), intent(inout) ::            j_fb(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    ::   grad_div_j_fb(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    ::   grad_rho_n_fb(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    :: grad_rho_np1_fb(nkx,nkr,nkO,3)
real(kind=8),    intent(in)    :: dt_inv, w2_inv(nkx,nkr,nkO)
integer         :: l,ikr,iO

!f2py intent(in) :: grad_div_j_fb,grad_rho_n_fb
!f2py intent(in) :: grad_rho_np1_fb,dt_inv,w2_inv
!f2py intent(in,out) :: j_fb
!f2py intent(hide) :: nkx,nkO,nkr

!$omp parallel default(shared) private(iO,l,ikr)
do iO=1,nkO
  do l=1,3
    !$omp do schedule(static)
    do ikr=1,nkr
      j_fb(:,ikr,iO,l) = j_fb(:,ikr,iO,l) & 
                       + ( grad_div_j_fb(:,ikr,iO,l) &
                           + (grad_rho_np1_fb(:,ikr,iO,l) &
                              -grad_rho_n_fb(:,ikr,iO,l) &
                             )*dt_inv & 
                         )*w2_inv(:,ikr,iO) 
    enddo
    !$omp end do
  enddo
enddo
!$omp end parallel
end subroutine

subroutine poiss_corr_stat(j_fb,grad_div_j_fb,&
                           grad_rho_n_fb,DT,w2_inv,&
                           nkx,nkr,nkO)
implicit none
integer, intent(in)            :: nkx,nkr,nkO
complex(kind=8), intent(inout) ::            j_fb(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    ::   grad_div_j_fb(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    ::   grad_rho_n_fb(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    :: DT(nkx)
real(kind=8),    intent(in)    :: w2_inv(nkx,nkr,nkO)
integer         :: l,ikr,iO

!f2py intent(in) :: grad_div_j_fb,grad_rho_n_fb
!f2py intent(in) :: DT,w2_inv
!f2py intent(in,out) ::j_fb
!f2py intent(hide) :: nkx,nkO,nkr

!$omp parallel default(shared) private(iO,l,ikr)
do iO=1,nkO
  do l=1,3
    !$omp do schedule(static)
    do ikr=1,nkr
      j_fb(:,ikr,iO,l) = j_fb(:,ikr,iO,l) & 
                       + ( grad_div_j_fb(:,ikr,iO,l)   &
                           + DT*grad_rho_n_fb(:,ikr,iO,l) & 
                          )*w2_inv(:,ikr,iO) 
    enddo
    !$omp end do
  enddo
enddo
!$omp end parallel
end subroutine

subroutine field_drift(EG_fb,kx,beta0,dt,nkx,nkr,nkO)
implicit none
integer, intent(in)            :: nkx,nkr,nkO
complex(kind=8), intent(inout) :: EG_fb(nkx,nkr,nkO,6)
real(kind=8),    intent(in)    :: kx(nkx),beta0,dt

complex(kind=8) :: ii=(0.0d0,1.0d0), fact,PROPGTR(nkx)
integer         :: l,ikr,iO

!f2py intent(in) :: kx,beta0,dt
!f2py intent(in,out) ::EG_fb
!f2py intent(hide) :: nkx,nkO,nkr

fact = -0.5*ii*dt*beta0
PROPGTR = EXP(fact*kx)

!$omp parallel default(shared) private(iO,l,ikr)
do l=1,6
  do iO=1,nkO
    !$omp do schedule(static)
    do ikr=1,nkr
      EG_fb(:,ikr,iO,l) = EG_fb(:,ikr,iO,l)*PROPGTR
    enddo
    !$omp end do
  enddo
enddo
!$omp end parallel
end subroutine

subroutine omp_mult_vec(vec_fb,A,nkx,nkr,nkO)
implicit none
integer, intent(in)            :: nkx,nkr,nkO
complex(kind=8), intent(inout) :: vec_fb(nkx,nkr,nkO,3)
real(kind=8),    intent(in)    :: A(nkx,nkr,nkO)
integer         :: l,ikr,iO

!f2py intent(in) :: A
!f2py intent(in,out) :: vec_fb
!f2py intent(hide) :: nkx,nkO,nkr

!$omp parallel default(shared) private(iO,l,ikr)
do l=1,3
  do iO=1,nkO
    !$omp do schedule(static)
    do ikr=1,nkr
      vec_fb(:,ikr,iO,l) = vec_fb(:,ikr,iO,l)*A(:,ikr,iO)
    enddo
    !$omp end do
  enddo
enddo
!$omp end parallel
end subroutine

subroutine omp_mult_scl(scl_fb,A,nkx,nkr,nkO)
implicit none
integer, intent(in)            :: nkx,nkr,nkO
complex(kind=8), intent(inout) :: scl_fb(nkx,nkr,nkO)
real(kind=8),    intent(in)    :: A(nkx,nkr,nkO)
integer         :: ikr,iO

!f2py intent(in) :: A
!f2py intent(in,out) :: scl_fb
!f2py intent(hide) :: nkx,nkO,nkr

!$omp parallel default(shared) private(iO,ikr)
do iO=1,nkO
  !$omp do schedule(static)
  do ikr=1,nkr
    scl_fb(:,ikr,iO) = scl_fb(:,ikr,iO)*A(:,ikr,iO)
  enddo
  !$omp end do
enddo
!$omp end parallel
end subroutine

subroutine omp_add_vec(vec_fb,A,nkx,nkr,nkO)
implicit none
integer, intent(in)            :: nkx,nkr,nkO
complex(kind=8), intent(inout) :: vec_fb(nkx,nkr,nkO,3)
complex(kind=8), intent(in)    :: A(nkx,nkr,nkO,3)
integer         :: l,ikr,iO

!f2py intent(in) :: A
!f2py intent(in,out) :: vec_fb
!f2py intent(hide) :: nkx,nkO,nkr

!$omp parallel default(shared) private(iO,l,ikr)
do l=1,3
  do iO=1,nkO
    !$omp do schedule(static)
    do ikr=1,nkr
      vec_fb(:,ikr,iO,l) = vec_fb(:,ikr,iO,l) + A(:,ikr,iO,l)
    enddo
    !$omp end do
  enddo
enddo
!$omp end parallel
end subroutine

subroutine omp_add_scl(scl_fb,A,nkx,nkr,nkO)
implicit none
integer, intent(in)            :: nkx,nkr,nkO
complex(kind=8), intent(inout) :: scl_fb(nkx,nkr,nkO)
complex(kind=8), intent(in)    :: A(nkx,nkr,nkO)
integer         :: ikr,iO

!f2py intent(in) :: A
!f2py intent(in,out) :: scl_fb
!f2py intent(hide) :: nkx,nkO,nkr

!$omp parallel default(shared) private(iO,ikr)
do iO=1,nkO
  !$omp do schedule(static)
  do ikr=1,nkr
    scl_fb(:,ikr,iO) = scl_fb(:,ikr,iO) + A(:,ikr,iO)
  enddo
  !$omp end do
enddo
!$omp end parallel
end subroutine


