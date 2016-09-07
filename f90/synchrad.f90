subroutine spectrum_comp_incoh_wgh(traject,comp,dt, phi_min,phi_max,NstepPhi,theta_min,theta_max,NstepTheta,&
    omega_min,omega_max,NstepOmega,Nout,out_step, spect,energy,num_nodes,num_part)

implicit none
integer, intent(in) :: num_nodes,num_part,comp, Nout,out_step
integer, intent(in) :: NstepOmega, NstepPhi,NstepTheta
integer             :: i,j,np,k,tt, outindx, outstep
real (kind=8), dimension(0:num_nodes,8,num_part), intent(in) :: traject
real (kind=8) , intent(in) :: dt, phi_min,phi_max, omega_min,omega_max, theta_min,theta_max
real (kind=8), dimension(1:NstepOmega+1,1:NstepTheta+1, 1:NstepPhi,1:Nout), intent(out) :: spect
real (kind=8), dimension(0:num_nodes/out_step) , intent(out) :: energy
real (kind=8) :: d_phi, d_o, d_theta ,phi, theta, omega, pi = 3.141592653589793, &
                 SinPhi, CosPhi, SinTheta, CosTheta
real (kind=8), dimension(0:num_nodes,num_part) :: velocX, velocY, velocZ, coordX, coordY, & 
                                                  coordZ, accelX, accelY, accelZ
real (kind=8), dimension(num_part) :: wgh
real (kind=8)    :: C1, C2, C4, C3, dV, dt2p
complex (kind=8) :: ii=(0.0, 1.0), integral, integral1, integral2, integral3

!f2py intent(in) :: traject,comp,dt,phi_min,phi_max,NstepPhi,omega_min,omega_max,NstepOmega
!f2py intent(in) :: theta_min, theta_max, NstepTheta, Nout,out_step
!f2py intent(out) :: spect,energy
!f2py intent(hide) :: num_nodes,num_part

dt2p = (dt/2.0/pi)

outstep = INT(CEILING(float(num_nodes)/float(Nout)))

phi = phi_min
theta = theta_min
omega = omega_min

if(NstepTheta>0) then
  d_theta = (theta_max-theta_min)/NstepTheta
else
  d_theta = 1.0
endif
if(NstepPhi>0) then
  d_phi = (phi_max-phi_min)/NstepPhi
else
  d_phi = 1.
endif
if(NstepOmega>0) then
  d_o=(omega_max-omega_min)/NstepOmega
else
  d_o=1.
endif

dV = d_theta*d_phi*d_o

spect = 0.0
energy = 0.0
integral1 = 0.0
integral2 = 0.0
integral3 = 0.0
integral = 0.

accelX(0:num_nodes,:) = 0.0
accelY(0:num_nodes,:) = 0.0
accelZ(0:num_nodes,:) = 0.0

coordX(0:num_nodes,:) = traject(0:num_nodes,1,:)
coordY(0:num_nodes,:) = traject(0:num_nodes,2,:)
coordZ(0:num_nodes,:) = traject(0:num_nodes,3,:)
velocX(0:num_nodes,:) = traject(0:num_nodes,4,:)/traject(0:num_nodes,7,:)
velocY(0:num_nodes,:) = traject(0:num_nodes,5,:)/traject(0:num_nodes,7,:)
velocZ(0:num_nodes,:) = traject(0:num_nodes,6,:)/traject(0:num_nodes,7,:)

accelX(1:num_nodes-1,:) = (velocX(2:num_nodes,:)-velocX(0:num_nodes-2,:))/dt/2.0
accelY(1:num_nodes-1,:) = (velocY(2:num_nodes,:)-velocY(0:num_nodes-2,:))/dt/2.0
accelZ(1:num_nodes-1,:) = (velocZ(2:num_nodes,:)-velocZ(0:num_nodes-2,:))/dt/2.0

wgh(:) = traject(num_nodes,8,:)

do np=1,num_part
 do i=1,NstepPhi
  SinPhi = sin(phi)
  CosPhi = cos(phi)
  do j=1,NstepTheta+1
   SinTheta = sin(theta)
   CosTheta = cos(theta)
   do k=1,NstepOmega+1
     if(comp==3) then
       integral1 = 0.0
       integral2 = 0.0
       integral3 = 0.0
     else
       integral = 0.
     endif
     outindx = 1
     do tt=0,num_nodes
       C2 = 1.0-(velocZ(tt,np)*SinTheta*CosPhi+velocY(tt,np)*SinTheta*SinPhi+velocX(tt,np)*CosTheta)
       if(omega <= 1.0/(2.0*dt*C2)) then
         C1 = accelZ(tt,np)*SinTheta*CosPhi+accelY(tt,np)*SinTheta*SinPhi+accelX(tt,np)*CosTheta
         C3 = 2.0*pi*(tt*dt - (coordZ(tt,np)*SinTheta*CosPhi+coordY(tt,np)*SinTheta*SinPhi+coordX(tt,np)*CosTheta))
         if(comp==0) then
           C4 = ((CosTheta-velocX(tt,np))*C1 - accelX(tt,np)*C2)/C2**2
           integral = integral + C4*exp(ii*omega*C3)*dt2p
         elseif(comp==1) then
           C4 = ((SinTheta*SinPhi-velocY(tt,np))*C1 - accelY(tt,np)*C2)/C2**2
           integral = integral + C4*exp(ii*omega*C3)*dt2p
         elseif(comp==2) then
           C4 = ((SinTheta*CosPhi-velocZ(tt,np))*C1 - accelZ(tt,np)*C2)/C2**2
           integral = integral + C4*exp(ii*omega*C3)*dt2p
         elseif(comp==3) then
           C4 = ((CosTheta-velocX(tt,np))*C1 - accelX(tt,np)*C2)/C2**2
           integral1 = integral1 + C4*exp(ii*omega*C3)*dt2p
           C4 = ((SinTheta*SinPhi-velocY(tt,np))*C1 - accelY(tt,np)*C2)/C2**2
           integral2 = integral2 + C4*exp(ii*omega*C3)*dt2p
           C4 = ((SinTheta*CosPhi-velocZ(tt,np))*C1 - accelZ(tt,np)*C2)/C2**2
           integral3 = integral3 + C4*exp(ii*omega*C3)*dt2p
         endif
       endif
       if((out_step>0).AND.(MOD(tt,out_step).eq.0)) then
         if(comp==3) then
           energy(tt/out_step) = energy(tt/out_step) + wgh(np)*(abs(integral1)**2 + abs(integral2)**2 + abs(integral3)**2) &
             *SinTheta*dV
         else
           energy(tt/out_step) = energy(tt/out_step) + wgh(np)*abs(integral)**2*SinTheta*dV
         endif
       endif
       if((tt.eq.outstep*outindx).OR.(tt.eq.num_nodes)) then
         if(comp==3) then
           spect(k,j,i,outindx) = spect(k,j,i,outindx)+ wgh(np)*(abs(integral1)**2 + abs(integral2)**2 + abs(integral3)**2)
         else
           spect(k,j,i,outindx) = spect(k,j,i,outindx)+ wgh(np)*abs(integral)**2
         endif
!         if((Nout>1).AND. outindx<Nout) then
!            spect(k,j,i,outindx+1) = spect(k,j,i,outindx)
         outindx = outindx+1
!         endif
       endif
     enddo
     omega = omega+d_o
   enddo
   omega = omega_min
   theta = theta + d_theta
  enddo
  theta = theta_min
  phi = phi+d_phi
 enddo
 phi = phi_min
enddo

end subroutine

subroutine spectrum_comp_coh(traject,comp,dt, phi_min,phi_max,NstepPhi,theta_min,theta_max,NstepTheta,&
    omega_min,omega_max,NstepOmega,Nout,out_step, spect,energy,num_nodes,num_part)

implicit none
integer, intent(in) :: num_nodes,num_part,comp, Nout,out_step
integer, intent(in) :: NstepOmega, NstepPhi,NstepTheta
integer             :: i,j,np,k,tt, outindx, outstep
real (kind=8), dimension(0:num_nodes,8,num_part), intent(in) :: traject
real (kind=8) , intent(in) :: dt, phi_min,phi_max, omega_min,omega_max, theta_min,theta_max
complex (kind=8), dimension(1:NstepOmega+1,1:NstepTheta+1, 1:NstepPhi,1:Nout), intent(out) :: spect
real (kind=8), dimension(0:num_nodes/out_step) , intent(out) :: energy
real (kind=8) :: d_phi, d_o, d_theta ,phi, theta, omega, pi = 3.141592653589793, &
                 SinPhi, CosPhi, SinTheta, CosTheta
real (kind=8), dimension(0:num_nodes,num_part) :: velocX, velocY, velocZ, coordX, coordY, & 
                                                  coordZ, accelX, accelY, accelZ
real (kind=8), dimension(num_part) :: wgh
real (kind=8)    :: C1, C2, C4, C3, dV, dt2p
complex (kind=8) :: ii=(0.0, 1.0), integral, integral1, integral2, integral3

!f2py intent(in) :: traject,comp,dt,phi_min,phi_max,NstepPhi,omega_min,omega_max,NstepOmega
!f2py intent(in) :: theta_min, theta_max, NstepTheta, Nout,out_step
!f2py intent(out) :: spect,energy
!f2py intent(hide) :: num_nodes,num_part

dt2p = (dt/2.0/pi)

outstep = INT(CEILING(float(num_nodes)/float(Nout)))

phi = phi_min
theta = theta_min
omega = omega_min

if(NstepTheta>0) then
  d_theta = (theta_max-theta_min)/NstepTheta
else
  d_theta = 1.0
endif
if(NstepPhi>0) then
  d_phi = (phi_max-phi_min)/NstepPhi
else
  d_phi = 1.
endif
if(NstepOmega>0) then
  d_o=(omega_max-omega_min)/NstepOmega
else
  d_o=1.
endif

dV = d_theta*d_phi*d_o

spect = 0.0
energy = 0.0
integral1 = 0.0
integral2 = 0.0
integral3 = 0.0
integral = 0.

accelX(0:num_nodes,:) = 0.0
accelY(0:num_nodes,:) = 0.0
accelZ(0:num_nodes,:) = 0.0

coordX(0:num_nodes,:) = traject(0:num_nodes,1,:)
coordY(0:num_nodes,:) = traject(0:num_nodes,2,:)
coordZ(0:num_nodes,:) = traject(0:num_nodes,3,:)
velocX(0:num_nodes,:) = traject(0:num_nodes,4,:)/traject(0:num_nodes,7,:)
velocY(0:num_nodes,:) = traject(0:num_nodes,5,:)/traject(0:num_nodes,7,:)
velocZ(0:num_nodes,:) = traject(0:num_nodes,6,:)/traject(0:num_nodes,7,:)

accelX(1:num_nodes-1,:) = (velocX(2:num_nodes,:)-velocX(0:num_nodes-2,:))/dt/2.0
accelY(1:num_nodes-1,:) = (velocY(2:num_nodes,:)-velocY(0:num_nodes-2,:))/dt/2.0
accelZ(1:num_nodes-1,:) = (velocZ(2:num_nodes,:)-velocZ(0:num_nodes-2,:))/dt/2.0

wgh(:) = traject(num_nodes,8,:)

do np=1,num_part
 do i=1,NstepPhi
  SinPhi = sin(phi)
  CosPhi = cos(phi)
  do j=1,NstepTheta+1
   SinTheta = sin(theta)
   CosTheta = cos(theta)
   do k=1,NstepOmega+1
     if(comp==3) then
       integral1 = 0.0
       integral2 = 0.0
       integral3 = 0.0
     else
       integral = 0.
     endif
     outindx = 1
     do tt=0,num_nodes
       C2 = 1.0-(velocZ(tt,np)*SinTheta*CosPhi+velocY(tt,np)*SinTheta*SinPhi+velocX(tt,np)*CosTheta)
       if(omega <= 1.0/(2.0*dt*C2)) then
         C1 = accelZ(tt,np)*SinTheta*CosPhi+accelY(tt,np)*SinTheta*SinPhi+accelX(tt,np)*CosTheta
         C3 = 2.0*pi*(tt*dt - (coordZ(tt,np)*SinTheta*CosPhi+coordY(tt,np)*SinTheta*SinPhi+coordX(tt,np)*CosTheta))
         if(comp==0) then
           C4 = ((CosTheta-velocX(tt,np))*C1 - accelX(tt,np)*C2)/C2**2
           integral = integral + C4*exp(ii*omega*C3)*dt2p
         elseif(comp==1) then
           C4 = ((SinTheta*SinPhi-velocY(tt,np))*C1 - accelY(tt,np)*C2)/C2**2
           integral = integral + C4*exp(ii*omega*C3)*dt2p
         elseif(comp==2) then
           C4 = ((SinTheta*CosPhi-velocZ(tt,np))*C1 - accelZ(tt,np)*C2)/C2**2
           integral = integral + C4*exp(ii*omega*C3)*dt2p
         elseif(comp==3) then
           C4 = ((CosTheta-velocX(tt,np))*C1 - accelX(tt,np)*C2)/C2**2
           integral1 = integral1 + C4*exp(ii*omega*C3)*dt2p
           C4 = ((SinTheta*SinPhi-velocY(tt,np))*C1 - accelY(tt,np)*C2)/C2**2
           integral2 = integral2 + C4*exp(ii*omega*C3)*dt2p
           C4 = ((SinTheta*CosPhi-velocZ(tt,np))*C1 - accelZ(tt,np)*C2)/C2**2
           integral3 = integral3 + C4*exp(ii*omega*C3)*dt2p
         endif
       endif
       if((out_step>0).AND.(MOD(tt,out_step).eq.0)) then
         if(comp==3) then
           energy(tt/out_step) = energy(tt/out_step) + wgh(np)*(abs(integral1)**2 + &
             abs(integral2)**2 + abs(integral3)**2)*SinTheta*dV
         else
           energy(tt/out_step) = energy(tt/out_step) + wgh(np)*abs(integral)**2*SinTheta*dV
         endif
       endif
       if((tt.eq.outstep*outindx).OR.(tt.eq.num_nodes)) then
         if(comp==3) then
           spect(k,j,i,outindx) = spect(k,j,i,outindx)+ integral1 + integral2 + integral3
         else
           spect(k,j,i,outindx) = spect(k,j,i,outindx)+ integral
         endif
         outindx = outindx+1
       endif
     enddo
     omega = omega+d_o
   enddo
   omega = omega_min
   theta = theta + d_theta
  enddo
  theta = theta_min
  phi = phi+d_phi
 enddo
 phi = phi_min
enddo

end subroutine

subroutine spectrum_comp_both(traject,comp,dt, phi_min,phi_max,NstepPhi,theta_min,theta_max,NstepTheta,&
    omega_min,omega_max,NstepOmega,Nout,out_step, spect_c,spect_i,energy,num_nodes,num_part)

implicit none
integer, intent(in) :: num_nodes,num_part,comp, Nout,out_step
integer, intent(in) :: NstepOmega, NstepPhi,NstepTheta
integer             :: i,j,np,k,tt, outindx, outstep
real (kind=8), dimension(0:num_nodes,8,num_part), intent(in) :: traject
real (kind=8) , intent(in) :: dt, phi_min,phi_max, omega_min,omega_max, theta_min,theta_max
complex (kind=8), dimension(1:NstepOmega+1,1:NstepTheta+1, 1:NstepPhi,1:Nout), intent(out) ::spect_c
real (kind=8), dimension(1:NstepOmega+1,1:NstepTheta+1, 1:NstepPhi,1:Nout), intent(out) :: spect_i
real (kind=8), dimension(0:num_nodes/out_step) , intent(out) :: energy
real (kind=8) :: d_phi, d_o, d_theta ,phi, theta, omega, pi = 3.141592653589793, &
                 SinPhi, CosPhi, SinTheta, CosTheta
real (kind=8), dimension(0:num_nodes,num_part) :: velocX, velocY, velocZ, coordX, coordY, & 
                                                  coordZ, accelX, accelY, accelZ
real (kind=8), dimension(num_part) :: wgh
real (kind=8)    :: C1, C2, C4, C3, dV, dt2p
complex (kind=8) :: ii=(0.0, 1.0), integral, integral1, integral2, integral3

!f2py intent(in) :: traject,comp,dt,phi_min,phi_max,NstepPhi,omega_min,omega_max,NstepOmega
!f2py intent(in) :: theta_min, theta_max, NstepTheta, Nout,out_step
!f2py intent(out) :: spect_c,spect_i,energy
!f2py intent(hide) :: num_nodes,num_part

dt2p = (dt/2.0/pi)

outstep = INT(CEILING(float(num_nodes)/float(Nout)))

phi = phi_min
theta = theta_min
omega = omega_min

if(NstepTheta>0) then
  d_theta = (theta_max-theta_min)/NstepTheta
else
  d_theta = 1.0
endif
if(NstepPhi>0) then
  d_phi = (phi_max-phi_min)/NstepPhi
else
  d_phi = 1.
endif
if(NstepOmega>0) then
  d_o=(omega_max-omega_min)/NstepOmega
else
  d_o=1.
endif

dV = d_theta*d_phi*d_o

spect_c = 0.0
spect_i = 0.0
energy = 0.0
integral1 = 0.0
integral2 = 0.0
integral3 = 0.0
integral = 0.

accelX(0:num_nodes,:) = 0.0
accelY(0:num_nodes,:) = 0.0
accelZ(0:num_nodes,:) = 0.0

coordX(0:num_nodes,:) = traject(0:num_nodes,1,:)
coordY(0:num_nodes,:) = traject(0:num_nodes,2,:)
coordZ(0:num_nodes,:) = traject(0:num_nodes,3,:)
velocX(0:num_nodes,:) = traject(0:num_nodes,4,:)/traject(0:num_nodes,7,:)
velocY(0:num_nodes,:) = traject(0:num_nodes,5,:)/traject(0:num_nodes,7,:)
velocZ(0:num_nodes,:) = traject(0:num_nodes,6,:)/traject(0:num_nodes,7,:)

accelX(1:num_nodes-1,:) = (velocX(2:num_nodes,:)-velocX(0:num_nodes-2,:))/dt/2.0
accelY(1:num_nodes-1,:) = (velocY(2:num_nodes,:)-velocY(0:num_nodes-2,:))/dt/2.0
accelZ(1:num_nodes-1,:) = (velocZ(2:num_nodes,:)-velocZ(0:num_nodes-2,:))/dt/2.0

wgh(:) = traject(num_nodes,8,:)


do np=1,num_part
 do i=1,NstepPhi
  SinPhi = sin(phi)
  CosPhi = cos(phi)
  do j=1,NstepTheta+1
   SinTheta = sin(theta)
   CosTheta = cos(theta)
   do k=1,NstepOmega+1
     if(comp==3) then
       integral1 = 0.0
       integral2 = 0.0
       integral3 = 0.0
     else
       integral = 0.
     endif
     outindx = 1
     do tt=0,num_nodes
       C2 = 1.0-(velocZ(tt,np)*SinTheta*CosPhi+velocY(tt,np)*SinTheta*SinPhi+velocX(tt,np)*CosTheta)
       if(omega <= 1.0/(2.0*dt*C2)) then
         C1 = accelZ(tt,np)*SinTheta*CosPhi+accelY(tt,np)*SinTheta*SinPhi+accelX(tt,np)*CosTheta
         C3 = 2.0*pi*(tt*dt - (coordZ(tt,np)*SinTheta*CosPhi+coordY(tt,np)*SinTheta*SinPhi+coordX(tt,np)*CosTheta))
         if(comp==0) then
           C4 = ((CosTheta-velocX(tt,np))*C1 - accelX(tt,np)*C2)/C2**2
           integral = integral + C4*exp(ii*omega*C3)*dt2p
         elseif(comp==1) then
           C4 = ((SinTheta*SinPhi-velocY(tt,np))*C1 - accelY(tt,np)*C2)/C2**2
           integral = integral + C4*exp(ii*omega*C3)*dt2p
         elseif(comp==2) then
           C4 = ((SinTheta*CosPhi-velocZ(tt,np))*C1 - accelZ(tt,np)*C2)/C2**2
           integral = integral + C4*exp(ii*omega*C3)*dt2p
         elseif(comp==3) then
           C4 = ((CosTheta-velocX(tt,np))*C1 - accelX(tt,np)*C2)/C2**2
           integral1 = integral1 + C4*exp(ii*omega*C3)*dt2p
           C4 = ((SinTheta*SinPhi-velocY(tt,np))*C1 - accelY(tt,np)*C2)/C2**2
           integral2 = integral2 + C4*exp(ii*omega*C3)*dt2p
           C4 = ((SinTheta*CosPhi-velocZ(tt,np))*C1 - accelZ(tt,np)*C2)/C2**2
           integral3 = integral3 + C4*exp(ii*omega*C3)*dt2p
         endif
       endif
       if((out_step>0).AND.(MOD(tt,out_step).eq.0)) then
         if(comp==3) then
           energy(tt/out_step) = energy(tt/out_step) + wgh(np)*(abs(integral1)**2 + &
                     abs(integral2)**2 + abs(integral3)**2)*SinTheta*dV
         else
           energy(tt/out_step) = energy(tt/out_step) + wgh(np)*abs(integral)**2*SinTheta*dV
         endif
       endif
       if((tt.eq.outstep*outindx).OR.(tt.eq.num_nodes)) then
         if(comp==3) then
           spect_c(k,j,i,outindx) = spect_c(k,j,i,outindx)+ integral1 + integral2 + integral3
           spect_i(k,j,i,outindx) = spect_i(k,j,i,outindx)+ wgh(np)*(abs(integral1)**2 + abs(integral2)**2 &
            + abs(integral3)**2)
         else
           spect_c(k,j,i,outindx) = spect_c(k,j,i,outindx)+ integral
           spect_i(k,j,i,outindx) = spect_i(k,j,i,outindx)+ wgh(np)*abs(integral)**2
         endif
         outindx = outindx+1
       endif
     enddo
     omega = omega+d_o
   enddo
   omega = omega_min
   theta = theta + d_theta
  enddo
  theta = theta_min
  phi = phi+d_phi
 enddo
 phi = phi_min
enddo

end subroutine
