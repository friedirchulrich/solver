module spannungstensor_mod

USE Variablen
USE derivatives_x_mod
USE derivatives_y_mod
USE derivatives_z_mod
USE neumann_mod

contains

subroutine tensor
integer::i,j,k,h

! delta_u berechnen

CALL Dphi_Dx_p(u,delta_u,qxi_i,qxi_nr,qxi_r)
! delta_v berechnen

CALL Dphi_Dy_p(v,delta_v,qyi_i,qyi_nr,qyi_r)

! delta_w berechnen

CALL Dphi_Dz_p(w,delta_w,qzi_i,qzi_nr,qzi_r)



sync all 

!tau_xx(:,:,:)[:,:,:]
tau_xx(1:m_local,1:n_local,1:o_local)=(2./3.)*&
mu(1:m_local,1:n_local,1:o_local)&
*(2.*delta_u(1:m_local,1:n_local,1:o_local)&
-delta_v(1:m_local,1:n_local,1:o_local)&
-delta_w(1:m_local,1:n_local,1:o_local))
!tau_yy(:,:,:)[:,:,:]
tau_yy(1:m_local,1:n_local,1:o_local)=(2./3.)*&
mu(1:m_local,1:n_local,1:o_local)&
*(-delta_u(1:m_local,1:n_local,1:o_local)&
+2.*delta_v(1:m_local,1:n_local,1:o_local)&
-delta_w(1:m_local,1:n_local,1:o_local))
!tau_zz(:,:,:)[:,:,:]
tau_zz(1:m_local,1:n_local,1:o_local)=(2./3.)*&
mu(1:m_local,1:n_local,1:o_local)&
*(-delta_u(1:m_local,1:n_local,1:o_local)&
-delta_v(1:m_local,1:n_local,1:o_local)&
+2.*delta_w(1:m_local,1:n_local,1:o_local))

!derivaitve_z=0
 !tau_zz(1:m_local,1:n_local,1:o_local)=0.0
! delta_u berechnen

CALL Dphi_Dx_p(v,delta_v,qxi_i,qxi_nr,qxi_r)
! delta_v berechnen
CALL Dphi_Dy_p(u,delta_u,qyi_i,qyi_nr,qyi_r)


!print *, "delta_v",  delta_v(12,1,20), delta_v(12,2,20),  delta_v(12,3,20),  delta_v(12,3,20)
!print *, "delta_u",  delta_u(12,1,20), delta_u(12,2,20),  delta_u(12,3,20),  delta_u(12,3,20)

!tau_xy(:,:,:)[:,:,:]
tau_xy(1:m_local,1:n_local,1:o_local)=&
mu(1:m_local,1:n_local,1:o_local)&
*(delta_u(1:m_local,1:n_local,1:o_local)&
+delta_v(1:m_local,1:n_local,1:o_local))



! delta_u berechnen
CALL Dphi_Dx_p(w,delta_w,qxi_i,qxi_nr,qxi_r)
!!ERROR??
! delta_w berechnen
CALL Dphi_Dz_p(u,delta_u,qzi_i,qzi_nr,qzi_r)


!tau_xz(:,:,:)[:,:,:]
tau_xz(1:m_local,1:n_local,1:o_local)=&
mu(1:m_local,1:n_local,1:o_local)&
*(delta_u(1:m_local,1:n_local,1:o_local)&
+delta_w(1:m_local,1:n_local,1:o_local))

!derivaitve_z=0
!  tau_xz(1:m_local,1:n_local,1:o_local)=0.0
! delta_w berechnen
CALL Dphi_Dy_p(w,delta_w,qyi_i,qyi_nr,qyi_r)

! delta_v berechnen
CALL Dphi_Dz_p(v,delta_v,qzi_i,qzi_nr,qzi_r)

!tau_zy(:,:,:)[:,:,:]
tau_zy(1:m_local,1:n_local,1:o_local)=&
mu(1:m_local,1:n_local,1:o_local)&
*(delta_v(1:m_local,1:n_local,1:o_local)&
+delta_w(1:m_local,1:n_local,1:o_local))



!derivaitve_z=0
!  tau_zy(1:m_local,1:n_local,1:o_local)=0.0
!print *,"tau", u(1,1:40,10)

!test
!tau_xx(1:m_local,1:n_local,1:o_local)=0.0
!tau_yy(1:m_local,1:n_local,1:o_local)=0.0
!tau_zz(1:m_local,1:n_local,1:o_local)=0.0
!tau_xz(1:m_local,1:n_local,1:o_local)=0.0
!tau_xy(1:m_local,1:n_local,1:o_local)=0.0
!tau_zy(1:m_local,1:n_local,1:o_local)=0.0
end subroutine

end module

