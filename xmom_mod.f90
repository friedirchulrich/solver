module xmom_mod
! Berechnet den X-Impuls

USE Variablen
USE derivatives_x_mod
USE derivatives_y_mod
USE derivatives_z_mod
USE neumann_mod
contains

subroutine xmom
integer::i,j,h,d
!delta_u, delta_v ,h, delta_w 
! werden als subroutinen genutzt

! dphi Hilfsvariable
! rho_u*u+p-tau_xx
!p(1:m_local,1:n_local,1:o_local)=101307.0


!CALL Dphi_Dx_p_neu(p,delta_v,qxi_i,qxi_nr,qxi_r)


!delta_v(1:m_local,1:n_local,1:o_local)=0.0

!rho_u(1:m_local,1:n_local,1:o_local)=432324.1

dphi(1:m_local,1:n_local,1:o_local)=&
rho_u(1:m_local,1:n_local,1:o_local)*&
u(1:m_local,1:n_local,1:o_local)+&
p(1:m_local,1:n_local,1:o_local)-&
tau_xx(1:m_local,1:n_local,1:o_local)

!ableiten
CALL Dphi_Dx_p(dphi,delta_u,qxi_i,qxi_nr,qxi_r)


!delta_u(1:m_local,1:n_local,1:o_local)=&
!delta_u(1:m_local,1:n_local,1:o_local)+&
!delta_v(1:m_local,1:n_local,1:o_local)


! rho_v*u-tau_xy
dphi(1:m_local,1:n_local,1:o_local)=&
rho_v(1:m_local,1:n_local,1:o_local)*&
u(1:m_local,1:n_local,1:o_local)-&
tau_xy(1:m_local,1:n_local,1:o_local)

!ableiten
CALL Dphi_Dy_p(dphi,delta_v,qyi_i,qyi_nr,qyi_r)

!j=20
!h=20
!do i=1,n_local
!do h=1,n_local
!do i=1,m_local
!	if(delta_v(i,j,h).NE.0)then
!print *, tau_xy(i,j,h),delta_v(i,j,h),i,j,h
!endif
!end do
!end do
!end do


! rho_w*u-tau_xz
dphi(1:m_local,1:n_local,1:o_local)=&
rho_w(1:m_local,1:n_local,1:o_local)*&
u(1:m_local,1:n_local,1:o_local)-&
tau_xz(1:m_local,1:n_local,1:o_local)

!ableiten
CALL Dphi_Dz_p(dphi,delta_w,qzi_i,qzi_nr,qzi_r)

!!bc neumann
!call bcneu

!do I=1,m_local
!print *, delta_u(i,20,20),delta_v(i,20,20),delta_w(i,20,20)
!end do

!rho_u_neu=rho_u-dt*(delta_u+delta_v+delta_w)
rho_u_neu(1:m_local,1:n_local,1:o_local)=&
rho_u(1:m_local,1:n_local,1:o_local)-&
dt*(delta_u(1:m_local,1:n_local,1:o_local)+&
delta_v(1:m_local,1:n_local,1:o_local)+&
delta_w(1:m_local,1:n_local,1:o_local))


end subroutine


end module

