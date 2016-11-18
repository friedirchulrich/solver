module energy_mod
! Berechnet die Energiegleichung

USE Variablen
USE derivatives_x_mod
USE derivatives_y_mod
USE derivatives_z_mod
USE neumann_mod

contains

subroutine energy

! delta_x Term
CALL Dphi_Dx_p(T,dphi,qxi_i,qxi_nr,qxi_r)
CALL Dphi_Dx_p_neu(p,delta_v,qxi_i,qxi_nr,qxi_r)

E_plus_p(1:m_local,1:n_local,1:o_local)=&
u(1:m_local,1:n_local,1:o_local)*&
(rho_E(1:m_local,1:n_local,1:o_local)&
!p(1:m_local,1:n_local,1:o_local))-&
)-tau_xx(1:m_local,1:n_local,1:o_local)*&
u(1:m_local,1:n_local,1:o_local)-&
tau_xy(1:m_local,1:n_local,1:o_local)*&
v(1:m_local,1:n_local,1:o_local)-&
tau_xz(1:m_local,1:n_local,1:o_local)*&
w(1:m_local,1:n_local,1:o_local)-&
k(1:m_local,1:n_local,1:o_local)*&
dphi(1:m_local,1:n_local,1:o_local)


CALL Dphi_Dx_p(E_plus_p,delta_u,qxi_i,qxi_nr,qxi_r)

delta_u(1:m_local,1:n_local,1:o_local)=&
delta_u(1:m_local,1:n_local,1:o_local)+&
delta_v(1:m_local,1:n_local,1:o_local)
! delta_y Term
CALL Dphi_Dy_p(T,dphi,qyi_i,qyi_nr,qyi_r)

E_plus_p(1:m_local,1:n_local,1:o_local)=&
v(1:m_local,1:n_local,1:o_local)*&
(rho_E(1:m_local,1:n_local,1:o_local)+&
p(1:m_local,1:n_local,1:o_local))-&
tau_xy(1:m_local,1:n_local,1:o_local)*&
u(1:m_local,1:n_local,1:o_local)-&
tau_yy(1:m_local,1:n_local,1:o_local)*&
v(1:m_local,1:n_local,1:o_local)-&
tau_zy(1:m_local,1:n_local,1:o_local)*&
w(1:m_local,1:n_local,1:o_local)-&
k(1:m_local,1:n_local,1:o_local)*&
dphi(1:m_local,1:n_local,1:o_local)


CALL Dphi_Dy_p(E_plus_p,delta_v,qyi_i,qyi_nr,qyi_r)

! delta_z Term
CALL Dphi_Dz_p(T,dphi,qzi_i,qzi_nr,qzi_r)

E_plus_p(1:m_local,1:n_local,1:o_local)=&
w(1:m_local,1:n_local,1:o_local)*&
(rho_E(1:m_local,1:n_local,1:o_local)+&
p(1:m_local,1:n_local,1:o_local))-&
tau_xz(1:m_local,1:n_local,1:o_local)*&
u(1:m_local,1:n_local,1:o_local)-&
tau_zy(1:m_local,1:n_local,1:o_local)*&
v(1:m_local,1:n_local,1:o_local)-&
tau_zz(1:m_local,1:n_local,1:o_local)*&
w(1:m_local,1:n_local,1:o_local)-&
k(1:m_local,1:n_local,1:o_local)*&
dphi(1:m_local,1:n_local,1:o_local)


CALL Dphi_Dz_p(E_plus_p,delta_w,qzi_i,qzi_nr,qzi_r)

!!bc neumann
call bcneu

!rho_u_neu=rho_u-dt*(delta_u+delta_v+delta_w)
rho_E_neu(1:m_local,1:n_local,1:o_local)=&
rho_E(1:m_local,1:n_local,1:o_local)-&
dt*(delta_u(1:m_local,1:n_local,1:o_local)+&
delta_v(1:m_local,1:n_local,1:o_local)+&
delta_w(1:m_local,1:n_local,1:o_local))
end subroutine


end module

