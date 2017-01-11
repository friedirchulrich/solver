module ymom_mod
! Berechnet den Y-Impuls

USE Variablen
USE derivatives_x_mod
USE derivatives_y_mod
USE derivatives_z_mod
USE neumann_mod
contains

subroutine ymom
integer::i,j,h
!delta_u, delta_v , delta_w 
! werden als subroutinen genutzt

! dphi Hilfsvariable
! rho_v*u-tau_xy
dphi(1:m_local,1:n_local,1:o_local)=&
rho_v(1:m_local,1:n_local,1:o_local)*&
u(1:m_local,1:n_local,1:o_local)+&
tau_xy(1:m_local,1:n_local,1:o_local)

!ableiten
CALL Dphi_Dx_p(dphi,delta_u,qxi_i,qxi_nr,qxi_r)


! rho_v*v+p-tau_yy
dphi(-1:m_local+2,-1:n_local+2,-1:o_local+2)=&
rho_v(-1:m_local+2,-1:n_local+2,-1:o_local+2)*&
v(-1:m_local+2,-1:n_local+2,-1:o_local+2)+&
p(-1:m_local+2,-1:n_local+2,-1:o_local+2)-&
tau_yy(-1:m_local+2,-1:n_local+2,-1:o_local+2)

!dphi=rho_v*v+p-tau_yy
!ableiten
!CALL Dphi_Dx_p(rho_w,delta_v,qxi_i,qxi_nr,qxi_r)

CALL Dphi_Dy_p(dphi,delta_v,qyi_i,qyi_nr,qyi_r)
!CALL Dphi_Dz_p(rho_w,delta_v,qzi_i,qzi_nr,qzi_r)


! rho_w*v-tau_yz
dphi(1:m_local,1:n_local,1:o_local)=&
rho_w(1:m_local,1:n_local,1:o_local)*&
v(1:m_local,1:n_local,1:o_local)-&
tau_zy(1:m_local,1:n_local,1:o_local)

!ableiten
CALL Dphi_Dz_p(dphi,delta_w,qzi_i,qzi_nr,qzi_r)

!rho_u_neu=rho_u-dt*(delta_u+delta_v+delta_w)
rho_v_neu(1:m_local,1:n_local,1:o_local)=&
rho_v(1:m_local,1:n_local,1:o_local)-&
dt*(delta_u(1:m_local,1:n_local,1:o_local)+&
delta_v(1:m_local,1:n_local,1:o_local)+&
delta_w(1:m_local,1:n_local,1:o_local))

!rho_v_neu=rho_v-dt*(delta_u+delta_v+delta_w)

end subroutine



end module

