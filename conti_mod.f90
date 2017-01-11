module conti_mod
! Konti-Gleichung

USE Variablen
USE derivatives_x_mod
USE derivatives_y_mod
USE derivatives_z_mod
USE neumann_mod
contains

subroutine conti

IMPLICIT NONE
integer::i,j,k
real::rho_0
	!Ableitung rho_u nach x
rho_0=1.293
	CALL Dphi_Dx_p(rho_u,delta_u,qxi_i,qxi_nr,qxi_r)! !Dphi_Dx_p(phi,sphi,qi_i,qi_nr,qi_r,rel)



!	dphi(1:m_local,1:n_local)=solve(1:m_local,1:n_local)

	! Ableitung rho_v nach y
	CALL Dphi_Dy_p(rho_v,delta_v,qyi_i,qyi_nr,qyi_r)

	! Ableitung rho_w nach z
	CALL Dphi_Dz_p(rho_w,delta_w,qzi_i,qzi_nr,qzi_r)
	! neues rho Berechnen


!bc_2
!call bcneu


	rho_neu(1:m_local,1:n_local,1:o_local)=rho(1:m_local,1:n_local,1:o_local)&
-dt*(delta_u(1:m_local,1:n_local,1:o_local)+delta_v(1:m_local,1:n_local,1:o_local)&
+delta_w(1:m_local,1:n_local,1:o_local))


rho_neu(1:m_local,1:n_local,1:o_local)=rho_0
	CALL SYNC_NEIGHBOURS_P	! rho aktuell

	RETURN

end subroutine


end module

