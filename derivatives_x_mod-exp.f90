module derivatives_x_mod
use wang_johnson_mod
USE Variablen
USE sync_mod


CONTAINS

subroutine Dphi_Dx_p(phi,sphi,qi_i,qi_nr,qi_r) !Dphi_Dx_p(rho_u,solve,qxi_i,qxi_nr,qxi_r,u_0*rho_0)

!sphi steht für solution (gem. ist die Ableitung) der Größe phi

	INTEGER::h,i
	REAL::phi(-1:m_local+2,-1:n_local+2,-1:o_local+2),sphi(1:m_local,1:n_local,1:o_local)
	REAL::qi_i,qi_r,qi_nr
	Integer::image(3)
!	CALL SYNC_NEIGHBOURS_P

!liest und speichert die image Koordinaten
image=this_image(rho)

sphi(1,1:n_local,1:o_local)=(phi(1,1:n_local,1:o_local)-phi(2,1:n_local,1:o_local))

	DO i=2,n_local-1
		DO h=1,m_local
			sphi(h,i,1:o_local)=(phi(h+1,i,1:o_local)-2.*phi(h,i,1:o_local)+phi(h-1,i,1:o_local))
		ENDDO
	ENDDO
sphi(n_local,1:n_local,1:o_local)=(phi(N_LOCAL-1,1:n_local,1:o_local)-phi(n_local,1:n_local,1:o_local))
	
	RETURN


end subroutine

!1. Ableitung
subroutine Dphi_Dx_p_neu(phi,sphi,qi_i,qi_nr,qi_r) !Dphi_Dx_p(rho_u,solve,qxi_i,qxi_nr,qxi_r,u_0*rho_0)

!sphi steht für solution (gem. ist die Ableitung) der Größe phi

	INTEGER::h,i
	REAL::phi(-1:m_local+2,-1:n_local+2,-1:o_local+2),sphi(1:m_local,1:n_local,1:o_local)
	REAL::qi_i,qi_r,qi_nr
	Integer::image(3)
!	CALL SYNC_NEIGHBOURS_P

!liest und speichert die image Koordinaten
image=this_image(rho)

sphi(1,1:n_local,1:o_local)=(phi(1,1:n_local,1:o_local)-phi(2,1:n_local,1:o_local))

	DO i=2,n_local-1
		DO h=1,m_local
			sphi(h,i,1:o_local)=(phi(1,1:n_local,1:o_local)-phi(2,1:n_local,1:o_local))
		ENDDO
	ENDDO
sphi(n_local,1:n_local,1:o_local)=0
	
	RETURN


end subroutine


! 2. Ableitung


end module
