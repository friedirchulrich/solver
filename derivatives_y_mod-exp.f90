module derivatives_y_mod

USE Variablen
USE sync_mod
USE wang_johnson_mod

CONTAINS


! 1.Ableitung nach y mit periodischen RB
SUBROUTINE Dphi_Dy_p(phi,sphi,qi_i,qi_nr,qi_r) !Dphi_Dy_p(rho_u,solve,qxi_i,qxi_nr,qxi_r)

	IMPLICIT NONE

	INTEGER::h,i,j
	REAL::phi(-1:m_local+2,-1:n_local+2,-1:o_local+2),sphi(1:m_local,1:n_local,1:o_local)
	REAL::qi_i,qi_r,qi_nr
	INTEGER::image(3)
	CALL SYNC_NEIGHBOURS_P	


Do j=1,m_local
	DO i=1,o_local
		sphi(j,1,i)=0.5*(phi(j,2,i)-phi(j,n_local,i))/dy
		DO h=2,n_local-1
			sphi(j,h,i)=0.5*(phi(j,h+1,i)-phi(j,h-1,i))/dy
		ENDDO
		sphi(j,n_local,i)=0.5*(phi(j,1,i)-phi(j,n_local-1,i))/dy
	ENDDO
ENDDO

	CALL SYNC_NEIGHBOURS_P

	RETURN

END SUBROUTINE

end module
