module derivatives_z_mod

USE Variablen
USE sync_mod
USE wang_johnson_mod

CONTAINS

subroutine Dphi_Dz_p(phi,sphi,qi_i,qi_nr,qi_r)

	IMPLICIT NONE

	INTEGER::h,i,j,k
	REAL::phi(-1:m_local+2,-1:n_local+2,-1:o_local+2)
	REAL,intent(out)::sphi(1:m_local,1:n_local,1:o_local)
	REAL::qi_i,qi_r,qi_nr,rel
	INTEGER::image(3)

Do j=1,m_local
	DO i=1,n_local

		sphi(j,i,1)=0.5*(phi(j,i,2)-phi(j,i,o_local))/dz
		DO h=2,o_local-1
			sphi(j,i,h)=0.5*(phi(j,i,h+1)-phi(j,i,h-1))/dz

		ENDDO
		sphi(j,i,o_local)=0.5*(phi(j,i,1)-phi(j,i,o_local-1))/dz

	ENDDO
ENDDO
 sync all

	CALL SYNC_NEIGHBOURS_P



	RETURN
end subroutine

end module
