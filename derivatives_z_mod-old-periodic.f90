module derivatives_z_mod

USE Variablen
USE sync_mod
USE wang_johnson_mod

CONTAINS

subroutine Dphi_Dz_p(phi,sphi,qi_i,qi_nr,qi_r)

	IMPLICIT NONE

	INTEGER::h,i,j,k
	REAL::phi(-1:m_local+2,-1:n_local+2,-1:o_local+2),sphi(1:m_local,1:n_local,1:o_local)
	REAL::qi_i,qi_r,qi_nr,rel
	INTEGER::image(3)
rel=phi(1,1,1)
DO h=1,m_local
DO i=1,n_local
DO j=1,o_local
rel=max(rel,phi(h,i,j))
end do
end do
end do
if(rel==0)then
rel =1
end if

	CALL SYNC_NEIGHBOURS_P
	
!liest und speichert die image Koordinaten
image=this_image(rho)
Do I=1,m_local

 h=1
	y_x(1,1:n_local,1)=qi_i*((phi(i,1:n_local,h+2)-phi(i,1:n_local,o_local-2))&
+28.*(phi(i,1:n_local,h+1)-phi(i,1:n_local,o_local-1)))/rel
h=2
	y_x(2,1:n_local,1)=qi_i*((phi(i,1:n_local,h+2)-phi(i,1:n_local,o_local-1))&
+28.*(phi(i,1:n_local,h+1)-phi(i,1:n_local,h-1)))/rel


	DO h=3,o_local-2
		!Entdimensionierung, Relativierung
		y_x(h,1:n_local,1)=qi_i*((phi(i,1:n_local,h+2)-phi(i,1:n_local,h-2))&
+28.*(phi(i,1:n_local,h+1)-phi(i,1:n_local,h-1)))/rel
	ENDDO
h=o_local-1
	y_x(o_local-1,1:n_local,1)=qi_i*((phi(i,1:n_local,1)-phi(i,1:n_local,h-2))&
+28.*(phi(i,1:n_local,h+1)-phi(i,1:n_local,h-1)))/rel
h=o_local
	y_x(o_local,1:n_local,1)=qi_i*((phi(i,1:n_local,2)-phi(i,1:n_local,h-2))&
+28.*(phi(i,1:n_local,1)-phi(i,1:n_local,h-1)))/rel

	CALL WANG_JOHNSSON_P_z(image(1))
do j=1,o_local
	sphi(i,1:n_local,j)=x_x(j,1:n_local,1)*rel
end do
	CALL SYNC_NEIGHBOURS_P
	
END DO


!error
!	sphi(1:m_local,1:n_local,1:o_local)=0


	RETURN
end subroutine

end module
