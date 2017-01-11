module derivatives_y2_mod

USE Variablen
USE sync_mod
USE wang_johnson_mod

CONTAINS


! 1.Ableitung nach y mit periodischen RB
SUBROUTINE Dphi_Dy_2(phi,sphi,qi_i,qi_nr,qi_r) !Dphi_Dy_p(rho_u,solve,qxi_i,qxi_nr,qxi_r)

	IMPLICIT NONE

	INTEGER::h,i,j
	REAL*8::phi(-1:m_local+2,-1:n_local+2,-1:o_local+2),sphi(1:m_local,1:n_local,1:o_local)
	REAL*8::qi_i,qi_r,qi_nr,rel
	INTEGER::image(3)
	CALL SYNC_NEIGHBOURS_P	

!liest und speichert die image Koordinaten
image=this_image(rho)

rel=1


Do h=1,o_local
    !Gitterzellen, Mitte, Nebenrand, Randrechts, Randlinks
	IF(R==1.AND.S==1)THEN
	
	

i=1
	y_y(1:m_local,1,1)=qi_i*((phi(1:m_local,i+2,h)-phi(1:m_local,n_local-1,h))+28.&
*(phi(1:m_local,i+1,h)-phi(1:m_local,n_local,h)))/rel

i=2
y_y(1:m_local,2,1)=qi_i*((phi(1:m_local,i+2,h)-phi(1:m_local,n_local,h))+28.&
*(phi(1:m_local,i+1,h)-phi(1:m_local,i-1,h)))/rel


       !Mitte
    	DO i=3,n_local-2
        y_y(1:m_local,i,1)=qi_i*((phi(1:m_local,i+2,h)-phi(1:m_local,i-2,h))+28.&
*(phi(1:m_local,i+1,h)-phi(1:m_local,i-1,h)))/rel
    	ENDDO
i=n_local-1
y_y(1:m_local,n_local-1,1)=qi_i*((phi(1:m_local,1,h)-phi(1:m_local,i-2,h))+28.&
*(phi(1:m_local,i+1,h)-phi(1:m_local,i-1,h)))/rel

i=n_local
y_y(1:m_local,n_local,1)=qi_i*((phi(1:m_local,2,h)-phi(1:m_local,i-2,h))+28.&
*(phi(1:m_local,1,h)-phi(1:m_local,i-1,h)))/rel


	ELSEIF(n_pos==S)THEN
	
i=n_local-1
y_y(1:m_local,n_local-1,1)=qi_i*((phi(1:m_local,1,h)-phi(1:m_local,i-2,h))+28.&
*(phi(1:m_local,i+1,h)-phi(1:m_local,i-1,h)))/rel

i=n_local
y_y(1:m_local,n_local,1)=qi_i*((phi(1:m_local,2,h)-phi(1:m_local,i-2,h))+28.&
*(phi(1:m_local,1,h)-phi(1:m_local,i-1,h)))/rel
	
	!Mitte
    	DO i=1,n_local-2

        y_y(1:m_local,i,1)=qi_i*((phi(1:m_local,i+2,h)-phi(1:m_local,i-2,h))+28.&
*(phi(1:m_local,i+1,h)-phi(1:m_local,i-1,h)))/rel
	
    	ENDDO

	ELSEIF(n_pos==1)THEN
	
	i=1
	y_y(1:m_local,1,1)=qi_i*((phi(1:m_local,i+2,h)-phi(1:m_local,n_local-1,h))+28.&
*(phi(1:m_local,i+1,h)-phi(1:m_local,n_local,h)))/rel

i=2
y_y(1:m_local,2,1)=qi_i*((phi(1:m_local,i+2,h)-phi(1:m_local,n_local,h))+28.&
*(phi(1:m_local,i+1,h)-phi(1:m_local,i-1,h)))/rel


       !Mitte
       DO i=3,n_local
       y_y(1:m_local,i,1)=qi_i*((phi(1:m_local,i+2,h)-phi(1:m_local,i-2,h))+28.&
*(phi(1:m_local,i+1,h)-phi(1:m_local,i-1,h)))/rel
	
    	ENDDO

	ELSE

	!Mitte
    	DO i=1,n_local
	
        y_y(1:m_local,i,1)=qi_i*((phi(1:m_local,i+2,h)-phi(1:m_local,i-2,h))+28.&
*(phi(1:m_local,i+1,h)-phi(1:m_local,i-1,h)))/rel	
   	ENDDO
	ENDIF


		
	CALL WANG_JOHNSSON_P_z(image(3))
	sphi(1:m_local,1:n_local,h)=x_y(1:m_local,1:n_local,1)*rel	!Dimensionierung	
	
	CALL SYNC_NEIGHBOURS_P
END DO
	RETURN

END SUBROUTINE

end module
