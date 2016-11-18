module derivatives_y_mod

USE Variablen
USE sync_mod
USE wang_johnson_mod
USE derivatives_y2_mod

CONTAINS


! 1.Ableitung nach y mit periodischen RB
SUBROUTINE Dphi_Dy_p(phi,sphi,qi_i,qi_nr,qi_r) !Dphi_Dy_p(rho_u,solve,qxi_i,qxi_nr,qxi_r)

	IMPLICIT NONE

	INTEGER::h,i,j,exi
	REAL::phi(-1:m_local+2,-1:n_local+2,-1:o_local+2),sphi(1:m_local,1:n_local,1:o_local)
	REAL::qi_i,qi_r,qi_nr,rel
	INTEGER::image(3)
	CALL SYNC_NEIGHBOURS_P	

!liest und speichert die image Koordinaten
image=this_image(rho)

rel=1.0
!rel=phi(1,1,1)
!DO h=1,m_local
!DO i=1,n_local
!DO j=1,o_local
!rel=max(rel,phi(h,i,j))
!end do
!end do
!end do
!if(rel==0)then
!rel =1
!end if
exi=1
if(exi==99)then
Do h=1,o_local
    !Gitterzellen, Mitte, Nebenrand, Randrechts, Randlinks
	IF(R==1.AND.S==1)THEN
	
	!RandUnten
        y_y(1:m_local,1,1)=qi_r *(-5.*phi(1:m_local,1,h)+4.*phi(1:m_local,2,h)+phi(1:m_local,3,h))/rel
	
       !RandOben
       y_y(1:m_local,n_local,1)=0.0
!qi_r *(-1.*phi(1:m_local,n_local-2,h)-4.&
!*phi(1:m_local,n_local-1,h)+5.*phi(1:m_local,n_local,h))/rel
       !Nebenrand_Unten
       y_y(1:m_local,2,1)=qi_nr*(3.*(phi(1:m_local,3,h)-phi(1:m_local,1,h)))/rel
       !Nebenrand_Oben
       y_y(1:m_local,n_local-1,1)=qi_nr*(3.*(phi(1:m_local,n_local,h)-phi(1:m_local,n_local-2,h)))/rel
       !Mitte
    	DO i=3,n_local-2
        y_y(1:m_local,i,1)=qi_i*((phi(1:m_local,i+2,h)-phi(1:m_local,i-2,h))+28.&
*(phi(1:m_local,i+1,h)-phi(1:m_local,i-1,h)))/rel
    	ENDDO

	ELSEIF(n_pos==S)THEN
	
 	!Randoben
        y_y(1:m_local,n_local,1)=0.0
!qi_r *(-1.*phi(1:m_local,n_local-3,h)-4.&
!*phi(1:m_local,n_local-2,h)+5.*phi(1:m_local,n_local-1,h))/rel
	!Nebenrand-Oben
        y_y(1:m_local,n_local-1,1)=qi_nr*(3.*(phi(1:m_local,n_local,h)-phi(1:m_local,n_local-2,h)))/rel	
	!Mitte
    	DO i=1,n_local-2

        y_y(1:m_local,i,1)=qi_i*((phi(1:m_local,i+2,h)-phi(1:m_local,i-2,h))+28.&
*(phi(1:m_local,i+1,h)-phi(1:m_local,i-1,h)))/rel
	
    	ENDDO

	ELSEIF(n_pos==1)THEN
	
	!Randunten
        y_y(1:m_local,1,1)=qi_r *(-5.*phi(1:m_local,1,h)+4.*phi(1:m_local,2,h)+phi(1:m_local,3,h))/rel	
       !Nebenrand_Unten
        y_y(1:m_local,2,1)=qi_nr*(3.*(phi(1:m_local,3,h)-phi(1:m_local,1,h)))/rel	


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


!	n_local=n_local-1
	CALL WANG_JOHNSSON_P_y(image(3))
!	n_local=n_local+1
!	sphi(1:m_local,1:n_local,h)=x_y(1:m_local,1:n_local,1)*rel	!Dimensionierung	
	sphi(1:m_local,n_local,h)=0.0
	CALL SYNC_NEIGHBOURS_P
END DO

endif
CALL Dphi_Dy_2(phi,sphi,qi_i,qi_nr,qi_r)
!	 sphi(1:m_local,1:n_local,1:o_local)=0.0
	RETURN

END SUBROUTINE

end module
