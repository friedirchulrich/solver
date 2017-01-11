module derivatives_x_mod
use wang_johnson_mod
USE Variablen
USE sync_mod


CONTAINS

!1. Ableitung
subroutine Dphi_Dx_p(phi,sphi,qi_i,qi_nr,qi_r) !Dphi_Dx_p(rho_u,solve,qxi_i,qxi_nr,qxi_r,u_0*rho_0)

!sphi steht für solution (gem. ist die Ableitung) der Größe phi

	INTEGER::h,i,j
	REAL*8::phi(-1:m_local+2,-1:n_local+2,-1:o_local+2),sphi(1:m_local,1:n_local,1:o_local)
	REAL*8::qi_i,qi_r,qi_nr,rel
	Integer::image(3)
	CALL SYNC_NEIGHBOURS_P

!liest und speichert die image Koordinaten
image=this_image(rho)
rel=1.0

!y_x(1:m_local,1:n_local,1:2)=0.0
DO I=1,o_local

    !Gitterzellen, Mitte, Nebenrand, Randrechts, Randlinks
	IF(R==1.AND.S==1)THEN
	
	!RandLinks
        !y_x(1,1:n_local,1)=qi_nr*(-37.*phi(1,1:n_local,i)+8.*phi(2,1:n_local,i)&
!+36*phi(3,1:n_local,i)-8.*phi(4,1:n_local,i)+phi(5,1:n_local,i)) 
	y_x(1,1:n_local,1)=qi_r *((phi(3,1:n_local,i)-phi(1,1:n_local,i))+4.*(phi(2,1:n_local,i)-phi(1,1:n_local,i)))
	
       !RandRechts
       y_x(m_local,1:n_local,1)=qi_r *((phi(m_local,1:n_local,i)-phi(m_local-2,1:n_local,i))&
+4.*(phi(m_local,1:n_local,i)-phi(m_local-1,1:n_local,i)))

       !Nebenrand_Links
       y_x(2,1:n_local,1)=qi_nr*(3.*(phi(3,1:n_local,i)-phi(1,1:n_local,i)))/rel
	
       !Nebenrand_Rechts
       y_x(m_local-1,1:n_local,1)=qi_nr*(3.*(phi(m_local,1:n_local,i)-phi(m_local-2,1:n_local,i)))/rel

       !Mitte
    	DO h=3,m_local-2
        y_x(h,1:n_local,1)=qi_i*((phi(h+2,1:n_local,i)-phi(h-2,1:n_local,i))+28.*(phi(h+1,1:n_local,i)-phi(h-1,1:n_local,i)))/rel
    	ENDDO

	ELSEIF(m_pos==R)THEN
	
 	!RandRechts
        y_x(m_local,1:n_local,1)=qi_r *(-1.*phi(m_local-2,1:n_local,i)-4.&
*phi(m_local-1,1:n_local,i)+5.*phi(m_local,1:n_local,i))/rel
	
	!Nebenrand-Rechts
        y_x(m_local-1,1:n_local,1)=qi_nr*(3.*(phi(m_local,1:n_local,i)-phi(m_local-2,1:n_local,i)))/rel
	
	!Mitte
    	DO h=3,n_local-2

        y_x(h,1:n_local,1)=qi_i*((phi(h+2,1:n_local,i)-phi(h-2,1:n_local,i))+28.*(phi(h+1,1:n_local,i)-phi(h-1,1:n_local,i)))/rel
    	ENDDO

	ELSEIF(m_pos==1)THEN
	
	!RandLinks
	y_x(1,1:n_local,1)=qi_r *(-5.*phi(1,1:n_local,i)+4.*phi(2,1:n_local,i)+phi(3,1:n_local,i))/rel
 !       y_x(1,1:n_local,1)=qi_i*(-37.*phi(1,1:n_local,i)+8.*phi(2,1:n_local,i)&
!+36*phi(3,1:n_local,i)-8.*phi(4,1:n_local,i)+phi(5,1:n_local,i))
        !Nebenrand_Links
       y_x(2,1:n_local,1)=qi_nr*(3.*(phi(3,1:n_local,i)-phi(1,1:n_local,i)))/rel	

        !Mitte
        DO h=3,m_local

        y_x(h,1:n_local,1)=qi_i*((phi(h+2,1:n_local,i)-phi(h-2,1:n_local,i))+28.&
*(phi(h+1,1:n_local,i)-phi(h-1,1:n_local,i)))/rel	
    	ENDDO

	ELSE

	!Mitte
    	DO h=1,m_local
	
        y_x(h,1:n_local,1)=qi_i*((phi(h+2,1:n_local,i)-phi(h-2,1:n_local,i))+28.*&
(phi(h+1,1:n_local,i)-phi(h-1,1:n_local,i)))/rel
	
   	ENDDO
	ENDIF

!	print *, y_x(1:m_local,20,1)
!	print *, "wang-Johnson"


	CALL WANG_JOHNSSON_P_x(image(3))
	sphi(1:m_local,1:n_local,i)=x_x(1:m_local,1:n_local,1)*rel!Dimensionierung	
!	IF(n_pos==R)THEN
!	sphi(m_local,1:n_local)=0.0	!Dimensionierung
!	ENDIF


	CALL SYNC_NEIGHBOURS_P

End Do


!	sphi(1:m_local,1:n_local,1:o_local)=0.0
	RETURN


end subroutine


! 1. Ableitung mit von-Neumann Bedingung
subroutine Dphi_Dx_p_neu(phi,sphi,qi_i,qi_nr,qi_r) !Dphi_Dx_p(rho_u,solve,qxi_i,qxi_nr,qxi_r,u_0*rho_0)

!sphi steht für solution (gem. ist die Ableitung) der Größe phi

	INTEGER::h,i,j
	REAL*8::phi(-1:m_local+2,-1:n_local+2,-1:o_local+2),sphi(1:m_local,1:n_local,1:o_local)
	REAL*8::qi_i,qi_r,qi_nr,rel
	Integer::image(3)
	CALL SYNC_NEIGHBOURS_P

!liest und speichert die image Koordinaten
image=this_image(rho)
rel=1.0
DO I=1,o_local

    !Gitterzellen, Mitte, Nebenrand, Randrechts, Randlinks
	IF(R==1.AND.S==1)THEN
	
	!RandLinks
        y_x(1,1:n_local,1)=qi_r *(-5.*phi(1,1:n_local,i)+4.*phi(2,1:n_local,i)+phi(3,1:n_local,i))/rel
	
       !RandRechts
       y_x(m_local,1:n_local,1)=0.0
!qi_r *(-1.*phi(m_local-3,1:n_local,i)-4.*&
!phi(m_local-2,1:n_local,i)+5.*phi(m_local-1,1:n_local,i))/rel

       !Nebenrand_Links
       y_x(2,1:n_local,1)=qi_nr*(3.*(phi(3,1:n_local,i)-phi(1,1:n_local,i)))/rel

       !Nebenrand_Rechts
       y_x(m_local-1,1:n_local,1)=qi_nr*(3.*(phi(m_local,1:n_local,i)-phi(m_local-2,1:n_local,i)))/rel

       !Mitte
    	DO h=3,m_local-2
        y_x(h,1:n_local,1)=qi_i*((phi(h+2,1:n_local,i)-phi(h-2,1:n_local,i))+28.*(phi(h+1,1:n_local,i)-phi(h-1,1:n_local,i)))/rel
    	ENDDO

	ELSEIF(m_pos==R)THEN
	
 	!RandRechts
        y_x(m_local,1:n_local,1)=0.0
!qi_r *(-1.*phi(m_local-3,1:n_local,i)-4.&
!*phi(m_local-2,1:n_local,i)+5.*phi(m_local-1,1:n_local,i))/rel
	
	!Nebenrand-Rechts
        y_x(m_local-1,1:n_local,1)=qi_nr*(3.*(phi(m_local,1:n_local,i)-phi(m_local-2,1:n_local,i)))/rel
	
	!Mitte
    	DO h=3,m_local-2

        y_x(h,1:n_local,1)=qi_i*((phi(h+2,1:n_local,i)-phi(h-2,1:n_local,i))+28.*(phi(h+1,1:n_local,i)-phi(h-1,1:n_local,i)))/rel
    	ENDDO

	ELSEIF(m_pos==1)THEN
	
	!RandLinks
        y_x(1,1:n_local,1)=qi_r *(-5.*phi(1,1:n_local,i)+4.*phi(2,1:n_local,i)+phi(3,1:n_local,i))/rel	
        !Nebenrand_Links
       y_x(2,1:n_local,1)=qi_nr*(3.*(phi(3,1:n_local,i)-phi(1,1:n_local,i)))/rel	

        !Mitte
        DO h=3,m_local

        y_x(h,1:n_local,1)=qi_i*((phi(h+2,1:n_local,i)-phi(h-2,1:n_local,i))+28.&
*(phi(h+1,1:n_local,i)-phi(h-1,1:n_local,i)))/rel	
    	ENDDO

	ELSE

	!Mitte
    	DO h=1,m_local
	
        y_x(h,1:n_local,1)=qi_i*((phi(h+2,1:n_local,i)-phi(h-2,1:n_local,i))+28.*&
(phi(h+1,1:n_local,i)-phi(h-1,1:n_local,i)))/rel
	
   	ENDDO
	ENDIF

!m_local=m_local-1	
	CALL WANG_JOHNSSON_P_x_2(image(3))
!m_local=m_local+1
	sphi(1:m_local-1,1:n_local,i)=x_x(1:m_local-1,1:n_local,1)*rel!Dimensionierung
	sphi(m_local,1:n_local,i)=0	
!	IF(n_pos==R)THEN
!	sphi(m_local,1:n_local)=0.0	!Dimensionierung
!	ENDIF
	CALL SYNC_NEIGHBOURS_P

End Do
	
	RETURN


end subroutine


end module
