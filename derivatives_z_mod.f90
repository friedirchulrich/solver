module derivatives_z_mod

USE Variablen
USE sync_mod
USE wang_johnson_mod
USE derivatives_y_mod
USE derivatives_y2_mod
CONTAINS

subroutine Dphi_Dz_p(phi,sphi,qi_i,qi_nr,qi_r)

	IMPLICIT NONE

	INTEGER::h,i,j,k,image(3)
	REAL*8::phi(-1:m_local+2,-1:n_local+2,-1:o_local+2),sphi(1:m_local,1:n_local,1:o_local)
	real*8::saver(-1:m_local+2,-1:n_local+2,-1:o_local+2)
	REAL*8::qi_i,qi_r,qi_nr,rel
! save phi
saver(-1:m_local+2,-1:n_local+2,-1:o_local+2)=&
phi(-1:m_local+2,-1:n_local+2,-1:o_local+2)
! transformation der Koordinaten
! x ---> x
! y ---> z
! z ---> y
sphi(1:m_local,1:n_local,1:o_local)=phi(1:m_local,1:n_local,1:o_local)
!sphi(1:m_local,1:n_local,1)=phi(1:m_local,1:n_local,o_local-1)
!sphi(1:m_local,1:n_local,2)=phi(1:m_local,1:n_local,o_local)
do I=1,m_local
	do j=1,n_local
		do h=1,o_local
			phi(i,h,j)=sphi(i,j,h)		
		end do
	end do
end do


CALL Dphi_Dy_2(phi,sphi,qi_i,qi_nr,qi_r)


phi(1:m_local,1:n_local,1:o_local)=sphi(1:m_local,1:n_local,1:o_local)
sync all
! Rücktransformation der Koordinaten
! x ---> x
! y ---> z
! z ---> y

do I=1,m_local
	do j=1,n_local
		do h=1,o_local
			sphi(i,h,j)=phi(i,j,h)	
		end do
	end do
end do
!set derivative to zero for testing
! sphi(1:m_local,1:n_local,1:o_local)=0.0
 
phi(-1:m_local+2,-1:n_local+2,-1:o_local+2)=&
saver(-1:m_local+2,-1:n_local+2,-1:o_local+2)
sync all
	RETURN
end subroutine

end module
