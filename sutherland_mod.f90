module sutherland_mod

USE Variablen

contains

subroutine sutherland
real:: n
integer:: i,j,h
! in dieser Subroutine wird die viskosität µ berechnet mit
! Hilfe des Sutherland Gesetzes
! und die Wärmeleitungskoefiziernt k


!Konstanten festlegen in INI_MOD
!cp=1004.5
!mu0=1.725
!Suth=110.3
!T_inf=275 !K
! Pr=0.72
! mu0=C1

mu0=1.458*10.**(-6)
Suth=110.4
 mu(1:m_local,1:n_local,1:o_local)=mu0*&
(T(1:m_local,1:n_local,1:o_local))**(3./2.)&
/(T(1:m_local,1:n_local,1:o_local)+Suth)

do I=1,m_local
	do j=1,n_local
		do h=1,o_local

			if((T(i,j,h)<110.3).AND.(T(i,j,h)>40))then
                      mu(i,j,h)=6.8070e-8 *T(i,j,h)
			endif
			if((T(i,j,h)<40))then
                      mu(i,j,h)=6.8070e-8 *T(i,j,h)
			endif

		enddo
	enddo
enddo
!print *, mu(1:10,1:10,1:10)
!read *, mu0
mu(1:m_local,1:n_local,1:o_local)=mu0
!cp=1004.5 cpo =: n
n=1004.5
Pr=0.72
k(1:m_local,1:n_local,1:o_local)=&
mu(1:m_local,1:n_local,1:o_local)*n/Pr


end subroutine

end module
