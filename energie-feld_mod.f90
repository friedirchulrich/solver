module EN_FELD_MOD
! in dieser Subroutine wird das Geschwindigkeitsfeld aktualisiert
USE Variablen

contains

subroutine EN_FELD
integer::h,i,j
real::hilf

E(1:m_local,1:n_local,1:o_local)=rho_E(1:m_local,1:n_local,1:o_local)&
/rho(1:m_local,1:n_local,1:o_local)



!hilf=1.293*(0.5*1+273*717.5)

!DO h=1,m_local
!DO i=1,n_local
!DO j=1,o_local


!if(E(h,i,j).NE.hilf)then
!print *, rho_E(h,i,j),rho(h,i,j),hilf*1.293
!endif
!end do
!end do
!end do

end subroutine

end module
