module TEMP_FELD_MOD
! in dieser Subroutine wird das Temperaturfeld aktualisiert

USE Variablen

contains

subroutine TEMP_FELD
implicit none
real*8:: konstant
integer:: h,i,j

! 1/cv
konstant=0.00139470013
! evtl E/rho?
! T= 1/c_v*(E/rho-0.5*(u^2+v^2+w^2))

T(1:m_local,1:n_local,1:o_local)=(E(1:m_local,1:n_local,1:o_local)&
-0.5*(rho_u(1:m_local,1:n_local,1:o_local)**2&
+rho_v(1:m_local,1:n_local,1:o_local)**2+rho_w(1:m_local,1:n_local,1:o_local)**2)&
/rho(1:m_local,1:n_local,1:o_local))/rho(1:m_local,1:n_local,1:o_local)*konstant


!do I=1,m_local

!print *, T(i,20,20),(T(i,20,20)-273.0)/273.0

!end do


end subroutine

end module
