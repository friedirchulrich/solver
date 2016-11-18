module TEMP_FELD_MOD
! in dieser Subroutine wird das Temperaturfeld aktualisiert

USE Variablen

contains

subroutine TEMP_FELD
implicit none
real:: konstant

konstant=717.5
! evtl E/rho?
! T= 1/c_v*(E-0.5*(u^2+v^2+w^2))

T(1:m_local,1:n_local,1:o_local)=1./konstant*(E(1:m_local,1:n_local,1:o_local)&
/rho(1:m_local,1:n_local,1:o_local)-&
0.5*(u(1:m_local,1:n_local,1:o_local)**2.&
+v(1:m_local,1:n_local,1:o_local)**2.+w(1:m_local,1:n_local,1:o_local)**2.))



end subroutine

end module
