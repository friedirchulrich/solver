module DRUCK_FELD_MOD
! in dieser Subroutine wird das Druckfeld aktualisiert
USE Variablen

contains

subroutine DRUCK_FELD

! p=rho*R*T (ideales Gas-Gesetz)
p(1:m_local,1:n_local,1:o_local)=&
(rho(1:m_local,1:n_local,1:o_local)*287.0)*T(1:m_local,1:n_local,1:o_local)



end subroutine

end module
