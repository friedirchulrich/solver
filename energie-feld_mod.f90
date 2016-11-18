module EN_FELD_MOD
! in dieser Subroutine wird das Geschwindigkeitsfeld aktualisiert
USE Variablen

contains

subroutine EN_FELD


E(1:m_local,1:n_local,1:o_local)=rho_E(1:m_local,1:n_local,1:o_local)&
/rho(1:m_local,1:n_local,1:o_local)
end subroutine

end module
