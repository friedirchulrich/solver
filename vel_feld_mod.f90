module VEL_FELD_MOD
! in dieser Subroutine wird das Geschwindigkeitsfeld aktualisiert
USE Variablen

contains

subroutine VEL_FELD
implicit none
!u_i=(rho*u_i)/rho
u(1:m_local,1:n_local,1:o_local)=rho_u(1:m_local,1:n_local,1:o_local)&
/rho(1:m_local,1:n_local,1:o_local)
	
v(1:m_local,1:n_local,1:o_local)=rho_v(1:m_local,1:n_local,1:o_local)&
/rho(1:m_local,1:n_local,1:o_local)

w(1:m_local,1:n_local,1:o_local)=rho_w(1:m_local,1:n_local,1:o_local)&
/rho(1:m_local,1:n_local,1:o_local)


end subroutine

end module
