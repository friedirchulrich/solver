module solver_mod
USE Variablen
USE xmom_mod
USE ymom_mod
USE zmom_mod
USE energy_mod
USE check_conti_mod
USE conti_mod
USE sync_mod
USE sutherland_mod
USE spannungstensor_mod

contains

subroutine solver
integer:: i,j,d


	!wichtig
	CALL SUTHERLAND
!print *, "sutherland"
	CALL tensor
!print *, "tensor"
!	CALL GC_EXCHANGE_P 


		!! Lösung mit kompakten Diff.			
!		CALL CHECK_CONTI

		CALL CONTI 
!print *, "conti"
		CALL XMOM
!print *, "xmom"
		CALL YMOM
!print *, "ymom"
		CALL ZMOM
!print *, "zmom"
		CALL ENERGY
!print *, "energy"

!Überschreiben der neuen Zustandsvariablen
! in die Alten
rho(1:m_local,1:n_local,1:o_local)=&
rho_neu(1:m_local,1:n_local,1:o_local)

rho_u(1:m_local,1:n_local,1:o_local)=&
rho_u_neu(1:m_local,1:n_local,1:o_local)

rho_v(1:m_local,1:n_local,1:o_local)=&
rho_v_neu(1:m_local,1:n_local,1:o_local)

rho_w(1:m_local,1:n_local,1:o_local)=&
rho_w_neu(1:m_local,1:n_local,1:o_local)

rho_E(1:m_local,1:n_local,1:o_local)=&
rho_E_neu(1:m_local,1:n_local,1:o_local)

end subroutine

end module
