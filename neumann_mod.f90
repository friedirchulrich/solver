module neumann_mod
USE Variablen
contains

subroutine bcneu

!bc_2
!delta_u(m_local:m_local,2:n_local,1:o_local)=0
!delta_v(m_local:m_local,2:n_local,1:o_local)=0
!delta_w(m_local:m_local,2:n_local,1:o_local)=0

!delta_u(1:1,2:n_local,1:o_local)=0
!delta_v(1:1,2:n_local,1:o_local)=0
!delta_w(1:1,2:n_local,1:o_local)=0


end subroutine


end module
