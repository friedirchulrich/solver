module bc_mod
USE Variablen
USE sync_mod

CONTAINS

subroutine bc(dx,dy)
integer ::i,j
real*8::p_0, gamma,rho_0,dummy,dx,dy
p_0=101325.0
rho_0=1.293
gamma=1.4

i=3
if(i==1)then
! j=2
!if(j==1)then

! Platte x=1

u(1,1:n_local,-1:o_local+2)=u_bc(1,1:n_local,-1:o_local+2)
v(1,1:n_local,-1:o_local+2)=v_bc(1,1:n_local,-1:o_local+2)
w(1,1:n_local,-1:o_local+2)=0.0

T(1,-1:n_local+2,-1:o_local+2)=T_bc(1,-1:n_local+2,-1:o_local+2)

rho(1,-1:n_local+2,-1:o_local+2)=rho_bc(1,-1:n_local+2,-1:o_local+2)

rho_u(1,1:n_local,-1:o_local+2)=u(1,1:n_local,-1:o_local+2)&
*rho(1,1:n_local,-1:o_local+2)
rho_v(1,1:n_local,-1:o_local+2)=v(1,1:n_local,-1:o_local+2)&
*rho(1,1:n_local,-1:o_local+2)
rho_w(1,1:n_local,-1:o_local+2)=0.0

!p(1,2:n_local,-1:o_local+2)=&
!(300.*p(1+1,2:n_local,-1:o_local+2)&
!-300.*p(1+2,2:n_local,-1:o_local+2)+200.*p(1+3,2:n_local,-1:o_local+2)&
!-75.*p(1+4,2:n_local,-1:o_local+2)+ 12.*p(1+5,2:n_local,-1:o_local+2))/137.

rho_E(1,-1:n_local+2,-1:o_local+2)=rho(1,-1:n_local+2,-1:o_local+2)**2*(&
0.5*(u(1,-1:n_local+2,-1:o_local+2)&
**2+v(1,-1:n_local+2,-1:o_local+2)**2&
+w(1,-1:n_local+2,-1:o_local+2)**2.)+T(1,-1:n_local+2,-1:o_local+2)*717.5)

!do i=1,n_local
!do j=1,o_local
!p(1,i,j)=(300.*p(2,i,j)-300.*p(3,i,j)+200.*p(4,i,j)&
!-75.*p(5,i,j)+ 12.*p(6,i,j))/137.
!end do
!end do
!p(1,-1:n_local+2,-1:o_local+2)=101325.0
!-i*(i-(n_local))/n_local*
!	v(-1:1,-1:n_local+2,-1:o_local+2)=0.00
!	w(-1:1,-1:n_local+2,-1:o_local+2)=0.00
!rho(-1:1,-1:n_local+2,-1:o_local+2)=rho_0
!	T(-1:1,2:n_local-1,-1:o_local+2)=273*(1+(gamma-1)/2*(Ma)**2)**(gamma/(gamma-1))
!     rho(-1:1,2:n_local-1,-1:o_local+2)=p(-1:1,2:n_local-1,-1:o_local+2)/(287.0*&
!T(-1:1,2:n_local-1,-1:o_local+2))
!	p(1:1,-1:n_local+2,-1:o_local+2)=p(2:2,-1:n_local+2,-1:o_local+2)

!p(-1:1,2:n_local-1,-1:o_local+2)=101325.0*1.001

! Platte x=M

p(m_local,2:n_local,-1:o_local+2)=&
(300.*p(m_local-1,2:n_local,-1:o_local+2)&
-300.*p(m_local-2,2:n_local,-1:o_local+2)+200.*p(m_local-3,2:n_local,-1:o_local+2)&
-75.*p(m_local-4,2:n_local,-1:o_local+2)+ 12.*p(m_local-5,2:n_local,-1:o_local+2))/137.

!rho(m_local,2:n_local,-1:o_local+2)=p(m_local,2:n_local,-1:o_local+2)/(287.0*&
!T(m_local,2:n_local,-1:o_local+2))

!rho_E(m_local,-1:n_local+2,-1:o_local+2)=rho(m_local,-1:n_local+2,-1:o_local+2)**2*(&
!0.5*(u(m_local,-1:n_local+2,-1:o_local+2)&
!**2+v(m_local,-1:n_local+2,-1:o_local+2)**2&
!+w(m_local,-1:n_local+2,-1:o_local+2)**2.)+T(m_local,-1:n_local+2,-1:o_local+2)*717.5)


! freestream
!call char_bc(rho(1:m_local,1:n_local,1:o_local)&
!,rho_v(1:m_local,1:n_local,1:o_local),&
!rho_u(1:m_local,1:n_local,1:o_local),&
!rho_w(1:m_local,1:n_local,1:o_local),&
!E(1:m_local,1:n_local,1:o_local),&
!T(1:m_local,1:n_local,1:o_local),&
!p(1:m_local,1:n_local,1)&
!,0.3,717.5,gamma,278.0,dy,dx,&         
!m_local,n_local,1,o_local)



! Platte y=0
	rho_u(2:m_local,-1:1,-1:o_local+2)=0.0
	rho_v(2:m_local,-1:1,-1:o_local+2)=0.0
	rho_w(2:m_local,-1:1,-1:o_local+2)=0.0
	u(2:m_local,-1:1,-1:o_local+2)=0.0
	v(2:m_local,-1:1,-1:o_local+2)=0.0
	w(2:m_local,-1:1,-1:o_local+2)=0.0
	T(2:m_local,-1:1,-1:o_local+2)=T_bc(1,1,1)
	rho(2:m_local-1,-1:1,-1:o_local+2)=rho_bc(1,1,1)

p(m_local,1,-1:o_local+2)=(300.*p(m_local-1,1,-1:o_local+2)&
-300.*p(m_local-2,1,-1:o_local+2)+200.*p(m_local-3,1,-1:o_local+2)&
-75.*p(m_local-4,1,-1:o_local+2)+ 12.*p(m_local-5,1,-1:o_local+2))/137.

rho(m_local,1,-1:o_local+2)=p(m_local,1,-1:o_local+2)/(287.0*&
T(m_local,1,-1:o_local+2))


rho_E(2:m_local,-1:n_local+2,-1:o_local+2)=rho(2:m_local,-1:n_local+2,-1:o_local+2)**2*(&
0.5*(u(2:m_local,-1:n_local+2,-1:o_local+2)&
**2+v(2:m_local,-1:n_local+2,-1:o_local+2)**2&
+w(2:m_local,-1:n_local+2,-1:o_local+2)**2.)+T(2:m_local,-1:n_local+2,-1:o_local+2)*717.5)



! Platte y=n
! 	rho_u(2:m_local,n_local:n_local+2,-1:o_local+2)=u_bc(1,n_local,1)&
!*rho(2:m_local,n_local:n_local+2,-1:o_local+2)
!u(2:m_local,n_local:n_local+2,-1:o_local+2)=u_bc(1,n_local,1)
!	rho_v(2:m_local,n_local,1:o_local)=v_bc(2:m_local,n_local,1:o_local)&
!*rho(2:m_local,n_local,1:o_local)
!rho_v(1:m_local,n_local,1:o_local)=(0.75*rho_v(1:m_local,n_local-1,1:o_local)+&
!0.5*rho_v(1:m_local,n_local-2,1:o_local)+0.25*rho_v(1:m_local,n_local-3,1:o_local))/1.5
!rho_v(2:m_local+2,n_local:n_local+2,-1:o_local+2)=0.00
!	v(2:m_local+2,n_local:n_local+2,-1:o_local+2)=v_bc(1,n_local,1)
!rho_w(2:m_local+2,n_local:n_local+2,-1:o_local+2)=0.00
!T(2:m_local+2,n_local:n_local+2,-1:o_local+2)=T_bc(1,n_local,1)
!	rho(2:m_local+2,n_local:n_local+2,-1:o_local+2)=rho_bc(1,n_local,1)
	
! Platte z=0
!	rho_u(1:m_local+2,-1:n_local+2,-1:1)=0.0
!	rho_v(1:m_local,1:n_local,1)=0.0
!	rho_w(1:m_local+2,-1:n_local+2,-1:1)=0.00
!	T(1:m_local+2,-1:n_local+2,-1:1)=T_bc(1,1,1)
!	rho(1:m_local+2,-1:n_local+2,-1:1)=rho_bc(1,1,1)
!	p(-1:m_local+2,-1:n_local+2,-1:1)=101325.0*1.2
!w(-1:m_local+2,-1:n_local+2,-1:1)=0.00
! Platte z=O_local
!	rho_u(1:m_local+2,-1:n_local+2,o_local)=0.0
!	rho_v(1:m_local,1:n_local,o_local)=0.0
!	rho_w(1:m_local+2,-1:n_local+2,o_local)=0.00
!	T(1:m_local+2,-1:n_local+2,o_local)=T_bc(1,1,1)
!	rho(1:m_local+2,-1:n_local+2,o_local)=rho_bc(1,1,1)
!	p(-1:m_local+2,-1:n_local+2,o_local:o_local+2)=101325.0-101325.0*0.1
!rho_w(-1:m_local+2,-1:n_local+2,o_local:o_local+2)=0.00
!w(1:m_local+2,1:n_local+2,-1:o_local+2)=0.0

endif

if(i==2)then
!v(1,1:n_local,-1:o_local+2)=v(m_local,1:n_local,-1:o_local+2)
u(1,1:n_local,-1:o_local+2)=u(m_local,1:n_local,-1:o_local+2)

!w(1,1:n_local,-1:o_local+2)=w(m_local,1:n_local,-1:o_local+2)
!u(1:m_local,1:n_local,-1:o_local+2)=1.0
!rho_u(1,1:n_local,-1:o_local+2)=u(1,1:n_local,-1:o_local+2)*&
!rho(1,1:n_local,-1:o_local+2)
!rho_E(1,-1:n_local+2,-1:o_local+2)=rho(1,-1:n_local+2,-1:o_local+2)**2*(&
!0.5*(u(1,-1:n_local+2,-1:o_local+2)&
!**2+v(1,-1:n_local+2,-1:o_local+2)**2&
!+w(1,-1:n_local+2,-1:o_local+2)**2)+T(1,-1:n_local+2,-1:o_local+2)*717.5)

!p(1:m_local+2,1:n_local+2,-1:o_local+2)=1.293*287.0*273.0
v(1:m_local,1,-1:o_local+2)=v(1:m_local,n_local,-1:o_local+2)
!T(1:m_local+2,1:n_local+2,-1:o_local+2)=273.0
!v(1:m_local,1,-1:o_local+2)=v(1:m_local,n_local,-1:o_local+2)
!rho_u(1:m_local,1,-1:o_local+2)=u(1:m_local,1,-1:o_local+2)*&
!rho(1:m_local,1,-1:o_local+2)
!rho_E(1:m_local,1,-1:o_local+2)=rho(1:m_local,1,-1:o_local+2)**2*(&
!0.5*(u(1:m_local,1,-1:o_local+2)&
!**2+v(1:m_local,1,-1:o_local+2)**2&
!+w(1:m_local,1,-1:o_local+2)**2.)+T(1:m_local,1,-1:o_local+2)*717.5)
endif

if(i==3)then
!v(1,1:n_local,-1:o_local+2)=v(m_local,1:n_local,-1:o_local+2)
u(1,1:n_local,-1:o_local+2)=u(m_local,1:n_local,-1:o_local+2)

!u(1:m_local,1,-1:o_local+2)=0.0
!u(1:m_local,1,-1:o_local+2)=1
!u(1:m_local,n_local,-1:o_local+2)=0.0


!p(1:m_local+2,1:n_local+2,-1:o_local+2)=1.293*287.0*273.0
v(1:m_local,1,-1:o_local+2)=v(1:m_local,n_local,-1:o_local+2)
!v(1,1:n_local,-1:o_local+2)=m_local
!v(m_local,1:n_local,-1:o_local+2)=1

!T(1:40,1:40,-1:o_local+2)=273.0



!w(1:m_local,1:n_local,-1:o_local+2)=0
!v(1,1:n_local,-1:o_local+2)=1
!v(1:m_local,1:n_local,-1:o_local+2)=m_local
!T(1:m_local+2,1:n_local+2,-1:o_local+2)=273.0
!v(1:m_local,1,-1:o_local+2)=v(1:m_local,n_local,-1:o_local+2)
rho_u(1:m_local,1,-1:o_local+2)=u(1:m_local,1,-1:o_local+2)*&
rho(1:m_local,1,-1:o_local+2)
rho_v(1:m_local,1,-1:o_local+2)=v(1:m_local,1,-1:o_local+2)*&
rho(1:m_local,1,-1:o_local+2)
rho_w(1:m_local,1,-1:o_local+2)=w(1:m_local,1,-1:o_local+2)*&
rho(1:m_local,1,-1:o_local+2)

!rho_E(1:m_local,1,-1:o_local+2)=rho(1:m_local,1,-1:o_local+2)**2*(&
!0.5*(u(1:m_local,1,-1:o_local+2)&
!**2+v(1:m_local,1,-1:o_local+2)**2&
!+w(1:m_local,1,-1:o_local+2)**2.)+T(1:m_local,1,-1:o_local+2)*717.5)
endif
SYNC ALL



end subroutine

end module
