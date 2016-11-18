module test_mod
USE Variablen
USE derivatives_x_mod
USE derivatives_y_mod
USE derivatives_z_mod

contains

subroutine test(dx,dy,dz)
real::dx,dy,dz,cosi(40)
integer::i,j,h
print *, " Subroutine Test is running..."
print *, dz
do I=1,m_local
	do j=1,n_local
		do h=1,o_local
			rho_u(i,j,h)=sin(2*3.141592654/(39.0)*(h-1))
			cosi(h)=cos(2*3.141592654/(39.0)*(h-1))		
		end do
	end do
end do



!CALL Dphi_Dx_p_neu(rho_u,delta_u,qxi_i,qxi_nr,qxi_r)
CALL Dphi_Dz_p(rho_u,delta_u,qzi_i,qzi_nr,qzi_r)
!delta_u(1:m_local,1,1)=delta_u(1:m_local,1,1)

do i=1,o_local
print *,rho_u(1,1,i),"Numerik:", delta_u(1,1,i)/6958.67773 , "Analytik:",cosi(i)
end do

print *, "miep"
end subroutine test

end module
