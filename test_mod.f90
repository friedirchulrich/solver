module test_mod
USE Variablen
USE derivatives_x_mod
USE derivatives_y_mod
USE derivatives_z_mod

contains

subroutine test(dx,dy,dz)
real*8::dx,dy,dz,cosi(40)
INTEGER :: i,j,h,l,k,f,d
	real*8:: tester,x,y,z,length_x,length_y,length_z
	CHARACTER(21)::dateiname_rho_u,dateiname_rho_v,dateiname_rho_r,dateiname_energ,dateiname_vel_u,dateiname_vel_v
	CHARACTER(21)::dateiname_press,dateiname_tempr,dateiname_konti,dateiname_ve_du,dateiname_ve_dv,dateiname_poiss,dateiname_resul
print *, " Subroutine Test is running..."
!print *, dz
do I=1,m_local
	do j=1,n_local
		do h=1,o_local
			!rho_u(i,j,h)=sin(2*3.141592654/(39.0)*(h-1))
			!cosi(h)=cos(2*3.141592654/(39.0)*(h-1))
T(i,j,h)=313259.65625000000		
		end do
	end do
end do

T(13,2,20)=313259.70000000

!CALL Dphi_Dx_p_neu(rho_u,delta_u,qxi_i,qxi_nr,qxi_r)

!CALL Dphi_Dz_p(T,dphi1,qzi_i,qzi_nr,qzi_r)
!CALL Dphi_Dx_p(T,dphi1,qxi_i,qxi_nr,qxi_r)
CALL Dphi_Dy_p(T,dphi1,qyi_i,qyi_nr,qyi_r)
!delta_u(1:m_local,1,1)=delta_u(1:m_local,1,1)
j=1
i=13
h=20
!do h=1,o_local
  do j=1,m_local
!	do i=1,n_local
		
!			if(dphi(i,j,h).NE.0)then
print *, dphi1(i,j,h),T(i,j,h),i,j,h
!endif		
		end do
!	end do
!end do


length_x=  1 
length_y= 1
length_z=  1

	IF(m_pos==1.AND.n_pos==1.AND.o_pos==1)THEN

			dateiname_resul='result.dat'
OPEN(unit=33,file=dateiname_resul)
		
print *, "printing file..."
tester=999.9
	
	!write results.dat header

	WRITE(33,FMT='(A16)') 'TITLE = "Result"'
	WRITE(33,FMT='(A11,A3,A1,A3,A1,A3,A1,A5,A1,A3,A1,A3,A1,A3,A1,A3,A1,A3,A1,A3)')'VARIABLES =','"x"',',','"y"',',',&
&'"z"',',','"rho"',',','"u"',',','"v"',',','"w"',',','"E"',',','"T"'&
,',','"p"'


		WRITE(33,FMT='(A8,I4,A5,I4,A4,I4,A4,I4)') 'ZONE T="',M,'", I=',M,', J=',N,', K=',O
DO f=1,T_p
	Do d=o_local,1,-1		
		DO l=S,1,-1
			DO j=n_local,1,-1
				DO h=1,R
					DO i=m_local,1,-1						
!	
		
x=(i+m_local*(R-1))*length_x
y=(j+n_local*(S-1))*length_y
z=(d+o_local*(T_p-1))*length_z


						WRITE(32,FMT='(A1)',ADVANCE='NO')' '
						WRITE(33,997)x,' ',y,' ',z,' ',dphi1(i,j,d)[h,l,f],' ',&
u(i,j,d)[h,l,f],' ',v(i,j,d)[h,l,f],' ',w(i,j,d)[h,l,f],&
' ',E(i,j,d)[h,l,f],' ',T(i,j,d)[h,l,f],' ',p(i,j,d)[h,l,f]
!print *, d,j,i,u(i,j,d)[h,l,f]
					ENDDO;
				ENDDO


			ENDDO

		ENDDO
	ENDDO
ENDDO



                CLOSE(33)
	ENDIF


999	FORMAT(F10.5)
991	FORMAT(E22.16)
997	FORMAT(F22.10,A1,F22.10,A1,F22.10,A1,F22.10,A1,F22.10,A1,F22.10,A1,F22.10,A1,F22.10,A1,F22.10,A1,F22.10)



print *, "miep"
end subroutine test

end module
