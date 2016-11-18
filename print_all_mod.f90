module print_all_mod
USE Variablen

contains

!! Ausgabe des Gesamtfeldes
SUBROUTINE PRINT_ALL

IMPLICIT NONE
	INTEGER :: i,j,h,l,k,f,d
	real:: test,x,y,z,length_x,length_y,length_z
	CHARACTER(21)::dateiname_rho_u,dateiname_rho_v,dateiname_rho_r,dateiname_energ,dateiname_vel_u,dateiname_vel_v
	CHARACTER(21)::dateiname_press,dateiname_tempr,dateiname_konti,dateiname_ve_du,dateiname_ve_dv,dateiname_poiss,dateiname_resul

length_x=  1 
length_y= 1
length_z=  1

	IF(m_pos==1.AND.n_pos==1.AND.o_pos==1)THEN

			dateiname_resul='result.dat'
OPEN(unit=33,file=dateiname_resul)
		

test=999.9
	
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
						WRITE(33,997)x,' ',y,' ',z,' ',rho(i,j,d)[h,l,f],' ',&
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


END SUBROUTINE


end module
