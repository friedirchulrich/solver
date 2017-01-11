program main
!**********************************
!* Navier-Stokes and Euler 3D 
!* Co-Array Fortran Solver
!* Friedrich Ulrich 2016
!* aufbauend auf
!* inin
!**********************************

! load the required modules
USE inputfile
USE variablen
USE INI_MOD
USE energy_tot_mod
USE Solver_mod
USE vel_feld_mod
USE druck_feld_mod
USE temp_feld_mod
USE PRINT_ALL_MOD
USE bc_mod
USE en_feld_mod
use TEST_MOD
! variables

	! Variable zur Modus-Selektion (Kanalströmung, Taylor-Green etc...)
	INTEGER :: Mode
!	REAL::length_x=  0.001 ,length_y= 1.110446458202763e-04 ,length_z=  0.00000001
!	REAL*8::length_x=3e-2 ,length_y=0.05e-3,length_z= 1e-2
!REAL::length_x=  0.1 ,length_y= 8.72e-5,length_z=  0.1
REAL*8::length_x=  1 ,length_y= 1,length_z=  1

	!! Laufvariablen
	INTEGER :: h,i,j,it,itmax,d
	! Abruchkrit
	REAL*8::energy_temp,epsilon_energy,cfl,Re,simtime,olddt,v_0,w_0,u_help


!program start

!intro
!choose image one to print only one text
IF (this_image()==1) THEN

	print *, "********************************"
	Print *, "*     _                 _      *"
	print *, "*   ___\ A E R 2 0 1 6 /___    *"
	print *, "*  |CAF - NSG - 3D - SOLVER|   *"
	Print *, "*                              *"
	print *, "********************************"

END IF

! read input file
! Die Variable T_p ist die Anzahl der Punkte T in z Richtung
! diese muss aber mit T-p bezeichnet werden, da T 
! für die Temperatur vorbehalten ist


call inputfileread(M,N,O,R,S,T_p,i,itmax,cfl,Re,Ma,mode)


!Parallel-System erstellen
	myrank=this_image()
	numprocs=num_images()
	m_local=M/R
	n_local=N/S
	o_local=O/T_p
	m_kop=m_local/S
	n_kop=n_local/R
	o_kop=o_local/T_p

!	nlocal=200


	! data 1:m_local, 2 ghostcells
	ALLOCATE(rho(-1:m_local+2,-1:n_local+2,-1:m_local+2)[1:R,1:S,1:*]) 
	ALLOCATE(rho_bc(-1:m_local+2,-1:n_local+2,-1:m_local+2)[1:R,1:S,1:*])
	ALLOCATE(rho_u(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(rho_v(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(rho_w(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
        ALLOCATE(rho_E(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(rho_neu(-1:m_local+2,-1:n_local+2,-1:m_local+2)[1:R,1:S,1:*]) 
	ALLOCATE(rho_u_neu(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(rho_v_neu(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(rho_w_neu(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
        ALLOCATE(rho_E_neu(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
        ALLOCATE(E_plus_p(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])

	ALLOCATE(p(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(E(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(Epu(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
        ALLOCATE(Epv(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
        ALLOCATE(Epw(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])

	ALLOCATE(u(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(u_bc(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(v_bc(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])

	ALLOCATE(v(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(w(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(T(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(T_bc(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(vis(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(mu(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(k(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(dphi(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(dphi1(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(dphi4(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(dphi2(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(dphi3(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])

	ALLOCATE(solve(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
   	ALLOCATE(E_k(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
   	ALLOCATE(a(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(uinit(1:nlocal,1:nlocal,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(rinit(1:nlocal,1:nlocal,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(vinit(1:nlocal,1:nlocal,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(tinit(1:nlocal,1:nlocal,1:o_local)[1:R,1:S,1:*])
	!!
	ALLOCATE(du(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(dv(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(drudx(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(drvdy(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(dudx(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(dvdx(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(dvdy(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(dudy(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(poisson(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(DDX_DDY(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])

	!Gleichungslöser in x

	ALLOCATE(b_x(1:m_local),c_x(1:m_local))
	ALLOCATE(a_x(1:m_local+1)[1:R,1:S,1:*],f_x(1:m_local+1)[1:R,1:S,1:*],g_x(1:m_local+1)[1:R,1:S,1:*])
	ALLOCATE(y_x(1:m_local+1,1:n_local,1:2)[1:R,1:S,1:*],x_x(0:m_local,1:n_local,1:2)[1:R,1:S,1:*])!????????????????????????
	ALLOCATE(a_x_s(1:R)[1:R,1:S,1:*],f_x_s(1:R)[1:R,1:S,1:*],g_x_s(1:R)[1:R,1:S,1:*])
	ALLOCATE(y_x_s(1:R,1:n_kop,1:2)[1:R,1:S,1:*],x_x_s(1:R,1:n_kop,1:2)[1:R,1:S,1:*])
	ALLOCATE(alpha_x(1:n_local,1:2)[1:R,1:S,1:*],beta_x(1:n_local)[1:R,1:S,1:*],z_x(1:m_local))	

	!Gleichungslöser in y
	ALLOCATE(b_y(1:n_local),c_y(1:n_local))
	ALLOCATE(a_y(1:n_local+1)[1:R,1:S,1:*],f_y(1:n_local+1)[1:R,1:S,1:*],g_y(1:n_local+1)[1:R,1:S,1:*])
	ALLOCATE(y_y(1:m_local,1:n_local+1,1:2)[1:R,1:S,1:*],x_y(1:m_local,0:n_local,1:2)[1:R,1:S,1:*])!????????????????????????
	ALLOCATE(a_y_s(1:S)[1:R,1:S,1:*],f_y_s(1:S)[1:R,1:S,1:*],g_y_s(1:S)[1:R,1:S,1:*])
	ALLOCATE(y_y_s(1:m_kop,1:S,1:2)[1:R,1:S,1:*],x_y_s(1:m_kop,1:S,1:2)[1:R,1:S,1:*])
	ALLOCATE(alpha_y(1:m_local,1:2)[1:R,1:S,1:*],beta_y(1:m_local)[1:R,1:S,1:*],z_y(1:n_local))

	!Gleichungslöser i z
	ALLOCATE(b_z(1:o_local),c_z(1:o_local))
	ALLOCATE(a_z(1:o_local+1)[1:R,1:S,1:*],f_z(1:o_local+1)[1:R,1:S,1:*],g_z(1:o_local+1)[1:R,1:S,1:*])
	ALLOCATE(y_z(1:m_local,1:n_local+1,1:2)[1:R,1:S,1:*],x_z(1:m_local,0:n_local,1:2)[1:R,1:S,1:*])!????????????????????????
	ALLOCATE(a_z_s(1:S)[1:R,1:S,1:*],f_z_s(1:S)[1:R,1:S,1:*],g_z_s(1:S)[1:R,1:S,1:*])
	ALLOCATE(y_z_s(1:m_kop,1:S,1:2)[1:R,1:S,1:*],x_z_s(1:m_kop,1:S,1:2)[1:R,1:S,1:*])
	ALLOCATE(alpha_z(1:m_local,1:2)[1:R,1:S,1:*],beta_z(1:m_local)[1:R,1:S,1:*],z_z(1:o_local))
	
	!! Spannungstensor
	ALLOCATE(tau_xx(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(tau_yy(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(tau_zz(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(tau_xy(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(tau_xz(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(tau_zy(-1:m_local+2,-1:n_local+2,-1:o_local+2)[1:R,1:S,1:*])
	ALLOCATE(delta_u(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(delta_v(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])
	ALLOCATE(delta_w(1:m_local,1:n_local,1:o_local)[1:R,1:S,1:*])

    !! Koordinaten des ausführenden Image bestimmen

	DO h=1,R
		DO i=1,S
			DO j=1,T_p
				IF(image_index(rho,[h,i,j])==myrank)THEN
					m_pos=h
					n_pos=i
					o_pos=j
					EXIT
				ENDIF
			ENDDO
		ENDDO
	ENDDO
    print *, 'myrank: ',myrank, 'm_pos:',m_pos,'n_pos: ',n_pos,'o_pos:'&
 ,o_pos,'procCount:',numprocs


	!! Synchronisierungs-Variablen (sync_n_up,sync_m_up) belegen
	IF(R==1)THEN
		sync_m_up=myrank
		sync_m_down=myrank
	ELSE
		IF(m_pos<R)sync_m_up=image_index(rho,[m_pos+1,n_pos,o_pos])
		IF(m_pos>1)sync_m_down=image_index(rho,[m_pos-1,n_pos,o_pos])
	ENDIF
	IF(S==1)THEN
		sync_n_up=myrank
		sync_n_down=myrank
	ELSE
		IF(n_pos<S)sync_n_up=image_index(rho,[m_pos,n_pos+1,o_pos])
		IF(n_pos>1)sync_n_down=image_index(rho,[m_pos,n_pos-1,o_pos])
	ENDIF
	IF(T_p==1)THEN
		sync_o_up=myrank
		sync_o_down=myrank
	ELSE
		IF(o_pos<T_p)sync_o_up=image_index(rho,[m_pos,n_pos,o_pos+1])
		IF(o_pos>1)sync_o_down=image_index(rho,[m_pos,n_pos,o_pos-1])
	ENDIF

	ALLOCATE(sync_m_array(1:R))
	DO i=1,R
		sync_m_array(i)=image_index(rho,[i,n_pos,o_pos])
	ENDDO
	ALLOCATE(sync_n_array(1:S))
	DO i=1,S
		sync_n_array(i)=image_index(rho,[m_pos,i,o_pos])
	ENDDO	
	ALLOCATE(sync_o_array(1:T_p))
	DO i=1,T_p
		sync_o_array(i)=image_index(rho,[m_pos,n_pos,i])
	ENDDO
	image_m_1=image_index(rho,[1,n_pos,o_pos])
	image_n_1=image_index(rho,[m_pos,1,o_pos])
	image_o_1=image_index(rho,[m_pos,n_pos,1])
	image_m_max=image_index(rho,[R,n_pos,o_pos])
	image_n_max=image_index(rho,[m_pos,S,o_pos])
	image_o_max=image_index(rho,[m_pos,n_pos,T_p])

	!Fall erstellen	
  	dx = length_x / (M-1)
   	dy = length_y / (N-1)
	dz = length_z / (O-1) 
!

! u_0=Ma*sqrt(gamma*R*T_0)
!	u_0=Ma*sqrt(1.4*287*273)
!u_0=0.11495E+03	
!v_0=0.1
!w_0=0.1
u_0=1	
v_0=1
w_0=0.1
print *, " dx,dy,dz:",dx,dy,dz
   	dt = cfl*min(dx/abs(u_0),dy/abs(u_0),dz/abs(u_0)) !v_0=u_0
!dt=1e-9
print *, "delta t ", dt
olddt=dt
simtime=0.0
mu_0=1.458*10.**(-6.)
	Re=u_0*length_x*mu_0
print *,"dt,Re, u0:", dt,Re, u_0
	print *,'dt::::::::::::::',dt,'dx/abs(u_0):',dx/abs(u_0),'dy/abs(u_0):',dy/abs(u_0),&
'dz/abs(u_0):',dz/abs(u_0)
	epsilon_energy=1e-06

	IF(m_pos==1.AND.n_pos==1.AND.o_pos==1)PRINT *,'dx: ',dx,' dy: ',dy,' dz: ',dz,' dt: ',dt,'Re: ',Re
	!kompakte Differenzen
	CALL INITIAL_COEFF
!	!Anfangswerte	
print *, "reading initial condition"
	  CALL INITIAL_CONDITION(mode,dx,dy)
	SYNC ALL
! uncomment falls ausgabe der Anfangswerte gewünscht werden.
!	start=.true. !variable für PRINT-Routinen, true->Ausgabe Anfangswerte, false->Endwerte
!	CALL PRINT_ALL

!	CALL CALC_ENERGY
	!energy_temp=total_energy
	IF(m_pos==1.AND.n_pos==1.AND.o_pos==1)PRINT *,'energy: ',total_energy
	SYNC ALL

!	CALL CALC_ENERGY_KIN
        IF(m_pos==1.AND.n_pos==1.AND.o_pos==1)PRINT *,'kin.energy: ',E_k_total
	SYNC ALL


! starting the calculation
	IF(m_pos==1.AND.n_pos==1.AND.o_pos==1)PRINT *,'Berechnung startet!'

it=0

! CALL BC
	! Geschwindigkeitsfelder werden aktualisiert
!	CALL VEL_FELD
!SYNC ALL
! CALL BC

	!Energiefelder werden aktualisiert
!CALL EN_FELD
!SYNC ALL
! CALL BC
	!Temperatur aus Energie- und GGeschwindigkeitsfeld 	
!	CALL TEMP_FELD
!SYNC ALL
! CALL BC
	! Druck aus dem idealen Gas-Gesetz
!	CALL DRUCK_FELD

CALL BC(dx,dy)
simtime=0
DO WHILE(it<itmax)

	it=it+1





	! Energie, Implsgleichungen etc...	
	CALL Solver


 CALL BC(dx,dy)


print *, "solver endend"
	! Geschwindigkeitsfelder werden aktualisiert
	CALL VEL_FELD
SYNC ALL
CALL BC(dx,dy)
	!Energiefelder werden aktualisiert
CALL EN_FELD
SYNC ALL
CALL BC(dx,dy)
!do i=1,m_local
!do j=1,n_local
!do d=1,o_local

!if(T(i,j,d).NE.273 )then
!print *, T(i,j,d),i,j,d
!endif
!end do
!end do
!end do

	!Temperatur aus Energie- und GGeschwindigkeitsfeld 	
CALL TEMP_FELD
SYNC ALL
CALL BC(dx,dy)


! Druck aus dem idealen Gas-Gesetz
CALL DRUCK_FELD


CALL BC(dx,dy)

!rho_E(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)**2&
!*(0.5*(u(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
!**2+v(-1:m_local+2,-1:n_local+2,-1:o_local+2)**2&
!+w(-1:m_local+2,-1:n_local+2,-1:o_local+2)**2.)+T(-1:m_local+2,-1:n_local+2,-1:o_local!+2)*717.5)

!rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)=p(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
!/287.0/T(-1:m_local+2,-1:n_local+2,-1:o_local+2)

!rho_u(-1:m_local+2,-1:n_local+2,-1:o_local+2)=u(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
!*rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)

!rho_v(-1:m_local+2,-1:n_local+2,-1:o_local+2)=v(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
!*rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)

!rho_w(-1:m_local+2,-1:n_local+2,-1:o_local+2)=w(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
!*rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)




		! Abbruchkrit
		! energy_temp=total_energy
	!	CALL CALC_ENERGY_KIN
		CALL CALC_ENERGY
		IF(m_pos==1.AND.n_pos==1.AND.o_pos==1)THEN
			PRINT *,'iterations: ',it,'total_energy: ',total_energy
		ENDIF
!if(it==10)then
!cfl=0.1
!endif

!if(it==250)then
!cfl=0.01
!endif

if(mod(it,1)==0)then
print *, "new cfl-time"

do i=1,m_local
do j=1,n_local
do h=1,o_local
u_0=max(u_0,u(i,j,h))
enddo
enddo
enddo
 print *, "max u:",u_0
do i=1,m_local
do j=1,n_local
do h=1,o_local
v_0=max(v_0,v(i,j,h))
enddo
enddo
enddo
 print *, "max v:",v_0
do i=1,m_local
do j=1,n_local
do h=1,o_local
w_0=max(w_0,w(i,j,h))
enddo
enddo
enddo
  print *, "max w:",w_0
u_help=max(u_0,v_0,w_0)
dt = cfl*min(dx/abs(u_help),dy/abs(u_help),dz/abs(u_help))
print *, "new dt:", dt
print *, "Abweichung dt:",1.0-(olddt-dt)/olddt
print *, "**********************************"
end if


simtime=simtime+dt
ENDDO


IF(m_pos==1.AND.n_pos==1.AND.o_pos==1)THEN
			PRINT *,"writing result.dat file..."
		ENDIF

 CALL PRINT_ALL

IF(m_pos==1.AND.n_pos==1.AND.o_pos==1)THEN
			PRINT *,"result.dat file complete"
		ENDIF
!test routine für ableitungen
! call test(dx,dy,dz)

! delocalisieren

	DEALLOCATE(rho,rho_u,rho_v,rho_w,rho_E,E,p,Epu,Epv,u,v,T,E_k,vis,mu,a)

	DEALLOCATE(dphi,dphi1,dphi2,dphi3,dphi4,solve,sync_m_array,sync_n_array,du,dv,uinit,rinit,vinit,tinit)!,arr1,arr2,arr3,arr4,arr5,arr6
	DEALLOCATE(a_x,b_x,c_x,y_x,f_x,g_x,x_x,a_x_s,f_x_s,g_x_s,y_x_s,x_x_s,z_x,alpha_x,beta_x)
	DEALLOCATE(a_y,b_y,c_y,y_y,f_y,g_y,x_y,a_y_s,f_y_s,g_y_s,y_y_s,x_y_s,z_y,alpha_y,beta_y)

IF(m_pos==1.AND.n_pos==1.AND.o_pos==1)THEN
			print *, "end of program"
print *,"Simulationtime:",simtime
!print *,"CFL-Zahl:",cfl

		ENDIF

end program
