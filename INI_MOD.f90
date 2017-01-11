Module INI_MOD
USE Variablen
USE derivatives_x_mod
USE derivatives_y_mod
contains
!Anfangsbedingungen einfügen
SUBROUTINE INITIAL_CONDITION(mode,dx,dy)


IMPLICIT NONE

integer::mode
	INTEGER::h,i,j,dymax,dummy,d,l,f
	REAL*8::pi,i_n,n_pi,m_pi,h_m,xh,yi
	REAL*8::	v_init,x,y,z,dx,dy
	REAL::  u_init,ybl(10000),eta(10000),rhobl(10000),tbl(10000),ubl(10000),rubl(10000)
	REAL*8::rhobl1(10000),tbl1(10000),ubl1(10000),rubl1(10000),x_c(10000)
	REAL*8:: ybl1(10000),sc,deta,v_ref,x0
REAL*8,allocatable::ublrec(:,:), tblrec(:,:),rhoblrec(:,:),rublrec(:,:),rvblrec(:,:),rublrecsol(:,:)
character::reader,car

ALLOCATE(ublrec(N,M+10))
ALLOCATE(tblrec(N,M+10))
ALLOCATE(rhoblrec(N,M+10))
ALLOCATE(rublrec(N,M+10))
ALLOCATE(rvblrec(N,M+10))
ALLOCATE(rublrecsol(N,M+10)) 

	rho_0=1.293
	p_0=101325.0
	kappa=5./3.
	gamma=1.4
	u_0=1.
	
	rho(1:m_local,1:n_local,1:o_local)=rho_0

	E_0=rho_0*u_0**2+p_0/(kappa-1.)
	vis=0.001

! Sutherland Variablen
cp=1004.5
mu0=1.725*(10**(-5))
Suth=110.3
T_inf=273 !K
Pr=0.72


!Input-File-Realisierung

IF(mode==1)then
	! Sinus-Profil
	pi=3.14159265358979324
	n_pi=2.*pi/N

	v_init=1
	u_init=1
	DO j=1,o_local
		DO i=1,n_local
			DO h=1,m_local
				i_n=real((n_pos-1)*n_local+i-1)
				u(h,i,j)=u_init+0.1*(1.+cos(pi+i_n*n_pi))
				p(h,i,j)=p_0
				E(h,i,j)=0.5*rho_0*(u(h,i,j)**2+v(h,i,j)**2)+p(h,i,j)/(kappa-1.)
	
			ENDDO
		ENDDO	
	ENDDO
	v(1:m_local,1:n_local,1:o_local)=v_init
end if

IF(mode==2)then
	!Taylor green
	pi=3.14159265358979324
	m_pi=2.*pi/M
	n_pi=2.*pi/N
DO j=1,o_local
	DO i=1,n_local
		DO h=1,m_local
			h_m=real((m_pos-1)*m_local+h-1)
			i_n=real((n_pos-1)*n_local+i-1)
			xh=-pi+h_m*m_pi
			yi=-pi+i_n*n_pi
			!Taylor Green Vortex
			u(h,i,j)=u_0*sin(xh)*cos(yi)
			v(h,i,j)=-u_0*cos(xh)*sin(yi)
			w(h,i,j)=0			
			p(h,i,j)=p_0+(0.25*rho_0*u_0**2*(cos(2*(xh))+cos(2*(yi))))
			E(h,i,j)=0.5*rho_0*(u(h,i,j)**2+v(h,i,j)**2)+p(h,i,j)/(kappa-1.)
			rho(h,i,j)=rho_0
			rho_u(h,i,j)=rho_0*u(h,i,j)
			rho_v(h,i,j)=rho_0*v(h,i,j)
			rho_w(h,i,j)=rho_0*w(h,i,j)
			T(h,i,j)=273
		ENDDO
	ENDDO
ENDDO	
	

end if

!Kanal
if(mode==3) then
	IF(R==1.AND.S==1)THEN

            u(1:m_local,1,1:o_local)=10.0
	    v(1:m_local,1,1:o_local)=0.0
	    p(1:m_local,1,1:o_local)=p_0


            u(1:m_local,n_local,1:o_local)=10.0
            v(1:m_local,n_local,1:o_local)=0.0
            u(1:m_local,2:n_local-1,1:o_local)=10.0
!            u(1:2,30:40)=0.0
            v(1:m_local,2:n_local-1,1:o_local)=0.0
            p(1:m_local,1:n_local,1:o_local)=p_0
            E(1:m_local,1:n_local,1:o_local)=0.5*rho_0*(u(1:m_local,1:n_local,1:o_local)&
**2+v(1:m_local,1:n_local,1:o_local)**2)+p(1:m_local,1:n_local,1:o_local)/(kappa-1.)
            rho_u(1:m_local,1:n_local,1:o_local)=rho_0*u(1:m_local,1:n_local,1:o_local)
            rho_v(1:m_local,1:n_local,1:o_local)=rho_0*v(1:m_local,1:n_local,1:o_local)

	ELSE
		IF(n_pos==S)THEN

                u(1:m_local,n_local,1:o_local)=0.0
                v(1:m_local,n_local,1:o_local)=0.0
                u(1:m_local,1:n_local-1,1:o_local)=1.0
				v(1:m_local,1:n_local-1,1:o_local)=0.0
				p(1:m_local,1:n_local,1:o_local)=p_0
				E(1:m_local,1:n_local,1:o_local)=0.5*rho_0*(&
u(1:m_local,1:n_local,1:o_local)**2+v(1:m_local,1:n_local,1:o_local)**2)+&
p(1:m_local,1:n_local,1:o_local)/(kappa-1.)
				rho_u(1:m_local,1:n_local,1:o_local)=rho_0*&
u(1:m_local,1:n_local,1:o_local)
				rho_v(1:m_local,1:n_local,1:o_local)=rho_0*&
v(1:m_local,1:n_local,1:o_local)

		ELSEIF(n_pos==1)THEN

                u(1:m_local,1,1:o_local)=0.0
                v(1:m_local,1,1:o_local)=0.0
                u(1:m_local,2:n_local,1:o_local)=1.0
                v(1:m_local,2:n_local,1:o_local)=0.0
                p(1:m_local,1:n_local,1:o_local)=p_0
                E(1:m_local,1:n_local,1:o_local)=0.5*rho_0*&
(u(1:m_local,1:n_local,1:o_local)**2+v(1:m_local,1:n_local,1:o_local)**2)&
+p(1:m_local,1:n_local,1:o_local)/(kappa-1.)
                rho_u(1:m_local,1:n_local,1:o_local)=rho_0*u(1:m_local,1:n_local,1:o_local)
                rho_v(1:m_local,1:n_local,1:o_local)=rho_0*v(1:m_local,1:n_local,1:o_local)


		ELSE

		u(1:m_local,1:n_local,1:o_local)=1.0
		v(1:m_local,1:n_local,1:o_local)=0.0
		p(1:m_local,1:n_local,1:o_local)=p_0
		E(1:m_local,1:n_local,1:o_local)=0.5*rho_0*(u(1:m_local,1:n_local,1:o_local)&
**2+v(1:m_local,1:n_local,1:o_local)**2)+p(1:m_local,1:n_local,1:o_local)/(kappa-1.)
		rho_u(1:m_local,1:n_local,1:o_local)=rho_0*u(1:m_local,1:n_local,1:o_local)
		rho_v(1:m_local,1:n_local,1:o_local)=rho_0*v(1:m_local,1:n_local,1:o_local)

		ENDIF
	ENDIF	
end if

if(mode==4) then

rho_0=1.29
 	w(-1:m_local+2,-1:n_local+2,-1:o_local+2)=0.0
 	v(-1:m_local+2,-1:n_local+2,-1:o_local+2)=0.0
	u(-1:m_local+2,-1:n_local+2,-1:o_local+2)=0.1

do i=-1,o_local+2
do j=-1,n_local+2
do h=-1,m_local+2
!	p(h,j,i)=-(((1+(gamma-1)/2*(Ma)**2)**(gamma/(gamma-1)))*rho_0*287.0*273-rho_0*287.0*273)/m_local*h&
!+(1+(gamma-1)/2*(Ma)**2)**(gamma/(gamma-1))*rho_0*287.0*273
p(h,j,i)=rho_0*287.0*273
end do
end do
end do

do i=1,m_local
do h=1,o_local
do j=1,n_local
	u(i,j,h)=-(j-1)*(j-(n_local))/n_local*150
end do
end do
end do

          
	rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho_0
	T(1:m_local,1:n_local,1:o_local)=p(1:m_local,1:n_local,1:o_local)/(rho_0*287.0)
	!p= rho*R*T
	!p(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho_0*287.0*273
p(1:m_local,1:n_local,1:o_local)=rho(1:m_local,1:n_local,1:o_local)&
*101325.0/273/1.29*T(1:m_local,1:n_local,1:o_local)
	E(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)*(&
0.5*(u(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
**2+v(-1:m_local+2,-1:n_local+2,-1:o_local+2)**2&
+w(-1:m_local+2,-1:n_local+2,-1:o_local+2)**2)+T(-1:m_local+2,-1:n_local+2,-1:o_local+2)*717.5)

	rho_E(-1:m_local+2,-1:n_local+2,-1:o_local+2)=E(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*rho_0
	

	rho_u(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho_0*u(-1:m_local+2,-1:n_local+2,-1:o_local+2)
	rho_v(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho_0*v(-1:m_local+2,-1:n_local+2,-1:o_local+2)
	rho_w(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho_0*w(-1:m_local+2,-1:n_local+2,-1:o_local+2)
end if

if(mode==5) then

!Inputfile

 open(99,file='Ma5_chris.dat',form='formatted')
!jump though header

  read(99,'(A)') reader
  read(99,'(A)') reader
  read(99,'(A)') reader

dymax=241
print *, "test!"
!einlesen des Ähnlichkeitsprofiles

do i=1,dymax
       read(99,'(8(e12.5))') x,ybl(i),z,eta(i),rhobl(i),tbl(i),ubl(i)
enddo
! Berechnung von rho*u der Grenzschricht

!openfile
print *, "close file!"
        close(99)




do i=1,dymax
       rubl(i)=ubl(i)*rhobl(i)
enddo


v_ref=0.11495E+03
x0=300790.*1.70098338E-05/0.11495E+03/1.29
print *, "x0",x0
deta=0.20000e-02
mu0=1.725e-5


do j=1,M+10



        ybl1(1)=0.
        x_c(j)=x0+float(j-10)*dx
        sc=sqrt(2.*rho_0*mu0*x_c(j)/v_ref)


        do i=2,dymax
         ybl1(i)=ybl1(i-1)+sc*deta/rhobl(i)
        enddo


        call lagrange_intpol(ybl1,dy,dymax,ubl,ubl1,tbl,tbl1,rhobl,&
rhobl1,rubl,rubl1,N)

        do i=1,INT(N)
         ublrec(i,j)  =ubl1(i)
        tblrec(i,j)  =tbl1(i)
         rhoblrec(i,j)=rhobl1(i)
         rublrec(i,j) =rubl1(i)

        enddo

enddo






do i=1,n_local
do h=1,o_local
do j=1,m_local
	rho(j,i,h)=rhoblrec(i,j+10)
	T(j,i,h)=tblrec(i,j+10)
	u(j,i,h)= ublrec(i,j+10)
end do
end do
end do


!do i=n_local-int(n_local/2)+1,n_local
!do h=1,o_local
!do j=1,m_local
!	rho(j,i,h)=rhoblrec(n_local-i+1,j+10)
!	T(j,i,h)=tblrec(n_local-i+1,j+10)
!	u(j,i,h)= ublrec(n_local-i+1,j+10)

!	rho(j,i,h)=rhoblrec(n_local/2,j+10)
!	T(j,i,h)=tblrec(n_local/2,j+10)
!	u(j,i,h)= ublrec(n_local/2,j+10)

!end do
!end do
!end do


rho_u(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*u(-1:m_local+2,-1:n_local+2,-1:o_local+2)
	
!rho_u(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho_u(-1:m_local+2,-1:n_local+2,-1:o_local+2)/148.259689
! Konti 2D
!d rho_u /dx=- d rho_v / dy
CALL Dphi_Dx_p(rho_u,delta_u,qxi_i,qxi_nr,qxi_r)!Dphi_Dx_p(rho_u,solve,qxi_i,qxi_nr,qxi_r,u_0*rho_0)


!do i=1,n_local
! print *, i,rho_u(1,i,1),delta_u(1,i,1),rho_u(30,i,1),delta_u(30,i,1)
!end do



dummy=1
do dummy=1,o_local

do i=1,n_local
do j=1,m_local
 rublrec(i,j)=delta_u(j,i,dummy)
end do
end do

call nin(M,N,rublrec,rvblrec,0.,-dy,1)

!do i=1,n_local/4
do i=1,n_local
h=dummy
do j=1,m_local
rho_v(j,i,h)= rvblrec(i,j)
end do

end do

!do i=n_local/2-int(n_local/4)+1,n_local/2
!do i=n_local-int(n_local/2)+1,n_local
!h=dummy
!do j=1,m_local
! rho_v(j,i,h)=rvblrec(n_local/2-i+1,j)
! rho_v(j,i,h)= rvblrec(n_local/2,j)
!end do
!end do

!do j=1,m_local
!do i=n_local/2,n_local
!rho_v(j,i,h)= -rho_v(j,i-n_local/2+1,h)
!end do
!end do

end do
!rho_v(1:M,N,1:O)=0.0

!print *, "temp",rho_v(1,1:40,1),"temp2",rho_v(30,1:40,1)

w(-1:m_local+2,-1:n_local+2,-1:o_local+2)=0.0

v(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho_v(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
/rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)

u_bc(1,1:n_local,-1:o_local+2)=u(1,1:n_local,-1:o_local+2)
v_bc(1:m_local,1:n_local,-1:o_local+2)=v(1:m_local,1:n_local,-1:o_local+2)
T_bc(1,1:n_local,-1:o_local+2)=T(1,1:n_local,-1:o_local+2)
rho_bc(1,1:n_local,-1:o_local+2)=rho(1,1:n_local,-1:o_local+2)

p(1:m_local,1:n_local,1:o_local)=101325.0
!rho(1:m_local,1:n_local,1:o_local)&
!*101325.0/273./1.29*T(1:m_local,1:n_local,1:o_local)

E(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)*(&
0.5*(u(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
**2+v(-1:m_local+2,-1:n_local+2,-1:o_local+2)**2&
+w(-1:m_local+2,-1:n_local+2,-1:o_local+2)**2)+T(-1:m_local+2,-1:n_local+2,-1:o_local+2)*717.5)

rho_E(-1:m_local+2,-1:n_local+2,-1:o_local+2)=E(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)
	
rho_w(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*w(-1:m_local+2,-1:n_local+2,-1:o_local+2)


print *, "ende!"
!print *, "u", u(1,1:40,1)
end if


if(mode==6) then

!readoldsolution

 open(99,file='oldresult.dat',form='formatted')
!jump though header

print *, "open old file!"
  read(99,'(A)') reader
  read(99,'(A)') reader
  read(99,'(A)') reader

	
	Do d=o_local,1,-1				
			DO j=n_local,1,-1
	do i=m_local,1,-1			
!read(99,997) x,y,z,rho(i,j,d),&
!u(i,j,d),v(i,j,d),w(i,j,d),&
!E(i,j,d),T(i,j,d),p(i,j,d)

read(99,997) x,car,y,car,z,car,rho(i,j,d),car,&
u(i,j,d),car,v(i,j,d),car,w(i,j,d),car,&
E(i,j,d),car,T(i,j,d),car,p(i,j,d)


				


	

		ENDDO
	ENDDO
ENDDO

!do i=1,m_local
!do j=1,n_local
!do d=1,o_local

!rho(i,j,d)=rho(i,j,40)

!u(i,j,d)=u(i,j,40)
!v(i,j,d)=v(i,j,40)
!w(i,j,d)=w(i,j,40)
!E(i,j,d)=E(i,j,40)
!T(i,j,d)=T(i,j,40)
!p(i,j,d)=p(i,j,40)

!end do
!end do
!end do
u_bc(1,1:n_local,-1:o_local+2)=u(1,1:n_local,-1:o_local+2)
v_bc(1:m_local,1:n_local,-1:o_local+2)=v(1:m_local,1:n_local,-1:o_local+2)
T_bc(1,1:n_local,-1:o_local+2)=T(1,1:n_local,-1:o_local+2)
rho_bc(1,1:n_local,-1:o_local+2)=rho(1,1:n_local,-1:o_local+2)

rho_E(-1:m_local+2,-1:n_local+2,-1:o_local+2)=E(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)

rho_u(-1:m_local+2,-1:n_local+2,-1:o_local+2)=u(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)
	
rho_v(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*v(-1:m_local+2,-1:n_local+2,-1:o_local+2)
	
rho_w(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*w(-1:m_local+2,-1:n_local+2,-1:o_local+2)


!997	FORMAT(F22.10,F22.10,F22.10,F22.10,F22.10,F22.10,F22.10,F22.10,F22.10,F22.10)
997	FORMAT(F22.10,A1,F22.10,A1,F22.10,A1,F22.10,A1,F22.10,A1,F22.10,A1,F22.10,A1,F22.10,A1,F22.10,A1,F22.10)

print *, "close file!"
        close(99)
end if

if(mode==7) then
u(1:m_local,1:n_local,1:o_local)=1.0
v(1:m_local,1:n_local,1:o_local)=1.0
w(1:m_local,1:n_local,1:o_local)=0.0
T(1:m_local,1:n_local,1:o_local)=273.0
rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho_0

p(1:m_local,1:n_local,1:o_local)=rho_0*287.0*273.0

E(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)*(&
0.5*(u(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
**2+v(-1:m_local+2,-1:n_local+2,-1:o_local+2)**2&
+w(-1:m_local+2,-1:n_local+2,-1:o_local+2)**2)+T(-1:m_local+2,-1:n_local+2,-1:o_local+2)*717.5)

rho_E(-1:m_local+2,-1:n_local+2,-1:o_local+2)=E(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)

rho_u(-1:m_local+2,-1:n_local+2,-1:o_local+2)=u(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)
	
rho_v(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*v(-1:m_local+2,-1:n_local+2,-1:o_local+2)
	
rho_w(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*w(-1:m_local+2,-1:n_local+2,-1:o_local+2)
endif

if(mode==8) then
do i=1,n_local
!u(1:m_local,i,1:o_local)=1.0*i
v(i,1:n_local,1:o_local)=1.0
enddo
!v(1:m_local,1:n_local,1:o_local)=1.0
u(1:m_local,1:n_local,1:o_local)=3.0
w(1:m_local,1:n_local,1:o_local)=0.0
T(1:m_local,1:n_local,1:o_local)=273.0
rho(1:m_local,1:n_local,1:o_local)=rho_0

p(1:m_local,1:n_local,1:o_local)=rho_0*287.0*273.0


E(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)*(&
0.5*(u(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
**2+v(-1:m_local+2,-1:n_local+2,-1:o_local+2)**2&
+w(-1:m_local+2,-1:n_local+2,-1:o_local+2)**2)+T(-1:m_local+2,-1:n_local+2,-1:o_local+2)*717.5)

rho_E(-1:m_local+2,-1:n_local+2,-1:o_local+2)=E(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)

rho_u(-1:m_local+2,-1:n_local+2,-1:o_local+2)=u(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)
	
rho_v(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*v(-1:m_local+2,-1:n_local+2,-1:o_local+2)
	
rho_w(-1:m_local+2,-1:n_local+2,-1:o_local+2)=rho(-1:m_local+2,-1:n_local+2,-1:o_local+2)&
*w(-1:m_local+2,-1:n_local+2,-1:o_local+2)
endif

DEALLOCATE(ublrec, tblrec,rhoblrec,rublrec,rvblrec,rublrecsol)
RETURN
END SUBROUTINE

!Koeffizienten für kompakte differenzen belegen
SUBROUTINE INITIAL_COEFF

	IMPLICIT NONE
	!1.Ableitung
	!innere Terme O6
	alpha_i=1.
	beta_i=3.
	gamma_i=1.
	q_i=12.
	a_i=-1.
	b_i=-28.
	c_i=0.
	d_i=28.
	e_i=1.
	!Terme neben Rand O4
	alpha_nr=1.
	beta_nr=4.
	gamma_nr=1.
	q_nr=1.
	b_nr=-3.
	c_nr=0.
	d_nr=3.
	!Terme am Rand O3
	beta_r=1.
	gamma_r=2.
	q_r=2.
	c_r=-5.
	d_r=4.
	e_r=1.
	
	!1/(dx*q)
	qxi_i=1./(dx*q_i)
	qxi_nr=1./(dx*q_nr)
	qxi_r=1./(dx*q_r)
	qyi_i=1./(dy*q_i)
	qyi_nr=1./(dy*q_nr)
	qyi_r=1./(dy*q_r)
	qzi_i=1./(dz*q_i)
	qzi_nr=1./(dz*q_nr)
	qzi_r=1./(dz*q_r)

	!2.Ableitung
	!innere Terme O6
	alpha2_i=2.
	beta2_i=11.
	gamma2_i=2.
	q2_i=4.
	a2_i=3.
	b2_i=48.
	c2_i=-102.
	d2_i=48.
	e2_i=3.
	!Terme neben Rand O4
	alpha2_nr=1.
	beta2_nr=10.
	gamma2_nr=1.
	q2_nr=1.
	b2_nr=12.
	c2_nr=-24.
	d2_nr=12.
	!Terme am Rand O3
	beta2_r=1.
	gamma2_r=11.
	q2_r=13.
	c2_r=-27.
	d2_r=15.
	e2_r=-1.
	
	!1/(dx^2*q)
	qxi2_i=1./(dx**2*q_i)
	qxi2_nr=1./(dx**2*q_nr)
	qxi2_r=1./(dx**2*q_r)
	qyi2_i=1./(dy**2*q_i)
	qyi2_nr=1./(dy**2*q_nr)
	qyi2_r=1./(dy**2*q_r)
	qzi2_i=1./(dz**2*q_i)
	qzi2_nr=1./(dz**2*q_nr)
	qzi2_r=1./(dz**2*q_r)	

	RETURN

END SUBROUTINE
 
end module
