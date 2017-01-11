module INPUTFILE


contains

subroutine inputfileread(M,N,O,R,S,T_p,i,itmax,cfl,Re,Ma,mode)

!variables
CHARACTER(5)::inputfilename
CHARACTER(10)::inputcommand
real*8, intent(out):: Re,Ma,cfl
integer, intent(out):: M,N,O,R,S,T_p,i,itmax,Mode
logical:: Mset,Nset,Oset,Rset,Sset,Tset,itmaxset,cflset,modeset,Reset,MAset
!subroutine program

!define file name
inputfilename='input'

! set all set variables to false
! if a variable is defined these variables
! will be set to true, if thez remain
! false, the program will sent a warning to the
! user and choose a default value
Mset=.FALSE.
Nset=.FALSE.
Rset=.FALSE.
Sset=.FALSE.
itmaxset=.FALSE.
cflset=.FALSE.
modeset=.FALSE.
!open input file
OPEN(unit=21,file=inputfilename)
Print *, "Input File was found and opened correctly!"


!read file
i=0
do while( i<1000) 
read(21,*) inputcommand

! define input commands
! with the following syntax:
! IF(inputcommand.EQ.'commandname' then
! ...read special parameter
! print value to screen for checking
! end if

! Number of Points in x-Direction
	IF(inputcommand.EQ.'M         ') then
	Print *, " M ="
	read(21,*) M
	print *, M
	Mset=.TRUE.
	END IF

! Number of Points in y-Direction
	IF(inputcommand.EQ.'N         ') then
	Print *, " N ="
	read(21,*) N
	print *, N
	Nset=.TRUE.
	END IF

! Number of Points in z-Direction
	IF(inputcommand.EQ.'O         ') then
	Print *, " O ="
	read(21,*) O
	print *, O
	Oset=.TRUE.
	END IF


! Number of Processes in x-Direction
	IF(inputcommand.EQ.'R         ') then
	Print *, " R ="
	read(21,*) R
	print *, R
	Rset=.TRUE.
	END IF


! Number of Processes in y-Direction
	IF(inputcommand.EQ.'S         ') then
	Print *, " S ="
	read(21,*) S
	print *, S
	Sset=.TRUE.
	END IF

! Number of Processes in z-Direction
	IF(inputcommand.EQ.'T         ') then
	Print *, " T ="
	read(21,*) T_p
	print *, T_p
	Tset=.TRUE.
	END IF

! Number of Iterations
	IF(inputcommand.EQ.'ITMAX         ') then
	Print *, " ITMAX ="
	read(21,*) itmax
	print *, itmax
	itmaxset=.TRUE.
	END IF

! CFL-Number
	IF(inputcommand.EQ.'CFL         ') then
	Print *, " CFL ="
	read(21,*) cfl
	print *, cfl
	cflset=.TRUE.
	END IF

! CFL-Number
	IF(inputcommand.EQ.'RE         ') then
	Print *, " RE ="
	read(21,*) Re
	print *, Re
	Reset=.TRUE.
	END IF

! Mach-Number
	IF(inputcommand.EQ.'MA         ') then
	Print *, " MA ="
	read(21,*) Ma
	print *, Ma
	MAset=.TRUE.
	END IF

! Select Initial Profile
	IF(inputcommand.EQ.'MODE         ') then
	Print *, " MODE ="
	read(21,*) mode
	print *, mode
		IF(mode==1)THEN
		print *, "Sinus-Profil"
		END IF

		IF(mode==2)THEN
		print *, " Taylor-Green"
		END IF

		IF(mode==3)THEN
		print *, "Kanal"
		END IF

		IF(mode==4)THEN
		print *, "Navier-Stokes"
		END IF

		IF(mode==5)THEN
		print *, "Kanalströmung"
		END IF

	modeset=.TRUE.
	END IF


	IF(inputcommand.EQ.'END        ') then
	Print *, " End of file"
	i=1000 
	END IF
i=i+1
end do
CLOSE(21)

! sent out warning and choose default backup value
! if certain values are not defined by the user in the
! input file

IF(Mset.eqv..FALSE.) THEN
print *, "*******WARNING**********"
print *, " Variable M is not defined in the input file!"
print *, " We therefore chose:"
print *, " M = 10"
M=10
END IF

IF(Nset.eqv..FALSE.) THEN
print *, "*******WARNING**********"
print *, " Variable N is not defined in the input file!"
print *, " We therefore chose:"
print *, " N = 10"
N=10
END IF

IF(Oset.eqv..FALSE.) THEN
print *, "*******WARNING**********"
print *, " Variable O is not defined in the input file!"
print *, " We therefore chose:"
print *, " O = 10"
O=10
END IF

IF(Rset.eqv..FALSE.) THEN
print *, "*******WARNING**********"
print *, " Variable R is not defined in the input file!"
print *, " We therefore chose:"
print *, " R= 10"
R=10
END IF

IF(Sset.eqv..FALSE.) THEN
print *, "*******WARNING**********"
print *, " Variable S is not defined in the input file!"
print *, " We therefore chose:"
print *, " S = 10"
S=10
END IF

IF(Tset.eqv..FALSE.) THEN
print *, "*******WARNING**********"
print *, " Variable T is not defined in the input file!"
print *, " We therefore chose:"
print *, " T = 10"
T_p=10
END IF

IF(itmaxset.eqv..FALSE.) THEN
print *, "*******WARNING**********"
print *, " Variable ITMAX is not defined in the input file!"
print *, " We therefore chose:"
print *, " ITMAX = 1000"
itmax=1000
END IF

IF(cflset.eqv..FALSE.) THEN
print *, "*******WARNING**********"
print *, " Variable CFL is not defined in the input file!"
print *, " We therefore chose:"
print *, " CFL = 0.001 "
cfl=0.001 
END IF

IF(reset.eqv..FALSE.) THEN
print *, "*******WARNING**********"
print *, " Variable RE is not defined in the input file!"
print *, " We therefore chose:"
print *, " RE = 0.001 "
Re=0.001 
END IF

IF(modeset.eqv..FALSE.) THEN
print *, "*******WARNING**********"
print *, " Variable MODE is not defined in the input file!"
print *, " We therefore chose:"
print *, " MODE = 2 (Taylor Green) "
mode=2 
END IF

IF(MAset.eqv..FALSE.) THEN
print *, "*******WARNING**********"
print *, " Variable M is not defined in the input file!"
print *, " We therefore chose:"
print *, " MA = 0.1"
MA=0.1
END IF

end subroutine

end module
