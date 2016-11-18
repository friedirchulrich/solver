module sync_mod
USE Variablen


contains

!Ghost-cells aktualisieren, periodisch
SUBROUTINE GC_EXCHANGE_P

	IMPLICIT NONE

	!!synchronisieren benachbarter Images
	CALL SYNC_NEIGHBOURS_P

	!!Datentausch
	!! m_pos > 1 bedeutet, dass das image einen linkes Nachbarimage besitzt
IF(m_pos>1)THEN
	rho(-1:0,1:n_local,1:o_local)=rho(m_local-1:m_local,1:n_local,1:o_local)[m_pos-1,n_pos,o_pos]

	rho_u(-1:0,1:n_local,1:o_local)=&
	rho_u(m_local-1:m_local,1:n_local,1:o_local)[m_pos-1,n_pos,o_pos]

	rho_v(-1:0,1:n_local,1:o_local)=&
	rho_v(m_local-1:m_local,1:n_local,1:o_local)[m_pos-1,n_pos,o_pos]

	rho_w(-1:0,1:n_local,1:o_local)=&
	rho_w(m_local-1:m_local,1:n_local,1:o_local)[m_pos-1,n_pos,o_pos]

	E(-1:0,1:n_local,1:o_local)=E(m_local-1:m_local,1:n_local,1:o_local)[m_pos-1,n_pos,o_pos]

	p(-1:0,1:n_local,1:o_local)=p(m_local-1:m_local,1:n_local,1:o_local)[m_pos-1,n_pos,o_pos]
		
	u(-1:0,1:n_local,1:o_local)=u(m_local-1:m_local,1:n_local,1:o_local)[m_pos-1,n_pos,o_pos]

	v(-1:0,1:n_local,1:o_local)=v(m_local-1:m_local,1:n_local,1:o_local)[m_pos-1,n_pos,o_pos]

	w(-1:0,1:n_local,1:o_local)=w(m_local-1:m_local,1:n_local,1:o_local)[m_pos-1,n_pos,o_pos]
	
ELSE  
	! m_pos==1 Das Linke Ende entspricht dem Rand
	! Bei einer Periodischen Randbedingung wird es dem rechten Rand gleichgesetzt
	rho(-1:0,1:n_local,1:o_local)=rho(m_local-1:m_local,1:n_local,1:o_local)[R,n_pos,o_pos]

	rho_u(-1:0,1:n_local,1:o_local)=&
	rho_u(m_local-1:m_local,1:n_local,1:o_local)[R,n_pos,o_pos]

	rho_v(-1:0,1:n_local,1:o_local)=&
	rho_v(m_local-1:m_local,1:n_local,1:o_local)[R,n_pos,o_pos]

	rho_w(-1:0,1:n_local,1:o_local)=&
	rho_w(m_local-1:m_local,1:n_local,1:o_local)[R,n_pos,o_pos]

	E(-1:0,1:n_local,1:o_local)=E(m_local-1:m_local,1:n_local,1:o_local)[R,n_pos,o_pos]

	p(-1:0,1:n_local,1:o_local)=p(m_local-1:m_local,1:n_local,1:o_local)[R,n_pos,o_pos]
		
	u(-1:0,1:n_local,1:o_local)=u(m_local-1:m_local,1:n_local,1:o_local)[R,n_pos,o_pos]

	v(-1:0,1:n_local,1:o_local)=v(m_local-1:m_local,1:n_local,1:o_local)[R,n_pos,o_pos]

	w(-1:0,1:n_local,1:o_local)=w(m_local-1:m_local,1:n_local,1:o_local)[R,n_pos,o_pos]

	ENDIF

	!! Alle Zellen mit m_pos kleiner als R besitzen eine rechte Nachbarzelle
IF(m_pos<R)THEN

	rho(m_local+1:m_local+2,1:n_local,1:o_local)=rho(1:2,1:n_local,1:o_local)[m_pos+1,n_pos,o_pos]

	rho_u(m_local+1:m_local+2,1:n_local,1:o_local)=&
	rho_u(1:2,1:n_local,1:o_local)[m_pos+1,n_pos,o_pos]

	rho_v(m_local+1:m_local+2,1:n_local,1:o_local)=&
	rho_v(1:2,1:n_local,1:o_local)[m_pos+1,n_pos,o_pos]

	rho_w(m_local+1:m_local+2,1:n_local,1:o_local)=&
	rho_w(1:2,1:n_local,1:o_local)[m_pos+1,n_pos,o_pos]

	E(m_local+1:m_local+2,1:n_local,1:o_local)=E(1:2,1:n_local,1:o_local)[m_pos+1,n_pos,o_pos]

	p(m_local+1:m_local+2,1:n_local,1:o_local)=p(1:2,1:n_local,1:o_local)[m_pos+1,n_pos,o_pos]
		
	u(m_local+1:m_local+2,1:n_local,1:o_local)=u(1:2,1:n_local,1:o_local)[m_pos+1,n_pos,o_pos]

	v(m_local+1:m_local+2,1:n_local,1:o_local)=v(1:2,1:n_local,1:o_local)[m_pos+1,n_pos,o_pos]

	w(m_local+1:m_local+2,1:n_local,1:o_local)=w(1:2,1:n_local,1:o_local)[m_pos+1,n_pos,o_pos]


	ELSE !!m_pos=R Rechter Rand, dieser wird mit dem linken Rand verbunden

	rho(m_local+1:m_local+2,1:n_local,1:o_local)=rho(1:2,1:n_local,1:o_local)[1,n_pos,o_pos]

	rho_u(m_local+1:m_local+2,1:n_local,1:o_local)=&
	rho_u(1:2,1:n_local,1:o_local)[1,n_pos,o_pos]

	rho_v(m_local+1:m_local+2,1:n_local,1:o_local)=&
	rho_v(1:2,1:n_local,1:o_local)[1,n_pos,o_pos]

	rho_w(m_local+1:m_local+2,1:n_local,1:o_local)=&
	rho_w(1:2,1:n_local,1:o_local)[1,n_pos,o_pos]

	E(m_local+1:m_local+2,1:n_local,1:o_local)=E(1:2,1:n_local,1:o_local)[1,n_pos,o_pos]

	p(m_local+1:m_local+2,1:n_local,1:o_local)=p(1:2,1:n_local,1:o_local)[1,n_pos,o_pos]
		
	u(m_local+1:m_local+2,1:n_local,1:o_local)=u(1:2,1:n_local,1:o_local)[1,n_pos,o_pos]

	v(m_local+1:m_local+2,1:n_local,1:o_local)=v(1:2,1:n_local,1:o_local)[1,n_pos,o_pos]

	w(m_local+1:m_local+2,1:n_local,1:o_local)=w(1:2,1:n_local,1:o_local)[1,n_pos,o_pos]	

	
	ENDIF


	!! Ein Image besitzt ein Image unter sich, Unterer Rand wird 		syncronisiert
IF(n_pos>1)THEN

	rho(1:m_local,-1:0,1:o_local)=rho(1:m_local,n_local-1:n_local,1:o_local)[m_pos,n_pos-1,o_pos]

	rho_u(1:m_local,-1:0,1:o_local)=&
	rho_u(1:m_local,n_local-1:n_local,1:o_local)[m_pos,n_pos-1,o_pos]

	rho_v(1:m_local,-1:0,1:o_local)=&
	rho_v(1:m_local,n_local-1:n_local,1:o_local)[m_pos,n_pos-1,o_pos]

	rho_w(1:m_local,-1:0,1:o_local)=&
	rho_w(1:m_local,n_local-1:n_local,1:o_local)[m_pos,n_pos-1,o_pos]

	E(1:m_local,-1:0,1:o_local)=E(1:m_local,n_local-1:n_local,1:o_local)[m_pos,n_pos-1,o_pos]

	p(1:m_local,-1:0,1:o_local)=p(1:m_local,n_local-1:n_local,1:o_local)[m_pos,n_pos-1,o_pos]
		
	u(1:m_local,-1:0,1:o_local)=u(1:m_local,n_local-1:n_local,1:o_local)[m_pos,n_pos-1,o_pos]

	v(1:m_local,-1:0,1:o_local)=v(1:m_local,n_local-1:n_local,1:o_local)[m_pos,n_pos-1,o_pos]

	w(1:m_local,-1:0,1:o_local)=w(1:m_local,n_local-1:n_local,1:o_local)[m_pos,n_pos-1,o_pos]


	!! n_pos==1 Wenn es kein Image unter dem Image gibt, gelten die 
	!  Randbedingungen
	ELSE

	rho(1:m_local,-1:0,1:o_local)=rho(1:m_local,n_local-1:n_local,1:o_local)[m_pos,S,o_pos]

	rho_u(1:m_local,-1:0,1:o_local)=&
	rho_u(1:m_local,n_local-1:n_local,1:o_local)[m_pos,S,o_pos]

	rho_v(1:m_local,-1:0,1:o_local)=&
	rho_v(1:m_local,n_local-1:n_local,1:o_local)[m_pos,S,o_pos]

	rho_w(1:m_local,-1:0,1:o_local)=&
	rho_w(1:m_local,n_local-1:n_local,1:o_local)[m_pos,S,o_pos]

	E(1:m_local,-1:0,1:o_local)=E(1:m_local,n_local-1:n_local,1:o_local)[m_pos,S,o_pos]

	p(1:m_local,-1:0,1:o_local)=p(1:m_local,n_local-1:n_local,1:o_local)[m_pos,S,o_pos]
		
	u(1:m_local,-1:0,1:o_local)=u(1:m_local,n_local-1:n_local,1:o_local)[m_pos,S,o_pos]

	v(1:m_local,-1:0,1:o_local)=v(1:m_local,n_local-1:n_local,1:o_local)[m_pos,S,o_pos]

	w(1:m_local,-1:0,1:o_local)=w(1:m_local,n_local-1:n_local,1:o_local)[m_pos,S,o_pos]

	ENDIF
	
	!! Jedes Image besitzt ein Image über sich
IF(n_pos<S)THEN

	rho(1:m_local,n_local+1:n_local+2,1:o_local)=&
rho(1:m_local,1:2,1:o_local)[m_pos,n_pos+1,o_pos]

	rho_u(1:m_local,n_local+1:n_local+2,1:o_local)=&
	rho_u(1:m_local,1:2,1:o_local)[m_pos,n_pos+1,o_pos]

	rho_v(1:m_local,n_local+1:n_local+2,1:o_local)=&
	rho_v(1:m_local,1:2,1:o_local)[m_pos,n_pos+1,o_pos]

	rho_w(1:m_local,n_local+1:n_local+2,1:o_local)=&
	rho_w(1:m_local,1:2,1:o_local)[m_pos,n_pos+1,o_pos]

	E(1:m_local,n_local+1:n_local+2,1:o_local)=&
E(1:m_local,1:2,1:o_local)[m_pos,n_pos+1,o_pos]

	p(1:m_local,n_local+1:n_local+2,1:o_local)=&
p(1:m_local,1:2,1:o_local)[m_pos,n_pos+1,o_pos]
		
	u(1:m_local,n_local+1:n_local+2,1:o_local)=&
u(1:m_local,1:2,1:o_local)[m_pos,n_pos+1,o_pos]

	v(1:m_local,n_local+1:n_local+2,1:o_local)=&
v(1:m_local,1:2,1:o_local)[m_pos,n_pos+1,o_pos]

	w(1:m_local,n_local+1:n_local+2,1:o_local)=&
w(1:m_local,1:2,1:o_local)[m_pos,n_pos+1,o_pos]


!! n_pos==S, d.h. es handelt sich um den unteren Rand	
ELSE

	rho(1:m_local,n_local+1:n_local+2,1:o_local)=&
rho(1:m_local,1:2,1:o_local)[m_pos,1,o_pos]

	rho_u(1:m_local,n_local+1:n_local+2,1:o_local)=&
	rho_u(1:m_local,1:2,1:o_local)[m_pos,1,o_pos]

	rho_v(1:m_local,n_local+1:n_local+2,1:o_local)=&
	rho_v(1:m_local,1:2,1:o_local)[m_pos,1,o_pos]

	rho_w(1:m_local,n_local+1:n_local+2,1:o_local)=&
	rho_w(1:m_local,1:2,1:o_local)[m_pos,1,o_pos]

	E(1:m_local,n_local+1:n_local+2,1:o_local)=&
E(1:m_local,1:2,1:o_local)[m_pos,1,o_pos]

	p(1:m_local,n_local+1:n_local+2,1:o_local)=&
p(1:m_local,1:2,1:o_local)[m_pos,1,o_pos]
		
	u(1:m_local,n_local+1:n_local+2,1:o_local)=&
u(1:m_local,1:2,1:o_local)[m_pos,1,o_pos]

	v(1:m_local,n_local+1:n_local+2,1:o_local)=&
v(1:m_local,1:2,1:o_local)[m_pos,1,o_pos]

	w(1:m_local,n_local+1:n_local+2,1:o_local)=&
w(1:m_local,1:2,1:o_local)[m_pos,1,o_pos]

	ENDIF

!! vorder und Hinterer Rand

!! vorderer Rand
!! Ein Image besitzt ein Image VOR sich
IF(o_pos>1)THEN

	rho(1:m_local,1:n_local,-1:0)=rho(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,o_pos-1]

	rho_u(1:m_local,1:n_local,-1:0)=&
	rho_u(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,o_pos-1]

	rho_v(1:m_local,1:n_local,-1:0)=&
	rho_v(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,o_pos-1]

	rho_w(1:m_local,1:n_local,-1:0)=&
	rho_w(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,o_pos-1]

	E(1:m_local,1:n_local,-1:0)=E(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,o_pos-1]

	p(1:m_local,1:n_local,-1:0)=p(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,o_pos-1]
		
	u(1:m_local,1:n_local,-1:0)=u(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,o_pos-1]

	v(1:m_local,1:n_local,-1:0)=v(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,o_pos-1]

	w(1:m_local,1:n_local,-1:0)=w(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,o_pos-1]


	!! n_pos==1 Wenn es kein Image unter dem Image gibt, gelten die 
	!  Randbedingungen
	ELSE

	rho(1:m_local,1:n_local,-1:0)=rho(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,T_p]

	rho_u(1:m_local,1:n_local,-1:0)=&
	rho_u(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,T_p]

	rho_v(1:m_local,1:n_local,-1:0)=&
	rho_v(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,T_p]

	rho_w(1:m_local,1:n_local,-1:0)=&
	rho_w(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,T_p]

	E(1:m_local,1:n_local,-1:0)=E(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,T_p]

	p(1:m_local,1:n_local,-1:0)=p(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,T_p]
		
	u(1:m_local,1:n_local,-1:0)=u(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,T_p]

	v(1:m_local,1:n_local,-1:0)=v(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,T_p]

	w(1:m_local,1:n_local,-1:0)=w(1:m_local,1:n_local,o_local-1:o_local)[m_pos,n_pos,T_p]

	ENDIF
	
	!! Jedes Image besitzt ein Image HINTER sich
IF(o_pos<T_p)THEN

	rho(1:m_local,1:n_local,o_local+1:o_local+2)=&
rho(1:m_local,1:n_local,1:2)[m_pos,n_pos,o_pos+1]

	rho_u(1:m_local,1:n_local,o_local+1:o_local+2)=&
	rho_u(1:m_local,1:n_local,1:2)[m_pos,n_pos,o_pos+1]

	rho_v(1:m_local,1:n_local,o_local+1:o_local+2)=&
	rho_v(1:m_local,1:n_local,1:2)[m_pos,n_pos,o_pos+1]

	rho_w(1:m_local,1:n_local,o_local+1:o_local+2)=&
	rho_w(1:m_local,1:n_local,1:2)[m_pos,n_pos,o_pos+1]

	E(1:m_local,1:n_local,o_local+1:o_local+2)=&
E(1:m_local,1:n_local,1:2)[m_pos,n_pos,o_pos+1]

	p(1:m_local,1:n_local,o_local+1:o_local+2)=&
p(1:m_local,1:n_local,1:2)[m_pos,n_pos,o_pos+1]
		
	u(1:m_local,1:n_local,o_local+1:o_local+2)=&
u(1:m_local,1:n_local,1:2)[m_pos,n_pos,o_pos+1]

	v(1:m_local,1:n_local,o_local+1:o_local+2)=&
v(1:m_local,1:n_local,1:2)[m_pos,n_pos,o_pos+1]

	w(1:m_local,n_local+1:n_local+2,1:o_local)=&
w(1:m_local,1:n_local,1:2)[m_pos,n_pos,o_pos+1]


!! o_pos==T_p, d.h. es handelt sich um den hinteren Rand	
ELSE

	rho(1:m_local,1:n_local,o_local+1:o_local+2)=&
rho(1:m_local,1:n_local,1:2)[m_pos,n_pos,1]

	rho_u(1:m_local,1:n_local,o_local+1:o_local+2)=&
	rho_u(1:m_local,1:n_local,1:2)[m_pos,n_pos,1]

	rho_v(1:m_local,1:n_local,o_local+1:o_local+2)=&
	rho_v(1:m_local,1:n_local,1:2)[m_pos,n_pos,1]

	rho_w(1:m_local,1:n_local,o_local+1:o_local+2)=&
	rho_w(1:m_local,1:n_local,1:2)[m_pos,n_pos,1]

	E(1:m_local,1:n_local,o_local+1:o_local+2)=&
E(1:m_local,1:n_local,1:2)[m_pos,n_pos,1]

	p(1:m_local,1:n_local,o_local+1:o_local+2)=&
p(1:m_local,1:n_local,1:2)[m_pos,n_pos,1]
		
	u(1:m_local,1:n_local,o_local+1:o_local+2)=&
u(1:m_local,1:n_local,1:2)[m_pos,n_pos,1]

	v(1:m_local,1:n_local,o_local+1:o_local+2)=&
v(1:m_local,1:n_local,1:2)[m_pos,n_pos,1]

	w(1:m_local,n_local+1:n_local+2,1:o_local)=&
w(1:m_local,1:n_local,1:2)[m_pos,n_pos,1]


	ENDIF


	!!synchronisieren benachbarter Images
	CALL SYNC_NEIGHBOURS_P

	RETURN

END SUBROUTINE

!Synchronisieren benachbarter Images, periodisch, mit workaround für Bug bei g95

SUBROUTINE SYNC_NEIGHBOURS_P

	IMPLICIT NONE

	SYNC ALL

! sollte obsolute sein, da open array fortran verwendet wird... mal sehen


!	IF(m_pos==1)THEN
!		IF(n_pos==1)THEN
!			SYNC IMAGES([sync_m_up,sync_n_up])
!		ELSEIF(n_pos==S)THEN
!			SYNC IMAGES([sync_m_up,sync_n_down])
!		ELSE
!			SYNC IMAGES([sync_m_up,sync_n_up,sync_n_down])
!		ENDIF
!	ELSEIF(m_pos==R)THEN
!		IF(n_pos==1)THEN
!			SYNC IMAGES([sync_m_down,sync_n_up])
!		ELSEIF(n_pos==S)THEN
!			SYNC IMAGES([sync_m_down,sync_n_down])
!		ELSE
!			SYNC IMAGES([sync_m_down,sync_n_up,sync_n_down])
!		ENDIF		
!	ELSE
!		IF(n_pos==1)THEN
!			SYNC IMAGES([sync_m_up,sync_m_down,sync_n_up])
!		ELSEIF(n_pos==S)THEN
!			SYNC IMAGES([sync_m_up,sync_m_down,sync_n_down])
!		ELSE
!			SYNC IMAGES([sync_m_up,sync_m_down,sync_n_up,sync_n_down])
!		ENDIF
!	ENDIF	
!	
! Workaround um Ringsynchronisation zu vermeiden, die in g95 einen Dead-Lock verursacht. 
!	IF(m_pos==1)SYNC IMAGES(image_m_max)
!	IF(m_pos==R)SYNC IMAGES(image_m_1)
!	IF(n_pos==1)SYNC IMAGES(image_n_max)
!	IF(n_pos==S)SYNC IMAGES(image_n_1)

	RETURN

END SUBROUTINE
end module
