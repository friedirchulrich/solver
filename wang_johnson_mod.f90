Module wang_johnson_mod
use variablen

contains

SUBROUTINE WANG_JOHNSSON_P_x(o_rank)

	IMPLICIT NONE

	INTEGER::i, o_rank,h,j
	REAL*8::hi,zi
	REAL*8::t_n,s_n
	INTEGER::n_anfang,n_ende

	!!SYNC
	SYNC IMAGES(sync_m_array)
	!!	
    !alpha_i=1.
    !beta_i=3.
    !gamma_i=1.
    !q_i=12.
    !a_i=-1.
    !b_i=-28.
    !c_i=0.
    !d_i=28.
    !e_i=1.
	!Belegung Matrix-Elemente
	IF(R==1.AND.S==1)THEN
	a_x(3:m_local-2)=beta_i
	a_x(1)=1.
	a_x(2)=4.
	a_x(m_local-1)=4.
	a_x(m_local)=1.
	b_x(1)=0.
	b_x(2:m_local-1)=alpha_i
	b_x(m_local)=2.
	c_x(1)=2.
	c_x(2:m_local-1)=gamma_i
	c_x(m_local)=0.


!	!Vektoren für Sherman-Morrison-Formel
!	a_x(1)=a_x(1)+s_n
!	a_x(m_local)=a_x(m_local)+t_n
!	y_x(1,1:n_local,2)=1.
!	y_x(2:m_local-1,1:n_local,2)=0.
!	y_x(m_local,1:n_local,2)=-1.
!	z_x(1)=s_n
!	z_x(2:m_local-1)=0.
!	z_x(m_local)=-t_n

!	print *,'b_x(1:m_local): ',b_x(1:m_local),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	print *,'a_x(1:m_local): ',a_x(1:m_local),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	print *,'c_x(1:m_local): ',c_x(1:m_local),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	print *,'z_x(1:m_local): ',z_x(1:m_local),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	print *,'y_x(1:m_local,i:n_local,2): ',y_x(1:m_local,1:n_local,2)
	DO i=2,m_local
		hi=b_x(i)/a_x(i-1)
		a_x(i)=a_x(i)-hi*c_x(i-1)
		y_x(i,1:n_local,1)=y_x(i,1:n_local,1)-hi*y_x(i-1,1:n_local,1)
	ENDDO



	x_x(m_local,1:n_local,1)=y_x(m_local,1:n_local,1)/a_x(m_local)

	Do i=m_local-1,1,-1
		x_x(i,1:n_local,1:2)=(y_x(i,1:n_local,1:2)-c_x(i)*x_x(i+1,1:n_local,1:2))/a_x(i)
	ENDDO

	ELSE
	IF(m_pos==1)THEN
		a_x(1)=1.
		a_x(2)=4.
		a_x(3:m_local)=beta_i
		b_x(1)=0.
		b_x(2:m_local)=alpha_i
		c_x(1)=2.
		c_x(2:m_local)=gamma_i
	ELSEIF(m_pos==R)THEN
		a_x(1:m_local-2)=beta_i
		a_x(m_local-1)=4.
		a_x(m_local)=1.
		b_x(1:m_local-1)=alpha_i
		b_x(m_local)=2.
		c_x(1:m_local-1)=gamma_i
		c_x(m_local)=0.
	ELSE
		a_x(1:m_local)=beta_i
		b_x(1:m_local)=alpha_i
		c_x(1:m_local)=gamma_i
	ENDIF	
!Belegung Matrix-Elemente
!	a_x(1:m_local)=beta_i
!	IF(m_pos==1)THEN
!		b_x(1)=0.
!		b_x(2:m_local)=alpha_i
!		c_x(1:m_local)=gamma_i
!		s_n=gamma_i
!	ELSEIF(m_pos==R)THEN
!		b_x(1:m_local)=alpha_i
!		c_x(1:m_local-1)=gamma_i
!		c_x(m_local)=0.
!		t_n=alpha_i
!	ELSE
!		b_x(1:m_local)=alpha_i
!		c_x(1:m_local)=gamma_i
!	ENDIF	
	
!	!Vektoren für Sherman-Morrison-Formel
!	IF(m_pos==1)THEN
!		a_x(1)=a_x(1)+s_n
!		y_x(1,1:n_local,2)=1.
!		y_x(2:m_local,1:n_local,2)=0.
!		z_x(1)=s_n
!		z_x(2:m_local)=0.
!	ELSEIF(m_pos==R)THEN
!		a_x(m_local)=a_x(m_local)+t_n
!		y_x(1:m_local-1,1:n_local,2)=0.
!		y_x(m_local,1:n_local,2)=-1.
!		z_x(1:m_local-1)=0.
!		z_x(m_local)=-t_n
!	ELSE
!		y_x(1:m_local,1:n_local,2)=0.
!		z_x(1:m_local)=0.
!	ENDIF	

	!!Schritt1
	IF(m_pos>1)f_x(1)=b_x(1)
	DO i=2,m_local
		hi=b_x(i)/a_x(i-1)
		a_x(i)=a_x(i)-hi*c_x(i-1)
		y_x(i,1:n_local,1:2)=y_x(i,1:n_local,1:2)-hi*y_x(i-1,1:n_local,1:2)
		IF(m_pos>1)f_x(i)=-hi*f_x(i-1)
	ENDDO
	
	!SYNC ALL
	!IF(m_pos==2.AND.n_pos==1)THEN
	!	PRINT *,' '
	!	PRINT *,'m_pos2'
	!	DO i=1,m_local
	!		PRINT *,y_x(i,1,1)
	!	ENDDO
	!ENDIF

	!!Schritt2
	g_x(m_local-1)=c_x(m_local-1)
	DO i=m_local-1,2,-1
		hi=c_x(i-1)/a_x(i)
		g_x(i-1)=-hi*g_x(i)
		y_x(i-1,1:n_local,1:2)=y_x(i-1,1:n_local,1:2)-hi*y_x(i,1:n_local,1:2)
		IF(m_pos>1)f_x(i-1)=f_x(i-1)-hi*f_x(i)
	ENDDO

	!!SYNC
	IF(m_pos==R)THEN
		SYNC IMAGES(sync_m_down)
	ELSEIF(m_pos==1)THEN
		SYNC IMAGES(sync_m_up)
	ELSE
		SYNC IMAGES([sync_m_up,sync_m_down])
	ENDIF
	!! Transfer
	IF(m_pos<R)THEN
		a_x(m_local+1)=a_x(1)[m_pos+1,n_pos,o_rank]
		f_x(m_local+1)=f_x(1)[m_pos+1,n_pos,o_rank]
		g_x(m_local+1)=g_x(1)[m_pos+1,n_pos,o_rank]
		y_x(m_local+1,1:n_local,1:2)=y_x(1,1:n_local,1:2)[m_pos+1,n_pos,o_rank]
	ENDIF
	!!SYNC
	IF(m_pos==R)THEN
		SYNC IMAGES(sync_m_down)
	ELSEIF(m_pos==1)THEN
		SYNC IMAGES(sync_m_up)
	ELSE
		SYNC IMAGES([sync_m_up,sync_m_down])
	ENDIF
	!!
	IF(m_pos<R)THEN
		hi=c_x(m_local)/a_x(m_local+1)
		a_x(m_local)=a_x(m_local)-hi*f_x(m_local+1)
		g_x(m_local)=-hi*g_x(m_local+1)
		y_x(m_local,1:n_local,1:2)=y_x(m_local,1:n_local,1:2)-hi*y_x(m_local+1,1:n_local,1:2)
	ENDIF

!	SYNC ALL
!	IF(m_pos==1.AND.n_pos==1)THEN
!		DO h=1,n_local
!			!PRINT *,'m_pos1'
!			!DO i=1,m_local			
!				PRINT *,y_x(m_local,h,1)
!			!ENDDO
!		ENDDO
!	ENDIF
	

	!!Schritt3: Kopplungssystem
	!!SYNC
	SYNC IMAGES(sync_m_array)
	!!
	!!Transfer
	n_anfang=(m_pos-1)*n_kop+1
	n_ende=m_pos*n_kop
	DO i=1,R
		IF(i/=m_pos)THEN
			a_x_s(i)=a_x(m_local)[i,n_pos,o_rank]
			f_x_s(i)=f_x(m_local)[i,n_pos,o_rank]
			g_x_s(i)=g_x(m_local)[i,n_pos,o_rank]
			y_x_s(i,1:n_kop,1:2)=y_x(m_local,n_anfang:n_ende,1:2)[i,n_pos,o_rank]
		ELSE
			a_x_s(m_pos)=a_x(m_local)
			f_x_s(m_pos)=f_x(m_local)
			g_x_s(m_pos)=g_x(m_local)
			y_x_s(m_pos,1:n_kop,1:2)=y_x(m_local,n_anfang:n_ende,1:2)
		ENDIF
	ENDDO
	!!SYNC
	SYNC IMAGES(sync_m_array)

!	IF(n_pos==1.AND.m_pos==2)THEN
!		DO i=1,n_kop
!			PRINT *,y_x_s(1,i,1)
!		ENDDO
!	ENDIF

	!!
	!!Gauss
	DO i=1,R-1
		zi=f_x_s(i+1)/a_x_s(i)			
		a_x_s(i+1)=a_x_s(i+1)-zi*g_x_s(i)
		y_x_s(i+1,1:n_kop,1:2)=y_x_s(i+1,1:n_kop,1:2)-zi*y_x_s(i,1:n_kop,1:2)
	ENDDO
	x_x_s(R,1:n_kop,1:2)=y_x_s(R,1:n_kop,1:2)/a_x_s(R)
	DO i=R-1,1,-1
		x_x_s(i,1:n_kop,1:2)=(y_x_s(i,1:n_kop,1:2)-g_x_s(i)*x_x_s(i+1,1:n_kop,1:2))/a_x_s(i)
	ENDDO

!	IF(n_pos==1.AND.m_pos==1)THEN
!		DO h=1,n_kop
!			DO i=1,R			
!				PRINT *,x_x_s(i,h,1)
!			ENDDO
!		ENDDO
!	ENDIF

	!!SYNC
	SYNC IMAGES(sync_m_array)
	!!
	!! Transfer
	DO i=1,R
		IF(i/=m_pos)THEN
			IF(i/=1)x_x(0,n_anfang:n_ende,1:2)[i,n_pos,o_rank]=x_x_s(i-1,1:n_kop,1:2)
			x_x(m_local,n_anfang:n_ende,1:2)[i,n_pos,o_rank]=x_x_s(i,1:n_kop,1:2)
		ELSE
			IF(m_pos/=1)x_x(0,n_anfang:n_ende,1:2)=x_x_s(m_pos-1,1:n_kop,1:2)
			x_x(m_local,n_anfang:n_ende,1:2)=x_x_s(m_pos,1:n_kop,1:2)
		ENDIF
	ENDDO

	!!SYNC
	SYNC IMAGES(sync_m_array)
	!!
	
	!!Schritt 4
	IF(m_pos==1)THEN
		DO i=1,m_local-1
			x_x(i,1:n_local,1:2)=(y_x(i,1:n_local,1:2)-g_x(i)*x_x(m_local,1:n_local,1:2))/a_x(i)
		ENDDO
	ELSE
		DO i=1,m_local-1
			x_x(i,1:n_local,1:2)=&
			(y_x(i,1:n_local,1:2)-f_x(i)*x_x(0,1:n_local,1:2)-g_x(i)*x_x(m_local,1:n_local,1:2))/a_x(i)
		ENDDO
	ENDIF

	!IF(m_pos==1.AND.n_pos==1)THEN
	!	DO h=1,n_local
	!		PRINT *,'x:'
	!		DO i=1,m_local
	!			PRINT *,x_x(i,h,2)
	!		ENDDO
	!	ENDDO
	!ENDIF	

ENDIF

!	print *,'x_x(m_local,1:n_local,2)',x_x(m_local,1:n_local,2),'m_pos: ',m_pos,'n_pos: ',n_pos	

!	!!Sherman-Morrison-Formel
!	alpha_x(1:n_local,1:2)=0.
!
!	DO i=1,m_local
!		alpha_x(1:n_local,1:2)=alpha_x(1:n_local,1:2)+z_x(i)*x_x(i,1:n_local,1:2)
!	ENDDO
!!	print *,'alpha_x(1:n_local,1:2): ',alpha_x(1:n_local,1),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	!!SYNC
!	IF(m_pos==1)THEN
!		SYNC IMAGES(sync_m_array)
!	ELSE
!		SYNC IMAGES(image_m_1)
!	ENDIF
!	!!
!	IF(m_pos==1)THEN
!		DO i=2,R
!!	print *,'KHODAAAAAAAAAAAAAAAAAAAAAAAAAAA'
!			alpha_x(1:n_local,1:2)=alpha_x(1:n_local,1:2)+alpha_x(1:n_local,1:2)[i,n_pos]
!		ENDDO
!		beta_x(1:n_local)=alpha_x(1:n_local,1)/(1.-alpha_x(1:n_local,2))
!		!IF(n_pos==2)THEN
!		!	DO h=1,n_local
!		!		PRINT *,beta_x(h)
!		!	ENDDO
!		!ENDIF
!	ENDIF
!	!!SYNC
!	IF(m_pos==1)THEN
!		SYNC IMAGES(sync_m_array)
!	ELSE
!		SYNC IMAGES(image_m_1)
!	ENDIF
!	!!
!	IF(m_pos>1)beta_x(1:n_local)=beta_x(1:n_local)[1,n_pos]	
!	!!SYNC
!	IF(m_pos==1)THEN
!		SYNC IMAGES(sync_m_array)
!	ELSE
!		SYNC IMAGES(image_m_1)
!	ENDIF
!	!!
!	DO i=1,n_local
!		x_x(1:m_local,i,1)=x_x(1:m_local,i,1)+beta_x(i)*x_x(1:m_local,i,2)
!		!IF(m_pos==1.AND.n_pos==1)THEN
!		!	PRINT *,'x: '
!		!	DO h=1,m_local
!		!		PRINT *,x_x(h,i,1)
!		!	ENDDO
!		!ENDIF
!	ENDDO	
	!!SYNC
	SYNC IMAGES(sync_m_array)

	RETURN

END SUBROUTINE

SUBROUTINE WANG_JOHNSSON_P_x_2(o_rank)

	IMPLICIT NONE

	INTEGER::i, o_rank,h,j
	REAL*8::hi,zi
	REAL*8::t_n,s_n
	INTEGER::n_anfang,n_ende

	!!SYNC
	SYNC IMAGES(sync_m_array)
	!!	
    !alpha_i=1.
    !beta_i=3.
    !gamma_i=1.
    !q_i=12.
    !a_i=-1.
    !b_i=-28.
    !c_i=0.
    !d_i=28.
    !e_i=1.
	!Belegung Matrix-Elemente
	IF(R==1.AND.S==1)THEN
	a_x(3:m_local-2)=beta_i
	a_x(1)=1.
	a_x(2)=4.
	a_x(m_local-1)=4.
	a_x(m_local)=1.
	b_x(1)=0.
	b_x(2:m_local-1)=alpha_i
	b_x(m_local)=0.
	c_x(1)=2.
	c_x(2:m_local-1)=gamma_i
	c_x(m_local)=0.


!	!Vektoren für Sherman-Morrison-Formel
!	a_x(1)=a_x(1)+s_n
!	a_x(m_local)=a_x(m_local)+t_n
!	y_x(1,1:n_local,2)=1.
!	y_x(2:m_local-1,1:n_local,2)=0.
!	y_x(m_local,1:n_local,2)=-1.
!	z_x(1)=s_n
!	z_x(2:m_local-1)=0.
!	z_x(m_local)=-t_n

!	print *,'b_x(1:m_local): ',b_x(1:m_local),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	print *,'a_x(1:m_local): ',a_x(1:m_local),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	print *,'c_x(1:m_local): ',c_x(1:m_local),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	print *,'z_x(1:m_local): ',z_x(1:m_local),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	print *,'y_x(1:m_local,i:n_local,2): ',y_x(1:m_local,1:n_local,2)
	DO i=2,m_local
		hi=b_x(i)/a_x(i-1)
		a_x(i)=a_x(i)-hi*c_x(i-1)
		y_x(i,1:n_local,1)=y_x(i,1:n_local,1)-hi*y_x(i-1,1:n_local,1)
	ENDDO



	x_x(m_local,1:n_local,1)=y_x(m_local,1:n_local,1)/a_x(m_local)

	Do i=m_local-1,1,-1
		x_x(i,1:n_local,1:2)=(y_x(i,1:n_local,1:2)-c_x(i)*x_x(i+1,1:n_local,1:2))/a_x(i)
	ENDDO

	ELSE
	IF(m_pos==1)THEN
		a_x(1)=1.
		a_x(2)=4.
		a_x(3:m_local)=beta_i
		b_x(1)=0.
		b_x(2:m_local)=alpha_i
		c_x(1)=2.
		c_x(2:m_local)=gamma_i
	ELSEIF(m_pos==R)THEN
		a_x(1:m_local-2)=beta_i
		a_x(m_local-1)=4.
		a_x(m_local)=1.
		b_x(1:m_local-1)=alpha_i
		b_x(m_local)=2.
		c_x(1:m_local-1)=gamma_i
		c_x(m_local)=0.
	ELSE
		a_x(1:m_local)=beta_i
		b_x(1:m_local)=alpha_i
		c_x(1:m_local)=gamma_i
	ENDIF	
!Belegung Matrix-Elemente
!	a_x(1:m_local)=beta_i
!	IF(m_pos==1)THEN
!		b_x(1)=0.
!		b_x(2:m_local)=alpha_i
!		c_x(1:m_local)=gamma_i
!		s_n=gamma_i
!	ELSEIF(m_pos==R)THEN
!		b_x(1:m_local)=alpha_i
!		c_x(1:m_local-1)=gamma_i
!		c_x(m_local)=0.
!		t_n=alpha_i
!	ELSE
!		b_x(1:m_local)=alpha_i
!		c_x(1:m_local)=gamma_i
!	ENDIF	
	
!	!Vektoren für Sherman-Morrison-Formel
!	IF(m_pos==1)THEN
!		a_x(1)=a_x(1)+s_n
!		y_x(1,1:n_local,2)=1.
!		y_x(2:m_local,1:n_local,2)=0.
!		z_x(1)=s_n
!		z_x(2:m_local)=0.
!	ELSEIF(m_pos==R)THEN
!		a_x(m_local)=a_x(m_local)+t_n
!		y_x(1:m_local-1,1:n_local,2)=0.
!		y_x(m_local,1:n_local,2)=-1.
!		z_x(1:m_local-1)=0.
!		z_x(m_local)=-t_n
!	ELSE
!		y_x(1:m_local,1:n_local,2)=0.
!		z_x(1:m_local)=0.
!	ENDIF	

	!!Schritt1
	IF(m_pos>1)f_x(1)=b_x(1)
	DO i=2,m_local
		hi=b_x(i)/a_x(i-1)
		a_x(i)=a_x(i)-hi*c_x(i-1)
		y_x(i,1:n_local,1:2)=y_x(i,1:n_local,1:2)-hi*y_x(i-1,1:n_local,1:2)
		IF(m_pos>1)f_x(i)=-hi*f_x(i-1)
	ENDDO
	
	!SYNC ALL
	!IF(m_pos==2.AND.n_pos==1)THEN
	!	PRINT *,' '
	!	PRINT *,'m_pos2'
	!	DO i=1,m_local
	!		PRINT *,y_x(i,1,1)
	!	ENDDO
	!ENDIF

	!!Schritt2
	g_x(m_local-1)=c_x(m_local-1)
	DO i=m_local-1,2,-1
		hi=c_x(i-1)/a_x(i)
		g_x(i-1)=-hi*g_x(i)
		y_x(i-1,1:n_local,1:2)=y_x(i-1,1:n_local,1:2)-hi*y_x(i,1:n_local,1:2)
		IF(m_pos>1)f_x(i-1)=f_x(i-1)-hi*f_x(i)
	ENDDO

	!!SYNC
	IF(m_pos==R)THEN
		SYNC IMAGES(sync_m_down)
	ELSEIF(m_pos==1)THEN
		SYNC IMAGES(sync_m_up)
	ELSE
		SYNC IMAGES([sync_m_up,sync_m_down])
	ENDIF
	!! Transfer
	IF(m_pos<R)THEN
		a_x(m_local+1)=a_x(1)[m_pos+1,n_pos,o_rank]
		f_x(m_local+1)=f_x(1)[m_pos+1,n_pos,o_rank]
		g_x(m_local+1)=g_x(1)[m_pos+1,n_pos,o_rank]
		y_x(m_local+1,1:n_local,1:2)=y_x(1,1:n_local,1:2)[m_pos+1,n_pos,o_rank]
	ENDIF
	!!SYNC
	IF(m_pos==R)THEN
		SYNC IMAGES(sync_m_down)
	ELSEIF(m_pos==1)THEN
		SYNC IMAGES(sync_m_up)
	ELSE
		SYNC IMAGES([sync_m_up,sync_m_down])
	ENDIF
	!!
	IF(m_pos<R)THEN
		hi=c_x(m_local)/a_x(m_local+1)
		a_x(m_local)=a_x(m_local)-hi*f_x(m_local+1)
		g_x(m_local)=-hi*g_x(m_local+1)
		y_x(m_local,1:n_local,1:2)=y_x(m_local,1:n_local,1:2)-hi*y_x(m_local+1,1:n_local,1:2)
	ENDIF

!	SYNC ALL
!	IF(m_pos==1.AND.n_pos==1)THEN
!		DO h=1,n_local
!			!PRINT *,'m_pos1'
!			!DO i=1,m_local			
!				PRINT *,y_x(m_local,h,1)
!			!ENDDO
!		ENDDO
!	ENDIF
	

	!!Schritt3: Kopplungssystem
	!!SYNC
	SYNC IMAGES(sync_m_array)
	!!
	!!Transfer
	n_anfang=(m_pos-1)*n_kop+1
	n_ende=m_pos*n_kop
	DO i=1,R
		IF(i/=m_pos)THEN
			a_x_s(i)=a_x(m_local)[i,n_pos,o_rank]
			f_x_s(i)=f_x(m_local)[i,n_pos,o_rank]
			g_x_s(i)=g_x(m_local)[i,n_pos,o_rank]
			y_x_s(i,1:n_kop,1:2)=y_x(m_local,n_anfang:n_ende,1:2)[i,n_pos,o_rank]
		ELSE
			a_x_s(m_pos)=a_x(m_local)
			f_x_s(m_pos)=f_x(m_local)
			g_x_s(m_pos)=g_x(m_local)
			y_x_s(m_pos,1:n_kop,1:2)=y_x(m_local,n_anfang:n_ende,1:2)
		ENDIF
	ENDDO
	!!SYNC
	SYNC IMAGES(sync_m_array)

!	IF(n_pos==1.AND.m_pos==2)THEN
!		DO i=1,n_kop
!			PRINT *,y_x_s(1,i,1)
!		ENDDO
!	ENDIF

	!!
	!!Gauss
	DO i=1,R-1
		zi=f_x_s(i+1)/a_x_s(i)			
		a_x_s(i+1)=a_x_s(i+1)-zi*g_x_s(i)
		y_x_s(i+1,1:n_kop,1:2)=y_x_s(i+1,1:n_kop,1:2)-zi*y_x_s(i,1:n_kop,1:2)
	ENDDO
	x_x_s(R,1:n_kop,1:2)=y_x_s(R,1:n_kop,1:2)/a_x_s(R)
	DO i=R-1,1,-1
		x_x_s(i,1:n_kop,1:2)=(y_x_s(i,1:n_kop,1:2)-g_x_s(i)*x_x_s(i+1,1:n_kop,1:2))/a_x_s(i)
	ENDDO

!	IF(n_pos==1.AND.m_pos==1)THEN
!		DO h=1,n_kop
!			DO i=1,R			
!				PRINT *,x_x_s(i,h,1)
!			ENDDO
!		ENDDO
!	ENDIF

	!!SYNC
	SYNC IMAGES(sync_m_array)
	!!
	!! Transfer
	DO i=1,R
		IF(i/=m_pos)THEN
			IF(i/=1)x_x(0,n_anfang:n_ende,1:2)[i,n_pos,o_rank]=x_x_s(i-1,1:n_kop,1:2)
			x_x(m_local,n_anfang:n_ende,1:2)[i,n_pos,o_rank]=x_x_s(i,1:n_kop,1:2)
		ELSE
			IF(m_pos/=1)x_x(0,n_anfang:n_ende,1:2)=x_x_s(m_pos-1,1:n_kop,1:2)
			x_x(m_local,n_anfang:n_ende,1:2)=x_x_s(m_pos,1:n_kop,1:2)
		ENDIF
	ENDDO

	!!SYNC
	SYNC IMAGES(sync_m_array)
	!!
	
	!!Schritt 4
	IF(m_pos==1)THEN
		DO i=1,m_local-1
			x_x(i,1:n_local,1:2)=(y_x(i,1:n_local,1:2)-g_x(i)*x_x(m_local,1:n_local,1:2))/a_x(i)
		ENDDO
	ELSE
		DO i=1,m_local-1
			x_x(i,1:n_local,1:2)=&
			(y_x(i,1:n_local,1:2)-f_x(i)*x_x(0,1:n_local,1:2)-g_x(i)*x_x(m_local,1:n_local,1:2))/a_x(i)
		ENDDO
	ENDIF

	!IF(m_pos==1.AND.n_pos==1)THEN
	!	DO h=1,n_local
	!		PRINT *,'x:'
	!		DO i=1,m_local
	!			PRINT *,x_x(i,h,2)
	!		ENDDO
	!	ENDDO
	!ENDIF	

ENDIF

!	print *,'x_x(m_local,1:n_local,2)',x_x(m_local,1:n_local,2),'m_pos: ',m_pos,'n_pos: ',n_pos	

!	!!Sherman-Morrison-Formel
!	alpha_x(1:n_local,1:2)=0.
!
!	DO i=1,m_local
!		alpha_x(1:n_local,1:2)=alpha_x(1:n_local,1:2)+z_x(i)*x_x(i,1:n_local,1:2)
!	ENDDO
!!	print *,'alpha_x(1:n_local,1:2): ',alpha_x(1:n_local,1),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	!!SYNC
!	IF(m_pos==1)THEN
!		SYNC IMAGES(sync_m_array)
!	ELSE
!		SYNC IMAGES(image_m_1)
!	ENDIF
!	!!
!	IF(m_pos==1)THEN
!		DO i=2,R
!!	print *,'KHODAAAAAAAAAAAAAAAAAAAAAAAAAAA'
!			alpha_x(1:n_local,1:2)=alpha_x(1:n_local,1:2)+alpha_x(1:n_local,1:2)[i,n_pos]
!		ENDDO
!		beta_x(1:n_local)=alpha_x(1:n_local,1)/(1.-alpha_x(1:n_local,2))
!		!IF(n_pos==2)THEN
!		!	DO h=1,n_local
!		!		PRINT *,beta_x(h)
!		!	ENDDO
!		!ENDIF
!	ENDIF
!	!!SYNC
!	IF(m_pos==1)THEN
!		SYNC IMAGES(sync_m_array)
!	ELSE
!		SYNC IMAGES(image_m_1)
!	ENDIF
!	!!
!	IF(m_pos>1)beta_x(1:n_local)=beta_x(1:n_local)[1,n_pos]	
!	!!SYNC
!	IF(m_pos==1)THEN
!		SYNC IMAGES(sync_m_array)
!	ELSE
!		SYNC IMAGES(image_m_1)
!	ENDIF
!	!!
!	DO i=1,n_local
!		x_x(1:m_local,i,1)=x_x(1:m_local,i,1)+beta_x(i)*x_x(1:m_local,i,2)
!		!IF(m_pos==1.AND.n_pos==1)THEN
!		!	PRINT *,'x: '
!		!	DO h=1,m_local
!		!		PRINT *,x_x(h,i,1)
!		!	ENDDO
!		!ENDIF
!	ENDDO	
	!!SYNC
	SYNC IMAGES(sync_m_array)

	RETURN

END SUBROUTINE

SUBROUTINE WANG_JOHNSSON_P_y(o_rank)

	IMPLICIT NONE

	INTEGER::i,o_rank
	REAL*8::hi,zi
	REAL*8::t_n,s_n
	INTEGER::m_anfang,m_ende

	!!SYNC
	SYNC IMAGES(sync_n_array)
	!!	
	IF(R==1.AND.S==1)THEN
    !!
!	alpha_i=1.
!	beta_i=3.
!	gamma_i=1.
        a_y(3:n_local-2)=beta_i
        a_y(1)=1.
        a_y(2)=4.
        a_y(n_local-1)=4.
        a_y(n_local)=1.
        b_y(1)=0.
        b_y(2:n_local-1)=alpha_i
        b_y(n_local)=0.
        c_y(1)=2.
        c_y(2:n_local-1)=gamma_i
        c_y(n_local)=0.

   !     	print *,'a_y(1:n_local): ',a_y(1:n_local),'m_pos: ',m_pos,'n_pos: ',n_pos
   !     	print *,'b_y(1:n_local): ',b_y(1:n_local),'m_pos: ',m_pos,'n_pos: ',n_pos
   !     	print *,'c_y(1:n_local): ',c_y(1:n_local),'m_pos: ',m_pos,'n_pos: ',n_pos
   !     !	print *,'z_y(1:n_local): ',z_y(1:n_local),'m_pos: ',m_pos,'n_pos: ',n_pos
   !     	print *,'y_y(1:m_local,1:n_local,1): ',y_y(1:m_local,1:n_local,1)

        DO i=2,n_local
            hi=b_y(i)/a_y(i-1)
            a_y(i)=a_y(i)-hi*c_y(i-1)
            y_y(1:m_local,i,1)=y_y(1:m_local,i,1)-hi*y_y(1:m_local,i-1,1)
        ENDDO



            x_y(1:m_local,n_local,1)=y_y(1:m_local,n_local,1)/a_y(n_local)

        Do i=n_local-1,1,-1
            x_y(1:m_local,i,1)=(y_y(1:m_local,i,1)-c_y(i)*x_y(1:m_local,i+1,1))/a_y(i)
        ENDDO

ELSE
!Belegung Matrix-Elemente
!        a_y(3:n_local-2)=beta_i
!         a_y(1)=1.
!         a_y(2)=4.
!        a_y(n_local-1)=4.
!        a_y(n_local)=1.
!         b_y(1)=0.
!         b_y(2:n_local-1)=alpha_i
!         b_y(n_local)=2.
!         c_y(1)=2.
!         c_y(2:n_local-1)=gamma_i
!        c_y(n_local)=0.

!	a_y(1:n_local)=beta_i
	IF(n_pos==1)THEN
		a_y(1)=1.
		a_y(2)=4.
		a_y(3:n_local)=beta_i
		b_y(1)=0.
		b_y(2:n_local)=alpha_i
		c_y(1)=2.
		c_y(2:n_local)=gamma_i
	ELSEIF(n_pos==R)THEN
		a_y(1:n_local-2)=beta_i
		a_y(n_local-1)=4.
		a_y(n_local)=1.
		b_y(1:n_local-1)=alpha_i
		b_y(n_local)=2.
		c_y(1:n_local-1)=gamma_i
		c_y(n_local)=0.
	ELSE
		a_y(1:n_local)=beta_i
		b_y(1:n_local)=alpha_i
		c_y(1:n_local)=gamma_i
	ENDIF	

!        	print *,'a_y(1:n_local): ',a_y(1:n_local),'m_pos: ',m_pos,'n_pos: ',n_pos
!        	print *,'b_y(1:n_local): ',b_y(1:n_local),'m_pos: ',m_pos,'n_pos: ',n_pos
!        	print *,'c_y(1:n_local): ',c_y(1:n_local),'m_pos: ',m_pos,'n_pos: ',n_pos
!        !	print *,'z_y(1:n_local): ',z_y(1:n_local),'m_pos: ',m_pos,'n_pos: ',n_pos
!        	print *,'y_y(1:m_local,1:n_local,1): ',y_y(1:m_local,1:n_local,1)
	!Vektoren für Sherman-Morrison-Formel
!	IF(n_pos==1)THEN
!		a_y(1)=a_y(1)+s_n
!		y_y(1:m_local,1,2)=1.
!		y_y(1:m_local,2:n_local,2)=0.
!		z_y(1)=s_n
!		z_y(2:n_local)=0.
!	ELSEIF(n_pos==R)THEN
!		a_y(n_local)=a_y(n_local)+t_n
!		y_y(1:m_local,1:n_local-1,2)=0.
!		y_y(1:m_local,n_local,2)=-1.
!		z_y(1:n_local-1)=0.
!		z_y(n_local)=-t_n
!	ELSE
!		y_y(1:m_local,1:n_local,2)=0.
!		z_y(1:n_local)=0.
!	ENDIF	

	!!Schritt1
	IF(n_pos>1)f_y(1)=b_y(1)
	DO i=2,n_local
		hi=b_y(i)/a_y(i-1)
		a_y(i)=a_y(i)-hi*c_y(i-1)
		y_y(1:m_local,i,1)=y_y(1:m_local,i,1)-hi*y_y(1:m_local,i-1,1)
		IF(n_pos>1)f_y(i)=-hi*f_y(i-1)
	ENDDO

	!!Schritt2
	g_y(n_local-1)=c_y(n_local-1)
	DO i=n_local-1,2,-1
		hi=c_y(i-1)/a_y(i)
		g_y(i-1)=-hi*g_y(i)
		y_y(1:m_local,i-1,1)=y_y(1:m_local,i-1,1)-hi*y_y(1:m_local,i,1)
		IF(n_pos>1)f_y(i-1)=f_y(i-1)-hi*f_y(i)
	ENDDO

	!!SYNC
	IF(n_pos==S)THEN
		SYNC IMAGES(sync_n_down)
	ELSEIF(n_pos==1)THEN
		SYNC IMAGES(sync_n_up)
	ELSE
		SYNC IMAGES([sync_n_up,sync_n_down])
	ENDIF
	!! Transfer
	IF(n_pos<S)THEN
		a_y(n_local+1)=a_y(1)[m_pos,n_pos+1,o_rank]
		f_y(n_local+1)=f_y(1)[m_pos,n_pos+1,o_rank]
		g_y(n_local+1)=g_y(1)[m_pos,n_pos+1,o_rank]
		y_y(1:m_local,n_local+1,1)=y_y(1:m_local,1,1)[m_pos,n_pos+1,o_rank]
	ENDIF
	!!SYNC
	IF(n_pos==S)THEN
		SYNC IMAGES(sync_n_down)
	ELSEIF(n_pos==1)THEN
		SYNC IMAGES(sync_n_up)
	ELSE
		SYNC IMAGES([sync_n_up,sync_n_down])
	ENDIF
	!!
	IF(n_pos<S)THEN
		hi=c_y(n_local)/a_y(n_local+1)
		a_y(n_local)=a_y(n_local)-hi*f_y(n_local+1)
		g_y(n_local)=-hi*g_y(n_local+1)
		y_y(1:m_local,n_local,1)=y_y(1:m_local,n_local,1)-hi*y_y(1:m_local,n_local+1,1)
	ENDIF

	!!Schritt3: Kopplungssystem
	!!SYNC
	SYNC IMAGES(sync_n_array)
	!!
	!!Transfer
	m_anfang=(n_pos-1)*m_kop+1
	m_ende=n_pos*m_kop
	DO i=1,S
		IF(i/=n_pos)THEN
			a_y_s(i)=a_y(n_local)[m_pos,i,o_rank]
			f_y_s(i)=f_y(n_local)[m_pos,i,o_rank]
			g_y_s(i)=g_y(n_local)[m_pos,i,o_rank]
			y_y_s(1:m_kop,i,1)=y_y(m_anfang:m_ende,n_local,1)[m_pos,i,o_rank]
		ELSE
			a_y_s(n_pos)=a_y(n_local)
			f_y_s(n_pos)=f_y(n_local)
			g_y_s(n_pos)=g_y(n_local)
			y_y_s(1:m_kop,n_pos,1)=y_y(m_anfang:m_ende,n_local,1)
		ENDIF
	ENDDO
	!!SYNC
	SYNC IMAGES(sync_n_array)
	!!
	!!Gauss
	DO i=1,S-1
		zi=f_y_s(i+1)/a_y_s(i)			
		a_y_s(i+1)=a_y_s(i+1)-zi*g_y_s(i)
		y_y_s(1:m_kop,i+1,1)=y_y_s(1:m_kop,i+1,1)-zi*y_y_s(1:m_kop,i,1)
	ENDDO
	x_y_s(1:m_kop,S,1)=y_y_s(1:m_kop,S,1)/a_y_s(S)
	DO i=S-1,1,-1
		x_y_s(1:m_kop,i,1)=(y_y_s(1:m_kop,i,1)-g_y_s(i)*x_y_s(1:m_kop,i+1,1))/a_y_s(i)
	ENDDO
	!!SYNC
	SYNC IMAGES(sync_n_array)
	!!
	!!
	!! Transfer
	DO i=1,S
		IF(i/=n_pos)THEN
			IF(i/=1)x_y(m_anfang:m_ende,0,1)[m_pos,i,o_rank]=x_y_s(1:m_kop,i-1,1)
			x_y(m_anfang:m_ende,n_local,1)[m_pos,i,o_rank]=x_y_s(1:m_kop,i,1)
		ELSE
			IF(n_pos/=1)x_y(m_anfang:m_ende,0,1)=x_y_s(1:m_kop,n_pos-1,1)
			x_y(m_anfang:m_ende,n_local,1)=x_y_s(1:m_kop,n_pos,1)
		ENDIF
	ENDDO
	!!SYNC
	SYNC IMAGES(sync_n_array)
	!!	
	!!Schritt 4
	IF(n_pos==1)THEN
		DO i=1,n_local-1
			x_y(1:m_local,i,1)=(y_y(1:m_local,i,1)-g_y(i)*x_y(1:m_local,n_local,1))/a_y(i)
		ENDDO
	ELSE
		DO i=1,n_local-1
			x_y(1:m_local,i,1)=&
			(y_y(1:m_local,i,1)-f_y(i)*x_y(1:m_local,0,1)-g_y(i)*x_y(1:m_local,n_local,1))/a_y(i)
		ENDDO
	ENDIF


ENDIF
!	!!Sherman-Morrison-Formel
!	alpha_y(1:m_local,1:2)=0.
!	DO i=1,n_local
!		alpha_y(1:m_local,1:2)=alpha_y(1:m_local,1:2)+z_y(i)*x_y(1:m_local,i,1:2)
!	ENDDO
!	!!SYNC
!	IF(n_pos==1)THEN
!		SYNC IMAGES(sync_n_array)
!	ELSE
!		SYNC IMAGES(image_n_1)
!	ENDIF
!	!!
!	IF(n_pos==1)THEN
!		DO i=2,S
!			alpha_y(1:m_local,1:2)=alpha_y(1:m_local,1:2)+alpha_y(1:m_local,1:2)[m_pos,i]
!		ENDDO
!		beta_y(1:m_local)=alpha_y(1:m_local,1)/(1.-alpha_y(1:m_local,2))
!	ENDIF
!	!!SYNC
!	IF(n_pos==1)THEN
!		SYNC IMAGES(sync_n_array) ! jauts do it fo th axamo: 
!	ELSE
!		SYNC IMAGES(image_n_1)
!	ENDIF
!	!!
!	IF(n_pos>1)beta_y(1:m_local)=beta_y(1:m_local)[m_pos,1]	
!	!!SYNC
!	IF(n_pos==1)THEN
!		SYNC IMAGES(sync_n_array)
!	ELSE
!		SYNC IMAGES(image_n_1)
!	ENDIF
!	!!
!	DO i=1,m_local
!		x_y(i,1:n_local,1)=x_y(i,1:n_local,1)+beta_y(i)*x_y(i,1:n_local,2)
!	ENDDO	
	!!SYNC
	SYNC IMAGES(sync_n_array)

	RETURN

END SUBROUTINE


!Lösung eines zyklischen Tridiagonal-Systems im Parallel-Fall mit Wang-Johnsson und Sherman-Morrison-Formel
!in y
SUBROUTINE WANG_JOHNSSON_P_z(o_rank)

	IMPLICIT NONE

	INTEGER::i,o_rank
	REAL*8::hi,zi
	REAL*8::t_n,s_n
	INTEGER::m_anfang,m_ende

	!!SYNC
	SYNC IMAGES(sync_n_array)
	!!	

	!Belegung Matrix-Elemente
	a_y(1:n_local)=beta_i
	IF(n_pos==1)THEN
		b_y(1)=0.
		b_y(2:n_local)=alpha_i
		c_y(1:n_local)=gamma_i
		s_n=gamma_i
	ELSEIF(n_pos==R)THEN
		b_y(1:n_local)=alpha_i
		c_y(1:n_local-1)=gamma_i
		c_y(n_local)=0.
		t_n=alpha_i
	ELSE
		b_y(1:n_local)=alpha_i
		c_y(1:n_local)=gamma_i
	ENDIF	
	
	!Vektoren für Sherman-Morrison-Formel
	IF(n_pos==1)THEN
		a_y(1)=a_y(1)+s_n
		y_y(1:m_local,1,2)=1.
		y_y(1:m_local,2:n_local,2)=0.
		z_y(1)=s_n
		z_y(2:n_local)=0.
	ELSEIF(n_pos==R)THEN
		a_y(n_local)=a_y(n_local)+t_n
		y_y(1:m_local,1:n_local-1,2)=0.
		y_y(1:m_local,n_local,2)=-1.
		z_y(1:n_local-1)=0.
		z_y(n_local)=-t_n
	ELSE
		y_y(1:m_local,1:n_local,2)=0.
		z_y(1:n_local)=0.
	ENDIF	

	!!Schritt1
	IF(n_pos>1)f_y(1)=b_y(1)
	DO i=2,n_local
		hi=b_y(i)/a_y(i-1)
		a_y(i)=a_y(i)-hi*c_y(i-1)
		y_y(1:m_local,i,1:2)=y_y(1:m_local,i,1:2)-hi*y_y(1:m_local,i-1,1:2)
		IF(n_pos>1)f_y(i)=-hi*f_y(i-1)
	ENDDO

	!!Schritt2
	g_y(n_local-1)=c_y(n_local-1)
	DO i=n_local-1,2,-1
		hi=c_y(i-1)/a_y(i)
		g_y(i-1)=-hi*g_y(i)
		y_y(1:m_local,i-1,1:2)=y_y(1:m_local,i-1,1:2)-hi*y_y(1:m_local,i,1:2)
		IF(n_pos>1)f_y(i-1)=f_y(i-1)-hi*f_y(i)
	ENDDO

	!!SYNC
	IF(n_pos==S)THEN
		SYNC IMAGES(sync_n_down)
	ELSEIF(n_pos==1)THEN
		SYNC IMAGES(sync_n_up)
	ELSE
		SYNC IMAGES([sync_n_up,sync_n_down])
	ENDIF
	!! Transfer
	IF(n_pos<S)THEN
		a_y(n_local+1)=a_y(1)[m_pos,n_pos+1,o_rank]
		f_y(n_local+1)=f_y(1)[m_pos,n_pos+1,o_rank]
		g_y(n_local+1)=g_y(1)[m_pos,n_pos+1,o_rank]
		y_y(1:m_local,n_local+1,1:2)=y_y(1:m_local,1,1:2)[m_pos,n_pos+1,o_rank]
	ENDIF
	!!SYNC
	IF(n_pos==S)THEN
		SYNC IMAGES(sync_n_down)
	ELSEIF(n_pos==1)THEN
		SYNC IMAGES(sync_n_up)
	ELSE
		SYNC IMAGES([sync_n_up,sync_n_down])
	ENDIF
	!!
	IF(n_pos<S)THEN
		hi=c_y(n_local)/a_y(n_local+1)
		a_y(n_local)=a_y(n_local)-hi*f_y(n_local+1)
		g_y(n_local)=-hi*g_y(n_local+1)
		y_y(1:m_local,n_local,1:2)=y_y(1:m_local,n_local,1:2)-hi*y_y(1:m_local,n_local+1,1:2)
	ENDIF

	!!Schritt3: Kopplungssystem
	!!SYNC
	SYNC IMAGES(sync_n_array)
	!!
	!!Transfer
	m_anfang=(n_pos-1)*m_kop+1
	m_ende=n_pos*m_kop
	DO i=1,S
		IF(i/=n_pos)THEN
			a_y_s(i)=a_y(n_local)[m_pos,i,o_rank]
			f_y_s(i)=f_y(n_local)[m_pos,i,o_rank]
			g_y_s(i)=g_y(n_local)[m_pos,i,o_rank]
			y_y_s(1:m_kop,i,1:2)=y_y(m_anfang:m_ende,n_local,1:2)[m_pos,i,o_rank]
		ELSE
			a_y_s(n_pos)=a_y(n_local)
			f_y_s(n_pos)=f_y(n_local)
			g_y_s(n_pos)=g_y(n_local)
			y_y_s(1:m_kop,n_pos,1:2)=y_y(m_anfang:m_ende,n_local,1:2)
		ENDIF
	ENDDO
	!!SYNC
	SYNC IMAGES(sync_n_array)
	!!
	!!Gauss
	DO i=1,S-1
		zi=f_y_s(i+1)/a_y_s(i)			
		a_y_s(i+1)=a_y_s(i+1)-zi*g_y_s(i)
		y_y_s(1:m_kop,i+1,1:2)=y_y_s(1:m_kop,i+1,1:2)-zi*y_y_s(1:m_kop,i,1:2)
	ENDDO
	x_y_s(1:m_kop,S,1:2)=y_y_s(1:m_kop,S,1:2)/a_y_s(S)
	DO i=S-1,1,-1
		x_y_s(1:m_kop,i,1:2)=(y_y_s(1:m_kop,i,1:2)-g_y_s(i)*x_y_s(1:m_kop,i+1,1:2))/a_y_s(i)
	ENDDO
	!!SYNC
	SYNC IMAGES(sync_n_array)
	!!
	!! Transfer
	DO i=1,S
		IF(i/=n_pos)THEN
			IF(i/=1)x_y(m_anfang:m_ende,0,1:2)[m_pos,i,o_rank]=x_y_s(1:m_kop,i-1,1:2)
			x_y(m_anfang:m_ende,n_local,1:2)[m_pos,i,o_rank]=x_y_s(1:m_kop,i,1:2)
		ELSE
			IF(n_pos/=1)x_y(m_anfang:m_ende,0,1:2)=x_y_s(1:m_kop,n_pos-1,1:2)
			x_y(m_anfang:m_ende,n_local,1:2)=x_y_s(1:m_kop,n_pos,1:2)
		ENDIF
	ENDDO
	!!SYNC
	SYNC IMAGES(sync_n_array)
	!!	
	!!Schritt 4
	IF(n_pos==1)THEN
		DO i=1,n_local-1
			x_y(1:m_local,i,1:2)=(y_y(1:m_local,i,1:2)-g_y(i)*x_y(1:m_local,n_local,1:2))/a_y(i)
		ENDDO
	ELSE
		DO i=1,n_local-1
			x_y(1:m_local,i,1:2)=&
			(y_y(1:m_local,i,1:2)-f_y(i)*x_y(1:m_local,0,1:2)-g_y(i)*x_y(1:m_local,n_local,1:2))/a_y(i)
		ENDDO
	ENDIF

	!!Sherman-Morrison-Formel
	alpha_y(1:m_local,1:2)=0.
	DO i=1,n_local
		alpha_y(1:m_local,1:2)=alpha_y(1:m_local,1:2)+z_y(i)*x_y(1:m_local,i,1:2)
	ENDDO
	!!SYNC
	IF(n_pos==1)THEN
		SYNC IMAGES(sync_n_array)
	ELSE
		SYNC IMAGES(image_n_1)
	ENDIF
	!!
	IF(n_pos==1)THEN
		DO i=2,S
			alpha_y(1:m_local,1:2)=alpha_y(1:m_local,1:2)+alpha_y(1:m_local,1:2)[m_pos,i,o_rank]
		ENDDO
		beta_y(1:m_local)=alpha_y(1:m_local,1)/(1.-alpha_y(1:m_local,2))
	ENDIF
	!!SYNC
	IF(n_pos==1)THEN
		SYNC IMAGES(sync_n_array)
	ELSE
		SYNC IMAGES(image_n_1)
	ENDIF
	!!
	IF(n_pos>1)beta_y(1:m_local)=beta_y(1:m_local)[m_pos,1,o_rank]	
	!!SYNC
	IF(n_pos==1)THEN
		SYNC IMAGES(sync_n_array)
	ELSE
		SYNC IMAGES(image_n_1)
	ENDIF
	!!
	DO i=1,m_local
		x_y(i,1:n_local,1)=x_y(i,1:n_local,1)+beta_y(i)*x_y(i,1:n_local,2)
	ENDDO	
	!!SYNC
	SYNC IMAGES(sync_n_array)

	RETURN

END SUBROUTINE


SUBROUTINE WANG_JOHNSSON_P_x_neu(o_rank)

	IMPLICIT NONE

	INTEGER::i, o_rank,h,j
	REAL*8::hi,zi
	REAL*8::t_n,s_n
	INTEGER::n_anfang,n_ende

	!!SYNC
	SYNC IMAGES(sync_m_array)
	!!	
    !alpha_i=1.
    !beta_i=3.
    !gamma_i=1.
    !q_i=12.
    !a_i=-1.
    !b_i=-28.
    !c_i=0.
    !d_i=28.
    !e_i=1.
	!Belegung Matrix-Elemente
	IF(R==1.AND.S==1)THEN
	a_x(3:m_local-3)=beta_i
	a_x(1)=1.
	a_x(2)=4.
	a_x(m_local-2)=4.
	a_x(m_local-1)=1.
	b_x(1)=0.
	b_x(2:m_local-2)=alpha_i
	b_x(m_local-1)=2.
	c_x(1)=2.
	c_x(2:m_local-2)=gamma_i
	c_x(m_local-1)=0.




!	!Vektoren für Sherman-Morrison-Formel
!	a_x(1)=a_x(1)+s_n
!	a_x(m_local)=a_x(m_local)+t_n
!	y_x(1,1:n_local,2)=1.
!	y_x(2:m_local-1,1:n_local,2)=0.
!	y_x(m_local,1:n_local,2)=-1.
!	z_x(1)=s_n
!	z_x(2:m_local-1)=0.
!	z_x(m_local)=-t_n

!	print *,'b_x(1:m_local): ',b_x(1:m_local),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	print *,'a_x(1:m_local): ',a_x(1:m_local),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	print *,'c_x(1:m_local): ',c_x(1:m_local),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	print *,'z_x(1:m_local): ',z_x(1:m_local),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	print *,'y_x(1:m_local,i:n_local,2): ',y_x(1:m_local,1:n_local,2)
	DO i=2,m_local-1
		hi=b_x(i)/a_x(i-1)
		a_x(i)=a_x(i)-hi*c_x(i-1)
		y_x(i,1:n_local,1)=y_x(i,1:n_local,1)-hi*y_x(i-1,1:n_local,1)
	ENDDO



	x_x(m_local-1,1:n_local,1)=y_x(m_local,1:n_local,1)/a_x(m_local)

	Do i=m_local-2,1,-1
		x_x(i,1:n_local,1:2)=(y_x(i,1:n_local,1:2)-c_x(i)*x_x(i+1,1:n_local,1:2))/a_x(i)
	ENDDO

	ELSE
	IF(m_pos==1)THEN
		a_x(1)=1.
		a_x(2)=4.
		a_x(3:m_local)=beta_i
		b_x(1)=0.
		b_x(2:m_local)=alpha_i
		c_x(1)=2.
		c_x(2:m_local)=gamma_i
	ELSEIF(m_pos==R)THEN
		a_x(1:m_local-3)=beta_i
		a_x(m_local-2)=4.
		a_x(m_local-1)=1.
		b_x(1:m_local-2)=alpha_i
		b_x(m_local-1)=2.
		c_x(1:m_local-2)=gamma_i
		c_x(m_local-1)=0.
	ELSE
		a_x(1:m_local-1)=beta_i
		b_x(1:m_local-1)=alpha_i
		c_x(1:m_local-1)=gamma_i
	ENDIF	
!Belegung Matrix-Elemente
!	a_x(1:m_local)=beta_i
!	IF(m_pos==1)THEN
!		b_x(1)=0.
!		b_x(2:m_local)=alpha_i
!		c_x(1:m_local)=gamma_i
!		s_n=gamma_i
!	ELSEIF(m_pos==R)THEN
!		b_x(1:m_local)=alpha_i
!		c_x(1:m_local-1)=gamma_i
!		c_x(m_local)=0.
!		t_n=alpha_i
!	ELSE
!		b_x(1:m_local)=alpha_i
!		c_x(1:m_local)=gamma_i
!	ENDIF	
	
!	!Vektoren für Sherman-Morrison-Formel
!	IF(m_pos==1)THEN
!		a_x(1)=a_x(1)+s_n
!		y_x(1,1:n_local,2)=1.
!		y_x(2:m_local,1:n_local,2)=0.
!		z_x(1)=s_n
!		z_x(2:m_local)=0.
!	ELSEIF(m_pos==R)THEN
!		a_x(m_local)=a_x(m_local)+t_n
!		y_x(1:m_local-1,1:n_local,2)=0.
!		y_x(m_local,1:n_local,2)=-1.
!		z_x(1:m_local-1)=0.
!		z_x(m_local)=-t_n
!	ELSE
!		y_x(1:m_local,1:n_local,2)=0.
!		z_x(1:m_local)=0.
!	ENDIF	

	!!Schritt1
	IF(m_pos>1)f_x(1)=b_x(1)
	DO i=2,m_local-1
		hi=b_x(i)/a_x(i-1)
		a_x(i)=a_x(i)-hi*c_x(i-1)
		y_x(i,1:n_local,1:2)=y_x(i,1:n_local,1:2)-hi*y_x(i-1,1:n_local,1:2)
		IF(m_pos>1)f_x(i)=-hi*f_x(i-1)
	ENDDO
	
	!SYNC ALL
	!IF(m_pos==2.AND.n_pos==1)THEN
	!	PRINT *,' '
	!	PRINT *,'m_pos2'
	!	DO i=1,m_local
	!		PRINT *,y_x(i,1,1)
	!	ENDDO
	!ENDIF

	!!Schritt2
	g_x(m_local-2)=c_x(m_local-2)
	DO i=m_local-2,2,-1
		hi=c_x(i-1)/a_x(i)
		g_x(i-1)=-hi*g_x(i)
		y_x(i-1,1:n_local,1:2)=y_x(i-1,1:n_local,1:2)-hi*y_x(i,1:n_local,1:2)
		IF(m_pos>1)f_x(i-1)=f_x(i-1)-hi*f_x(i)
	ENDDO

	!!SYNC
	IF(m_pos==R)THEN
		SYNC IMAGES(sync_m_down)
	ELSEIF(m_pos==1)THEN
		SYNC IMAGES(sync_m_up)
	ELSE
		SYNC IMAGES([sync_m_up,sync_m_down])
	ENDIF
	!! Transfer
	IF(m_pos<R)THEN
		a_x(m_local+1)=a_x(1)[m_pos+1,n_pos,o_rank]
		f_x(m_local+1)=f_x(1)[m_pos+1,n_pos,o_rank]
		g_x(m_local+1)=g_x(1)[m_pos+1,n_pos,o_rank]
		y_x(m_local+1,1:n_local,1:2)=y_x(1,1:n_local,1:2)[m_pos+1,n_pos,o_rank]
	ENDIF
	!!SYNC
	IF(m_pos==R)THEN
		SYNC IMAGES(sync_m_down)
	ELSEIF(m_pos==1)THEN
		SYNC IMAGES(sync_m_up)
	ELSE
		SYNC IMAGES([sync_m_up,sync_m_down])
	ENDIF
	!!
	IF(m_pos<R)THEN
		hi=c_x(m_local)/a_x(m_local+1)
		a_x(m_local)=a_x(m_local)-hi*f_x(m_local+1)
		g_x(m_local)=-hi*g_x(m_local+1)
		y_x(m_local,1:n_local,1:2)=y_x(m_local,1:n_local,1:2)-hi*y_x(m_local+1,1:n_local,1:2)
	ENDIF

!	SYNC ALL
!	IF(m_pos==1.AND.n_pos==1)THEN
!		DO h=1,n_local
!			!PRINT *,'m_pos1'
!			!DO i=1,m_local			
!				PRINT *,y_x(m_local,h,1)
!			!ENDDO
!		ENDDO
!	ENDIF
	

	!!Schritt3: Kopplungssystem
	!!SYNC
	SYNC IMAGES(sync_m_array)
	!!
	!!Transfer
	n_anfang=(m_pos-1)*n_kop+1
	n_ende=m_pos*n_kop
	DO i=1,R
		IF(i/=m_pos)THEN
			a_x_s(i)=a_x(m_local-1)[i,n_pos,o_rank]
			f_x_s(i)=f_x(m_local-1)[i,n_pos,o_rank]
			g_x_s(i)=g_x(m_local-1)[i,n_pos,o_rank]
			y_x_s(i,1:n_kop,1:2)=y_x(m_local,n_anfang:n_ende,1:2)[i,n_pos,o_rank]
		ELSE
			a_x_s(m_pos)=a_x(m_local-1)
			f_x_s(m_pos)=f_x(m_local-1)
			g_x_s(m_pos)=g_x(m_local-1)
			y_x_s(m_pos,1:n_kop,1:2)=y_x(m_local-1,n_anfang:n_ende,1:2)
		ENDIF
	ENDDO
	!!SYNC
	SYNC IMAGES(sync_m_array)

!	IF(n_pos==1.AND.m_pos==2)THEN
!		DO i=1,n_kop
!			PRINT *,y_x_s(1,i,1)
!		ENDDO
!	ENDIF

	!!
	!!Gauss
	DO i=1,R-1
		zi=f_x_s(i+1)/a_x_s(i)			
		a_x_s(i+1)=a_x_s(i+1)-zi*g_x_s(i)
		y_x_s(i+1,1:n_kop,1:2)=y_x_s(i+1,1:n_kop,1:2)-zi*y_x_s(i,1:n_kop,1:2)
	ENDDO
	x_x_s(R,1:n_kop,1:2)=y_x_s(R,1:n_kop,1:2)/a_x_s(R)
	DO i=R-1,1,-1
		x_x_s(i,1:n_kop,1:2)=(y_x_s(i,1:n_kop,1:2)-g_x_s(i)*x_x_s(i+1,1:n_kop,1:2))/a_x_s(i)
	ENDDO

!	IF(n_pos==1.AND.m_pos==1)THEN
!		DO h=1,n_kop
!			DO i=1,R			
!				PRINT *,x_x_s(i,h,1)
!			ENDDO
!		ENDDO
!	ENDIF

	!!SYNC
	SYNC IMAGES(sync_m_array)
	!!
	!! Transfer
	DO i=1,R
		IF(i/=m_pos)THEN
			IF(i/=1)x_x(0,n_anfang:n_ende,1:2)[i,n_pos,o_rank]=x_x_s(i-1,1:n_kop,1:2)
			x_x(m_local-1,n_anfang:n_ende,1:2)[i,n_pos,o_rank]=x_x_s(i,1:n_kop,1:2)
		ELSE
			IF(m_pos/=1)x_x(0,n_anfang:n_ende,1:2)=x_x_s(m_pos-1,1:n_kop,1:2)
			x_x(m_local-1,n_anfang:n_ende,1:2)=x_x_s(m_pos,1:n_kop,1:2)
		ENDIF
	ENDDO

	!!SYNC
	SYNC IMAGES(sync_m_array)
	!!
	
	!!Schritt 4
	IF(m_pos==1)THEN
		DO i=1,m_local-2
			x_x(i,1:n_local,1:2)=(y_x(i,1:n_local,1:2)-g_x(i)*x_x(m_local-1,1:n_local,1:2))/a_x(i)
		ENDDO
	ELSE
		DO i=1,m_local-2
			x_x(i,1:n_local,1:2)=&
			(y_x(i,1:n_local,1:2)-f_x(i)*x_x(0,1:n_local,1:2)-g_x(i)*x_x(m_local-1,1:n_local,1:2))/a_x(i)
		ENDDO
	ENDIF

	!IF(m_pos==1.AND.n_pos==1)THEN
	!	DO h=1,n_local
	!		PRINT *,'x:'
	!		DO i=1,m_local
	!			PRINT *,x_x(i,h,2)
	!		ENDDO
	!	ENDDO
	!ENDIF	

ENDIF

!	print *,'x_x(m_local,1:n_local,2)',x_x(m_local,1:n_local,2),'m_pos: ',m_pos,'n_pos: ',n_pos	

!	!!Sherman-Morrison-Formel
!	alpha_x(1:n_local,1:2)=0.
!
!	DO i=1,m_local
!		alpha_x(1:n_local,1:2)=alpha_x(1:n_local,1:2)+z_x(i)*x_x(i,1:n_local,1:2)
!	ENDDO
!!	print *,'alpha_x(1:n_local,1:2): ',alpha_x(1:n_local,1),'m_pos: ',m_pos,'n_pos: ',n_pos	
!	!!SYNC
!	IF(m_pos==1)THEN
!		SYNC IMAGES(sync_m_array)
!	ELSE
!		SYNC IMAGES(image_m_1)
!	ENDIF
!	!!
!	IF(m_pos==1)THEN
!		DO i=2,R
!!	print *,'KHODAAAAAAAAAAAAAAAAAAAAAAAAAAA'
!			alpha_x(1:n_local,1:2)=alpha_x(1:n_local,1:2)+alpha_x(1:n_local,1:2)[i,n_pos]
!		ENDDO
!		beta_x(1:n_local)=alpha_x(1:n_local,1)/(1.-alpha_x(1:n_local,2))
!		!IF(n_pos==2)THEN
!		!	DO h=1,n_local
!		!		PRINT *,beta_x(h)
!		!	ENDDO
!		!ENDIF
!	ENDIF
!	!!SYNC
!	IF(m_pos==1)THEN
!		SYNC IMAGES(sync_m_array)
!	ELSE
!		SYNC IMAGES(image_m_1)
!	ENDIF
!	!!
!	IF(m_pos>1)beta_x(1:n_local)=beta_x(1:n_local)[1,n_pos]	
!	!!SYNC
!	IF(m_pos==1)THEN
!		SYNC IMAGES(sync_m_array)
!	ELSE
!		SYNC IMAGES(image_m_1)
!	ENDIF
!	!!
!	DO i=1,n_local
!		x_x(1:m_local,i,1)=x_x(1:m_local,i,1)+beta_x(i)*x_x(1:m_local,i,2)
!		!IF(m_pos==1.AND.n_pos==1)THEN
!		!	PRINT *,'x: '
!		!	DO h=1,m_local
!		!		PRINT *,x_x(h,i,1)
!		!	ENDDO
!		!ENDIF
!	ENDDO	
	!!SYNC
	SYNC IMAGES(sync_m_array)

	RETURN

END SUBROUTINE


end module
