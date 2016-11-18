module energy_tot_mod
USE variablen

contains

subroutine CALC_ENERGY

! Berechnet Gesamt Energie des Systems

	IMPLICIT NONE

	INTEGER::h,i,j	

	total_energy=0.
	DO j=1,o_local
		DO i=1,n_local
			DO h=1,m_local
				total_energy=total_energy+E(h,i,j)
			ENDDO
		ENDDO
	ENDDO

	print *,'total Energie[',myrank,']: ',total_energy[myrank]
	
	IF(myrank==1)THEN
		SYNC IMAGES(*)
	ELSE
		SYNC IMAGES(1)
	ENDIF
	IF(myrank==1)THEN	
		DO i=2,numprocs
			total_energy=total_energy+total_energy[i]
		ENDDO
		PRINT *,'Gesamtenergie:' ,total_energy
		DO i=2,numprocs
			total_energy[i]=total_energy
		ENDDO		
	ENDIF
	IF(myrank==1)THEN
		SYNC IMAGES(*)
	ELSE
		SYNC IMAGES(1)
	ENDIF

	RETURN

end subroutine

subroutine CALC_ENERGY_KIN

    IMPLICIT NONE

    INTEGER::h,i,j


    E_k(1:m_local,1:n_local,1:o_local)=0.5*(u(1:m_local,1:n_local,1:o_local)**2&
+v(1:m_local,1:n_local,1:o_local)**2+w(1:m_local,1:n_local,1:o_local)**2)

    E_k_total=0.
	DO j=1,o_local
        	DO i=1,n_local
            		DO h=1,m_local
                		 E_k_total=E_k_total+E_k(h,i,j)
            		ENDDO
        	ENDDO
	ENDDO

    IF(myrank==1)THEN
        SYNC IMAGES(*)
    ELSE
        SYNC IMAGES(1)
    ENDIF
    IF(myrank==1)THEN
        DO i=2,numprocs
            E_k_total=E_k_total+E_k_total[i]
        ENDDO
        !PRINT *,'Gesamtenergie:' ,total_energy
        DO i=2,numprocs
            E_k_total[i]=E_k_total
        ENDDO
    ENDIF
    IF(myrank==1)THEN
        SYNC IMAGES(*)
    ELSE
    SYNC IMAGES(1)
    ENDIF

    RETURN

end subroutine

end module
