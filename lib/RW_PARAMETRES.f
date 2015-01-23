CC **********************************************************************
CC			DATA : 22/05/2003
CC	GUARDA LA TEMPERATURA A 4 FONDARIES: iz=0, 1, izmoho, izmoho+1

      SUBROUTINE READP (TFPARA, rosed, roc, roalfa, RHOAST, RHOH2O, 
     + 		    TSURF, TBOTT, ZASTH, TISOTER, THDIFF, CONDUC, HSURF,
     +		    HEXP, PHEAT_m, IRheology_type, Qarray, enarray, 
     +		    Aarray, NELROW,NELCOL, NELZ, strainrate,
     +		    nitermax, tallmax, alfa, tinicial, Dtany, npasmax, 
     +	            NPASINT, ITSR, visco_cnst, Iremoval, Time_removal,
     +		    Zremoval, TISO_rem, Zctall, iremesh, dvis_allow, 
     +		    Dif_K, Te_flexure, hydro_model, Kerosdif, rain, 
     + 		    Krain, windazimut, relative_humidity,
     +		    evaporation, erodability, erodability_sed, 
     +		    K_river_cap, l_fluv_sedim, CXrain, CYrain)

C    READS PARAMETERS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER*2 hydro_model
      REAL Te_flexure, Kerosdif, rain, Krain, 
     +     windazimut, relative_humidity, evaporation, erodability,
     +	   erodability_sed, K_river_cap, l_fluv_sedim, CXrain, CYrain
      DIMENSION CONDUC(0:3),Qarray(3),enarray(3),Aarray(3)
      CHARACTER*80 TFPARA
      LOGICAL RESTART

      OPEN(8,FILE='parametres.in')
      READ(8,*)
         WRITE(6,"('==================== PARAMETERS (DEFINE THE MODEL)',
     + 		   ' =======================')")
      READ(8,"(A70)") TFPARA
         WRITE(6,"(A70)") TFPARA

      READ(8,*) RHOH2O, rosed, roc, RHOAST
         WRITE(6,"(1X,4F7.1,3X,'Density of water, sediments, crust and',
     +		' asthenosphere (Kg/m**3)')") RHOH2O, rosed, roc, RHOAST

      READ(8,*) roalfa
         WRITE(6,"(1X,1P,E10.3,3X,'Volumetric thermal expansion of ',
     +		'lithospheric mantle (1/K)')") roalfa

      READ(8,*) TSURF, TBOTT, ZASTH
         WRITE(6,"(1X,2F8.1,3X,'Surface (z=0) and bottom (z=',F7.2,
     +		' km) Temperature (K)')") TSURF, TBOTT, ZASTH/1000

      READ(8,*) TISOTER
         WRITE(6,"(1X,F10.1,3X,'Isotherm base lithosphere (GLit_ter)',
     +				' (K)')") TISOTER

      READ(8,*) THDIFF
         WRITE(6,"(1X,1P,E10.3,3X,'Thermal diffusivity (m**2/s)')") 
     +				THDIFF
       
      READ(8,*) CONDUC(0), CONDUC(1),CONDUC(2), CONDUC(3)
         WRITE(6,"(1X,4F7.2,3X,'Conductivity of sediments, crust, ', 
     +		'lithospheric mantle and asthenosphere (W/m*K)' )")
     +		CONDUC(0),CONDUC(1),CONDUC(2), CONDUC(3)
       
      READ(8,*) HSURF,HEXP
       WRITE(6,"(1X,1P,2E10.3,3X,'Surface and Exponent Crustal Heat ',
     +			'Production')") HSURF, HEXP

      READ(8,*) PHEAT_m
         WRITE(6,"(1X,1P,E10.3,3x,'Lithospheric Mantle Heat ',
     +			'Production (constant)')") PHEAT_m

      READ(8,*) IRheology_type
         WRITE(6,"(I5,3X,'Rheological parameters ')") IRheology_type
	 IF(IRheology_type==99) THEN
		READ(8,*) Qarray(1),Qarray(2),Qarray(3)
		READ(8,*) enarray(1),enarray(2),enarray(3)
		READ(8,*) Aarray(1),Aarray(2),Aarray(3)
		WRITE(6,"(8X,3E10.2,6X,'Power flow activation energy, '
     +		     'Qarray [J/mol]'/8X,3F6.1,10X,'Power law exponent,'
     +		     ' enarray'/8X,3E11.3,3X,'Pre-exponencial constant,'
     +		     ' Aarray[MPa-n s-1]')") Qarray(1), Qarray(2),
     +			Qarray(3), enarray(1), enarray(2), enarray(3),
     +			Aarray(1), Aarray(2), Aarray(3)
	 ENDIF

      READ(8,*)
         WRITE(6,"(' ============================',
     +    	' HOW TO FIND THE SOLUTION ===================')")

      READ(8,*) NELCOL, NELROW, NELZ
         WRITE(6,"(1X,3I6,3X,'n, m, NELZ')") NELCOL, NELROW, NELZ

      READ(8,*) strainrate
         WRITE(6,"(1X,1P,E10.3,3X,
     +		'reference strain rate (s-1)')") strainrate

      READ(8,*) nitermax
         WRITE(6,"(1X,I10,3X,'Maximum iterations within velocity',
     +			' solution in each timestep')") nitermax

       READ(8,*) tallmax
       WRITE(6,"(1X,F11.4,3X,'Acceptable rms fractional error',
     +		' (stops iteration)')") tallmax
     
       READ(8,*) alfa
       WRITE(6,"(1X,1P,E10.2,' Relacio entre la velocitat abans i',
     +		' despres de iterar (alfa)')") alfa
     
      READ(8,*) tinicial
         WRITE(6,53) tinicial
  53   FORMAT(' ',1P,E10.2,' BEGINNING OF CALCULATION')
      READ(8,*) Dtany
         WRITE(6,55) Dtany
  55   FORMAT(' ',1P,E10.2,' SIZE OF TIME STEPS (anys)')
      READ(8,*) npasmax
         WRITE(6,57) npasmax
  57   FORMAT(' ',I10,' NUMBER OF TIME STEPS')
      READ(8,*) NPASINT
         WRITE(6,58) NPASINT
  58   FORMAT(' ',I10,' NUMERO DE PASSOS INTERMITJOS QUE GUARDARE')
      READ(8,*)
         WRITE(6,68)
  68   FORMAT(' ================== PARTICULAR CONDITIONS',
     +      ' ========')
      READ(8,*) ITSR
	 IF(ITSR<=0) THEN
	      BACKSPACE (UNIT=8)
              READ(8,*) ITSR, visco_cnst
	      IF(ITSR==0) THEN
		   WRITE(6,"(1X,I10,' Constant Viscosity:',
     +  		     1P,G12.3,' Pa.s')") ITSR, visco_cnst
     	      ELSE
		   WRITE(6,"(1X,I10,' Constant lithosphere strength:',
     +  		     1P,G12.3,' N/m')") ITSR, visco_cnst
	      ENDIF
	  ELSE
	      WRITE(6,"(1X,I10,
     +  	       ' Viscosity strain rate dependent',G10.3)") ITSR
	 ENDIF
      READ(8,*) Iremoval
	 IF(Iremoval>0) THEN
	    BACKSPACE (UNIT=8)
	    READ(8,*) Iremoval, Time_removal, Rem_tmp, Zctall
	      WRITE(6,"(3X,I3,' Convective removal at',F6.1,' My')")
     +				Iremoval, Time_removal		
            IF((Iremoval>1.AND.Iremoval<6) .OR. (Iremoval==22)) THEN
	      Zremoval=Rem_tmp
	      WRITE(6,"(6X,'To lithosphere depth:',F6.1,' km')")
     +			Zremoval/1.D3	
	      IF(Iremoval==2.OR.Iremoval==3.OR.Iremoval==22) 
     +			WRITE(6,"(8X,'if lithosphere >',F6.1,' km ')")
     +				Zctall/1.D3
	      IF(Iremoval==4) WRITE(6,"(8X,'if crustal thickness >',
     +			F6.1,' km ')") Zctall/1.D3
	    ELSE
	      TISO_rem=Rem_tmp
              WRITE(6,"(6X,'To the isotherm T =',F8.2,' K')") TISO_rem
	      IF(Iremoval==1) THEN
	         WRITE(6,"(8X,'if lithosphere thickness >',F6.1,' km')")
     +				 Zctall/1.D3
	      ELSE	      
	         WRITE(6,"(8X,'if crustal thickness >',F6.1,' km')")
     +				Zctall/1.D3
	      ENDIF	      
	    ENDIF
	 ENDIF
      READ(8,*) iremesh
         IF(iremesh==1) 
     +		WRITE(6,"(' iremesh =',I3,' => Re-mesh')") iremesh
      READ(8,*) dvis_allow, Dif_K
      WRITE(6,"(1X,1P,2E11.3,' maximum gradient permited on variables',
     +  	  ' and constant diffusive filter')") dvis_allow, Dif_K

      READ(8,*)
         WRITE(6,"(' ============================ ELASTIC THICKNESS',
     +    	' AND SURFACE PROCESSES ===================')")
      READ(8,*) Te_flexure
         WRITE(6,"(1X,F12.1,3X,'elastic thickness [m]')") Te_flexure
      READ(8,*) hydro_model
         WRITE(6,"(' hydro_model=',I3)") hydro_model
      READ(8,*) Kerosdif
         WRITE(6,"(1X,F10.1,' DIFFUSIVE TRANSPORT EROSION COEFF.',
     +		 ' [m2/a]')") Kerosdif
      READ(8,*) rain
         WRITE(6,144) rain
 144   FORMAT(' ',F10.1,' BACKGROUND RUNOFF [l/m2/a]=[mm/a]')
      READ(8,*) Krain
         WRITE(6,146) Krain
 146   FORMAT(' ',F10.1,
     + 		' PROPORTIONALITY OF RUNOFF WITH ALTITUDE [l/m2/a/km]')
      READ(8,*) windazimut, relative_humidity
         WRITE(6,"(1X,2F12.3,3X,' windazimut AND relative_humidity')")
     +		 windazimut, relative_humidity    
      READ(8,*) evaporation
         WRITE(6,"(1X,F12.1,3X,'evaporation')") evaporation       
      READ(8,*) erodability, erodability_sed 
         WRITE(6,"(1X,2F12.1,3X,'erodability AND erodability_sed')")
     +		  erodability, erodability_sed       
      READ(8,*) K_river_cap, l_fluv_sedim
         WRITE(6,"(1X,2F12.1,3X,'K_river_cap AND l_fluv_sedim')")
     +		  K_river_cap, l_fluv_sedim
      READ(8,*) CXrain, CYrain
         WRITE(6,"(1X,2F12.1,3X,'CXrain AND CYrain'/)")
     +		  CXrain, CYrain

       CLOSE(8)
       
       RETURN
       END SUBROUTINE READP
       
C **********************************************************************
C     SUBRUTINA PER LLEGIR ELS RESULTATS D'UN FITXER
C
      SUBROUTINE RERESULT (fitxer, nporta, TITLE, Tsegons, Tanys, m, n,
     +            NELZ, nn, nnsd, AX, BY, Dz, TLITOS, scrust, GLit_ter,
     +            GLit, visco, u, v, epuntzz, Tmoho, Qsurface,
     +		  elevation, sediment, kpuntssd, D_Bodies, curve)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*170 TITLE
      CHARACTER*15 fitxer,Achar,Achar2
      DIMENSION elevation(nnsd),sediment(nnsd),scrust(nnsd),GLit(nnsd),
     +          GLit_ter(nnsd),visco(nnsd),u(nnsd),v(nnsd),
     +          epuntzz(nnsd),Qsurface(nnsd),Tmoho(nnsd),
     +		D_Bodies(kpuntssd,2),curve(0:kpuntssd,2)

      PRINT*,'   LLEGEIXO ELS RESULTATS DEL FITXER :',fitxer
      OPEN(nporta,FILE=fitxer)

 1004  FORMAT(' npasos:',I3,'  TIME = ',1P,E10.4,' (',0P,F9.3,' Manys)')
 1003  FORMAT(' TEMPERATURA BASE LITOSFERA (TLITOS):',F23.13,' Kelvins')
 1005  FORMAT(1P,4E19.12)
 1006  FORMAT(1P,6E13.6)

       READ(nporta,"(A170)") TITLE
       WRITE(6,*) TITLE

!      READ(nporta,"(' npasos:',I3,'  TIME = ',1P,E10.4,' (',0P,F9.3,
!    +             ' Manys)')") npasos,Tsegons,TEMPSMa
       READ(nporta,*) npasos,Tsegons,TEMPSMa
       WRITE(6,*) npasos
       WRITE(6,*) Tsegons
       WRITE(6,*) TEMPSMa
       Tanys=TEMPSMa*1.D6
       
!      READ(nporta,1003) TLITOS
       READ(nporta,*) TLITOS
       WRITE(6,*) TLITOS
       
!      READ(nporta,"('n =',I4,' m =',I4,4X,'h=',F16.8,' m   tau =',
!    +       F16.8,' m   Dz =',F18.13,' m')") n,m,h,tau,Dz
        READ(nporta,*) n,m,h,tau,Dz
        WRITE(6,*) n,h
        WRITE(6,*) m,tau
        WRITE(6,*) Dz

       nn=(n+1)*(m+1)
       AX=h*n
       BY=tau*m 
       READ(nporta,*)
       READ(nporta,1005) (elevation(kxy), kxy=1,nn)
       READ(nporta,*)
       READ(nporta,1005) (sediment(kxy), kxy=1,nn)
       READ(nporta,*)
       READ(nporta,1005) (scrust(kxy), kxy=1,nn)
       READ(nporta,*)
       READ(nporta,1005) (GLit_ter(kxy), kxy=1,nn)
       READ(nporta,*)
       READ(nporta,1005) (GLit(kxy), kxy=1,nn)
       READ(nporta,*)
       READ(nporta,1005) (visco(kxy), kxy=1,nn)
       READ(nporta,*)
       READ(nporta,1005) (u(kxy), kxy=1,nn)
       READ(nporta,*)
       READ(nporta,1005) (v(kxy), kxy=1,nn)
       READ(nporta,*)
       READ(nporta,1005) (epuntzz(kxy), kxy=1,nn)
       READ(nporta,*)
       READ(nporta,1005) (Tmoho(kxy), kxy=1,nn)
       READ(nporta,*)
       READ(nporta,1005) (Qsurface(kxy), kxy=1,nn)

      READ(nporta,*,END=11) Achar, kpunts1, Achar2, qFlit1
      WRITE(6,*) Achar, kpunts1, Achar2, qFlit1
      D_Bodies(1,1)=kpunts1
      D_Bodies(1,2)=qFlit1
      IF(kpunts1/=0) THEN
	      READ(nporta,*) 
	      READ(nporta,1006) ((D_Bodies(kp,i),i=1,2),kp=2,kpunts1+1)
      END IF
      
      READ(nporta,*,END=11) Achar, kpunts2, Achar2, qFlit2
      WRITE(6,*) Achar, kpunts2, Achar2, qFlit2
      D_Bodies(kpunts1+2,1)=kpunts2
      D_Bodies(kpunts1+2,2)=qFlit2
      IF(kpunts2/=0) THEN
		READ(nporta,*) 
		READ(nporta,1006) ((D_Bodies(kp,i),i=1,2),
     +       			kp=kpunts1+3,kpunts1+kpunts2+2)
      END IF

      READ(nporta,*,END=11) Achar, kpunts3, Achar2, qFlit3
      WRITE(6,*) Achar, kpunts3, Achar2, qFlit3
      D_Bodies(kpunts1+kpunts2+3,1)=kpunts3
      D_Bodies(kpunts1+kpunts2+3,2)=qFlit3
      IF(kpunts3/=0) THEN
		READ(nporta,*) 
		READ(nporta,1006) ((D_Bodies(kp,i),i=1,2),
     + 			kp=kpunts1+kpunts2+4,kpunts1+kpunts2+kpunts3+3)
      END IF

      READ(nporta,*,END=11) Achar, kpunts4, Achar2, qFlit4
      WRITE(6,*) Achar, kpunts4, Achar2, qFlit4
      D_Bodies(kpunts1+kpunts2+kpunts3+4,1)=kpunts4
      D_Bodies(kpunts1+kpunts2+kpunts3+4,2)=qFlit4
      IF(kpunts4/=0) THEN
		READ(nporta,*) 
		READ(nporta,1006) ((D_Bodies(kp,i),i=1,2),
     + 			kp=kpunts1+kpunts2+kpunts3+5,
     +			    kpunts1+kpunts2+kpunts3+kpunts4+4)
      END IF

      READ(nporta,*,END=11) Achar, kpunts5, Achar2, qFlit5
      WRITE(6,*) Achar, kpunts5, Achar2, qFlit5
      D_Bodies(kpunts1+kpunts2+kpunts3+kpunts4+5,1)=kpunts5
      D_Bodies(kpunts1+kpunts2+kpunts3+kpunts4+5,2)=qFlit5
      IF(kpunts5/=0) THEN
		READ(nporta,*) 
		READ(nporta,1006) ((D_Bodies(kp,i),i=1,2),
     + 			kp=kpunts1+kpunts2+kpunts3+kpunts4+6,
     +			   kpunts1+kpunts2+kpunts3+kpunts4+kpunts5+6)
      END IF

      READ(nporta,*,END=11) Achar, kpcurve
      curve(0,1)=kpcurve
      READ(nporta,*) 
      READ(nporta,1006) ((curve(kp,i),i=1,2),kp=1,kpcurve)

 11   CONTINUE  
      PRINT*,'   return reresults'   
      RETURN
      END SUBROUTINE RERESULT

C*****************************************************************
C     SUBRUTINA PER GUARDAR ELS RESULTATS EN UN FITXER
C   TOTES LES VARIABLES ENTREN AMB SISTEMA INTERNACIONAL
C             LA VELOCITAT LA GUARDO AMB mm/any

      SUBROUTINE WRRESULT (fitxer, TITLE, TEMPSMa, npasos,
     +                     m, n, NELZ, nn, FACTEMP, FACVEL, AX, BY, 
     +                     Dz, TLITOS, scrust, GLit_ter, GLit, visco,
     +                     u, v, epuntzz, Tmoho, Qsurface, elevation,
     +			   sediment, kpuntssd, D_Bodies, curve)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*170 TITLE
      CHARACTER*15 fitxer
      DIMENSION elevation(nn),sediment(nn),scrust(nn),GLit_ter(nn),
     +         GLit(nn),visco(nn),u(nn),v(nn),epuntzz(nn),Qsurface(nn),
     +         Tmoho(nn),D_Bodies(kpuntssd,2),curve(0:kpuntssd,2)

      OPEN(1,FILE=fitxer)

      h=AX/n
      tau=BY/m
      nn=(n+1)*(m+1)
      Tsegons=TEMPSMa*1.D6*FACTEMP

 1004  FORMAT(' npasos:',I3,'  TIME = ',1P,E10.4,' (',0P,F9.3,' Manys)')
 1003  FORMAT(' TEMPERATURA BASE LITOSFERA (TLITOS):',F23.13,' Kelvins')
 1005  FORMAT(1P,4E19.12)
 1006  FORMAT(1P,6E13.6)

       WRITE(1,"(A170)") TITLE
       
!      WRITE(1,1004) npasos,Tsegons,TEMPSMa
       WRITE(1,*) npasos,Tsegons,TEMPSMa
       
!      WRITE(1,1003) TLITOS
       WRITE(1,*) TLITOS
       
!      WRITE(1,"('n =',I4,' m =',I4,4X,'h=',F16.8,' m   tau =',
!    +       F16.8,' m   Dz =',F18.13,' m')") n,m,h,tau,Dz
       WRITE(1,*) n,m,h,tau,Dz
       
      WRITE(1,*) '  elevation, (m)  ( elevation(kxy), kxy=1,nn )'
      WRITE(1,1005) (elevation(kxy), kxy=1,nn)
      WRITE(1,*) '  sediments, (m)  ( sediment(kxy), kxy=1,nn )'
      WRITE(1,1005) (sediment(kxy), kxy=1,nn)
      WRITE(1,*) '   Gruix cortical, (m)  ( scrust(kxy), kxy=1,nn )'
      WRITE(1,1005) (scrust(kxy), kxy=1,nn)
      WRITE(1,*) '  Thermal Lithospheric thickness, (m)  ',
     + 			'( GLit_ter(kxy), kxy=1,nn )'
      WRITE(1,1005) (GLit_ter(kxy), kxy=1,nn)
      WRITE(1,*) '  Depth of conductivity channge (lith-asth) , (m)  ',
     + 			'( GLit(kxy), kxy=1,nn )'
      WRITE(1,1005) (GLit(kxy), kxy=1,nn)
      WRITE(1,*) '    viscositat  ( visco(kxy), kxy=1,nn )'
      WRITE(1,1005) (visco(kxy), kxy=1,nn)
      WRITE(1,*) '   velocitat x, (mm/any)  ( u(kxy), kxy=1,nn )'
      WRITE(1,1005) (u(kxy)*FACVEL, kxy=1,nn)
      WRITE(1,*) '   velocitat y, (mm/any)  ( v(kxy), kxy=1,nn )'
      WRITE(1,1005) (v(kxy)*FACVEL, kxy=1,nn)
      WRITE(1,*) '    vertical strain rate  ( epuntzz(kxy), kxy=1,nn )'
      WRITE(1,1005) (epuntzz(kxy), kxy=1,nn)
      WRITE(1,*) '   moho Temperature (K)  ( Tmoho(kxy), kxy=1,nn )'
      WRITE(1,1005) (Tmoho(kxy), kxy=1,nn)
      WRITE(1,*) '   surface heat flow, (W/m2) ',
     + 			'( Qsurface(kxy), kxy=1,nn )'
      WRITE(1,1005) (Qsurface(kxy), kxy=1,nn)

      kpunts1=D_Bodies(1,1)			!! BODY 1
      qFlit1=D_Bodies(1,2)
      IF(kpunts1/=0) THEN
         WRITE(1,"(' kpunts',I4,'  qFlit',F11.2)") kpunts1, qFlit1
         WRITE(1,*) 'BODY 1:  (x,y)(kp) kp=1,kpunts'
         WRITE(1,1006) ((D_Bodies(kp,i),i=1,2),kp=2,kpunts1+1)
      ENDIF
      kpunts2=D_Bodies(kpunts1+2,1)		!! BODY 2
      qFlit2=D_Bodies(kpunts1+2,2)
      WRITE(1,"(' kpunts',I4,'  qFlit',F11.2)") kpunts2, qFlit2
      IF(kpunts2/=0) THEN
         WRITE(1,*) 'BODY 2:  (x,y)(kp) kp=1,kpunts'
         WRITE(1,1006) 
     +		((D_Bodies(kp,i),i=1,2),kp=kpunts1+3,kpunts1+kpunts2+2)
      ENDIF
      kpunts3=D_Bodies(kpunts1+kpunts2+3,1)	!! BODY 3
      qFlit3=D_Bodies(kpunts1+kpunts2+3,2)
      k1=kpunts1+kpunts2+4
      k2=kpunts1+kpunts2+kpunts3+3
      WRITE(1,"(' kpunts',I4,'  qFlit',F11.2)") kpunts3, qFlit3
      IF(kpunts3/=0) THEN
         WRITE(1,*) 'BODY 3:  (x,y)(kp) kp=1,kpunts'
         WRITE(1,1006) 
     +		((D_Bodies(kp,i),i=1,2),kp=k1,k2)
      ENDIF
      kpunts4=D_Bodies(kpunts1+kpunts2+kpunts3+4,1)	!! BODY 4
      qFlit4=D_Bodies(kpunts1+kpunts2+kpunts3+4,2)
      k1=kpunts1+kpunts2+kpunts3+5
      k2=kpunts1+kpunts2+kpunts3+kpunts4+4
      WRITE(1,"(' kpunts',I4,'  qFlit',F11.2)") kpunts4, qFlit4
      IF(kpunts4/=0) THEN
         WRITE(1,*) 'BODY 4:  (x,y)(kp) kp=1,kpunts'
         WRITE(1,1006) 
     +		((D_Bodies(kp,i),i=1,2),kp=k1,k2)
      ENDIF
      kpunts5=D_Bodies(kpunts1+kpunts2+kpunts3+kpunts4+5,1)	!! BODY 5
      qFlit5=D_Bodies(kpunts1+kpunts2+kpunts3+kpunts4+5,2)
      k1=kpunts1+kpunts2+kpunts3+kpunts4+6
      k2=kpunts1+kpunts2+kpunts3+kpunts4+kpunts5+5
      WRITE(1,"(' kpunts',I4,'  qFlit',F11.2)") kpunts5, qFlit5
      IF(kpunts5/=0) THEN
         WRITE(1,*) 'BODY 5:  (x,y)(kp) kp=1,kpunts'
         WRITE(1,1006) 
     +		((D_Bodies(kp,i),i=1,2),kp=k1,k2)
      ENDIF

      kpcurve=curve(0,1)		!! curve
      IF(kpcurve/=0) THEN
         WRITE(1,"(' kpcurve',I5)") kpcurve
         WRITE(1,*) 'curve:  (x,y)(kp) kp=1,kpcurve'
         WRITE(1,1006) ((curve(kp,i),i=1,2),kp=1,kpcurve)
      ENDIF

      RETURN
      END SUBROUTINE WRRESULT
       
!! ********************************************************************
      SUBROUTINE vertical_strain_rate_2D (AX, BY, n, m,nn,u,v,epuntzz)

CC    Calculation of the vertical strain rate and control that 
CC		it should be lower than ep_limit.

CC  epuntzz [1/s]  Array returning the vertical strain rate.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ep_limit=5.D-15, zero=0.D-100)
      DIMENSION u(nn), v(nn), epuntzz(nn)

      Dx=AX/n
      Dy=BY/m
      epuntzz_min=1.D30
      epuntzz_max=-1.D30
     
       DO iy=1,m-1
         DO ix=1,n-1
	    	kxy=ix+1+iy*(n+1)
		ux=(u(kxy+1)-u(kxy-1))/(2.D0*Dx)
		vy=(v(kxy+n+1)-v(kxy-n-1))/(2.D0*Dy)
		epuntzz(kxy)=-(ux+vy)
          	epuntzz_max=MAX(epuntzz(kxy),epuntzz_max)
              	epuntzz_min=MIN(epuntzz(kxy),epuntzz_min)
	 END DO
       END DO

!!    VERTICAL STRAIN RATE TO THE BOUNDARIES
             DO ix=1,n-1
                 kxys=ix+1
                 kxyn=ix+1+m*(n+1)
		 IF(ABS(u(kxys))<zero .AND. ABS(v(kxys))<zero) THEN
		 	epuntzz(kxys)=0.D0
		  ELSE
			epuntzz(kxys)=epuntzz(kxys+n+1)
		 ENDIF
		 IF(ABS(u(kxyn))<zero .AND. ABS(v(kxyn))<zero) THEN
		 	epuntzz(kxyn)=0.D0
		  ELSE			
			epuntzz(kxyn)=epuntzz(kxyn-n-1)
		 ENDIF
	     END DO
             DO iy=0,m
                 kxyw=1+iy*(n+1)
                 kxye=n+1+iy*(n+1)
		 IF(ABS(u(kxyw))<zero .AND. ABS(v(kxyw))<zero) THEN		 
			epuntzz(kxyw)=0.D0
		 ELSE
			epuntzz(kxyw)=epuntzz(kxyw+1)
		 ENDIF
		 IF(ABS(u(kxye))<zero .AND. ABS(v(kxye))<zero) THEN		 
			epuntzz(kxye)=0.D0
		  ELSE
                 	epuntzz(kxye)=epuntzz(kxye-1)
		 ENDIF
	     END DO
CC	LIMIT TO THE VERTICAL STRAIN RATE
          DO kxy=1,nn
C  		limit the calculated vert.str.rate 'epuntzz'.
    		IF(epuntzz(kxy).GT.ep_limit) epuntzz(kxy)=ep_limit
            	IF(epuntzz(kxy).LT.(-1.D0*ep_limit)) 
     +			epuntzz(kxy)=-1.D0*ep_limit
          	epuntzz_max=MAX(epuntzz(kxy),epuntzz_max)
              	epuntzz_min=MIN(epuntzz(kxy),epuntzz_min)
	  END DO
 
      WRITE(6,61) epuntzz_min,epuntzz_max,ep_limit
 61   FORMAT(4X,'Vertical strain rate: minim:',1P,G12.4,
     +	  ' s-1,    maxim:',1P,G12.4,' s-1,  limit:',1P,G9.2,' s-1')

      RETURN
      END SUBROUTINE vertical_strain_rate_2D

!! ********************************************************************
      SUBROUTINE vertical_strain_rate (ix_main,iy_main, m,n,nn,Dx,Dy,
     +					   u, v, ep_limit, epuntzz_P)
      
!  Calculate the vertical strain rate at point (ix_main,iy_main) 
!	from the velocity field (u,v).
!	This vertical strain rate should be lower than ep_limit
!  Input: m, n, nn, Dx, Dy, u, v.	Output: epuntzz_P [1/s]
!  Dx, Dy : grid interval on x and y.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (zero=0.D-100)
      DIMENSION u(nn), v(nn)

      ddx=2.D0*Dx
      ddy=2.D0*Dy
      ix=ix_main
      iy=iy_main
      kxy=ix+1+iy*(n+1)
      !IF(ix==0. OR. ix==n .OR. iy==0 .OR. iy==m) THEN			! If the velocity is nul on the borders => 
	!   IF(ABS(u(kxy))<zero .AND. ABS(v(kxy))<zero) THEN		!	 => Fix the vertical strain rate to zero
	!	epuntzz_P=0.D0
	!   	RETURN
	!   END IF
      !END IF
      
      IF(ix==0) ix=1    
      IF(ix==n) ix=n-1
      IF(iy==0) iy=1
      IF(iy==m) iy=m-1
        kxy=ix+1+iy*(n+1)
        !ux=(u(kxy+1)-u(kxy-1))/ddx  		!! Derivada centrada
        !vy=(v(kxy+n+1)-v(kxy-n-1))/ddy  	!! Derivada centrada

	sign_x=SIGN(1.D0,u(kxy))						!! SIGN(A,B)=Valor de A amb signe de B 
	sign_y=SIGN(1.D0,v(kxy))
        ux_sup=(u(kxy+1)-u(kxy))/Dx		!! Derivada segons el fluxe del fluid
        ux_inf=(u(kxy)-u(kxy-1))/Dx	
	ux=ux_inf*MAX(sign_x,0D0)-ux_sup*MIN(sign_x,0D0)
        vy_sup=(v(kxy+n+1)-v(kxy))/Dy 
        vy_inf=(v(kxy)-v(kxy-n-1))/Dy	
	vy=vy_inf*MAX(sign_y,0D0)-vy_sup*MIN(sign_y,0D0)

        epuntzz_P=-(ux+vy)

        epuntzz_P=MIN(epuntzz_P,ep_limit)
        epuntzz_P=MAX(epuntzz_P,(-1.D0*ep_limit))

      RETURN
      END SUBROUTINE vertical_strain_rate

!! ********************************************************************
      Double Precision FUNCTION effective_strainrate (ix_main, iy_main,
     +					m, n, nn, Dx, Dy, u, v)
      
!  Calculate the Effective strain rate at point (ix,iy) 
!	from the velocity field (u,v)
!  Input: m, n, nn, Dx, Dy, u, v.	Output: epeffec
!  Dx, Dy : grid interval on x and y.
!  Effective strain rate = sqrt[1/2(exx**2+eyy**2+ezz**2)+exy**2]
!  Second Invariant E=sqrt[2*(exx**2+eyy**2+exy**2+exx*eyy)]

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (epeffec_inf=1.D-18, epeffec_sup=1.D-15)
      DIMENSION u(nn), v(nn)

      ddx=2.D0*Dx
      ddy=2.D0*Dy
      ix=ix_main
      iy=iy_main
      IF(ix==0) ix=1
      IF(ix==n) ix=n-1
      IF(iy==0) iy=1
      IF(iy==m) iy=m-1

        kxy=ix+1+iy*(n+1)
        !ux=(u(kxy+1)-u(kxy-1))/ddx	  		!! Derivada centrada
        !vx=(v(kxy+1)-v(kxy-1))/ddx  			!! Derivada centrada
        !vy=(v(kxy+n+1)-v(kxy-n-1))/ddy  		!! Derivada centrada
        !uy=(u(kxy+n+1)-u(kxy-n-1))/ddy  		!! Derivada centrada

	sign_x=SIGN(1.D0,u(kxy))			!! Derivada segons el fluxe del fluid		
	sign_y=SIGN(1.D0,v(kxy))			!! Derivada segons el fluxe del fluid
        ux_sup=(u(kxy+1)-u(kxy))/Dx			!! Derivada segons el fluxe del fluid
        ux_inf=(u(kxy)-u(kxy-1))/Dx	
	ux=ux_inf*MAX(sign_x,0D0)-ux_sup*MIN(sign_x,0D0)
        vy_sup=(v(kxy+n+1)-v(kxy))/Dy 
        vy_inf=(v(kxy)-v(kxy-n-1))/Dy	
	vy=vy_inf*MAX(sign_y,0D0)-vy_sup*MIN(sign_y,0D0)
	uy_sup=(u(kxy+n+1)-u(kxy))/Dy 
        uy_inf=(u(kxy)-u(kxy-n-1))/Dy 
	uy=uy_inf*MAX(sign_y,0D0)-uy_sup*MIN(sign_y,0D0)
        vx_sup=(v(kxy+1)-v(kxy))/Dx 
        vx_inf=(v(kxy)-v(kxy-1))/Dx
	vx=vx_inf*MAX(sign_x,0D0)-vx_sup*MIN(sign_x,0D0)	!! Derivada segons el fluxe del fluid	
	
        epuntzz=-(ux+vy)
        ep=((ux*ux+vy*vy+(epuntzz*epuntzz))/2.D0)+
     +                            ((1.D0/4.D0)*(uy+vx)*(uy+vx))
        epeffec=DSQRT(ep)
        epeffec=MAX(epeffec,epeffec_inf)
        epeffec=MIN(epeffec,epeffec_sup)
	
	effective_strainrate=epeffec

      END FUNCTION 

!! ********************************************************************
      Double Precision FUNCTION P_average (elev, sed, crust, hlitos,
     +					     z_compens, RHOH2O, rosed,
     +					     roc, rom, RHOAST, g)
!      SUBROUTINE P_averagelib (elev, sed, crust, hlitos,
!     +				z_compens, RHOH2O, rosed,
!     +				roc, rom, RHOAST, g, P_average)
      
!  Output: Average Pressure: Depth-averaged vertical stress over the plate
!	Depth integral of the vertical stress/thin sheet thickness
!		integral[from surface to thin sheet thickness][vertical stress *dz]

!  Input: 
!    elev:	Elevation
!    sed:	Sediment thickness
!    crust:	Crustal thickness
!    hlitos:	Lithospheric thickness (sediments+crust+lithospheric mantle)
!    z_compens:	Depth of compensation.
!    layers densities: RHOH2O: water, rosed: sediments, roc: crust, 
!		       rom: lithospheric mantle, RHOAST: asthenosphere
!    g [m/s2]:	Gravity acceleration at the surface of the planet

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

   !   thick_layer=z_compens
   !   IF(elev.GT.0.D0) thick_layer=z_compens+elev
      thick_layer=z_compens+elev
      hm=hlitos-crust-sed
      ha=z_compens+elev-hlitos		!     ha=thick_layer-hlitos
      water=ABS(elev)
      IF(elev.GT.0.D0) water=0.D0
      
      Fw=(g*RHOH2O*water*water)/2.D0
      Fsed=(g*RHOH2O*water*sed)+((g*rosed*sed*sed)/2.D0)
      Fc=(g*(RHOH2O*water+rosed*sed)*crust)+
     +			((g*roc*crust*crust)/2.D0)
      Fm=(g*(RHOH2O*water+rosed*sed+roc*crust)*hm)+
     +			(g*rom*hm*hm/2.D0)
      Fa=(g*(RHOH2O*water+rosed*sed+roc*crust+rom*hm)*ha)+
     +			(g*RHOAST*ha*ha/2.D0)
      Pressure_integr=Fw+Fsed+Fc+Fm+Fa
      P_average=Pressure_integr/thick_layer
   !    WRITE(6,"(5X,5G20.11)")  Fw,Fsed,Fc,Fm,Fa

   !   END SUBROUTINE P_averagelib
      END FUNCTION 
