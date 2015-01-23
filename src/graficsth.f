	PROGRAM graficsth
C __________ Data: 2013.    Ultima modifica: 17/12/2003 _______________
C !!!!!!! OJO Les dimensions, que siguin les mateixes que a uhuru !!!!!
CC   LIBRARY NEEDED TO THE COMPILATION :
!		lib/RW_PARAMETRES.f
!		/usr/users/ivone/lib/lib_uhuru/outin.f
C  FITXERS QUE LLEGEIX : fresult_fin, fresult_ini.
C  ELS GRAFICS SON DEL FITXER fresult_fin 
C  AQUESTS RESULTATS ELS COMPARO AMB ELS DEL FITXER
C  fresult_ini , PER TROBAR EL GRAU D'ENGRUIXIMENT (beta)
C  VOLTOT : volum total nomes tenint en compte els punts interiors.


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL Te_flexure, Kerosdif, rain, Krain, windazimut, 
     +	   relative_humidity, evaporation, erodability,
     +	   erodability_sed, K_river_cap, l_fluv_sedim, CXrain, CYrain
      INTEGER*2 hydro_model
      CHARACTER*170 TITLE,TITLERE,TITPB1,TITPB2,TITPB3
      CHARACTER*15 fresult_ini,fresult_fin,ENDAVANT
      PARAMETER (nnsd=50000,nzsd=1000,kpuntssd=500,
     +        	PI=3.1415926535897932D0,FACVEL=3.1536D10,VLLEG=10.0,
     +		ep_limit=5.D-15,g=9.8D0) 
      DIMENSION scrust(nnsd),GLit_ter(nnsd),GLit(nnsd),u(nnsd),v(nnsd),
     +	       elevation(nnsd),sediment(nnsd),Qsurface(nnsd),vis(nnsd),
     +         epuntzz(nnsd),Tmoho(nnsd),D_Bodies(kpuntssd,2),
     +	       CONDUC(0:3),x_pol(50),y_pol(50),curve(0:kpuntssd,2)
      DATA dedtmax /-1.D20/ dedtmin /1.D20/
     
      WRITE(6,"(/'------------ PROGRAMA graficsth -----------------')")

      CALL READP (TITLE, rosed,roc, roalfa, RHOAST, RHOH2O, TSURF,TBOTT,
     +            ZASTH, TISOTER, THDIFF, CONDUC, HSURF, HEXP, PHEAT_m,
     +            IRheology_type, Qarray, enarray, Aarray, 
     +            m, n, NELZ, vrmma, nitermax, tallmax, alfa, tinicial,
     +            Dtany, npasmax, NPASINT, ITSR, visco_cnst, 
     +            Iremoval, Time_removal, Zremoval, TISO_rem, Zctall, 
     +            iremesh, dvis_allow, Dif_K, Te_flexure, hydro_model, 
     +		  Kerosdif, rain, Krain, windazimut, relative_humidity,
     +		  evaporation, erodability, erodability_sed, 
     +		  K_river_cap, l_fluv_sedim, CXrain, CYrain)


!      nn=(n+1)*(m+1)
!      fresult_ini='res1'
C      fresult_fin='res2'
C        fresult_fin='resultats5'

      nporta2=18
      WRITE(6,"(2X,'ep_limit should be the same in UHURU Program:',
     +		1P,G9.2,' s-1'//
     +		2X,'WRITE THE RESULTS FILE')") ep_limit
!      fresult_fin='resultats8'
      READ(5,"(A15)") fresult_fin
      WRITE(6,"(/3X,'PROFUNDITAT DE LA ISOTERMA:',F7.1,' K'/
     +      3X,'ELS GRAFICS SON DELS RESULTATS DE FITXER : ',A15,/)")
     +      TISOTER, fresult_fin
    
      CALL RERESULT (fresult_fin, nporta2, TITLERE, Tsegons, Tanys,
     +              m, n, NELZ, nn, nnsd, AX, BY, Dz, TLITOS, scrust,
     +              GLit_ter, GLit, vis, u, v, epuntzz, Tmoho, Qsurface,
     +		    elevation, sediment, kpuntssd, D_Bodies, curve)

      PRINT*,'   after call results'   
C  PASSO LA VELOCITAT A m/s (INICIALMENT A mm/any)
      u=u/FACVEL
      v=v/FACVEL

      ZASTH=NELZ*Dz
      Dx=AX/n
      Dy=BY/m
      WRITE(6,62) AX/1D3,BY/1D3,ZASTH/1D3,Dx/1D3,Dy/1D3,Dz/1D3,n,m
 62   FORMAT(1X,'AX=',F9.3,' km     BY: ',F9.3,' km    ZASTH:',
     +     F9.3,'km ',/,1X,'Dx:',F9.3,' km     Dy: ',F9.3,' km    Dz:',
     +     F9.3,' km',/,1X,'n:',I3,'     m:',I3)
      
      XMAX=(Dx*(n-1))/1.D3
      XMIN=Dx/1.D3
      YMAX=(Dy*(m-1))/1.D3
      YMIN=Dy/1.D3
C ----------------------------------------------------------------------
C   CRUSTAL AND LITHOSPHERIC MANTLE THICNESS. CONCENTRACIO DE DEFORMACIO (CD)
C ----------------------------------------------------------------------
C  elevacio > 0 -> per sobre el nivell del mar. 

      OPEN(1,FILE='velocity.xy')
      OPEN(3,FILE='e_sedsL_Tm_vis_Q_epeff_epzz.xy')
      OPEN(25,FILE='GLit.xy')
      WRITE(3,"('# TIME =',F9.4,' My')") Tanys/1.D6
      scrmin=500.D5
      scrmax=-500.D5
      GLmin=500.D5
      GLmax=-500.D5

      Rows: DO iy=0,m
         Columns: DO ix=0,n  
             kxy=ix+1+iy*(n+1)        
             XCORD=(ix*Dx)/1.D3
             YCORD=(iy*Dy)/1.D3 
             rom=RHOAST*(1.D0+(roalfa/2.D0)*(TISOTER-Tmoho(kxy)))	!! temperature depending
	  !   rom=3250.0							!! constant value
	     ELEVAC=elevation(kxy)
	     sed=sediment(kxy)
	     crust=scrust(kxy)
	     hlitos=GLit_ter(kxy)
	     epeffec=1.D0
	     epeffec=effective_strainrate (ix,iy,m,n,nn,Dx,Dy,u,v)
!	     CALL vertical_strain_rate (ix, iy, m, n, nn, Dx, Dy,
!     +					    u, v, ep_limit, epuntzz_P)
             epuntzz_P=epuntzz(kxy)
             
CC ***  ds/dt  *****
C		  sx=(scrust(kxy+1)-scrust(kxy-1))/(2.D0*Dx)
C		  sy=(scrust(kxy+n+1)-scrust(kxy-n-1))/(2.D0*Dy)
CC		  dsdt=(scrust(kxy)*epuntzz)
C		  dsdt=(scrust(kxy)*epuntzz)-
C     +  		(u(kxy)*sx+v(kxy)*sy)
CC     		  dhmdt=-dsdt
C     		  dhmdt=((GLit_ter(kxy)-scrust(kxy))/scrust(kxy))*dsdt
C     		  dedt=((RHOAST-roc)*dsdt+(RHOAST-rom)*dhmdt)/RHOAST
C     		  IF(ELEVAC.LT.0.D0) dedt=dedt*(RHOAST/(RHOAST-RHOH2O))
C     		  dedtmmany=dedt*FACVEL
C                  	IF(dedtmmany.GT.dedtmax) dedtmax=dedtmmany
C                  	IF(dedtmmany.LT.dedtmin) dedtmin=dedtmmany
	     WRITE(3,"(1X,7F11.4,E13.5,F9.3,1P,2G15.7)") XCORD, YCORD,
     +  		ELEVAC/1.D3,sediment(kxy)/1.D3,scrust(kxy)/1.D3,
     +     		GLit_ter(kxy)/1.D3,Tmoho(kxy),vis(kxy),
     +			Qsurface(kxy)*1.D3,epeffec,epuntzz_P
             WRITE(25,"(1X,3F13.4)") XCORD,YCORD,GLit(kxy)/1.D3
             WRITE(1,"(1X,2F11.4,2F13.6)") 
     +			XCORD,YCORD,u(kxy)*FACVEL,v(kxy)*FACVEL
	    scrmin=MIN(scrmin,scrust(kxy))
	    scrmax=MAX(scrmax,scrust(kxy))
	    GLmin=MIN(GLmin,GLit_ter(kxy))
	    GLmax=MAX(GLmax,GLit_ter(kxy))

	 END DO Columns
      END DO Rows
!      scrmin=MINVAL(scrust)
!      scrmax=MAXVAL(scrust)
!      GLmin=MINVAL(GLit_ter)
!      GLmax=MAXVAL(GLit_ter)
	
      WRITE(6,"(3X,'GRUIX CORTICAL MINIM:  ',F9.3,' km',5X,
     +        'GRUIX CORTICAL MAXIM: ',F10.3,' km'/
     +        3X,'GRUIX LITOSFERIC MINIM:',F9.3,' km',5X,
     +        'GRUIX LITOSFERIC MAXIM:',F10.3,' km'/)")
     +       	     scrmin/1.D3, scrmax/1.D3, GLmin/1.D3, GLmax/1.D3
!     +        3X,'d(elevacio)/dt MINIM:',F8.4,' mm/any',3X,
!     +        'd(elevacio)/dt MAXIM:',F8.4,' mm/any'/)")      
!     +       	     dedtmin, dedtmax
 
      CLOSE(25)
      CLOSE(3)
      CLOSE(1)

! ------  Different bodies  --------------------------------------------------
      OPEN(1,FILE='corba.tmp')	
      kpunts1=D_Bodies(1,1)
      qFlit1=D_Bodies(1,2)
      IF(kpunts1/=0) THEN
	DO kp=2,kpunts1+1
          WRITE(1,"(1X,2F10.2)") D_Bodies(kp,1)/1.D3,D_Bodies(kp,2)/1.D3
        END DO
      ENDIF 
      kpunts2=D_Bodies(kpunts1+2,1)			!! BODY 2
      IF(kpunts2/=0) THEN
        qFlit2=D_Bodies(kpunts1+2,2)
        WRITE(1,"('>')")
        DO kp=kpunts1+3,kpunts1+kpunts2+2
          WRITE(1,"(1X,2F10.2)") D_Bodies(kp,1)/1.D3,D_Bodies(kp,2)/1.D3
        END DO
      ENDIF 
      kpunts3=D_Bodies(kpunts1+kpunts2+3,1)		!! BODY 3
      IF(kpunts3/=0) THEN
        qFlit3=D_Bodies(kpunts1+kpunts2+3,2)
        WRITE(1,"('>')")
        DO kp=kpunts1+kpunts2+4,kpunts1+kpunts2+kpunts3+3
          WRITE(1,"(1X,2F10.2)") D_Bodies(kp,1)/1.D3,D_Bodies(kp,2)/1.D3
        END DO
      ENDIF
      kpunts4=D_Bodies(kpunts1+kpunts2+kpunts3+4,1)		!! BODY 4
      IF(kpunts4/=0) THEN
        qFlit4=D_Bodies(kpunts1+kpunts2+kpunts3+4,2)
        WRITE(1,"('>')")
        DO kp=kpunts1+kpunts2+kpunts3+5,
     +		kpunts1+kpunts2+kpunts3+kpunts4+4
          WRITE(1,"(1X,2F10.2)") D_Bodies(kp,1)/1.D3,D_Bodies(kp,2)/1.D3
        END DO
      ENDIF
      kpunts5=D_Bodies(kpunts1+kpunts2+kpunts3+kpunts4+5,1)		!! BODY 5
      IF(kpunts5/=0) THEN
        qFlit5=D_Bodies(kpunts1+kpunts2+kpunts3+kpunts4+5,2)
        WRITE(1,"('>')")
        DO kp=kpunts1+kpunts2+kpunts3+kpunts4+6,
     +		kpunts1+kpunts2+kpunts3+kpunts4+kpunts5+5
          WRITE(1,"(1X,2F10.2)") D_Bodies(kp,1)/1.D3,D_Bodies(kp,2)/1.D3
        END DO
      ENDIF



      CLOSE(1)
      OPEN(1,FILE='Points.tmp')	
      kpcurve=curve(0,1)
      IF(kpcurve/=0) THEN
 !       WRITE(1,"('>')")
        DO kp=1,kpcurve
          WRITE(1,"(1X,2F10.2)") curve(kp,1)/1.D3,curve(kp,2)/1.D3
        END DO
      ENDIF

      CLOSE(1)

C ---------------------------------------------------------------------
!      call Concentration_Deformation (Tanys, n, m, nn, Dx, Dy, 
!     +     			elevation, sediment, scrust, kpuntssd)

C ---------------------------------------------------------------------
       OPEN(1,FILE='titol.tmp')
       IXTITOL=0
       IYTITOL=INT(YMAX)-250
       TManys=REAL(Tanys/1.D6)
       IZHEXP=HEXP/1.D3
              
      IF(qFlit1==0) qFlit1=1.0
      IF(qFlit2==0) qFlit2=1.0
      IF(qFlit3==0) qFlit3=1.0

      WRITE(1,"(6I5,2X,A150)") IXTITOL,IYTITOL,11,0,5,0,TITLERE
      
	iyt2=IYTITOL-YMAX/10		! IYTITOL-100,
      IF(qFlit1==1.AND.qFlit2==1.AND.qFlit3==1) THEN
       	   WRITE(1,"(6I5,' H=',F5.2,'exp(-z/',G9.3,')  @~m@~W/m@+3@+ ,',
     +		   '   T@-litos@-= ',F5.0,'@+o@+C')") IXTITOL,iyt2,
     +		     11,0,5,0,HSURF*1.D6,HEXP,TISOTER-273
	ELSE
       	   WRITE(1,"(6I5,' H=',F5.2,'exp(-z/',G9.3,')  @~m@~W/m@+3@+',
     + 	    ',   Q@-1@-=',F9.1,',  Q@-2@-=',F9.1,',  Q@-3@-=',F9.1)")
     +           IXTITOL,iyt2,11,0,5,0,
     +                 HSURF*1.D6,HEXP,qFlit1,qFlit2,qFlit3
      END IF
      
	iyt3=IYTITOL-YMAX/5		! IYTITOL-200,
      IF(Kerosdif==0.AND.rain==0.AND.Krain==0) THEN
       	   WRITE(1,"(6I5,' TIME = ',F6.2,' My')")
     +           IXTITOL,iyt3,11,0,5,0,TManys
       ELSE
		  WRITE(1,"(6I5,' TIME = ',F6.2,' My,   Kerosdif =',F6.1,
     +	   ' m@+2@+/a, rain =',F6.1,' mm/a, Krain =',F6.1)")
     +       IXTITOL,iyt3,11,0,5,0,TManys,Kerosdif,rain,Krain

 !      	   WRITE(1,"(6I5,' TIME = ',F5.1,' My,   Kerosdif =',F6.1,
 !    + 		   ' m@+2@+/a,  rain =',
 !    +		   F6.1,' mm/a,  Krain =',F6.1,' l/m@+2@+/a/km')")
 !    +           IXTITOL,IYTITOL-200,11,0,5,0,TManys,Kerosdif,rain,Krain
      END IF
      CLOSE(1)
       
      OPEN(1,FILE='TITOL.tmp')
      WRITE(1,1105) TITLERE,
     +               HSURF*1.D6,IZHEXP,
!     +		     roc,RHOAST,CONDUC(1),CONDUC(2),
     +               TISOTER-273,Kerosdif,rain
!     + qFlit1,qFlit2
 1105 FORMAT(A150,
     +       /' H=',F5.2,'exp(-z/',I2,') ,',
!     +      ' @~r@~@-c@-= ',F5.0,'kg/m@+3@+, @~r@~@-a@-= ',F5.0,
!     +      'kg/m@+3@+ , K@-c@-= ',F4.2,', K@-m@-= ',F4.2, 
     +    ' T@-L@-= ',F5.0,'@+o@+C, Kerosdif=',F6.1,' m@+2@+/a, rain=',
     +       F6.1,' mm/a,  Krain=')
!     +	    '   Q@-1@-=',F4.1,',   Q@-2@-=',F4.1)

       CLOSE(1)

 128   CONTINUE
C **********************************************************************

       OPEN(8,FILE='limits.d')
           ZX=120.0
           ZY=120.0
           XMAXVEL=(Dx*n)/1.D3
           YMAXVEL=(Dy*m)/1.D3
           write(8,"(6F10.3)") XMIN,XMAX,YMIN,YMAX, 
     +                 XMAXVEL,YMAXVEL
       CLOSE(8)

 901   CONTINUE
C        GOTO 9000
CC ********************************************************************
CC --------------------------------------------------------------------
C        CALL TEMPEFOND (m, n, NELZ, nnsd, nzsd, Dx, Dy, Dz, 
C     +       TEMPE, scrust)
CC *********************************************************************
CC ---------------------------------------------------------------------

          CALL PRINCIPAL_STRESS (m, n, nn, AX, BY, u, v, vis, 
     +            PI, FACVEL, scrust, GLit_ter, epuntzz) 
CC *********************************************************************
CC --------------------------------------------------------------------
CC  TROBA LA VELOCITAT MITJA PER SOBRE I PER SOTA DE LA FALLA:
C          CALL vel_EUAF (n, m, nn, AX, BY, u, v, FACVEL)
          
CC *********************************************************************

 9000  CONTINUE
       PRINT*,'   PROGRAMA graficsth FINALITZAT '
       STOP 
       END PROGRAM graficsth

C **********************************************************************
C *************  CALCUL DE LA CONCENTRACIO DE DEFORMACIO  **************

      SUBROUTINE Concentration_Deformation (Tanys, n, m, nn, Dx, Dy,
     +			          elevation, sediment, scrust, kpuntssd)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*170 TITLE0
      CHARACTER*15 fresult0
      DIMENSION elevation(nn),sediment(nn),scrust(nn),elevation0(nn),
     +          sediment0(nn),scrust0(nn),GLit_ter0(nn),GLit0(nn),
     +          u0(nn),v0(nn),vis0(nn),epuntzz0(nn),Qsurface0(nn),
     +          Tmoho0(nn),D_Bodies_0(kpuntssd,2),scrint(nn),sord(nn),
     +		D_var(nn),Delta_ord(nn),curve(0:kpuntssd,2)

      DATA VOLTOT/0.D0/ VOLPAR/0.D0/ APARVOL/0.D0/ Deftot/0.D0/
     +     APARDEF/0.D0/

      nnint=(n-1)*(m-1)           
      AELEM=Dx*Dy
      
      TP100CD=0.05D0			! % del Volum cortical
      kk=1
      DO ix=1,n-1			! Only internal points
         DO iy=1,m-1      
             kxy=ix+1+iy*(n+1)
             scrint(kk)=scrust(kxy)
             kk=kk+1
             VOLTOT=VOLTOT+scrust(kxy)*AELEM
	 END DO     
      END DO    
      VTP=VOLTOT*TP100CD  		
       
      CALL ORDRE_MAG (scrint, nnint, sord)
      DO IK=1,nnint   
             VOLPAR=VOLPAR+sord(IK)*AELEM
             APARVOL=APARVOL+AELEM
             IF(VOLPAR >= VTP) THEN
                  VSUF=sord(IK)*AELEM-(VOLPAR-VTP)
                  ASUF=VSUF/sord(IK)
                  VOLPAR=VOLPAR-sord(IK)*AELEM+VSUF
                  APARVOL=APARVOL-AELEM+ASUF
                  GOTO 10
             ENDIF     
             nkpuntA=IK
      END DO     
 10   CONTINUE 
      VVV=(VOLPAR*100.D0)/VOLTOT
      AAAVOL=(APARVOL*100.D0)/(nnint*AELEM)
         PRINT*,'  '
         PRINT*,' nkpuntA, DE LA AREA :',nkpuntA
         PRINT*,'% de area sota la que hi ha el',VVV,'% del volum ',
     +           'total:',AAAVOL,'%'


!   Calcul de la deformacio a apartir dels valors inicials
      nporta=17
      fresult0='resultats0'
      CALL RERESULT (fresult0, nporta, TITLE0, Tsegons0, Tanys0,
     +               m, n, NELZ, nn, nnsd, nzsd, AX, BY, Dz, TLITOS,
     +               scrust0, GLit_ter0, GLit0, vis0, u0, v0, epuntzz0,
     +               Tmoho0, Qsurface0, elevation0, sediment0, 
     +		     kpuntssd, D_Bodies_0, curve)


      TP100CD=0.05D0			! % of crustal deformation	
      D_total=0.D0
      D_var=0.D0
      Delta_ord=0.D0
      kk=1
      DO ix=1,n-1			! Only internal points
         DO iy=1,m-1      
             kxy=ix+1+iy*(n+1)
             D_var(kk)=ABS(scrust(kxy)-scrust0(kxy))  
             D_total=D_total+D_var(kk)
             kk=kk+1
         END DO    
      END DO    
      Delta_TP=D_total*TP100CD
      
      Def_elev=0.0
      Area_Delev=0.0
      IF(D_total.NE.0.D0) THEN
           CALL ORDRE_MAG (D_var, nnint, Delta_ord)
           Delta_par=0.D0
           APARDEF=0.D0
           DO IK=1,nnint
               Delta_par=Delta_par+Delta_ord(IK)
               APARDEF=APARDEF+AELEM
               IF(Delta_par >= Delta_TP) THEN  
                    EXCES=(Delta_par-Delta_TP)/Delta_ord(IK)
                    Delta_par=Delta_par-EXCES*Delta_ord(IK)
                    APARDEF=APARDEF-EXCES*AELEM
                    GOTO 13 
               ENDIF             
               nkpuntD=IK
           END DO     
 13        CONTINUE 
           Def_crust=(Delta_par*100.D0)/D_total
           Area_Dcrust=(APARDEF*100.D0)/(nnint*AELEM)
      ENDIF             
         PRINT*,'  '
         PRINT*,' nkpuntD, DE LA DEFORMACIO :',nkpuntD
         PRINT*,'% de area sota la que hi ha el',Def_crust,'% de la ',
     +           'deformacio total:',Area_Dcrust,'%'


      TP100CD=0.05D0			! elevation variations
      D_total=0.D0
      kk=1
      DO ix=1,n-1			! Only internal points
         DO iy=1,m-1      
             kxy=ix+1+iy*(n+1)
             D_var(kk)=ABS(elevation(kxy)-elevation0(kxy))  
             D_total=D_total+D_var(kk)
             kk=kk+1
         END DO    
      END DO    
      Delta_TP=D_total*TP100CD
      
      Def_elev=0.0
      Area_Delev=0.0
      IF(D_total.NE.0.D0) THEN
           CALL ORDRE_MAG (D_var, nnint, Delta_ord)
           Delta_par=0.D0
           APARDEF=0.D0
           DO IK=1,nnint
               Delta_par=Delta_par+Delta_ord(IK)
               APARDEF=APARDEF+AELEM
               IF(Delta_par >= Delta_TP) THEN  
                    EXCES=(Delta_par-Delta_TP)/Delta_ord(IK)
                    Delta_par=Delta_par-EXCES*Delta_ord(IK)
                    APARDEF=APARDEF-EXCES*AELEM
                    GOTO 23 
               ENDIF             
               nkpuntD=IK
           END DO     
 23        CONTINUE 
           Def_elev=(Delta_par*100.D0)/D_total
           Area_Delev=(APARDEF*100.D0)/(nnint*AELEM)
      ENDIF             

         PRINT*,' Punts DE LA DEFORMACIO :',nkpuntD
         PRINT*,'% de area sota la que hi ha el',Def_elev,'% de la ',
     +           'elevation variation total:',Area_Delev,'%'
     
!   Sediments created Volum
      D_sediment=0.D0
      DO ix=1,n-1			! Only internal points
         DO iy=1,m-1      
            kxy=ix+1+iy*(n+1)
            D_sediment=D_sediment+(sediment(kxy)-sediment0(kxy))
         END DO    
      END DO         
      D_sed_average=D_sediment/nnint
!  Time, %Volum, Area sota la que hi ha el %Volum,
!       %Deformacio, Area sota la que hi ha el %Deformacio   
       OPEN(1,FILE='t_V_A_Def')
       WRITE(1,"('# Time     %Volum  Area vol'
     + 	         '  %DCrustal Area crust '
     +           '  %Dtopo Area topo    Dsediments average (m)')")
            WRITE(1,"(1X,F5.1,3(F10.1,F10.5),F15.3)") 
     +                Tanys/1.D6,VVV,AAAVOL,Def_crust,Area_Dcrust,
     +			Def_elev,Area_Delev,D_sed_average
       CLOSE(1)  
      RETURN 
      END SUBROUTINE Concentration_Deformation

C ******************************************************************
C            SUBROUTINE PRINCIPAL_STRESS 
C CALCULA ELS ESFOROS xx yy xy i zz. 
C UN COP OBTINGUTS ELS ESFORSOS xx,yy,xy TROBO ELS PRINCIPALS
C EU,ED: ESFOROS PRINCIPALS (EU < ED)  NEGATIU=COMPRESSIU. 
C EZZ: ESFOROS VERTICAL
C NF: Normal Fault (EZZ<EU<ED).
C SS: Strike-Slip Fault (EU<EZZ<ED).
C TF: Thrust Fault (EU<ED<EZZ).

       SUBROUTINE PRINCIPAL_STRESS (m, n, nn, AX, BY, u, v, vis, 
     +             PI, FACVEL, scrust, GLit_ter, epuntzz)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION vis(nn),u(nn),v(nn),scrust(nn),GLit_ter(nn),epuntzz(nn),
     +          azimuth(nn),epunt1(nn),epunt2(nn),esc_srsuau(nn)
      CHARACTER*5 CHFALLAe,CHFALLAs
 
       Dx=AX/n
       Dy=BY/m
      WRITE(6,62) AX/1D3,BY/1D3,Dx/1D3,Dy/1D3,n,m
 62   FORMAT(/' SUBROUTINE PRINCIPAL_STRESS'//
     +       1X,'AX=',F9.3,' km     BY: ',F9.3,' km '/
     +       1X,'Dx:',F9.3,' km     Dy: ',F9.3,' km '/
     +       1X,'n:',I3,'     m:',I3/)
        
 87    FORMAT(4F10.3,2F8.1,2F10.3,3F8.1)
 88    FORMAT(4F10.3,2F8.1,2F10.3,4F8.1)

      OPEN(1,FILE='principal_stress.xy')
      OPEN(10,FILE='epunt12zz.xy')
C      OPEN(3,FILE='stress_azimuth.xy')

      PRINT*,' Per el fitxer: stress_azimuth.xy, he de calcular ',
     +    'els esforcos a tots els punts. Mirar que NINTX i NINTY ',
     +    'valguin 1'
   
CC NOMES BUSCO ELS TERMES CADA 'NINT' 
         NINTSTR=1
C         NINTX=2
C         NINTY=3
         NINTX=1
         NINTY=1
          
        PRINT*,'NINTX =',NINTX,'   NINTY =',NINTY 
         MGUA=INT((m-1)/NINTY)
         NGUA=INT((n-1)/NINTX)
         NNI=0 
        npstress=0
      DO 5 iyint=1,MGUA
         iy=NINTY*iyint
         Y=(iy*Dy)/1.D3
         DO 5 ixint=1,NGUA 
             ix=NINTX*ixint
              X=(ix*Dx)/1.D3
              npstress=npstress+1 
C              PRINT*,'punts stress:',npstress
              NNI=NNI+1
              kxy=ix+1+iy*(n+1)
                 ux=(u(kxy+1)-u(kxy-1))/(2.D0*Dx)
                 vx=(v(kxy+1)-v(kxy-1))/(2.D0*Dx)
                 vy=(v(kxy+n+1)-v(kxy-n-1))/(2.D0*Dy)
                 uy=(u(kxy+n+1)-u(kxy-n-1))/(2.D0*Dy)
		vorticity=vx-uy
 !              epuntzz(kxy)=-(ux+vy)
               EXX=2.D0*vis(kxy)*ux
               EYY=2.D0*vis(kxy)*vy
               EXY=vis(kxy)*(uy+vx)
               EZZ=-2.D0*vis(kxy)*(ux+vy)
         
CC   PRINCIPAL STRESS
               IF(EXY.EQ.0.D0) THEN
                   EU=EXX
                   ED=EYY
                   TETA=0.D0
                   TETAG=0.D0
                   GOTO 100
               ENDIF
               IF(EXX.EQ.EYY) THEN
                   EU=EXX+EXY
                   ED=EXX-EXY
                   TETA=PI/4.D0
                   TETAG=45.D0
               ELSE
                   A1=(EXX+EYY)/2.D0
                   B1=((EXX-EYY)*(EXX-EYY))/4.D0
                   B2=EXY*EXY
                   BS=DSQRT(B1+B2)
                   EU= A1+BS
                   ED=A1-BS
                   CC=(2.D0*EXY)/(EXX-EYY)
                   TETA=(1.D0/2.D0)*DATAN(CC)
                   TETAG=(TETA*180.D0)/PI
               ENDIF
100           CONTINUE
              IF(EU.GT.ED) THEN
                   EUTMP=EU
                   EU=ED
                   ED=EUTMP
              ENDIF 
C Quin angle correspon a cada 'principal stress'
              TERM1=(EU+ED)/2.D0
              TERM2=((EU-ED)/2.D0)*DCOS(2.D0*TETA)
              EXXPR=TERM1+TERM2
              ZERO1=ABS(ABS(EXX)-ABS(EXXPR))
              ZERO2=ABS(ABS(EYY)-ABS(EXXPR))
              IF(ZERO2.LT.ZERO1) THEN
                  TETA=TETA+PI/2.D0
                  TETAG=TETAG+90.D0
              ENDIF
              TERM2=((EU-ED)/2.D0)*DCOS(2.D0*TETA)
              EXXPR=TERM1+TERM2
              EYYPR=TERM1-TERM2
              EXYPR=-((EU-ED)/2.D0)*DSIN(2.D0*TETA)  

              epunt1(kxy)=EU/(2.D0*vis(kxy))
              epunt2(kxy)=ED/(2.D0*vis(kxy))
              ep1mep2=(epunt1(kxy)-epunt2(kxy))

CC  -------------------------------------------------------------
              IF(NNI.EQ.NINTSTR) THEN

CC  TIPUS DE FALLA A PARTIR DELS STRAIN RATES PRINCIPALS (TFALLA)
CC	 CHFALLAe=TF (FALLA INVERSA)   CHFALLAe=NF (FALLA NORMAL)
CC  	 CHFALLAe=SS  (STRIKE SLIP)
C                    IF(epunt1(kxy).LT.0.D0.AND.epunt2(kxy).LT.0.D0)
C     +                      CHFALLAe="TF"
C                    IF(epunt1(kxy).GT.0.D0.AND.epunt2(kxy).GT.0.D0)
C     +                      CHFALLAe="NF"
C                    IF(epunt1(kxy).LT.0.D0.AND.epunt2(kxy).GT.0.D0)
C     +                      CHFALLAe="SS"
C                    IF(epunt1(kxy).GT.0.D0.AND.epunt2(kxy).LT.0.D0)
C     +                      CHFALLAe="SS"
C      
CC  TIPUS DE FALLA A PARTIR DELS ESFORCOS (TFESF)
                    IF(EZZ.LT.EU) THEN
                    	IF(EU.GE.0.D0) THEN
	         	   CHFALLAs="NF"
	         	ELSE
	         	   CHFALLAs="NS"
	         	ENDIF  
                    ENDIF
	            IF(EZZ.GE.EU.AND.EZZ.LE.ED) CHFALLAs="SS"
                    IF(EZZ.GT.ED) THEN
                    	IF(ED.LE.0.D0) THEN
                           CHFALLAs="TF"
                        ELSE 
                           CHFALLAs="TS"
                        ENDIF    
                    ENDIF 

                   WRITE(1,112) X,Y,TETAG,EU,ED,EZZ,CHFALLAs
                   NNI=0
              ENDIF
CC  -------------------------------------------------------------
CC  Passo els valors a azimut i els poso en el primer o segon quadrant
CC     0 <= azimuth < 180
            azimuth(kxy)=90.D0-TETAG
            IF(azimuth(kxy).LT.0.D0) azimuth(kxy)=azimuth(kxy)+360.D0
 123        IF(azimuth(kxy).GT.180.D0) THEN
                   azimuth(kxy)=azimuth(kxy)-180.D0
                   GOTO 123
            ENDIF 
C            WRITE(3,110) X*1.D3,Y*1.D3,azimuth(kxy)
CC  -------------------------------------------------------------

           WRITE(10,"(1X,2F10.2,3G16.7,1F10.2,G16.7)") X,Y,epunt1(kxy),
     + 		epunt2(kxy),epuntzz(kxy),azimuth(kxy),vorticity
 5     CONTINUE
C       WRITE(10,*) '-30  -40 ',ITF
C       WRITE(10,*) '650 -40 ',INF
C       WRITE(10,*) '1270 -40 ',ISS
        PRINT*,'punts on he calculat el stress:',npstress
 110   FORMAT(1X,3F14.1) 
 111   FORMAT(1X,4E17.5)
 112   FORMAT(1X,3F10.3,3E17.5,3X,1A5)
 333   FORMAT(1X,2F10.2)
       CLOSE(1)
       CLOSE(10)
C       CLOSE(3)
CCCC ------  Comparar amb els azimuths d'una base de dades  -------
C      PRINT*,' '
C      PRINT*,'******************** CALL stress_azimuth ****************'
C       CALL stress_azimuth (n, m, nn, AX, BY, azimuth)   
CCCC --------------------------------------------------------------- 
CCC -- Trobar un escalar per l'strain rate i suavitzar-lo a la falla --
!       PRINT*,'*************** CALL suau_strainrate *******************'
!       CALL suau_strainrate (n, m, nn, AX, BY, epunt1, epunt2, epuntzz,
!     +                       esc_srsuau, u, v, FACVEL)
CCCC ---------------------------------------------------------------
CCCC     Coeficient de correlacio entre:
C  strain rate del model i strain rate sismic (a partir de terratremols)
C       PRINT*,'************** CALL correlation ************************'
C       CALL correlation (n, m, nn, AX, BY, esc_srsuau)   
CCCC ---------------------------------------------------------------
        IF(NINTX.NE.1.OR.NINTY.NE.1) WRITE(6,666)
 666    FORMAT(/'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',
     +            '!!!!!!!!!!!!!!!!'/
     +  '!!! NO CALCULA ELS ESFORCOS PRINCIPALS EN TOTS ELS PUNTS !!!'/
     +  '!!!  PER TANT LES DESVIACIONS TROBADES NO SERAN BONES    !!!'/
     +  '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',
     +  '!!!!!')
     
      RETURN
      END
CC ****************************************************************
CC ****************************************************************
CC  Compara amb els azimuths d'una base de dades 
 
       SUBROUTINE stress_azimuth (n, m, nn, AX, BY, azimuth)
 
       IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
       CHARACTER*70 DATAFIL
       CHARACTER*2 TFAULT,QUALITY
       PARAMETER(NMAX=20000,DQUALITY_A=5.D0,DQUALITY_B=4.D0,
     +       DQUALITY_C=3.D0,DQUALITY_D=2.D0,DQUALITY_E=1.D0)
       DIMENSION azimuth(nn)
      
C       DATAFIL='/model/ivone/DADES/stress/wsm_zone.xy'  
       DATAFIL='/model/ivone/DADES/Atlantic_Octubre99/wsm_zone.xy'  
       OPEN(3,FILE='azuhu_azdata.xy')
       Dx=AX/n
       Dy=BY/m
       Xmin=Dx
C       Xmax=(n-1)*Dx
C       Xmax=2250.D3
       Xmax=1764.D3
       Ymin=Dy
       Ymax=(m-1)*Dy  
       ixmin=Xmin/Dx
       ixmax=Xmax/Dx
       iymin=Ymin/Dy
       iymax=Ymax/Dy
       PRINT*,'Legeixo el fitxer:',DATAFIL
       Deltamax=0.D0
       Deltamin=1.D5
       Deltaaver=0.D0
       DeltaaverQ=0.D0
       DQ_sum=0.D0
       OPEN(UNIT=10,FILE=DATAFIL) 
       NDATA=0  
       NDATAIN=0
       DO 100 i=1,NMAX 
          READ(10,*,IOSTAT=IOS) XDATA,YDATA,AZDATA,TFAULT,QUALITY
            IF (IOS < 0) THEN
               WRITE(6,*) 'END OF FILE. Number of poits:',NDATA
               GO TO 120
           ELSE IF (IOS > 0) THEN
               WRITE (*, 90) NDATA
   90          FORMAT (/' ERROR: Encountered bad data after ',
     +                  I6,' successful READs.')
               STOP
           END IF
CC --  l'azimut al 1 i 2 quadrants:  
           IF(AZDATA.LT.0.D0) AZDATA=AZDATA+360.D0
 80        IF(AZDATA.GT.180.D0) THEN
                   AZDATA=AZDATA-180.D0
                   GOTO 80
            ENDIF
CC --------
           NDATA=NDATA+1
           ix=XDATA/Dx
           iy=YDATA/Dy  
         IF(ix.LE.ixmin.OR.ix.GE.ixmax.OR.iy.LE.iymin.OR.iy.GE.iymax)
C     +         .OR.QUALITY.EQ."E")
C     +         .OR.QUALITY.EQ."D".OR.QUALITY.EQ."C".OR.QUALITY.EQ."B")
     +         GOTO 100
           NDATAIN=NDATAIN+1
           IF(QUALITY.EQ."E") DQUALITY=DQUALITY_E
           IF(QUALITY.EQ."D") DQUALITY=DQUALITY_D
           IF(QUALITY.EQ."C") DQUALITY=DQUALITY_C
           IF(QUALITY.EQ."B") DQUALITY=DQUALITY_B
           IF(QUALITY.EQ."A") DQUALITY=DQUALITY_A
           kxy=ix+1+iy*(n+1)
           azimuth1=azimuth(kxy)
           Delta1=ABS(AZDATA-azimuth1)
           Delta2=180.D0-Delta1 
           Deltaaz1=MIN(Delta1,Delta2) 
           IF((ix+1).LE.ixmax.AND.(iy+1).LE.iymax) THEN
                 distan1=(XDATA-ix*Dx)**2+(YDATA-iy*Dy)**2
                 distan2=(XDATA-(ix+1)*Dx)**2+(YDATA-iy*Dy)**2
                 distan3=(XDATA-ix*Dx)**2+(YDATA-(iy+1)*Dy)**2
                 distan4=(XDATA-(ix+1)*Dx)**2+(YDATA-(iy+1)*Dy)**2
                 azimuth2=azimuth(kxy+1)
                 azimuth3=azimuth(kxy+n+1)
                 azimuth4=azimuth(kxy+n+2)
                 azim_mig=(distan2*distan3*distan4*azimuth1+
     +                 distan1*distan3*distan4*azimuth2+
     +                 distan1*distan2*distan4*azimuth3+
     +                 distan1*distan2*distan3*azimuth4)/
     +                (distan2*distan3*distan4+distan1*distan3*distan4+
     +                 distan1*distan2*distan4+distan1*distan2*distan3)
                 Del_mig1=ABS(AZDATA-azim_mig)
                 Del_mig2=180.D0-Del_mig1
                 Deltaaz=MIN(Del_mig1,Del_mig2) 
                 XPUNT=XDATA
                 YPUNT=YDATA
            ELSE
                 XPUNT=ix*Dx
                 YPUNT=iy*Dy
                 Deltaaz=Deltaaz1 
                 azim_mig=azimuth1
           ENDIF     
C           WRITE(6,63) azimuth(kxy),AZDATA,Deltaaz
C 63     FORMAT(1X,5F10.2)
           Deltamax=MAX(Deltamax,Deltaaz)
           Deltamin=MIN(Deltamin,Deltaaz)
           Deltaaver=Deltaaver+Deltaaz 
           DeltaaverQ=DeltaaverQ+Deltaaz*DQUALITY
           DQ_sum=DQ_sum+DQUALITY
           WRITE(3,33) XPUNT/1.D3,YPUNT/1.D3,azim_mig,
     +               XDATA/1.D3,YDATA/1.D3,AZDATA,TFAULT 
 33      FORMAT(1X,2F12.1,F8.2,2F12.1,F8.2,A4)
  100 CONTINUE
      WRITE(6,110)NMAX,NDATA
  110 FORMAT(/' Dataset contained more than ',
     &          I5,' data; first ',I5,' were read.')
  120 CONTINUE
      WRITE(6,130)NDATA,NDATAIN
  130 FORMAT(//6X,'Reading of data completed:',I7,' data points.'/
     +         6X,'Punts dins de la finestra escollida:',I7/)   
      Deltaaver=Deltaaver/NDATAIN
      DeltaaverQ=DeltaaverQ/DQ_sum

      CLOSE(10)
      WRITE(3,34) DeltaaverQ
 34      FORMAT(1X,1F10.2)      
      CLOSE(3)  
       WRITE(6,6) (ixmin*Dx)/1D3,(ixmax*Dx)/1D3,
     +                    (iymin*Dy)/1D3,(iymax*Dy)/1D3
 6    FORMAT(/2X,'Finestra on comparare les direccions d esforcos:'/
     +      8X,'Xmin =',F10.3,' km',5X,'Xmax =',F10.3,/
     +      8X,'Ymin =',F10.3,' km',5X,'Ymax =',F10.3/)
       
      WRITE(6,65) Deltamin,Deltamax,Deltaaver,DQUALITY_A,DQUALITY_B,
     +       DQUALITY_C,DQUALITY_D,DQUALITY_E,DeltaaverQ
 65   FORMAT(/3X,'DESVIACIO MINIMA:',F8.2,' graus',3X,
     +        'DESVIACIO MAXIMA:',F8.2,' graus'/
     +        8X,'DESVIACIO MITJA:',F8.2,' graus'//
     +        6X,'PESOS SEGONS LA QUALITAT DE LES DADES:'/
     +        9X,'QUALITAT A:',F5.2,'    QUALITAT B:',F5.2,/
     +        9X,'QUALITAT C:',F5.2,'    QUALITAT D:',F5.2,/
     +        9X,'QUALITAT E:',F5.2,/
     +        9X,'DESVIACIO MITJA AMB AQUESTS PESOS:',F8.2,' graus'/)

      RETURN
      END   
CC ****************************************************************
CC              suau_strainrate
CC  Per suavitzar l'strain rate al voltant de la falla

       SUBROUTINE suau_strainrate (n, m, nn, AX, BY, epunt1, epunt2,
     +       epuntzz, esc_srsuau, u, v, FACVEL)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
       CHARACTER*20 DATAFIL
       PARAMETER(NFMAX=200,TP0=0.3D0,TPDx=0.4D0,TP2Dx=0.3D0)
       DIMENSION epunt1(nn),epunt2(nn),epuntzz(nn),esc_strainr(nn),
     +           esc_srsuau(nn),u(nn),v(nn)
 
       WRITE(6,1) TP0,TPDx,TP2Dx
 1     FORMAT(1X,/'  SUBROUTINE suau_strainrate:'/
     +        '  Suavitzo el maxim strain rate, amb els valors:'/
     +              F5.2,' / ',F5.2,' / ',F5.2,' / '/)
       Dx=AX/n
       Dy=BY/m 
       TPtotal=TP0+TPDx+TP2Dx
       DIFTP=TPtotal-1.D0
       IF(ABS(DIFTP).GT.1.D-13) THEN
             WRITE(6,2) TPtotal 
 2           FORMAT(1X,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',
     +   '!!!!!!!!!!!!!!!!'/'La suma dels tants per cent no val 1: ',
     +    F6.4,' Aturo el programa'/'!!!!!!!!!!!!!!!!!!!!!!!!!!!!',
     +   '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') 
             STOP
       ENDIF
CC ---- Busco un escalar per l'strain rate  ---- 
       DO 3 ix=1,n-1
          DO 3 iy=1,m-1 
              kxy=ix+1+iy*(n+1)
C              esc_strainr(kxy)=MAX(ABS(epunt1(kxy)),ABS(epunt2(kxy)),
C     +                             ABS(epuntzz(kxy))) 
              esc_strainr(kxy)=(ABS(epunt1(kxy))+ABS(epunt2(kxy))+
     +                       ABS(epuntzz(kxy)))/2.D0  
              esc_srsuau(kxy)=esc_strainr(kxy)
 3     CONTINUE             
CCC -----  Llegeix d'un fitxer la falla (punts de debilitat) ------- 
       DATAFIL='punts_falla_xy.res'
       WRITE(6,6)DATAFIL
 6     FORMAT(3X,'Llegeix els punts del fitxer:',A20) 
        OPEN(16,FILE=DATAFIL)  
         NPFALLA=0  
         DO 100 nf=1,NFMAX
             READ(16,*,IOSTAT=IOS) xfalla,yfalla,qpfalla
             IF (IOS < 0) THEN
                   GO TO 120
             ELSE IF (IOS > 0) THEN
               WRITE (*, 90)NPFALLA 
   90          FORMAT (/' ERROR: Encountered bad data after ',
     +                  I6,' successful READs.'/)
               STOP
             ENDIF
             NPFALLA=NPFALLA+1
             ix=xfalla/Dx
             iy=yfalla/Dy
             kxy=ix+1+iy*(n+1)
             kxy1p=ix+1+(iy+1)*(n+1)
             kxy2p=ix+1+(iy+2)*(n+1)
             kxy1n=ix+1+(iy-1)*(n+1)
             kxy2n=ix+1+(iy-2)*(n+1)
             esc_srsuau(kxy)=TP0*esc_strainr(kxy)+
     +            TPDx*((esc_strainr(kxy1p)+esc_strainr(kxy1n))/2.D0)+
     +            TP2Dx*((esc_strainr(kxy2p)+esc_strainr(kxy2n))/2.D0)
             esc_srsuau(kxy1p)=TP0*esc_strainr(kxy1p)+
     +                 TPDx*esc_strainr(kxy)+TP2Dx*esc_strainr(kxy2p)
             esc_srsuau(kxy1n)=TP0*esc_strainr(kxy1n)+
     +                 TPDx*esc_strainr(kxy)+TP2Dx*esc_strainr(kxy2n) 
             esc_srsuau(kxy2p)=TP0*esc_strainr(kxy2p)+
     +                 TPDx*esc_strainr(kxy1p)+TP2Dx*esc_strainr(kxy) 
             esc_srsuau(kxy2n)=TP0*esc_strainr(kxy2n)+
     +                 TPDx*esc_strainr(kxy1n)+TP2Dx*esc_strainr(kxy)

  100    CONTINUE
         WRITE(6,110)NFMAX,NPFALLA
  110    FORMAT(/' Dataset contained more than ',
     &          I5,' data; first ',I5,' were read.'/)
  120    CONTINUE
         WRITE(6,130)NPFALLA
  130    FORMAT(3X,' Reading of points fault completed:',I7,
     +             ' data points.'///)
        CLOSE(16)
 77    FORMAT(1X,A3,2F10.3,2F11.4)
CCCC  -------------------------------------------------- 
        OPEN(16,FILE='escalar_strainrate.xy')
        DO 53 iy=1,m-1
            ykm=(iy*Dy)/1.D3
            DO 53 ix=1,n-1 
               xkm=(ix*Dx)/1.D3
               kxy=ix+1+iy*(n+1)
               WRITE(16,55) xkm,ykm,esc_strainr(kxy),esc_srsuau(kxy)
 53     CONTINUE    
 55    FORMAT(1X,2F12.1,2G15.6) 
       CLOSE(16)   
       RETURN
       END  
CCC ****************************************************************
CCC           SUBROUTINE  correlation
CCC    CORRELATION COEFFICIENT STRAIN RATES . DELS LOGARITMES.

       SUBROUTINE correlation (n, m, nn, AX, BY, esc_srsuau)  
       IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
       CHARACTER*70 FSMOOT,TBANDA
       DIMENSION esc_srsuau(nn),epuntlog(nn),epEQlog(nn) 
       DATA ep_mean/0.D0/ EQ_mean/0.D0/ sprod/0.D0/ sumEQ/0.D0/ 
     +      sumepunt/0.D0/     

C       FSMOOT='/model/ivone/DADES/Earthquakes/smooth_ISC_b127.xy'
       FSMOOT='/model/ivone/DADES/Atlantic_Octubre99/seismic_ep.xy'
C       FSMOOT='/model/ivone/DADES/Earthquakes/smooth_CNSS_b127.xy'
 160   FORMAT(/' LLEGEIXO EL FITXER AMB EL STRAIN RATE SISMIC:'/A70/)
       WRITE(6,160) FSMOOT
       Dx=AX/n
       Dy=BY/m    
       OPEN(16,FILE=FSMOOT)
 166   FORMAT(A70)       
       READ(16,166) TBANDA
       PRINT*,'EL SUAVITZAT DELS TERRATREMOLS HA ESTAT FET AMB UNA',
     +        ' GAUSSIANA DE BANDA:'
       WRITE(6,166) TBANDA
       DO 9 iy=1,m-1 
          ynode=iy*Dy
          DO 9 ix=1,n-1 
              xnode=ix*Dx
              kxy=ix+1+iy*(n+1)
              READ(16,*) xfile,yfile,epEQ
                epEQlog(kxy)=log10(epEQ) 
                epuntlog(kxy)=log10(esc_srsuau(kxy))
                  IF(xfile.NE.xnode.OR.yfile.NE.ynode) THEN
                     PRINT*,'EL fitxer:',FSMOOT,' no esta ben ordenat'
                     STOP
                  ENDIF 
                EQ_mean=EQ_mean+epEQlog(kxy)
                ep_mean=ep_mean+epuntlog(kxy)
 9     CONTINUE  
       CLOSE(16) 
       EQ_mean=EQ_mean/((n-1)*(m-1)) 
       ep_mean=ep_mean/((n-1)*(m-1)) 
      
       OPEN(16,FILE='epsismic_uhu.xy')
       DO 15 iy=1,m-1 
         DO 15 ix=1,n-1
            kxy=ix+1+iy*(n+1)
            sprod=sprod+(epEQlog(kxy)-EQ_mean)*(epuntlog(kxy)-ep_mean)
            sumEQ=sumEQ+((epEQlog(kxy)-EQ_mean)**2)
            sumepunt=sumepunt+((epuntlog(kxy)-ep_mean)**2) 
            WRITE(16,55) (ix*Dx)/1.D3,(iy*Dy)/1.D3,epEQlog(kxy),
     +                    epuntlog(kxy)
 15    CONTINUE
 55    FORMAT(1X,2F12.3,2G15.6) 
       cor_coef=sprod/(DSQRT(sumEQ*sumepunt))
       WRITE(16,56) cor_coef
 56    FORMAT(1X,F6.4)
       CLOSE(16) 
      
       WRITE(6,21) EQ_mean,ep_mean,cor_coef
 21    FORMAT(/'VALOR MIG DELS LOGARITMES DEL STRAIN RATE:'/
     +        4X,'- SEISMIC, APARTIR DELS TERRATREMOLS:',F10.4/
     +        4X,'- OBTINGUT APARTIR DE uhuru:         ',F10.4//
     +        1X,'COEFICIENT DE CORRELACIO:',F6.4)
 
       RETURN
       END
       
C ****************************************************************
C ****************************************************************
C            SUBROUTINE TEMPEFOND 
C  TEMPERATURA EN FOND·RIA EN EL PUNT (ixt,iyt) 

       SUBROUTINE TEMPEFOND (m, n, NELZ, nnsd, nzsd, Dx, Dy, Dz, 
     +          TEMPE, scrust)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TEMPE(nnsd,0:nzsd),scrust(nnsd) 

      PRINT*,'       SUBROUTINE TEMPEFOND'

       NMIG=n/2
       NQUART=n/4
       NVUIT=n/8
       MMIG=m/2
       MQUART=m/4
       MVUIT=m/8
       MMM=(M-2*MVUIT)/4
       ix1=48
       iy1=44
       ix2=49
       iy2=44
       ix3=48
       iy3=43
       ix4=49
       iy4=43
       ix5=48
       iy5=19   
       ix6=10
       iy6=10
C       ix1=MVUIT
C       iy1=NMIG
C       ix2=MMIG-MMM
C       iy2=NMIG
C       ix3=MMIG
C       iy3=NMIG
C       ix4=MMIG+MMM
C       iy4=NMIG
C       ix5=M-MVUIT
C       iy5=NMIG
C       ix6=MMIG+MQUART
C       iy6=NMIG
      
       kxy1=ix1+1+iy1*(n+1)   
       kxy2=ix2+1+iy2*(n+1)
       kxy3=ix3+1+iy3*(n+1)
       kxy4=ix4+1+iy4*(n+1)
       kxy5=ix5+1+iy5*(n+1)
       kxy6=ix6+1+iy6*(n+1)
      OPEN(11,FILE='punts.xy')
           WRITE(11,*) ix1*Dx,iy1*Dy
           WRITE(11,*) ix2*Dx,iy2*Dy
           WRITE(11,*) ix3*Dx,iy3*Dy
           WRITE(11,*) ix4*Dx,iy4*Dy
           WRITE(11,*) ix5*Dx,iy5*Dy
           WRITE(11,*) ix6*Dx,iy6*Dy
      CLOSE(11)
      OPEN(12,FILE='Temp1_z_ter')
      OPEN(11,FILE='TEMPE.z')
         DO 10 iz=0,NELZ 
            depth=(iz*Dz)/1.D3
            WRITE(11,11) depth,TEMPE(kxy1,iz), TEMPE(kxy2,iz),
     +     TEMPE(kxy3,iz),TEMPE(kxy4,iz),TEMPE(kxy5,iz),TEMPE(kxy6,iz)
           WRITE(12,22) TEMPE(kxy1,iz),depth  
 10      CONTINUE  
       PRINT*,' Gruix cortical de la columna: Temp1_z_ter',
     +            scrust(kxy1),' m'
 11   FORMAT(1X,7F11.3) 
 22   FORMAT(1X,2F11.4)
      CLOSE(11)
      CLOSE(12)
      RETURN
      END
C *********************************************************************
C     ORDENA EL VECTOR 'VIN' DEL VALOR MAJOR AL MENOR
C        VIN NOMES AMB VALÑORS POSITIUS

        SUBROUTINE ORDRE_MAG (VIN, NDIM, VOUT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION VIN(NDIM),VINTMP(NDIM),VOUT(NDIM)
      
       DO 5 KN=1,NDIM
            VINTMP(KN)=VIN(KN)
 5     CONTINUE          

       DO 9 IOUT=1,NDIM
           VMAX=0.D0
           DO 19 KN=1,NDIM
                 IF(VINTMP(KN).GE.VMAX) THEN
                        VMAX=VINTMP(KN)
                        KIT=KN
                 ENDIF  
 19         CONTINUE                
            VOUT(IOUT)= VMAX
            VINTMP(KIT)=0.D0
 9     CONTINUE         

      RETURN
      END
      
CC ****************************************************************
CC              vel_EUAF
CC  CALCULA LA VELOCITAT MITJA PER SOBRE I PER SOTA DE LA FALLA
CC  LA DIFERENCIA ENTRE LA VELOCITAT PER SOBRE I PER SOTA DELS PUNTS
CC    ON HI HAN DADES DEL NUVEL-1

       SUBROUTINE vel_EUAF (n, m, nn, AX, BY, u, v, FACVEL)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
       CHARACTER*70 DATAFIL, f_NUVEL
       PARAMETER(NFMAX=200,PI=3.14159265D0)
       DIMENSION u(nn),v(nn)
 
C       DATAFIL='/model/ivone/DADES/Golf/falla_Model3.xy'
C       f_NUVEL='/model/ivone/DADES/velocitats/NUVEL1_EUAF.xy'
      DATAFIL='/model/ivone/DADES/Atlantic_Octubre99/falla_Model3.xy'
      f_NUVEL='/model/ivone/DADES/Atlantic_Octubre99/NUVEL1_EUAF.xy'
       WRITE(6,1)DATAFIL
 
 1     FORMAT(1X,/'	    SUBROUTINE vel_EUAF:'/
     +     'CALCULA LA VELOCITAT MITJA PER SOBRE I PER SOTA DELS '/
     +     'PUNTS DEL FITXER:'/3X,A60/)

       Dx=AX/n
       Dy=BY/m 

CCC ------------------  Llegeix DATAFIL   --------------------------- 
        OPEN(16,FILE=DATAFIL)  
        OPEN(17,FILE='falla_vel.xy')  
         NPFALLA=0  
         DO 100 nf=1,NFMAX
             READ(16,*,IOSTAT=IOS) xfalla,yfalla
             IF (IOS < 0) THEN
                   GO TO 120
             ELSE IF (IOS > 0) THEN
               WRITE (*, 90)NPFALLA 
   90          FORMAT (/' ERROR: Encountered bad data after ',
     +                  I6,' successful READs.'/)
               STOP
             ENDIF
             NPFALLA=NPFALLA+1
             ix=xfalla/Dx
             iy=yfalla/Dy
             kxy=ix+1+iy*(n+1)

C        IF(xfalla.GT.1900.D3) GOTO 100
CC Velocitat mitja per sobre i per sota la falla:
             usup=0.D0
             vsup=0.D0
             uinf=0.D0
             vinf=0.D0
             nsup=0.D0
             ninf=0.D0
             ysupp=0.D0
             yinfp=0.D0
             DO 25 iyp=iy-10,m-1
                   kxyp=ix+1+iyp*(n+1)
             	   IF(iyp.GT.iy) THEN
                	usup=usup+u(kxyp)
                	vsup=vsup+v(kxyp)
                	ysupp=ysupp+iyp*Dy
                	nsup=nsup+1
                   ENDIF	
             	   IF(iyp.LT.iy) THEN
                	uinf=uinf+u(kxyp)
                	vinf=vinf+v(kxyp)
                	yinfp=yinfp+iyp*Dy
                	ninf=ninf+1
                   ENDIF
 25          CONTINUE
              usup=usup/nsup 
              vsup=vsup/nsup
              ysupp=ysupp/nsup
              uinf=uinf/ninf 
              vinf=vinf/ninf  
              yinfp=yinfp /ninf 	
                WRITE(17,77) "Sup",xfalla/1.D3,ysupp/1.D3,
     +                       usup*FACVEL,vsup*FACVEL
                WRITE(17,77) "Mig",xfalla/1.D3,yfalla/1.D3,
     +                      u(kxy)*FACVEL,v(kxy)*FACVEL
                WRITE(17,77) "Inf",xfalla/1.D3,yinfp/1.D3,
     +                       uinf*FACVEL,vinf*FACVEL

  100    CONTINUE
         WRITE(6,110)NFMAX,NPFALLA
  110    FORMAT(/' Dataset contained more than ',
     &          I5,' data; first ',I5,' were read.'/)
  120    CONTINUE
         WRITE(6,130)NPFALLA
  130    FORMAT(3X,' Reading of points fault completed:',I7,
     +             ' data points.'///)
        CLOSE(16)
        CLOSE(17)
 77    FORMAT(1X,A3,2F10.3,2F11.4)

CCC ------------------  Llegeix NUVEL-1   --------------------------- 
        OPEN(16,FILE=f_NUVEL)  
        OPEN(17,FILE='vel_pNUVEL.xy')  
        NPNUVEL=0
        Domega=0.D0
         DO 200 nf=1,NFMAX
             READ(16,*,IOSTAT=IOS) xn,yn,teta,vnuvel
             IF (IOS < 0) THEN
                   GO TO 220
             ELSE IF (IOS > 0) THEN
               WRITE (*, 90)NPFALLA 
               STOP
             ENDIF
             NPNUVEL=NPNUVEL+1
             ix1=INT(xn*1.D3/Dx)
             iy1=INT(yn*1.D3/Dy)
             u1=0.D0
             v1=0.D0
             u2=0.D0
             v2=0.D0
             ni=0
             DO 225 iy=iy1-10,iy1
                 kxy=ix1+1+iy*(n+1)
                 u1=u1+u(kxy)
		 v1=v1+v(kxy)
                 u2=u2+u(kxy+1)
		 v2=v2+v(kxy+1)
		 ni=ni+1
 225         CONTINUE
             u3=0.D0
             v3=0.D0
             u4=0.D0
             v4=0.D0
             ns=0
             DO 230 iy=iy1+1,iy1+10
                 kxy=ix1+1+iy*(n+1)
                 u3=u3+u(kxy)
		 v3=v3+v(kxy)
                 u4=u4+u(kxy+1)
		 v4=v4+v(kxy+1)
		 ns=ns+1
 230         CONTINUE
 	     u1=u1/ni
 	     v1=v1/ni
 	     u2=u2/ni
 	     v2=v2/ni
             u3=u3/ns
             v3=v3/ns
             u4=u4/ns
             v4=v4/ns
             x1=ix1*Dx
             y1=iy1*Dy
             d1=DSQRT((xn-x1)*(xn-x1)+(yn-y1)*(yn-y1)) 
             d2=DSQRT((xn-x1-Dx)*(xn-x1-Dx)+(yn-y1-Dy)*(yn-y1-Dy)) 
             d3=DSQRT((xn-x1)*(xn-x1)+(yn-y1-Dy)*(yn-y1-Dy)) 
             d4=DSQRT((xn-x1-Dx)*(xn-x1-Dx)+(yn-y1-Dy)*(yn-y1-Dy))
             uinf=(d2*u1+d1*u2)/(d1+d2)
             vinf=(d2*v1+d1*v2)/(d1+d2)
             usup=(d4*u3+d4*u3)/(d3+d4)
             vsup=(d4*v3+d4*v3)/(d3+d4)
             upN=usup-uinf
             vpN=vsup-vinf
             WRITE(17,277) xn,yn,
     +               usup*FACVEL,vsup*FACVEL,uinf*FACVEL,vinf*FACVEL
     		urel=uinf-usup
     		vrel=vinf-vsup
                tetamodel=DATAN(vrel/urel)
                tetamodel=(tetamodel*180.D0)/PI
                   if(urel.lt.0.D0) tetamodel=tetamodel+180.D0
     		Domega=Domega+ABS(teta-tetamodel)
  200    CONTINUE
         WRITE(6,110)NFMAX,NPNUVEL
  220    CONTINUE
	Domega=Domega/NPNUVEL 
	PRINT*,'file read:  ',f_NUVEL        	   
         WRITE(6,260) NPNUVEL,Domega
  260    FORMAT(3X,'number of data points :',I7,/
     +	'DESVIACIO MITJA DE LES VELOCITATS (Domega):',F10.3,' GRAUS'//)
        CLOSE(16)
        CLOSE(17)
 277    FORMAT(1X,2F11.3,4F12.6)
  
       RETURN
       END  
CCC ****************************************************************

