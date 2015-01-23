CC    SUBRUTINES  Velocitat_Slip.f
CC __________ Condicions de Contorn per el camp de velocitats _________
CC
CC        OJO EN QUIN PUNT COMENCA EN CADA ZONA (0,n),(0,m)
CC
CC *********************************************************************
CC                  SUBROUTINE fixCCSud
CC
CC              FIXO LES CONDICIONS DE CONTORN A iy=0
CC                iy=0   -  C.C. BASE - SUD   (ix=1,n-1)
CC *********************************************************************
CC POSICIO DEL POL DE ROTACIO (degrees) : (dlonpol,dlatpol)
CC POSICIO DEL POL DE ROTACIO (m): (xpol,ypol)
CC   Passo de [degrees] a [m] segons la projeccio Mercator, considerant
CC   l'origen (x,y)=(0,0) a la latitut-longitut (dlon0,dlat0)
CC VELOCITAT DE ROTACIO (deg/Ma): omegada
CC VELOCITAT DE ROTACIO (rad/s): omega=(omegada*PI)/(180.D6*FACTEMP)
CC DISTANCIA ENTRE EL POL I EL PUNT ON ES CALCULA V: RADI 
CC                          RADI=DSQRT((XCORD-xpol)**2+(YCORD-ypol)**2)

      SUBROUTINE fixCCSud (a,b,m,n,nincogn,nbanda,Dx,Dy,vr,
     +				dlonpol,dlatpol,omegada)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION a(nincogn,nbanda),b(nincogn)
      PARAMETER (PI=3.1415926535897932D0,FACVEL=3.1536D10,
     +       FACTEMP=3.1536D7,REarth=6.371D6,dlon0=2.D0,dlat0=43.D0) 

C  -----------  POSICIO (x,y) I VELOCITAT DEL POL  ------------------
           xpol=PI*REarth*(DCOS(dlatpol*PI/180.D0))*
     +           (dlonpol-dlon0)/180.D0
           ypol=PI*REarth*(dlatpol-dlat0)/180.D0
           omega=(omegada*PI)/(180.D6*FACTEMP) 
C  -----------------------------------------------------------------
         Li=2*(n+1)+3
         Ls=Li
         Ld=Li+1
         kp=Ld+2*(n+1)
         kn=Ld-2*(n+1)
	 ivora1=(503.D3)/Dx
	 ivora2=(1101.D3)/Dx
	 ix1=0
	 ix2=ivora1
	 ix3=n
	 ix4=n
C	 ix1=ivora1-3
C	 ix2=ivora1+3
C	 ix3=ivora2-3
C	 ix4=ivora2+3
	     TETA1=(68.D0*PI)/180.D0
	     TETA2=(90.D0*PI)/180.D0
	     ixm1=ix1
	     ixm2=ix2
	     ixm3=ivora2
	     ixf=n
	     vi=4.5D0/FACVEL
	     vf=5.2D0/FACVEL
CCC  ITBC= 1:VELOCITAT (INDENTER, u=0),	
CCC	   2:FREE SLIP o v constant (du/dy=0, v=cnst.)
CC	   3: du/dy=0, v=INDENTER 
CC	   4: FIXO u=cnst1, v=cnst2
CC	   5: VELOCITAT INDENTER: u i v
CC	   6: VELOCITAT INDENTER West, creix cap a l'Est
	ITBC=4
C ----------------  EQUATION 1 : --------------------------------
C   ****  C.C. DE VELOCITAT EN X, u  **** 
	IF(ITBC.EQ.1.OR.ITBC.EQ.4) THEN
             YCORD=0.D0
             u0=0.D0
             DO 3 ix=0,n
             	XCORD=ix*Dx
             	leq1=2*ix+1
             	a(leq1,Ld)=1.D0
             	b(leq1)=u0
C             	b(leq1)=-omega*(YCORD-ypol)
 3          CONTINUE
 	    GOTO 222
 	ENDIF
	IF(ITBC.EQ.5) THEN
	     ix1h=ix1
	     ix2h=ix2
	     ix3h=ivora2
	     ix4h=n-5
	     unul=0.D0
             u0=-0.634D-10
             DO 173 ix=0,ix1h
             	leq1=2*ix+1
              	a(leq1,Ld)=1.D0
              	b(leq1)=unul
 173          CONTINUE
             FACA=(PI/2.D0)/(ix1h-ix2h-1)
             FACB=-FACA*(ix2h+1)
             DO 175 ix=ix1h+1,ix2h
             	leq1=2*ix+1
              	a(leq1,Ld)=1.D0
              	b(leq1)=u0*(DCOS(FACA*ix+FACB))**2
 175          CONTINUE
             DO 179 ix=ix2h+1,ix3h
             	leq1=2*ix+1
              	a(leq1,Ld)=1.D0
              	b(leq1)=u0
 179          CONTINUE
             FACA=(PI/2.D0)/(ix4h+1-ix3h)
             FACB=-FACA*ix3h
             DO 182 ix=ix3h+1,ix4h
             	leq1=2*ix+1
              	a(leq1,Ld)=1.D0
              	b(leq1)=u0*(DCOS(FACA*ix+FACB))**2
 182          CONTINUE
             DO 184 ix=ix4h+1,n
             	leq1=2*ix+1
              	a(leq1,Ld)=1.D0
              	b(leq1)=unul
 184          CONTINUE
 	    GOTO 222
 	ENDIF
	IF(ITBC.EQ.6) THEN
             FACA=(PI/2.D0)/(ixm1-ixm2-1)
             FACB=-FACA*(ixm2+1)
 	     Pvel=(vf-vi)/(ixf-ixm2)
	     Pteta=(TETA2-TETA1)/(ixf-ixm3)
             DO 273 ix=0,ixm1
             	leq1=2*ix+1
              	a(leq1,Ld)=1.D0
              	b(leq1)=0.D0
 273         CONTINUE
             DO 275 ix=ixm1+1,ixm2
	     	velmod=vi*(DCOS(FACA*ix+FACB))**2
             	leq1=2*ix+1
              	a(leq1,Ld)=1.D0
              	b(leq1)=-velmod*DCOS(TETA1)
 275         CONTINUE
             DO 279 ix=ixm2+1,ixm3
	        velmod=Pvel*(ix-ixm2)+vi
             	leq1=2*ix+1
              	a(leq1,Ld)=1.D0
              	b(leq1)=-velmod*DCOS(TETA1)
 279         CONTINUE
             DO 282 ix=ixm3+1,n
	     	TETAx=Pteta*(ix-ixm3)+TETA1
		velmod=Pvel*(ix-ixm2)+vi
             	leq1=2*ix+1
              	a(leq1,Ld)=1.D0
              	b(leq1)=-velmod*DCOS(TETAx)
 282          CONTINUE
 	    GOTO 222
 	ENDIF

C   ****  C.C. du/dy=0  ****
	IF(ITBC.EQ.2.OR.ITBC.EQ.3) THEN
             DO 9 ix=0,n
             	leq1=2*ix+1
             	a(leq1,kp)=1.D0
             	a(leq1,Ld)=-1.D0
 9           CONTINUE
 	     GOTO 222
 	ENDIF
C   ****  C.C. D'ESFORS XY, TAUXY=0  ->  du/dy+dv/dx=0  ****
C       DIAGONAL=1.D0/Dy   NODIAGONAL=1.D0/(2.D0*Dx)
C         DO 12 ix=1,n-1
C             leq1=2*ix+1
C             a(leq1,kp)=1.D0
C             a(leq1,Ld)=-1.D0
C             a(leq1,Ld+3)=Dy/(2.D0*Dx)
C             a(leq1,Ld-1)=-Dy/(2.D0*Dx)
C12       CONTINUE
C ----------------  EQUATION 2 : --------------------------------
C   ****  C.C. DE VELOCITAT EN Y, v 
 222	CONTINUE
	IF(ITBC.EQ.2.OR.ITBC.EQ.4) THEN
C            YCORD=0.D0    
C             vconst=0.D0
             vconst=vr
             DO 63 ix=0,n
C	        vconst=0.951D-10+DBLE((3.D0*ix*0.634D-10)/n)
             	XCORD=ix*Dx
             	leq2=2*ix+2
             	a(leq2,Ld)=1.D0
             	b(leq2)=vconst
C             	b(leq2)=omega*(XCORD-xpol)
 63          CONTINUE
 	    GOTO 999
   	ENDIF
	IF(ITBC.EQ.6) THEN
             FACA=(PI/2.D0)/(ixm1-ixm2-1)
             FACB=-FACA*(ixm2+1)
 	     Pvel=(vf-vi)/(ixf-ixm2)
	     Pteta=(TETA2-TETA1)/(ixf-ixm3)
             DO 373 ix=0,ixm1
              	leq2=2*ix+2
              	a(leq2,Ld)=1.D0
              	b(leq2)=0.D0
 373         CONTINUE
             DO 375 ix=ixm1+1,ixm2
	     	velmod=vi*(DCOS(FACA*ix+FACB))**2
              	leq2=2*ix+2
              	a(leq2,Ld)=1.D0
              	b(leq2)=velmod*DSIN(TETA1)
 375         CONTINUE
             DO 379 ix=ixm2+1,ixm3
	        velmod=Pvel*(ix-ixm2)+vi
              	leq2=2*ix+2
              	a(leq2,Ld)=1.D0
              	b(leq2)=velmod*DSIN(TETA1)
 379         CONTINUE
             DO 382 ix=ixm3+1,n
	     	TETAx=Pteta*(ix-ixm3)+TETA1
		velmod=Pvel*(ix-ixm2)+vi
              	leq2=2*ix+2
              	a(leq2,Ld)=1.D0
              	b(leq2)=velmod*DSIN(TETAx)
 382          CONTINUE
 	    GOTO 999
 	ENDIF
             vnul=0.D0
             v0=vr
             DO 73 ix=0,ix1
              	leq2=2*ix+2
              	a(leq2,Ld)=1.D0
              	b(leq2)=vnul
 73          CONTINUE
             FACA=(PI/2.D0)/(ix1-ix2-1)
             FACB=-FACA*(ix2+1)
             DO 75 ix=ix1+1,ix2
              	leq2=2*ix+2
              	a(leq2,Ld)=1.D0
              	b(leq2)=v0*(DCOS(FACA*ix+FACB))**2
 75          CONTINUE
             DO 79 ix=ix2+1,ix3
              	leq2=2*ix+2
              	a(leq2,Ld)=1.D0
              	b(leq2)=v0
 79          CONTINUE
             FACA=(PI/2.D0)/(ix4+1-ix3)
             FACB=-FACA*ix3
             DO 82 ix=ix3+1,ix4
              	leq2=2*ix+2
              	a(leq2,Ld)=1.D0
              	b(leq2)=v0*(DCOS(FACA*ix+FACB))**2
 82          CONTINUE
             DO 84 ix=ix4+1,n
              	leq2=2*ix+2
              	a(leq2,Ld)=1.D0
              	b(leq2)=vnul
 84          CONTINUE
 	GOTO 999
CC ---
C  ****  C.C. D'ESFORS YY, TAUYY=0  ->  dv/dy=0  ****
C         DO 102 ix=1,n-1
C             leq2=2*ix+2
C             a(leq2,kp)=1.D0
C             a(leq2,Ld)=-1.D0
C102       CONTINUE

 999	CONTINUE
       RETURN
       END
C **********************************************************************
C                  SUBROUTINE fixCCNort
C
C              FIXO LES CONDICIONS DE CONTORN A iy=m
C               iy=m   -  C.C. TOP - NORT   (ix=1,n-1)
C **********************************************************************
C  TOTES LES VARIABLES JA ENTREN ADIMENSIONALS

       SUBROUTINE fixCCNort (a,b,m,n,nincogn,nbanda,Dx,Dy,vr)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION a(nincogn,nbanda),b(nincogn)

         Li=2*(n+1)+3
         Ls=Li
         Ld=Li+1
         kp=Ld+2*(n+1)
         kn=Ld-2*(n+1)
CCC  ITBC= 1:VELOCITAT,	2:FREE SLIP (du/dy=0, v=0)
	ITBC=1
CC ----------------  EQUATION 1 : --------------------------------
CC   ****  C.C. DE VELOCITAT EN X, u
	IF(ITBC.EQ.1) THEN
             u0=0.D0
             DO 3 ix=0,n
              	leq1=2*ix+1+2*m*(n+1)
              	a(leq1,Ld)=1.d0
             	b(leq1)=u0
 3           CONTINUE
	     GOTO 111 	
	ENDIF
CC   ****  C.C. du/dy=0
	IF(ITBC.EQ.2) THEN
             DO 9 ix=0,n
              	leq1=2*ix+1+2*m*(n+1)
              	a(leq1,Ld)=1.D0
              	a(leq1,kn)=-1.D0
 9           CONTINUE
	     GOTO 111 	
	ENDIF
CC   ****  C.C. D'ESFORS XY, TAUXY=0  ->  du/dy+dv/dx=0  ****
CC       DIAGONAL=1.D0/Dy   NODIAGONAL=1.D0/(2.D0*Dx)
C         DO 12 ix=1,n-1
C             leq1=2*ix+1+2*m*(n+1)
C             a(leq1,Ld)=1.D0
C             a(leq1,kn)=-1.D0
C             a(leq1,Ld+3)=Dy/(2.D0*Dx)
C             a(leq1,Ld-1)=-Dy/(2.D0*Dx)
C12       CONTINUE
C ----------------  EQUATION 2 : --------------------------------
C  ****  C.C. DE VELOCITAT EN Y, v
 111	CONTINUE
             v0=0.D0
             DO 75 ix=0,n
              	leq2=2*ix+2+2*m*(n+1)
                a(leq2,Ld)=1.d0
                b(leq2)=v0
 75          CONTINUE
C  ****  C.C. D'ESFORS YY, TAUYY=0  ->  dv/dy=0  ****
C         DO 102 ix=1,n-1
C             leq2=2*ix+2+2*m*(n+1)
C             a(leq2,Ld)=1.D0
C             a(leq2,kn)=-1.D0
C102       CONTINUE

       RETURN
       END
C **********************************************************************
C                  SUBROUTINE fixCCWest
C
C              FIXO LES CONDICIONS DE CONTORN A ix=0
C           ix=0   -  C.C. ESQUERRA - OEST   (iy=0,m)
C **********************************************************************
C  TOTES LES VARIABLES JA ENTREN ADIMENSIONALS

	SUBROUTINE fixCCWest (a,b,m,n,nincogn,nbanda,Dx,Dy,vr,
     +				dlonpol,dlatpol,omegada,vis,nn)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION a(nincogn,nbanda),b(nincogn),vis(nn) 
       PARAMETER (PI=3.1415926535897932D0,FACVEL=3.1536D10,
     +       FACTEMP=3.1536D7,REarth=6.371D6,dlon0=2.D0,dlat0=43.D0,
     +       iylim=21)

CC  -----------  POSICIO (x,y) I VELOCITAT DEL POL  ------------------
           xpol=PI*REarth*(DCOS(dlatpol*PI/180.D0))*
     +           (dlonpol-dlon0)/180.D0
           ypol=PI*REarth*(dlatpol-dlat0)/180.D0
           omega=(omegada*PI)/(180.D6*FACTEMP) 
CC  -----------------------------------------------------------------
         Li=2*(n+1)+3
         Ls=Li
         Ld=Li+1
         kp=Ld+2*(n+1)
         kn=Ld-2*(n+1)
         XCORD=0.D0
CCC  ITBC= 1:VELOCITAT,	2:FREE SLIP (u=0, dv/dx=0),  3: FREE
	ITBC=2
C ----------------  EQUATION 1 : --------------------------------
C  ****  C.C. DE VELOCITAT EN X, u 
CC *** Segons un pol de rotacio  
C        DO 3 iy=1,iylim
C             YCORD=iy*Dy
C             leq1=1+2*iy*(n+1)
C             a(leq1,Ld)=1.D0
C             b(leq1)=-omega*(YCORD-ypol)
C 3      CONTINUE
C  ****  u =  constant ***
	IF(ITBC.EQ.1.OR.ITBC.EQ.2) THEN
             uconst=0.D0
	     IF(ITBC.EQ.2) uconst=0.D0
             DO 3 iy=1,m-1
              	leq1=1+2*iy*(n+1)
              	a(leq1,Ld)=1.D0
              	b(leq1)=uconst
 3            CONTINUE
C 	      GOTO 111
	ENDIF
CC ****  C.C. STRAIN RATE xx, epuntxx=0  ->  du/dx=0  ****
CC ***     			   u(1,iy)-u(0,iy)=0  ****
	IF(ITBC.EQ.3) THEN
             DO 12 iy=1,m-1
        	leq1=1+2*iy*(n+1)
             	a(leq1,Ld)=1.D0
             	a(leq1,Ld+2)=-1.D0
 12          CONTINUE
	ENDIF

C ----------------  EQUACIO 2 : --------------------------------    
CC *** Segons el pol de rotacio
C         DO 73 iy=1,iylim
C             YCORD=iy*Dy
C             leq2=2+2*iy*(n+1)
C             a(leq2,Ld)=1.D0
C             b(leq2)=omega*(XCORD-xpol)
C 73        CONTINUE
C  ****  C.C. DE VELOCITAT EN Y, v
	IF(ITBC.EQ.1) THEN
             vconst=0.D0
             DO 65 iy=1,m-1
             	leq2=2+2*iy*(n+1)
                a(leq2,Ld)=1.D0
                b(leq2)=vconst
 65          CONTINUE
 	ENDIF
C  ****  C.C.  dv/dx=0
	IF(ITBC.EQ.2) THEN
             DO 87 iy=1,m-1
            	leq2=2+2*iy*(n+1)
               	a(leq2,Ld)=-1.d0
               	a(leq2,Ld+2)=1.D0
 87          CONTINUE
 	ENDIF
C  ****  C.C. TAUMXY=0
C         DO 89 iy=0,m
C               leq2=2+2*iy*(n+1)
C               a(leq2,kp-1)=1.D0/(2.D0*Dy)
C               a(leq2,kn-1)=-1.D0/(2.D0*Dy)
C               a(leq2,Ld+2)=1.D0/Dx
C               a(leq2,Ld)=-1.D0/Dx
C89      CONTINUE
CC ****  C.C. STRAIN RATE XY, epuntxy=0  ->  du/dy+dv/dx=0  ****
CC **** (u(0,iy+1)-u(0,iy-1))/(2*Dy)+(v(1,iy)-v(0,iy))/Dx=0 ****
CC ****     DIAGONAL=1.D0/Dx   NODIAGONAL=1.D0/(2.D0*Dy)    ****
	IF(ITBC.EQ.3) THEN
             DO 97 iy=1,m-1
             	leq2=2+2*iy*(n+1)
             	a(leq2,Ld)=1.D0
             	a(leq2,Ld+2)=-1.D0
             	a(leq2,kp-1)=-Dx/(2.D0*Dy)
             	a(leq2,kn-1)=Dx/(2.D0*Dy)
 97          CONTINUE
 	ENDIF

       RETURN
       END
C **********************************************************************
C                  SUBROUTINE fixCCEst
C
C              FIXO LES CONDICIONS DE CONTORN A ix=n
C           ix=n   -  C.C. ESQUERRA - EST   (iy=0,m)
C **********************************************************************
C  TOTES LES VARIABLES JA ENTREN ADIMENSIONALS

	SUBROUTINE fixCCEst (a,b,m,n,nincogn,nbanda,Dx,Dy,vr,
     +				dlonpol,dlatpol,omegada,IFAULT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION a(nincogn,nbanda),b(nincogn)   
       PARAMETER (PI=3.1415926535897932D0,FACVEL=3.1536D10,
     +       FACTEMP=3.1536D7,REarth=6.371D6,dlon0=2.D0,dlat0=43.D0)

        DATA  limcc /20/

CC  -----------  POSICIO (x,y) I VELOCITAT DEL POL  ------------------
           xpol=PI*REarth*(DCOS(dlatpol*PI/180.D0))*
     +           (dlonpol-dlon0)/180.D0
           ypol=PI*REarth*(dlatpol-dlat0)/180.D0
           omega=(omegada*PI)/(180.D6*FACTEMP) 
CC  -----------------------------------------------------------------
         Li=2*(n+1)+3
         Ls=Li
         Ld=Li+1
         kp=Ld+2*(n+1)
         kn=Ld-2*(n+1) 
         XCORD=n*Dx
CCC  BOUNDARY CONDITION:
CCC		ITBC= 1: VELOCITAT,	ITBC=2: FREE SLIP (u=0, dv/dx=0),  
CCC		ITBC=3: FREE (exx=du/dx=0, exy=0)
	ITBC=2
CC ----------------  EQUATION 1 : -------------------------------- 
CC *** Segons un pol de rotacio  
CC  ****  C.C. DE VELOCITAT EN X, u
	IF(ITBC.EQ.1.OR.ITBC.EQ.2) THEN
             uconst=0.D0
	     IF(ITBC.EQ.2) uconst=0.D0
             DO 5 iy=1,m-1
              	leq1=2*n+1+2*iy*(n+1)
              	a(leq1,Ld)=1.D0
              	b(leq1)=uconst
 5           CONTINUE
 	     GOTO 222
 	ENDIF
CC  ****  C.C. STRAIN RATE XX = 0  (TAUXX=TAUZZ) ->  du/dx=0  ****
	IF(ITBC.EQ.3) THEN
             DO 12 iy=1,m-1
       		leq1=2*n+1+2*iy*(n+1)
             	a(leq1,Ld)=1.D0
             	a(leq1,Ld-2)=-1.D0
 12          CONTINUE
 	     GOTO 222
 	ENDIF
CC ----------------  EQUATION 2 : -------------------------------- 
CC ****  C.C. DE VELOCITAT EN Y, v
 222    CONTINUE
	IF(ITBC.EQ.1) THEN
             v0=0.D0
             DO 75 iy=1,m-1
             	leq2=2*n+2+2*iy*(n+1)
               	a(leq2,Ld)=1.D0
               	b(leq2)=v0
 75          CONTINUE
 	     GOTO 999
 	ENDIF
CC  ****  C.C.  dv/dx=0
	IF(ITBC.EQ.2) THEN
             DO 79 iy=1,m-1
               	leq2=2*n+2+2*iy*(n+1)
               	a(leq2,Ld)=1.D0
               	a(leq2,Ld-2)=-1.D0
 79          CONTINUE
 	     GOTO 999
 	ENDIF
CC  ****  C.C. STRAIN RATE XY = 0  (TAUXY=0)->  du/dy+dv/dx=0  ****
CC       DIAGONAL=1.D0/Dx   NODIAGONAL=1.D0/(2.D0*Dy)
CC **    iy=0
C          leq2=2*n+2
C          a(leq2,Ld)=1.D0
C          a(leq2,Ld-2)=-1.D0
C          a(leq2,kp-1)=Dx/Dy
C          a(leq2,Ld-1)=-Dx/Dy
	IF(ITBC.EQ.3) THEN
             DO 87 iy=1,m-1
       		leq2=2*n+2+2*iy*(n+1)
             	a(leq2,Ld)=1.D0
             	a(leq2,Ld-2)=-1.D0
             	a(leq2,kp-1)=Dx/(2.D0*Dy)
             	a(leq2,kn-1)=-Dx/(2.D0*Dy)
 87          CONTINUE
 	     GOTO 999
 	ENDIF
CC **    iy=m
C          leq2=2*n+2+2*m*(n+1)
C          a(leq2,Ld)=1.D0
C          b(leq2)=0.D0


 999	CONTINUE
       RETURN
       END
