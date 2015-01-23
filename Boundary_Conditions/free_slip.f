CC    SUBRUTINES  free_slip.f
CC __________ Condicions de Contorn per el camp de velocitats _________
CC        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
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
     +       FACTEMP=3.1536D7,REarth=6.371D6,dlon0=-25.D0,dlat0=30.D0) 

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
         nquart=n/4
         nvuit=n/8
          ix1=nquart
          ix2=nquart+nvuit
          ix3=n-nquart-nvuit-1
          ix4=n-nquart-1
C ----------------  Equaci¢ 1 : --------------------------------
C   ****  C.C. DE VELOCITAT EN X, u  **** 
C        YCORD=0.D0
C        DO 3 ix=0,n
C             XCORD=ix*Dx
C             leq1=2*ix+1
C             a(leq1,Ld)=1.D0
C             b(leq1)=-omega*(YCORD-ypol)
C 3      CONTINUE
C   ****  C.C. du/dy=0  ****
         DO 9 ix=0,n
             leq1=2*ix+1
             a(leq1,kp)=1.D0
             a(leq1,Ld)=-1.D0
9        CONTINUE
C   ****  C.C. D'ESFORS XY, TAUXY=0  ->  du/dy+dv/dx=0  ****
C       DIAGONAL=1.D0/Dy   NODIAGONAL=1.D0/(2.D0*Dx)
C         DO 12 ix=1,n-1
C             leq1=2*ix+1
C             a(leq1,kp)=1.D0
C             a(leq1,Ld)=-1.D0
C             a(leq1,Ld+3)=Dy/(2.D0*Dx)
C             a(leq1,Ld-1)=-Dy/(2.D0*Dx)
C12       CONTINUE
C ----------------  Equaci¢ 2 : --------------------------------
C   ****  C.C. DE VELOCITAT EN Y, v 
C         YCORD=0.D0    
C         DO 73 ix=0,n
C             XCORD=ix*Dx
C             leq2=2*ix+2
C             a(leq2,Ld)=1.D0
C             b(leq2)=omega*(XCORD-xpol)
C 73        CONTINUE
C  ****   v = constant ***
          vconst=0.D0
          DO 93 ix=0,n
             leq2=2*ix+2
             a(leq2,Ld)=1.D0
             b(leq2)=vconst
 93        CONTINUE
C  ****  C.C. D'ESFORS YY, TAUYY=0  ->  dv/dy=0  ****
C         DO 102 ix=1,n-1
C             leq2=2*ix+2
C             a(leq2,kp)=1.D0
C             a(leq2,Ld)=-1.D0
C102       CONTINUE

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
CC ----------------  Equaci¢ 1 : --------------------------------
CC   ****  C.C. DE VELOCITAT EN X, u
C          u0=0.D0
C          DO 3 ix=0,n
C              leq1=2*ix+1+2*m*(n+1)
C              a(leq1,Ld)=1.d0
C              b(leq1)=u0
C 3        CONTINUE
CC   ****  C.C. du/dy=0
          DO 9 ix=0,n
              leq1=2*ix+1+2*m*(n+1)
              a(leq1,Ld)=1.D0
              a(leq1,kn)=-1.D0
 9        CONTINUE
CC   ****  C.C. D'ESFORS XY, TAUXY=0  ->  du/dy+dv/dx=0  ****
CC       DIAGONAL=1.D0/Dy   NODIAGONAL=1.D0/(2.D0*Dx)
C         DO 12 ix=1,n-1
C             leq1=2*ix+1+2*m*(n+1)
C             a(leq1,Ld)=1.D0
C             a(leq1,kn)=-1.D0
C             a(leq1,Ld+3)=Dy/(2.D0*Dx)
C             a(leq1,Ld-1)=-Dy/(2.D0*Dx)
C12       CONTINUE
C ----------------  Equaci¢ 2 : --------------------------------
C  ****  C.C. DE VELOCITAT EN Y, v
          v0=0.D0
          DO 75 ix=0,n
                leq2=2*ix+2+2*m*(n+1)
                a(leq2,Ld)=1.d0
                b(leq2)=v0
 75        CONTINUE
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
     +       FACTEMP=3.1536D7,REarth=6.371D6,dlon0=-25.D0,dlat0=30.D0,
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
         mmig=m/2
         mters=m/3
         mnou=m/9
         mquart=m/4
         mvuit=m/8
         XCORD=0.D0
C ----------------  Equaci¢ 1 : --------------------------------
C  ****  C.C. DE VELOCITAT EN X, u 
CC *** Segons un pol de rotacio  
C        DO 3 iy=1,iylim
C             YCORD=iy*Dy
C             leq1=1+2*iy*(n+1)
C             a(leq1,Ld)=1.D0
C             b(leq1)=-omega*(YCORD-ypol)
C 3      CONTINUE
C  ****  u =  constant ***
          uconst=0.D0
          DO 3 iy=1,m-1
              leq1=1+2*iy*(n+1)
              a(leq1,Ld)=1.D0
              b(leq1)=uconst
 3        CONTINUE
CC ****  C.C. STRAIN RATE xx, epuntxx=0  ->  du/dx=0  ****
CC ***     			   u(1,iy)-u(0,iy)=0  ****
C          DO 12 iy=iylim+1,m-1
C          DO 12 iy=1,m-1
C             leq1=1+2*iy*(n+1)
C             a(leq1,Ld)=1.D0
C             a(leq1,Ld+2)=-1.D0
C 12       CONTINUE
CC **** C.C. STRAIN RATE xx = STRAIN RATE zz, epuntxx=epuntzz
CC   Condicio d'incompresssibilitat: epuntzz=-(epuntxx+epuntyy)
C          DO 23 iy=1,m-1
C             leq1=1+2*iy*(n+1)
C             a(leq1,Ld)=1.D0
C             a(leq1,Ld+2)=-1.D0
C  ACABAR-HO
C 23       CONTINUE

C ----------------  Equaci¢ 2 : --------------------------------    
CC *** Segons el pol de rotacio
C         DO 73 iy=1,iylim
C             YCORD=iy*Dy
C             leq2=2+2*iy*(n+1)
C             a(leq2,Ld)=1.D0
C             b(leq2)=omega*(XCORD-xpol)
C 73        CONTINUE
C  ****  C.C. DE VELOCITAT EN Y, v
C          vconst=0.D0
C          vpetit=vconst/1.6D0
C          DO 65 iy=0,m
C                leq2=2+2*iy*(n+1)
C                a(leq2,Ld)=1.D0
C                b(leq2)=vconst
C 65       CONTINUE
C  ****  C.C.  dv/dx=0
          DO 87 iy=1,m-1
               leq2=2+2*iy*(n+1)
               a(leq2,Ld)=-1.d0
               a(leq2,Ld+2)=1.D0
 87       CONTINUE
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
C         DO 87 iy=iylim+1,m-1
C         DO 87 iy=1,m-1
C             leq2=2+2*iy*(n+1)
C             a(leq2,Ld)=1.D0
C             a(leq2,Ld+2)=-1.D0
C             a(leq2,kp-1)=-Dx/(2.D0*Dy)
C             a(leq2,kn-1)=Dx/(2.D0*Dy)
C 87       CONTINUE

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
     +       FACTEMP=3.1536D7,REarth=6.371D6,dlon0=-25.D0,dlat0=30.D0,
     +       CClimlat=36.0D0,iylim1=15,iylim2=25)

        DATA  limcc /20/

CC  -----------  POSICIO (x,y) I VELOCITAT DEL POL  ------------------
           xpol=PI*REarth*(DCOS(dlatpol*PI/180.D0))*
     +           (dlonpol-dlon0)/180.D0
           ypol=PI*REarth*(dlatpol-dlat0)/180.D0
           omega=(omegada*PI)/(180.D6*FACTEMP) 
CC  -----------------------------------------------------------------
C        CCrad=(PI*CClimlat)/180.D0
C        YCClim=(REarth*log(tan(PI/4.D0+CCrad/2.D0)))-ycord0  
C        limcc=YCClim/Dy
         Li=2*(n+1)+3
         Ls=Li
         Ld=Li+1
         kp=Ld+2*(n+1)
         kn=Ld-2*(n+1) 
         XCORD=n*Dx
CCC  Falla amb dues Branques -> GOTO 991
CC ----------------  Equaci¢ 1 : -------------------------------- 
CC *** Segons un pol de rotacio  
CC  ****  C.C. DE VELOCITAT EN X, u
          u0=0.D0
          DO 5 iy=1,m-1
              leq1=2*n+1+2*iy*(n+1)
              a(leq1,Ld)=1.D0
              b(leq1)=u0
 5        CONTINUE
CC  ****  C.C. STRAIN RATE XX = 0  (TAUXX=TAUZZ) ->  du/dx=0  ****
C         DO 12 iy=1,m-1
C             leq1=2*n+1+2*iy*(n+1)
C             a(leq1,Ld)=1.D0
C             a(leq1,Ld-2)=-1.D0
C 12       CONTINUE
CC      iy=m
C         leq1=2*n+1+2*m*(n+1)
C         a(leq1,Ld)=1.D0
C         b(leq1)=0.D0

CC ----------------  Equaci¢ 2 : -------------------------------- 
CC *** Segons el pol de rotacio
CC ****  C.C. DE VELOCITAT EN Y, v
C          v0=0.D0
C          DO 75 iy=1,m-1
C               leq2=2*n+2+2*iy*(n+1)
C               a(leq2,Ld)=1.D0
C               b(leq2)=v0
C 75       CONTINUE
CC  ****  C.C.  dv/dx=0
          DO 79 iy=1,m-1
               leq2=2*n+2+2*iy*(n+1)
               a(leq2,Ld)=1.D0
               a(leq2,Ld-2)=-1.D0
 79        CONTINUE
CC  ****  C.C. STRAIN RATE XY = 0  (TAUXY=0)->  du/dy+dv/dx=0  ****
CC       DIAGONAL=1.D0/Dx   NODIAGONAL=1.D0/(2.D0*Dy)
CC **    iy=0
C          leq2=2*n+2
C          a(leq2,Ld)=1.D0
C          a(leq2,Ld-2)=-1.D0
C          a(leq2,kp-1)=Dx/Dy
C          a(leq2,Ld-1)=-Dx/Dy
C         DO 87 iy=limcc+1,m-1
C             leq2=2*n+2+2*iy*(n+1)
C             a(leq2,Ld)=1.D0
C             a(leq2,Ld-2)=-1.D0
C             a(leq2,kp-1)=Dx/(2.D0*Dy)
C             a(leq2,kn-1)=-Dx/(2.D0*Dy)
C 87       CONTINUE
CC **    iy=m
C          leq2=2*n+2+2*m*(n+1)
C          a(leq2,Ld)=1.D0
C          b(leq2)=0.D0

 999   CONTINUE
       RETURN
       END
