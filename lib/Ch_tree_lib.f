
!	Compressio : negatiu.	Extensio : positiu.
!  Bexten, Bcompres : Brittle failure coefs 
!  Quc, Qlc, Qml : Power flow activation energies
!  QDornml : Dorn law activation energy
!  Auc, Alc, Aml : Pre-exponencial constant in upper,    
!                        lower crust and mantle in (s-1 Pa-n) 
!  en : Power law exponent
!  epsDorn : Dorn strain rate reference. Pre-exponencial constant
!              for olivine Dorn law creep.
!  sigmaD : Dorn stress reference
!  crustup, crustlow : Base upper and lower crust
!  strainrate : Strain rate (s-1)
!  stresslim :  mecanic stress limit (Pa) (-> espesor mecanic)
!  RGAS : Constant dels gasos (J/mol/K)
!  stressPD : stress de transicio Power law - Dorn law

!  iswitch = 1	=> keep in a file the stress envelope and the strength

      SUBROUTINE Ch_tree_lib (NELZ, Dz, T, strainrate, crustlow,
     +			       ZLITOS, Flitcomp, Flitext, ZMECANIC, 
     +			       IRheology_type, Qarray, enarray, Aarray,
     +			       TIPPAR, iswitch)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(0:NELZ),Qarray(3),enarray(3),Aarray(3)
      CHARACTER*40 TIPPAR
      CHARACTER*20 TIPPAR1
      PARAMETER (RGAS=8.314D0,Bexten=16.D3,Bcompres=40.D3,
     +           epsDorn=5.7D11,stressPD=200.D6,stresslim=20.D6,
     +		grav=9.8D0,ffmuf=0.85D0,yieldmax=1.3D9,
     +		stress_max=200.D6)		!! Ranalli et al. (J.Geodynamics, 2006)
     						!! Renshaw & Schulson (JGR, 2004)
      DATA sigmaD/8.5D9/ QDornml/100.D3/ 
CC ---------------------------------------------------------------------
CC     litosfera oceanica
C       DATA lithtype/2/
CC ***[1] Lynch&Morgan, 1987, escorça oceanica ******
C          PARAMETER (Quc=520.D3,Qlc=520.D3,Qml=523.D3,
C     +               enuc=3.D0,enlc=3.D0,enml=3.D0) 
C          DATA AucMPa/7.D4/ AlcMPa/7.D4/ AmlMPa/1.D3/ 
C          TIPPAR1='Lynch&Morgan 1987'
CC ---------------------------------------------------------------------
CC     litosfera continental
       DATA lithtype/1/
CC ---------------------------------------------------------------------
CC *********  A=[MPa**-n s-1]  DIFERENTS AUTORS ****************
!! ***[1] Lynch&Morgan, 1987 ******
      IF(IRheology_type.EQ.1) THEN
           Qarray=(/138.D3, 251.D3, 523.D3/)
	   enarray=(/3.0D0, 3.0D0, 3.0D0/)
	   Aarray=(/2.5D-8, 3.2D-3, 1.D3/)
           TIPPAR1='Lynch&Morgan 1987'
      END IF
!! ***[2] Braun&Beaumont, 1989 ******
      IF(IRheology_type.EQ.2) THEN
           Qarray=(/151.D3, 239.D3, 498.D3/)
	   enarray=(/1.8D0, 3.2D0, 4.5D0/)
	   Aarray=(/1.D-2, 3.D-2, 1.9D5/)
           TIPPAR1='Braun&Beaumont 1989'
      END IF
!! ***[3] Fadaie&Ranalli, 1990  ******
      IF(IRheology_type.EQ.3) THEN
           Qarray=(/219.D3, 268.D3, 535.D3/)
	   enarray=(/2.4D0, 3.3D0, 3.6D0/)
	   Aarray=(/1.3D-3, 3.2D-3, 3.2D4/)
           TIPPAR1='Fadaie&Ranalli 1990'
      END IF
!! ***[4] Buck, 1991 ******
      IF(IRheology_type.EQ.4) THEN
           Qarray=(/149.D3, 238.D3, 500.D3/)
	   enarray=(/2.9D0, 3.2D0, 3.D0/)
	   Aarray=(/1.7D-7, 8.9D-4, 1.D3/)
           TIPPAR1='Buck 1991'
      END IF
!! ***[5] Liu&Furlong, 1993  ******
      IF(IRheology_type.EQ.5) THEN
           Qarray=(/123.D3, 260.D3, 420.D3/)
	   enarray=(/3.D0, 3.4D0, 3.D0/)
	   Aarray=(/1.6D-9, 2.D-4, 1.9D3/)
           TIPPAR1='Liu&Furlong 1993'
      END IF
!! ***[6] Lowe&Ranalli, 1993  ******
      IF(IRheology_type.EQ.6) THEN
           Qarray=(/144.D3, 238.D3, 535.D3/)
	   enarray=(/3.2D0, 3.2D0, 3.5D0/)
	   Aarray=(/1.3D-9, 3.3D-4, 1.4D5/)
           TIPPAR1='Lowe&Ranalli 1993'
           sigmaD=15.D9
      END IF
!! ***[7] Mareschal, 1994  ******
      IF(IRheology_type.EQ.7) THEN
           Qarray=(/141.D3, 445.D3, 527.D3/)
	   enarray=(/1.9D0, 4.2D0, 3.D0/)
	   Aarray=(/2.D-4, 1.4D-4, 4.3D2/)
           TIPPAR1='Mareschal 1994'
           sigmaD=12.D9
      END IF
!! ***[8] Bassi, 1995 (dry) ******
      IF(IRheology_type.EQ.8) THEN
           Qarray=(/185.D3, 235.D3, 535.D3/)
	   enarray=(/2.8D0, 3.9D0, 3.6D0/)
	   Aarray=(/3.4D-6, 2.3D-6, 2.9D4/)
           TIPPAR1='Bassi 1995'
      END IF
!! ***[9] Bassi, 1991 (wet) ******
      IF(IRheology_type.EQ.9) THEN
           Qarray=(/150.D3, 238.D3, 445.D3/)
	   enarray=(/1.8D0, 3.2D0, 3.4D0/)
	   Aarray=(/2.9D-3, 3.3D-4, 1.4D4/)
           TIPPAR1='Bassi 1991 (wet)'
      END IF
!! ***[99] Rheilogical parameters from the input file ******
      IF(IRheology_type.EQ.99) THEN
           TIPPAR1='Input file'
      END IF

      AucMPa=Aarray(1)
      AlcMPa=Aarray(2)
      AmlMPa=Aarray(3)
      enuc=enarray(1)
      enlc=enarray(2)
      enml=enarray(3)
!   Conversio al S.I.
      Auc=AucMPa*(10.D0**(-6.D0*enuc))
      Alc=AlcMPa*(10.D0**(-6.D0*enlc))  
      Aml=AmlMPa*(10.D0**(-6.D0*enml))  
      
!! *** IRheology_type>=100 : PARAMETRES DEL 'shells'  A shells=[MPa s**1/n] ******
! Second invariant = 3(strain_rate**2)	==> A=(2*sqrt(3))**(n-1) * Ashells**(-n)
! Second invariant = strain_rate**2	==> A= 2**(n-1) * Ashells**(-n)
! On the crust (Bird & Kong, 1994) : log(Ashells)=12.2-2.8D-5*Q
      IF(IRheology_type>=100) THEN
      	   IF(IRheology_type==100) THEN
           	TIPPAR1=' Peter+Kirby'
           	Qarray=(/100.D3, 100.D3, 456.D3/)	! JGR Mediterrani
           ENDIF
      	   IF(IRheology_type==110) THEN
           	TIPPAR1='Peter+Kirby, soft'
                Qarray=(/100.D3, 100.D3, 420.D3/)	! Q minimes (soft)
           ENDIF
      	   IF(IRheology_type==120) THEN
           	TIPPAR1='Peter+Kirby, hard'
                Qarray=(/445.D3, 445.D3, 535.D3/)	! Q maximes  (hard)
           ENDIF
	   enarray=(/3.D0, 3.D0, 3.D0/)
	   enuc=enarray(1)
	   enlc=enarray(2)
	   enml=enarray(3)

	   !Ashell_c=2.3D9
	   Ashell_c=10.D0**(12.2-2.8D-5*Qarray(1))	! Bird & Kong (1994)
 	   Ashell_m=9.5D4
           Auc=((2.D0)**(enuc-1))/(Ashell_c**(enuc))
           Alc=Auc
           Aml=((2.D0)**(enml-1))/(Ashell_m**(enml))

           AucMPa=Auc*(10.D0**(6.D0*enuc))
           AlcMPa=Alc*(10.D0**(6.D0*enlc))
           AmlMPa=Aml*(10.D0**(6.D0*enml))
      ENDIF

      Quc=Qarray(1)
      Qlc=Qarray(2)
      Qml=Qarray(3)
CC =====================================================================

       crustup=(2.D0/3.D0)*crustlow
       crustup=MIN(crustup,40.D3)
C       crustup=crustlow
C       TIPPAR=TIPPAR1//'no Dorn law'
CCC ------  CALCUL DE LA  QDornml A PARTIR DELS PARAMETRES -----------
	   IF(lithtype.EQ.2) THEN
        	crustup=0
        	TIPPAR=TIPPAR1//' oceanica + Dorn law'
           ELSE
         	TIPPAR=TIPPAR1//'+ Dorn law          '	
           ENDIF	
       termln=log(stressPD*((Aml/strainrate)**(1.D0/enml)))
       TPLDL=Qml/(enml*RGAS*termln)
       QDcalc=RGAS*TPLDL*log(epsDorn/strainrate)*
     +         ((1.D0-stressPD/sigmaD)**(-2.D0))
       QDornml=QDcalc
CCC ------------------------------------------------------------------
      IF(iswitch == 1) THEN
	   WRITE(6,"(2X,'Rheological parameters:'/
     +		5X,'Upper crust:',9X,'A=',G13.5,' MPa-n s-1,   n=',
     +			F4.1,',   Q=',F8.2,' kJ mol-1'/
     +		5X,'Lower crust:',9X,'A=',G13.5,' MPa-n s-1,   n=',
     +			F4.1,',   Q=',F8.2,' kJ mol-1'/
     +		5X,'Lithospheric mantle: A=',G13.5,' MPa-n s-1,   n=',
     +			F4.1,',   Q=',F8.2,' kJ mol-1'/
     +		5X,'Dorn law:   SigmaD=',G14.5,' Pa,  QD=',F8.2,
     +			' kJ mol-1,  srD=',G14.5,' s-1'/)") 
     +		AucMPa, enuc, Quc/1D3, AlcMPa, enlc, Qlc/1D3,
     +		AmlMPa, enml, Qml/1D3, sigmaD, QDornml/1D3, epsDorn
     
      	   WRITE(6,"(2X,'Crustal thickness =',F6.2,' km',
     +       7X,'Upper crust =',F6.2,' km')") crustlow/1.D3,crustup/1.D3
           OPEN(20,FILE='stress_envelope.yield')
	   WRITE(20,"('# depth(km), comp.stress,',
     +	           ' ext.stress (MPa), Temperature(K)')")
      ENDIF
      espmecan=0.D0
      strengthc=0.D0
      strengthe=0.D0
      ycant=0.D0
      yeant=0.D0
      icont=0	
      NDEPTH=ZLITOS/Dz    
      Depth: DO iz=0,NDEPTH 
               z=iz*Dz
               icont=icont+1
               IF(z.LT.crustup) THEN
                   Q=Quc
                   epref=Auc
                   en=enuc
               ENDIF 
               IF(z.GE.crustup.AND.z.LT.crustlow) THEN 
                   Q=Qlc
                   epref=Alc
                   en=enlc
               ENDIF 
               IF(z.GE.crustlow) THEN
                   Q=Qml
                   epref=Aml
                   en=enml
               ENDIF     

!   	Brittle failure:
               brittleext = Bexten*z 
               brittlecomp = Bcompres*z 
!       Brittle failure, coeficient de friccio:
C	       IF(z.LE.crustlow) then	
C		   brittleext = ffmuf*grav*(roc-RHOH2O)*z
C	   	   brittlecomp = brittleext
C	       else
C		   brittleext = ffmuf*grav*
C     +			    (roc*crustlow+rom*(z-crustlow)-RHOH2O*z)
C	   	   brittlecomp = brittleext
C	       ENDIF

C   Ductil flow: Power law flow
               aux1=(strainrate/epref)**(1.D0/en)
               ductilpl=aux1*DEXP(Q/(en*RGAS*T(iz)))       

C   Dorn flow when stress is higher than stressPD
               ductilDl=sigmaD*(1.D0-
     +              DSQRT(RGAS*T(iz)/QDornml*log(epsDorn/strainrate)))
               IF(ductilDl.LT.0.D0) ductilDl=0.D0

C  Power and Dorn laws have to coincide at stressPD
               IF(ductilpl.GT.stressPD.AND.ductilDl.GT.0.AND.
     +             z.GE.crustlow) THEN
                      ductil=ductilDl
               ELSE
                      ductil=ductilpl
               ENDIF
           
               yieldcompres=-MIN(brittlecomp,ductil,stress_max)
               yieldextens=MIN(brittleext,ductil,stress_max)
CC ------------ Fixar un stress maxim---------------------------------
C	       IF(ABS(yieldcompres).GT.yieldmax) yieldcompres=-yieldmax
C	       IF(yieldextens.GT.yieldmax) yieldextens=yieldmax
CC--------------------------------------------------------------------
           
               IF((yieldextens-yieldcompres).GT.stresslim.OR.
     +                 z.LT.crustlow) espmecan=espmecan+Dz
     
               IF(icont.NE.1) THEN
                    ycmig=(yieldcompres+ycant)/2.D0
                    yemig=(yieldextens+yeant)/2.D0
                    strengthc=strengthc+ycmig*Dz
                    strengthe=strengthe+yemig*Dz
               ENDIF
               ycant=yieldcompres
               yeant=yieldextens
               IF(iswitch == 1) WRITE(20,"(1X,4F13.4)") 
     +	              z/1.D3,yieldcompres/1.D6,yieldextens/1.D6,T(iz)
      END DO Depth
       
      ZMECANIC=espmecan
      Flitcomp=strengthc
      Flitext=strengthe
      IF(iswitch == 1) THEN	    
            CLOSE(20)
            WRITE(6,"(4X,'Mecanic thickness (stress >',F5.1,
     +             ' MPa) :',F10.2,' m'/
     +             4X,'Strength,  compress: ',1P,E10.3,' N/m,',
     +             4X,'extens: ',1P,E10.3,' N/m')")
     +             stresslim/1.D6, ZMECANIC, Flitcomp, Flitext
            OPEN(20,FILE='Strength.yield')
               WRITE(20,"('# ',A40)") TIPPAR
               WRITE(20,"('# crust,  upper crust, Mecanic thickness,',
     +		    ' Strength comp, Strength ext,  strainrate')")
               WRITE(20,"(2F8.3,7X,F10.1,3X,3(1P,G15.4))")crustlow/1.D3,
     +               crustup/1.D3,ZMECANIC,Flitcomp,Flitext,strainrate
            CLOSE(20)
      ENDIF

      RETURN
      END SUBROUTINE Ch_tree_lib       
