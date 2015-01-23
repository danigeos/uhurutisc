	PROGRAM UHURU
!!   Data inici: 1996.    Ultima modifica: 01/08/2007

!!   MODEL WITH DIFFERENT LAYERS:  
!!	(sea), (sediment), crust, lithospheric mantle and asthenosphere.
!!
!!   Neotectonic study: steady state Temperature only on the lithosphere
!!		boundary conditions: surface Temperature and heat flow.
!!		Inputs (e,s) -> Q,L   o   (e,Q) -> s,L 
!!   Evolutive study: steady state Temperature from z=0 to ZASTH.
!!		boundary conditions: Temperature at T(z=0)=TSURF and T(z=ZASTH)=TBOTT
!!		Input s,L
!!   EXTERNAL LIBRARIES :
!!		thin_sheet.f
!!		lib/RW_PARAMETRES.f
!!		/ivone/lib/lib_uhuru/mathematics.f
!!		/ivone/lib/lib_uhuru/Ch_tree_lib.f
!!		/ivone/lib/lib_uhuru/outin.f
!!		/ivone/lib/lib_uhuru/Termica_1D_lib.f
!!		/src/tao3D/

!!   Input files: fort.11		Initial geometry: x, y, crust and lithosphere thickness (Grid.in)
!!		  parametres.in		Input parameters - open(8)
!!		  fort.4		Velocity Boundary conditions: ix,iy,type of BC,condition1,condition2 (BC.in)
!!		  fort.12 (optional)	Bodies: Zones with diferent strength (bodies.in)
!!		  fort.13 (optional)	Time and files with different velocity boundary conditions (Time_BCfiles.in)
!!		  fort.10 (optional)	Points which uhuru follow with time. (Points.in)
!!		  fort.44 (optional)	Contour of the area where the lithosphere will be removed
!!		  fort.54 (optional) => Re-write the 'Grid.in' file to 'Grid_outUhuru.in'

!!   Output files: resultats*		open(1)
!!		   Flit_Zmec_xy.res	open(14) Initial lithospheric strength, different bodies 
!!		   res_Szzaver_Flit*	open(30) Depth-averaged vertical stress over the plate, lithosphere strength	
!!			res_tisc*.st	Surface processes (erosion/sedimentation)
!!			res_tisc*.xyw
 
!!  Dt [s]	time increment, if the amount of deformation is too 
!!				large then Dt will be reduced to fit 'quotadef'
!!  AX [m]	domain length in x
!!  BY [m]	domain length in y
!!  n, m		Number of nodes (in x,y directions) minus 1: column index runs 
!!				from 0 (west) to n (east), 0 (south) to m (north) 
!!  nn		Number of nodes (n+1)*(m+1)

!!  nnsd > (n+1)(m+1)
!!  ix=0,n  /  iy=0,m  /  kxy=1,nn  /  kxy=ix+1+iy*(n+1)
!!  kvel=1,2*nn   /   iz=0,NELZ
!!  kpunts : Numero de punts dels que vaig seguint la seva deformacio, amb
!!		la Subrutina DEFORMA.
!!  nitermax	max. number of iterations for convergence between 
!!			viscosity and velocity (if exceeded => divergence).
!!			If 0 => viscosity doesn't depend on strain rate.
!!  g [m/s2]		Gravity acceleration at the surface of the study planet
!!  RHOH2O [kg/m3]	Density of sea water.
!!  rosed [kg/m3]	Density of sediments.	
!!  roc	[kg/m3]		Density of the crust.
!!  RHOAST [kg/m3]	Density of asthenosphere.	
!!  roalfa [1/K]	Volumetric thermal expansion coeff. of lith. mantle (3.5e-5 K-1).
!!  TISOTER [K]		Temperature at which the lithosphere has density RHOAST (1573 K)
!!					rom = RHOAST * (1+roalfa/2*(TISOTER-Tmoho))

!!  TSURF : Temperature at surface, z=0 (K).
!!  TBOTT : Temperature to the bottom of the model, z=ZASTH

!!  ZASTH [m]	Depth of the asthenosphere (e.g., 300e3 m)
 
!! Lithospheric mantle = GLit_ter - scrust - sediment

!!	ARRAYS:
!!  elevation [m]	Array of topography (positive above sea level).
!!  sediment [m]	Array of sediment thickness.	
!!  scrust [m]		Array of crustal thickness. 
!!  GLit_ter [m]	Array of lithospheric thickness = mantle + crust + sediment
!				             (defined with the isotherm TISOTER)
!!  GLit [m]		Array of lithospheric thickness
!				             (where the conductivity change)
!!  Tmoho [K]		Array of temperature distribution at Moho.
!!  Flitos [N/m]  	Array of the strength of the lithosphere (integral of the stress envelope)
!!				to Subroutine VELOCITYFIELD  send  Flitos/2*thick_layer  (visTer)
!!  visco [Pa*s]	Array of the therm of the viscosity	visco=Flitos/(2.0*thick_layer*strainrate)
 
!!  epuntzz [s-1]	Array of the vertical strain rate.
!!  TEMPE{1/2}(kxy,iz) [K]  Temperature

!!  D_Bodies:	Array with the different bodies (kpunts1+kpunts2+kpunts3+3)
!!		D_Bodies(1,1)=kpunts1, D_Bodies(1,2)=qFlit1			
!!		D_Bodies((i=2,kpunts1+1),1/2)= posicio1((kp=1,kpunts1),1) posicio1((kp=1,kpunts1),2)
!!		D_Bodies(kpunts1+2,1)=kpunts2, D_Bodies(kpunts1+2,2)=qFlit2
!!		D_Bodies((i=kpunts1+3,kpunts1+kpunts2+2),1/2)= posicio2((kp=1,kpunts2),1) posicio2((kp=1,kpunts2),2)
!!		D_Bodies(kpunts1+kpunts2+3,1)=kpunts3, D_Bodies(kpunts1+kpunts2+3,2)=qFlit3
!!		D_Bodies((i=kpunts1+kpunts2+4,kpunts1+kpunts2+kpunts3+3),1/2)= posicio3((kp=1,kpunts3),1) posicio3((kp=1,kpunts3),2)

!!  Dz : interval de discretitzacio en profunditat.
!!  h = AX/n = Dx
!!  tau = BY/m = Dy

!!  VELOCITY(m/s) = VELOCITY(mm/any)/FACVEL 

!!  RESTART : Variable logica
!!     SI RESTART=T -> TINC L'ESTAT INICIAL EN UN FITXER : estat_inicial
!!     SI RESTART=F -> CREO L'ESTAT INICIAL AMB LES SUBRUTINES :
!!   PROHEAT=HSURF*DEXP(-Z/HEXP)   HEXP=constant
!!
!!  dGLit_allow:maximum Lithosphere gradient permited = [d(GLit)/dx]/GLit   [e.g. (15km/30km)/70km = 7.14D-6 ] 
!!  Dif_K_L: 	constant diffusive filter. [Dvis=Dif_K*(d2zeta/dx2+d2zeta/dy2)]  ( Dif_K_L=1.D8  Dif_K_L<3.D8 (?))

!!	If rain=0 and Krain=0 -> no rivers		(rain~500)
!!	If Kerosdif=0 ->	no diffusive erosion.	(0 <= Kerosdif <=5000)
!!	No surface porcesses -> Kerosdif=0, rain=0 and Krain=0.

!!	iflexure=1 : Flexure calculation

!!	irel_ZL=0 ==> nivell de relaxacio fixe, ZLITOS_rel => GLit(kxy)=ZLITOS_rel+elevation(kxy)
!!	irel_ZL=1 ==> GLit(kxy)= depth of isotherm (TLITOS)
!!	irel_ZL=2 ==> GLit(kxy)=GLit_ter(kxy)+DGL_ter

!!	Ldepth_var=0 -> constant Lithospheric depth (Depth_lit) and no temporal variations
!!	iTemp_var=0  -> No Temperature temporal variations. TEMPE2=TEMPE1

!!	iremesh=0 -> No remeshing
!!	iremesh=1 -> Remeshing, when the indenter go norther than Y0+Dy

!! 	ifixv=1 => keep in a file the points inside body 1 (subrutine VISCOSITAT) => velocity fixed in these points
!!			in subrutine VELOCITYFIELD

!!		type *, system ('mv '//BCfile//' pippo')

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*170 TITLE 
      CHARACTER*90 TIT_Struct
      CHARACTER*40 TIPPAR
      CHARACTER*84 BCfile,Bodyfile,BCtime,BCfNew,file_P
      CHARACTER*80 TFPARA
      CHARACTER*55 zname
      CHARACTER*20 fresult_ini,file_result,file_eros,file_rot,file_Pav,
     +		   file_TK,file_ep
      LOGICAL RESTART
      PARAMETER (kpuntssd=500, nsBC=5, Dtmin=5.D3, FACVEL=3.1536D10,
     +		FACTEMP=3.1536D7, TMOHOlim=1496.D0, DBANDAKZ=20.D0,
     +		PI=3.1415926535897932D0, g=9.8D0, crustinf=10.D3,
     +          visinf=1.D17, Depth_lit=120.D3, irel_ZL=2, 
     +		ZLITOS_rel=130.D3, DGL_ter=20.D3, Ldepth_var=1, 
     +		iTemp_var=1, ifixv=1)
      DATA iw_erosion /0/, vissup /1.D25/
      REAL dt_eros, Kerosdif,K_river_cap, erodability, erodability_sed,
     +		l_fluv_sedim,rain, Krain, CXrain, CYrain, evaporation, 
     +		xmin, xmax, ymin, ymax, Dt_tisc, Te_flexure
      INTEGER*2 Nx, Ny
      PARAMETER (dt_eros=0.05)
      DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: elevation, scrust,
     +		elevation0, scrust0, hmantle0, sediment, GLit, GLit_ter, 
     +		u, v, epuntzz, visco, Flitos, Tmoho, Qsurface, erosion,
     +		fin_rotat, A_removal, B_removal
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: TEMPE1,TEMPE2
      DIMENSION CONDUC(0:3),Qarray(3),enarray(3),Aarray(3),
     +		D_Bodies(kpuntssd,2),BCMy(nsBC),BCfNew(nsBC),
     +		Points_t(0:kpuntssd,2)

      CONDUC=0.D0
      Qarray=0.D0
      enarray=0.D0
      Aarray=0.D0

      D_Bodies=0.D0
      BCMy=0.D0
      Points_t=0.D0

      BCfile="fort.4"		! file with the velocity boundary conditions
      Bodyfile="fort.12"	! file with the diferent bodies
      BCtime="fort.13"		! file with the time(My) when the BC change     
      file_P="fort.10"		! file with points     

      CALL READP (TFPARA, rosed, roc, roalfa, RHOAST, RHOH2O, TSURF,
     +            TBOTT, ZASTH, TISOTER, THDIFF, CONDUC, HSURF, HEXP,
     +            PHEAT_m, IRheology_type, Qarray, enarray, Aarray, 
     +		  m, n, NELZ, strainrate, nitermax, tallmax,
     +            alfa, tinicial, Dtany, npasmax, NPASINT, ITSR, 
     +            visco_cnst, Iremoval, Time_removal,Zremoval,TISO_rem, 
     +            Zctall, iremesh, dvis_allow, Dif_K, Te_flexure,
     +		  Kerosdif, rain, Krain, evaporation, erodability,
     +		  erodability_sed, K_river_cap, l_fluv_sedim,
     +		  CXrain, CYrain)
     
      dGLit_allow=dvis_allow
      Dif_K_L=Dif_K
      WRITE(6,"(3X,'Limit values: minimum crustal thickness',F8.3,' km'/
     +		17X,'Viscosity minimum',1P,G10.3,' Pa s,  maximum ',
     +		1P,G10.3,' Pa s')") crustinf/1.D3, visinf, vissup

  !! Control of the existence of some files file (file_P). if (Ifile_P=0 => file exist) else (no exist)
      OPEN(10,FILE=file_P,STATUS='OLD',ACTION='READ',IOSTAT=Ifile_P)
	 IF(Ifile_P==0) WRITE(6,"(3X,'File with points: ',A84)") file_P
      OPEN(12,FILE=Bodyfile,STATUS='OLD',ACTION='READ',IOSTAT=IBodyex)
	 IF(IBodyex==0) 
     +	     WRITE(6,"(3X,'File with different bodies: ',A84)") Bodyfile  
      OPEN(13,FILE=BCtime,STATUS='OLD',ACTION='READ',IOSTAT=IBCtime)
	 IF(IBCtime==0) THEN
      	     CALL TIME_BCFILE (BCtime, nsBC, BCMy, BCfNew, NBC)
	     BCfile=BCfNew(1)
	     WRITE(6,"(3X,'Velocity boundary conditions change in time'/
     +		  3X,'Initial velocity boundary conditions: ',A84,/
     +		  3X,'File with the time(My) when the BC change:',A84)")
     +		  BCfile, BCtime
	     IBC_next=2
	 ENDIF

      nresults=48
      INTERVWK=npasmax/(NPASINT+1)
      NPASGRUP=1
      NPASKEEP=INTERVWK*NPASGRUP
		
      nn=(n+1)*(m+1)
      WRITE(6,"(/3X,'Thermal conductivity transition from lithospheric',
     +	       ' mantle to astenosphere:'/7X,I4,' +',I4,' points, up +',
     +	       ' down [DBANDAKZ]')") INT(DBANDAKZ), INT(DBANDAKZ)
	
       npasos=0
C ---------------------------------------------------------------------
C -------------------------  ESTAT INICIAL  ---------------------------
C ---------------------------------------------------------------------
      IF(RESTART) THEN
		! No funciona !!!!
          PRINT*,'The initial state in a file. Is not working'
	  STOP
   !       fresult_ini='estat_inicial'
   !       nporta=1
   !       CALL RERESULT (fresult_ini, nporta, TITLE, Tsegons,Tanys,
   !  +                      m, n, NELZ, nn, nnsd, AX, BY, Dz,
   !  +                      TLITOS, scrust, GLit, visco, u, v, Tmoho,
   !  +			    Qsurface, elevation, sediment, kpuntssd,
   !  +                      D_Bodies,Points_t)
   !       u=u/FACVEL	! From mm/year to m/s
   !       v=v/FACVEL
   !       tinicial=Tanys          
   !       TEMPSMa=tinicial/1.D6
   !       ZASTH=NELZ*Dz
   !       Ndepth=NELZ+1
   !       PRINT*,'      TEMPS INICIAL =',TEMPSMa,'  Ma'
      ELSE
          WRITE(6,"(/'===== CONSTRUCTION OF THE INITIAL MODEL ======')")
          TEMPSMa=tinicial/1.D6
          Dtsegons=Dtany*FACTEMP
          Dz=ZASTH/NELZ
	  Dz3dc=Dz*100.D0
	  Dz3dc_int=DINT(Dz3dc)
	  Dz=Dz3dc_int/100.D0
	  ZASTH=NELZ*Dz
          WRITE(6,"(/'Number of vertical points, iz = 0,...,',I4,
     +		' (NELZ)'/5X,'Maximun depth (ZASTH):',F13.3,' m  ',
     +		'With one point every:',F10.3,' m')") NELZ, ZASTH, Dz
         Ndepth=NELZ+1

!! ----------------------  INPUT DATA  -------------------------
!!   INPUT_DATA=1 ==>  x, y, elevation, surface heat flow (only for neotectonic)
!!   INPUT_DATA=2 ==>  x, y, elevation, crustal thickness (only for neotectonic)
!!   INPUT_DATA=3 ==>  x, y, crustal and lithospheric thickness
          READ(11,"(A90)") TIT_Struct
          READ(11,*) n_file,m_file,Dx_file,Dy_file,INPUT_DATA
           IF(n_file/=n.OR.m_file/=m) WRITE(6,"(/'Number of ',
     +            'horizontal points in x and y directions: n =',
     +            I4,'  m =',I4)") n_file,m_file
          n=n_file
          m=m_file
          nn=(n+1)*(m+1)

          ALLOCATE(elevation(nn),STAT=istat1) 
          ALLOCATE(sediment(nn),STAT=istat2)
          ALLOCATE(scrust(nn),STAT=istat3) 
          ALLOCATE(scrust0(nn),STAT=istat3) 
          ALLOCATE(elevation0(nn),STAT=istat3) 
          ALLOCATE(hmantle0(nn),STAT=istat3) 
          ALLOCATE(GLit_ter(nn),STAT=istat4)       
          ALLOCATE(GLit(nn),STAT=istat5) 
          ALLOCATE(u(nn),STAT=istat6) 
          ALLOCATE(v(nn),STAT=istat7)
	  ALLOCATE(epuntzz(nn),STAT=istat8)
	  ALLOCATE(Tmoho(nn),STAT=istat9)
	  ALLOCATE(Qsurface(nn),STAT=istat10)
	  ALLOCATE(erosion(nn),STAT=istat12)
	  ALLOCATE(A_removal(nn),STAT=istat12)
	  ALLOCATE(B_removal(nn),STAT=istat12)
          ALLOCATE(TEMPE1(nn,0:NELZ),STAT=istat11) 
        elevation=0.D0
        sediment=0.D0
        scrust=0.D0
        scrust0=0.D0
        elevation0=0.D0
        hmantle0=0.D0
        GLit_ter=0.D0
        GLit=0.D0
        u=0.D0
        v=0.D0
        epuntzz=0.D0
        Tmoho=0.D0
        Qsurface=0.D0
        erosion=0.D0
        A_removal=0.D0
        B_removal=0.D0
        TEMPE1=0.D0
          IF(istat1>0.or.istat2>0.or.istat3>0.or.istat4>0.or.istat5>0. 
     +	   or.istat6>0.or.istat7>0.or.istat8>0.or.istat9>0.or.istat10>0.
     +	    or.istat11>0.or.istat12>0) STOP " Error allocation main" 
CCC =========== INPUT_DATA=1 : elevation, surface heat flow =============
!!  SUBROUTINE HEATFLOW_TOPO: WITH THE ELEVATION AND SURFACE HEAT FLOW
!!    AS A INPUT DATA, CALCULATE THE CRUSTAL AND THE DISTRIBUTION OF
!!    TEMPERATURE WITH DEPTH (LITHOSPHERIC THICKNESS). ASSUMING STEADY
!!    STATE AND LOCAL ISOSTASY.
          IF(INPUT_DATA.EQ.1) THEN
              WRITE (6,"(/'READ THE FILE WITH: ',4X,A90/
     +               4X,'(x, y, elevation, surface heat flow)'/
     +               5X,'where (x,y) is a regular grid. ',
     +               'Everything in S.I.')") TIT_Struct
	      TLITOS=TISOTER+15.D0
              CALL HEATFLOW_TOPO (m, n, nn, Dz, NELZ, AX, BY, 
     +                         roc, roalfa, RHOAST, RHOH2O, CONDUC,
     +                         HSURF, HEXP, PHEAT_m, elevation, 
     +                         scrust, GLit, GLit_ter, TEMPE1,
     +                         TMOHOlim, TSURF, TBOTT, TISOTER, TLITOS,
     +                         Qsurface, Tmoho)
          ENDIF
CCC =========== INPUT_DATA=2 : elevation, crustal thickness ==============
CCC  SUBROUTINE CRUST_TOPO: A PARTIR DEL GRUIX CORTICAL I LA TOPOGRAFIA
!!  SUBROUTINE CRUST_TOPO: WITH THE ELEVATION AND CRUSTAL THICKNESS 
!!    AS A INPUT DATA, CALCULATE THE DISTRIBUTION OF TEMPERATURE 
!!    WITH DEPTH (LITHOSPHERIC THICKNESS). ASSUMING STEADY
!!    STATE AND LOCAL ISOSTASY.
CC  !!!!  Aquesta opcio dona problemes !!!!!
          IF(INPUT_DATA.EQ.2)
     +        CALL CRUST_TOPO (m, n, nn, Ndepth, Dz, NELZ,
     +                        AX,BY,roc,RHOAST,RHOH2O,CONDUC,HSURF,
     +                        HEXP,TSURF,TBOTT,TISOTER,scrust,GLit, 
     +                        GLit_ter, TEMPE1, TLITOS, elevation)
CCC =========== INPUT_DATA=3 : crustal and lithospheric thickness ===========
CCC                   C.C. T superficial i base 
          IF(INPUT_DATA.EQ.3) THEN
              WRITE (6,"(/'READ THE FILE WITH: ',4X,A90/
     +               4X,'(x, y, crustal thick, lithospheric thick,',
     +               ' sediment thick)'/
     +               5X,'where (x,y) is a regular grid. ',
     +               'Everything in S.I.')") TIT_Struct

              CALL GRUIXOSINI (m, n, nn, AX, BY, scrust, GLit, sediment)
            	  Dx=AX/n
           	  Dy=BY/m
  		  zname='Input Crustal thickness (scrust)'
		  !CALL smooth_horizontal (m, n, nn, Dx, Dy, scrust)
                  ds_allow=dGLit_allow
		  iwr_filt=-1
      	          CALL fit_gradients (m, n, nn, Dx, Dy, Dif_K_L, 
     +	    		              ds_allow, scrust, zname, iwr_filt)
		  zname='Input Lithosphere thickness (GLit)'
		  !CALL smooth_horizontal (m, n, nn, Dx, Dy, GLit)
		  iwr_filt=-1
		  CALL fit_gradients (m, n, nn, Dx, Dy, Dif_K_L, 
     +	    			     dGLit_allow, GLit, zname, iwr_filt)

!	      IF(Ifile_P==0) 
!     +		    CALL DEFORMA_Points (Dtsegons,TEMPSMa, m,n, AX,BY,
!     +				   u,v,nn, kpuntssd, Points_t, file_P)
                 
!	      IF(IBodyex==0) 
!     +	   	    CALL DEFORMA_Bodies (Dtsegons, TEMPSMa, m,n, AX,BY,
!     +                           u, v,nn, kpuntssd, D_Bodies, Bodyfile)
	
              CALL TEMP_STEADY_DF (m, n, nn, NELZ, Dz, DBANDAKZ, CONDUC,
     +                       HSURF, HEXP, PHEAT_m, TSURF, TBOTT, scrust,
     +                       GLit_ter, GLit, sediment, TEMPE1, TMOHOlim,
     +                       TLITOS, TISOTER, Tmoho, Qsurface, AX, BY, 
     +			     dGLit_allow, Dif_K_L, irel_ZL, DGL_ter)

       	      CALL ISOSTASY (m, n, nn, RHOH2O, rosed, roc, RHOAST,
     +			 roalfa, elevation, sediment, scrust, GLit_ter,
     +			 Tmoho, TISOTER, Dx, Dy, kpuntssd, D_Bodies)
	       IF(Ldepth_var==0) THEN		!!! constant lithospheric depth 
	            PRINT*,'CONSTANT LITHOSPHERIC DEPTH'
	            GLit_ter=Depth_lit+elevation
	            GLmax=VAL_MAX(GLit_ter,nn)
	            GLmin=VAL_MIN(GLit_ter,nn)
	            WRITE(6,"('GLit_ter maxim:',F8.3,' km,     ',
     +		      'GLit_ter minim',F8.3,' km')") GLmax/1D3,GLmin/1D3
	            DO i=1,5						
	               CALL ISOSTASY (m,n,nn,RHOH2O,rosed,roc, RHOAST, 
     +			    roalfa,elevation,sediment,scrust,GLit_ter,
     +			    Tmoho,TISOTER,Dx,Dy,kpuntssd,D_Bodies)
	            GLit_ter=Depth_lit+elevation
	            GLmax=VAL_MAX(GLit_ter,nn)
	            GLmin=VAL_MIN(GLit_ter,nn) 
	            WRITE(6,"('GLit_ter maxim:',F8.3,' km,    ',
     +		     ' GLit_ter minim',F8.3,' km')") GLmax/1D3,GLmin/1D3
	            END DO	
	       END IF
 	       IF(irel_ZL==0) GLit=ZLITOS_rel+elevation		!!  Fixed relaxation depth

          ENDIF
CCC ================================================================
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	      IF(Ifile_P==0) 
     +		    CALL DEFORMA_Points (Dtsegons,TEMPSMa, m,n, AX,BY,
     +				   u,v,nn, kpuntssd, Points_t, file_P)
                 
	      IF(IBodyex==0) 
     +	   	    CALL DEFORMA_Bodies (Dtsegons, TEMPSMa, m,n, AX,BY,
     +                           u, v,nn, kpuntssd, D_Bodies, Bodyfile)
          Dx=AX/n
          Dy=BY/m
	  ix_write_T=aint(4000.D3/Dx)
	  IF(npasmax<2) ix_write_T=aint(2530.D3/Dx)
	  	  
          IF(Dx_file/=Dx.OR.Dy_file/=Dy) THEN
                WRITE (6,"(/4X,'The file (x,y,s,L) has some problem'/
     +                    6X,' Dx_file=',F8.3,' km,   Dx=',F8.3,' km'/
     +                    6X,' Dy_file=',F8.3,' km,   Dy=',F8.3,' km'/
     +                    2X,'STOP PROGRAM')")Dx_file,Dx,Dy_file,Dy
	        STOP
          ENDIF
          ALLOCATE(visco(nn),STAT=istat2) 
          ALLOCATE(Flitos(nn),STAT=istat3) 
              visco=0.D0
              Flitos=0.D0
              IF(istat2>0.OR.istat3>0.) 
     +		    STOP " Error allocation  visco or Flitos"
            iwrite=1
            file_ep='res_StrainRateEff0'
          Flitos=0.D0
          CALL VISCOSITAT (TEMPSMa, npasmax, m, n, NELZ, nn, AX, BY, Dz,
     +         	      u, v, epuntzz, elevation, sediment, scrust,
     +		      GLit_ter,TEMPE1, visco, nitermax, strainrate, 
     +                TISOTER, IRheology_type, Qarray, enarray, Aarray,
     +                TIPPAR,kpuntssd, D_Bodies, vissup, visinf, ITSR,
     +		      visco_cnst, Flitos, dvis_allow, Dif_K, iwrite,
     +		      file_ep, ifixv)

            file_Pav='res_Szzaver_Flit0'
      OPEN(13,FILE=BCfile,STATUS='OLD',ACTION='READ',IOSTAT=IBCfile)
      CLOSE (13)
      IF(IBCfile==0) THEN
          CALL velvis (TEMPSMa, Dtsegons, AX,BY, n,m,nn, nitermax,
     +			tallmax, alfa, g, Flitos, visco, vissup, visinf,
     +			elevation, sediment, scrust, GLit_ter, u, v,
     +                  epuntzz, RHOH2O, rosed, roc, RHOAST, roalfa, 
     +			TISOTER, Tmoho, ZASTH, BCfile, iwrite, 
     +			file_Pav, kpuntssd, D_Bodies, ifixv)
     
!	        srmax=VAL_MAX(epuntzz,nn)
!		IF(srmax>=1.D-18) THEN
! 		     zname='vertical strain rate (epuntzz)'
!	 	     !CALL smooth_horizontal (m, n, nn, Dx, Dy, epuntzz)
!		     iwr_filt=-1
!  !  	             CALL fit_gradients (m, n, nn, Dx, Dy, Dif_K,	!! Tibet
!  !   +	    		          dvis_allow, epuntzz, zname, iwr_filt)
!    	             CALL fit_gradients (m, n, nn, Dx, Dy, Dif_K_L,	!! Alps
!     +	    		          dGLit_allow, epuntzz, zname, iwr_filt)
!		END IF
      ENDIF
          DEALLOCATE(Flitos)  
          Dtany=Dtsegons/FACTEMP
	  OPEN(25,FILE='res_Tempe0')
	    DO ix=0,n
		 IF(ix==ix_write_T) THEN
		     WRITE(25,"(' # ',2F12.3,' x')") TEMPSMa,ix*Dx/1.D3
		     DO iy=0,m
		        kxy=ix+1+iy*(n+1)
			WRITE(25,"(1X,3F12.3)") 
     +				(iy*Dy,iz*Dz,TEMPE1(kxy,iz), iz=0,NELZ)
		     END DO
		 ENDIF		    
	    END DO
	  CLOSE(25)
          file_result='resultats'//CHAR(nresults)
          TITLE=TFPARA//TIPPAR
          !WRITE(6,*) TITLE
          CALL WRRESULT (file_result,TITLE, TEMPSMa, npasos,
     +                   m, n, NELZ, nn, FACTEMP, FACVEL, AX, BY, Dz,
     +                   TLITOS, scrust, GLit_ter, GLit, visco, u, v,
     +                   epuntzz, Tmoho, Qsurface, elevation, 
     +			 sediment, kpuntssd, D_Bodies, Points_t)

          WRITE(6,"(/'======== INITIAL MODEL CONSTRUCTED ========')")
      ENDIF
      WRITE(6,"(/2X,'Depth of the asthenosphere (ZASTH) :',F8.3,' km'/
     +    2X,'Vertical nodes distance: Dz =',F7.2,' m'/)") ZASTH/1.D3,Dz
      
      scrust0=scrust						!! keep initial crustal thickness
      elevation0=elevation					!! keep initial elevation
      DO kxy=1,nn						!! keep initial mantle weight
	hmantle0(kxy)=GLit_ter(kxy)-scrust(kxy)-sediment(kxy)				!! keep initial mantle weight
      END DO
      XP0=D_Bodies(2,2)
C --------------------------------------------------------------------
C -------------------  CONSTRUIT L`ESTAT INICIAL  --------------------
C --------------------------------------------------------------------
      IF(npasmax.EQ.0) GOTO 9876 
      WRITE (6,"(' ***************************************************'/
     +		 ' *                     TIME STEPS                  *'/
     +	    ' ***************************************************'/)")

      NumPoints=Points_t(0,1)
      ALLOCATE(fin_rotat(NumPoints),STAT=istat)
	fin_rotat=0.D0
      IF(istat1>0) STOP " Error allocation fin_rotat" 
      fin_rotat=0.D0

      Nx=n+1
      Ny=m+1	
      xmin=0.0
      xmax=REAL(AX)
      ymin=0.0
      ymax=REAL(BY)
      IF(Kerosdif/=0.0.OR.rain/=0.0.OR.Krain/=0.0)  iw_erosion=1	! =>  Model with erosion and sedimantation
      IF(Te_flexure/=0.0)  iflexure=1						! =>  Model with FLEXURE

      IF(iw_erosion==1 .OR. iflexure==1) THEN
      	  IF(iflexure==1) WRITE(6,"(3X,'Elastic thickness',F8.3,' km')")
     +				Te_flexure/1.D3
	  CALL init_tisc (Nx, Ny, xmin, xmax, ymin, ymax, dt_eros,
     +		Kerosdif, K_river_cap, erodability, erodability_sed,
     +		l_fluv_sedim, rain, Krain, CXrain, CYrain, evaporation, 
     +		Te_flexure)
      ENDIF

1000     CONTINUE 

      npasos=npasos+1
      Dtany=Dtsegons/FACTEMP
      TEMPSMa=TEMPSMa+Dtany/1.D6
      iwrite=0
      IF(npasos==NPASKEEP.OR.npasos==npasmax.OR.npasos==1) THEN
                iwrite=1
                nresults=nresults+1
                nres_dec1=AINT((nresults-48.D0)/10.D0)
                nres_unit=nresults-10*nres_dec1
                nres_dec=nres_dec1+48
      ENDIF

      XP1=D_Bodies(2,2)-XP0	!! re-mesh, move all the model -Dy
      IF(iremesh==1. AND. XP1>Dy) THEN
	    IF(D_Bodies(2,1)<0.0 .OR. D_Bodies(2,2)<0.0) THEN
		 WRITE(6,"(///3X,'THE REFERENCE RE-MESHING POINT IS ',
     +		      'OUT OF THE DOMAIN. THE RE-MESHING WILL NOT WORK',
     +		      /5X,'CHANGE THE REFERENCE RE-MESHING POINT.'///)")
     		 STOP
	    ENDIF
	    WRITE(6,"(3X,'XP1:',F7.2,' Remesh, move the model Dy:',
     +			F8.2)") XP1/1.D3, Dy/1.D3
	    CALL RE_MESH (m, n, nn, NELZ, kpuntssd, Dx, Dy, elevation,
     +			  scrust, sediment, GLit, GLit_ter, u, v,
     +			  epuntzz, visco, Tmoho, Qsurface,
     +			  erosion, TEMPE1, D_Bodies, Points_t)
            PRINT*,'RE_MESHED'
      ENDIF
      WRITE (6,335) npasos,Dtany,TEMPSMa
 335  FORMAT(/'=====================================================',
     +		'======================'/10X,'Time Step Number :',I4//
     +		5X,'Dt =',F12.2,' years,',6X,'Time =',F10.4,' My'/)

           IF(IBCtime==0.AND.IBC_next<=NBC) THEN	!!BC change in time
	       IF(TEMPSMa>=BCMy(IBC_next)) THEN
		  BCfile=BCfNew(IBC_next)
		  WRITE(6,"(3X,'New velocity boundary conditions,',
     +	                        ' file: ',A84/)") BCfile
		  IBC_next=IBC_next+1
	       ENDIF
	   ENDIF
	   	
	   IF(Ifile_P==0) THEN
      	        CALL DEFORMA_Points (Dtsegons,TEMPSMa, m,n, AX,BY, u,v,
     +				      nn, kpuntssd, Points_t, file_P)
		NumPoints=Points_t(0,1)
		CALL FINITE_ROTATION (Dtsegons, m, n, AX, BY, u, v, nn,
     +				 kpuntssd,Points_t,NumPoints,fin_rotat)
	   ENDIF

           IF(IBodyex==0) 
     +	   	CALL DEFORMA_Bodies (Dtsegons, TEMPSMa, m, n, AX, BY,
     +                    	u, v, nn, kpuntssd, D_Bodies, Bodyfile)
!!      Sediment thickness
           sedmax=VAL_MAX(sediment,nn)
	   IF(sedmax>1.D-100) THEN
	       IBC_thicken=0   !  0 => thickness=thickness_old  (boundaries)
	       IF(IBC_thicken==0) THEN
                   WRITE(6,"(4X,'...New Sediment thickness. BC: ',
     +                     'thickness(t+Dt)=thickness(t)')")
               ELSE
                   WRITE(6,"(4X,'...New Sediment thickness. BC: ',
     +            '[d(thickness)/dx](t+Dt)=[d(thickness)/dy](t+Dt)=0')")
	       ENDIF  
               CALL THICKEN (Dtsegons, n, m, nn, AX, BY, u, v, 
     +				epuntzz, sediment, IBC_thicken)
	   ENDIF

!!      Total erosion
           erosmax=VAL_MAX(erosion,nn)
	   IF(erosmax>1.D-100) THEN
	       IBC_thicken=0   !  0 => thickness=thickness_old  (boundaries)
	       IF(IBC_thicken==0) THEN
                   WRITE(6,"(4X,'...New erosion acumulation. BC: ',
     +                     'thickness(t+Dt)=thickness(t)')")
               ELSE
                   WRITE(6,"(4X,'...New erosion acumulation. BC: ',
     +            '[d(thickness)/dx](t+Dt)=[d(thickness)/dy](t+Dt)=0')")
	       ENDIF  
               CALL THICKEN (Dtsegons, n, m, nn, AX, BY, u, v, 
     +				epuntzz, erosion, IBC_thicken)
	   ENDIF

!!      Crustal thickness
	   IBC_thicken=0   !  1 => d(thickness)/dx=d(thickness)/dy=0 (boundaries)
	   IF(IBC_thicken==0) THEN
                WRITE(6,"(4X,'...New Crustal thickness.   BC: ',
     +                     'thickness(t+Dt)=thickness_old(t)')")
           ELSE
                WRITE(6,"(4X,'...New Crustal thickness.   BC: ',
     +           '[d(thickness)/dx][t+Dt]=[d(thickness)/dy][t+Dt]=0')")
	   ENDIF  
           CALL THICKEN (Dtsegons, n, m, nn, AX, BY, u, v, epuntzz,
     +				scrust, IBC_thicken)

		    crustmin=VAL_MIN(scrust,nn)		!! Limit the lower values of the crustal thickness
		    IF(crustmin<crustinf) THEN
			DO ixy=1,nn
			   sc=scrust(ixy)
			   scrust(ixy)=MAX(sc,crustinf)
			END DO
			cmin=VAL_MIN(scrust,nn)
			WRITE(6,"(5X,'After limited to ',F10.2,' m,  ',
     +			    'Crust minim:',F10.2,' m')") crustinf,cmin
		    ENDIF
	
 !          CALL ISOSTASY (m, n, nn, RHOH2O, rosed, roc, RHOAST,
 !    +			 roalfa, elevation, sediment, scrust, GLit_ter,
 !    +			 Tmoho, TISOTER, Dx, Dy, kpuntssd, D_Bodies)

           ALLOCATE(TEMPE2(nn,0:NELZ),STAT=istat) 
		TEMPE2=0.D0
           IF (istat>0) STOP " Error allocation TEMPE2(:,:)" 

	   IF(iTemp_var==0) THEN
		PRINT*,'  Temperature constant with time'
		TEMPE2=TEMPE1
	    ELSE
	    	Niremoved_old=Niremoved
	        Niremoved=0
		CALL TEMPEDT (TEMPSMa,Dtsegons,m,n,nn,Dz,NELZ,DBANDAKZ,
     +                 	CONDUC, THDIFF, HSURF, HEXP, PHEAT_m, TLITOS,
     +			elevation,sediment,scrust,GLit, TEMPE1,TEMPE2,
     +			epuntzz,u, v, AX,BY, GLit_ter, Tmoho, TISOTER,
     +			TMOHOlim,Qsurface, dGLit_allow, Dif_K_L,irel_ZL,
     +			DGL_ter, Iremoval, Time_removal, Zremoval, 
     +			TISO_rem, Zctall, Niremoved,A_removal,B_removal)
	   ENDIF
	   DEALLOCATE(TEMPE1)  
	   IF(Niremoved>0 .OR. Niremoved_old>0) THEN	! convective removal at some points of the grid or in the previous time step
		IF(Niremoved>0) WRITE(6,"(3X
     +		      'Lithosphere removed in',I8,' points')") Niremoved
     			!vissup=0.6D23
			!IRheology_type=2  !! Change rheological parameteers
			!visco_cnst=visco_cnst/2.D0	!!!!!!!! ojo, change lithosphere strength
			!ITSR=3				!! Change (viscosity) - (strain rate) relation
			!visco_cnst=1.D22
			!WRITE(6,"(/3X,'Constant viscosity everywhere:',
     +			!	G13.5,' Pa/s')") visco_cnst
		IF(iwrite/=1 .AND. Iremoval/=3) THEN		! After some convective removal, keep results in a file
		     iwrite=1
		     nresults=nresults+1
		     nres_dec1=AINT((nresults-48.D0)/10.D0)
		     nres_unit=nresults-10*nres_dec1
		     nres_dec=nres_dec1+48
	   	ENDIF
	   ENDIF
	    
	     IF(Ldepth_var==0) THEN		!!! constant lithospheric depth  
           		PRINT*,'CONSTANT LITHOSPHERIC DEPTH'	 
	        	GLit_ter=Depth_lit+elevation
             END IF

           CALL ISOSTASY (m, n, nn, RHOH2O, rosed, roc, RHOAST,
     +			 roalfa, elevation, sediment, scrust, GLit_ter,
     +			 Tmoho, TISOTER, Dx, Dy, kpuntssd, D_Bodies)
	     IF(irel_ZL==0) GLit=ZLITOS_rel+elevation		!!! Fixed relaxation depth

	     IF(Ldepth_var==0) THEN		!!! constant lithospheric depth  
	    		GLit_ter=Depth_lit+elevation					
              		GLmax=VAL_MAX(GLit_ter,nn)	
			GLmin=VAL_MIN(GLit_ter,nn)
	        	WRITE(6,"('GLit_ter maxim:',F8.3,' km,    ',
     +		     ' GLit_ter minim',F8.3,' km')") GLmax/1D3,GLmin/1D3
           END IF

           IF(iwrite==1) THEN
              IF(nres_dec==48) THEN
                file_Pav='res_Szzaver_Flit'//CHAR(nres_unit)
		file_ep='res_StrainRateEff'//CHAR(nres_unit)
              ELSE
                file_Pav=
     +		     'res_Szzaver_Flit'//CHAR(nres_dec)//CHAR(nres_unit)
		file_ep=
     +		    'res_StrainRateEff'//CHAR(nres_dec)//CHAR(nres_unit)
             ENDIF
           ENDIF
          ALLOCATE(Flitos(nn),STAT=istat3) 
		Flitos=0.D0
           CALL VISCOSITAT (TEMPSMa, npasmax, m, n, NELZ, nn, AX, BY,
     +			    Dz, u, v, epuntzz, elevation, sediment,
     +			    scrust, GLit_ter, TEMPE2, visco, nitermax,
     +			    strainrate, TISOTER, IRheology_type, Qarray,
     +			    enarray, Aarray, TIPPAR, kpuntssd, D_Bodies,
     +			    vissup, visinf, ITSR, visco_cnst, Flitos,
     +			    dvis_allow, Dif_K, iwrite, file_ep, ifixv)
           TITLE=TFPARA//TIPPAR

	   !IF(TEMPSMa>Time_removal) THEN	!!!!!!!! ojo, reduccio de la viscositat
	!	visco=visco/100.0
	!	PRINT*,'!!! OJO, reduccio factor 100 de la viscositat!!'
        !   END IF
           CALL velvis (TEMPSMa, Dtsegons, AX, BY, n, m, nn, nitermax,
     +			tallmax, alfa, g, Flitos, visco, vissup, visinf,
     +			elevation, sediment, scrust, GLit_ter, u, v,
     +			epuntzz, RHOH2O, rosed, roc, RHOAST, roalfa,
     +			TISOTER, Tmoho, ZASTH, BCfile, iwrite, file_Pav,
     +			kpuntssd, D_Bodies, ifixv)
!	        srmax=VAL_MAX(epuntzz,nn)
!		IF(srmax>=1.D-18) THEN
! 		     zname='vertical strain rate (epuntzz)'
!		     !CALL smooth_horizontal (m, n, nn, Dx, Dy, epuntzz)
!		     iwr_filt=-1
!  !  	             CALL fit_gradients (m, n, nn, Dx, Dy, Dif_K, 	!! Tibet
!  !   +	    		           dvis_allow, epuntzz, zname, iwr_filt)
!    	             CALL fit_gradients (m, n, nn, Dx, Dy, Dif_K_L,	!! Alps
!     +	    		           dGLit_allow, epuntzz, zname, iwr_filt)
!		END IF
      		Dtany=Dtsegons/FACTEMP
          DEALLOCATE(Flitos)  

	   IF(iw_erosion==1 .OR. iflexure==1) THEN
		nn_tao=nn-1
		ndrainage=nresults
		CALL call_uhuru_tisc (iw_erosion, iflexure, iwrite, 
     +			  Dtsegons, m,n, nn, nn_tao, rosed, roc,
     +			  roalfa, RHOAST,  elevation, elevation0, 
     +			  sediment, scrust, scrust0, erosion, 
     +			  ndrainage, Dx, Dy, hmantle0, GLit_ter, 
     +			  Tmoho, TISOTER, ZASTH)

	   ENDIF
     
        IF(Dtany.LT.Dtmin) THEN
            PRINT*,'       Dtany = ',Dtany,' anys'
            PRINT*,'INCREMENT DE TEMPS MASSA PETIT, ATURO EL PROGRAMA'
            GOTO 9875
        ENDIF   
C ---------------------------------------------------------------------
C PER TORNAR A FER UN INCREMENT DE TEMPS INICIALITZO LA TEMPERATURA
        ALLOCATE(TEMPE1(nn,0:NELZ),STAT=istat11)
	TEMPE1=TEMPE2	! re-initialize the temperature for the next time step
        DEALLOCATE(TEMPE2)
         
          IF(iwrite==1.AND.npasos/=npasmax) THEN
              IF(nres_dec==48) THEN
                file_result='resultats'//CHAR(nres_unit)
                file_eros='res_erosion'//CHAR(nres_unit)
                file_rot='res_rotation'//CHAR(nres_unit)
                file_TK='res_Tempe'//CHAR(nres_unit)
              ELSE
                file_result='resultats'//CHAR(nres_dec)//CHAR(nres_unit)
                file_eros='res_erosion'//CHAR(nres_dec)//CHAR(nres_unit)
                file_rot='res_rotation'//CHAR(nres_dec)//CHAR(nres_unit)
                file_TK='res_Tempe'//CHAR(nres_dec)//CHAR(nres_unit)
             ENDIF
             
             IF(iw_erosion==1) THEN  
                 OPEN(30,FILE=file_eros)
                   WRITE(30,"('# Time:',F6.2,' My.   x[m], y[m],',
     +			' erosion [m]')") TEMPSMa
      		   DO iy=0,m
		     DO ix=0,n
			kxy=ix+1+iy*(n+1)
                        WRITE(30,"(3F15.3)") ix*Dx,iy*Dy,erosion(kxy)
		     END DO
		   END DO
		 CLOSE (30)
             ENDIF
             IF(Ifile_P==0) THEN  
                 OPEN(30,FILE=file_rot)
                   WRITE(30,"('# Time:',F6.2,' My.   x[km], y[km],',
     +			' Finite rotation [degree]')") TEMPSMa
      		   DO kp=1,NumPoints
                        WRITE(30,"(3F15.3)") Points_t(kp,1)/1.D3,
     + 			    Points_t(kp,2)/1.D3,fin_rotat(kp)
		   END DO
		 CLOSE (30)
             ENDIF

             OPEN(25,FILE=file_TK)
             DO ix=0,n
		 IF(ix==ix_write_T) THEN
		     WRITE(25,"(' # ',2F12.3,' x')") TEMPSMa,ix*Dx/1.D3
		     DO iy=0,m
		        kxy=ix+1+iy*(n+1)
			WRITE(25,"(1X,3F12.3)") 
     +				(iy*Dy,iz*Dz,TEMPE1(kxy,iz), iz=0,NELZ)
		     END DO
		 ENDIF		    
             END DO
             CLOSE(25)

              CALL WRRESULT (file_result, TITLE, TEMPSMa, 
     +                         npasos, m, n, NELZ, nn, FACTEMP, FACVEL,
     +                         AX, BY, Dz, TLITOS, scrust, GLit_ter, 
     +                         GLit, visco, u, v, epuntzz, Tmoho,
     +                         Qsurface, elevation, sediment,
     +			       kpuntssd, D_Bodies, Points_t)
     
              IF(npasos==NPASKEEP.AND.npasos/=npasmax) THEN
                 	NPASGRUP=NPASGRUP+1
                  	NPASKEEP=INTERVWK*NPASGRUP
              ENDIF
          ENDIF

      IF(npasos.LT.npasmax) GOTO 1000 

 9875     CONTINUE

          IF(nres_dec==48) THEN
	     file_result='resultats'//CHAR(nres_unit)
	     file_eros='res_erosion'//CHAR(nres_unit)
             file_rot='res_rotation'//CHAR(nres_unit)
	     file_TK='res_Tempe'//CHAR(nres_unit)
          ELSE
	     file_result='resultats'//CHAR(nres_dec)//CHAR(nres_unit)
	     file_eros='res_erosion'//CHAR(nres_dec)//CHAR(nres_unit)
             file_rot='res_rotation'//CHAR(nres_dec)//CHAR(nres_unit)
             file_TK='res_Tempe'//CHAR(nres_dec)//CHAR(nres_unit)
          ENDIF
          IF(iw_erosion==1) THEN  
            OPEN(30,FILE=file_eros)
      		DO iy=0,m
		   DO ix=0,n
			kxy=ix+1+iy*(n+1)
                        WRITE(30,"(3F15.3)") ix*Dx,iy*Dy,erosion(kxy)
		   END DO
		END DO
	    CLOSE (30)
          ENDIF
          IF(Ifile_P==0) THEN  
              OPEN(30,FILE=file_rot)
                WRITE(30,"('# Time:',F6.2,' My.   x[km], y[km],',
     +			'  Finite rotation [degree]')") TEMPSMa
      		DO kp=1,NumPoints
                    WRITE(30,"(3F15.3)") Points_t(kp,1)/1.D3,
     + 			  Points_t(kp,2)/1.D3,fin_rotat(kp)
		END DO
	      CLOSE (30)
          ENDIF
          OPEN(25,FILE=file_TK)
	      DO ix=0,n
		 IF(ix==ix_write_T) THEN
		     WRITE(25,"(' # ',2F12.3,' (x)')") TEMPSMa,ix*Dx/1D3
		     DO iy=0,m
		        kxy=ix+1+iy*(n+1)
			WRITE(25,"(1X,3F12.3)") 
     +				(iy*Dy,iz*Dz,TEMPE1(kxy,iz), iz=0,NELZ)
		     END DO
		 ENDIF		    
	      END DO
          CLOSE(25)
          CALL WRRESULT (file_result, TITLE, TEMPSMa, npasos,
     +                   m, n, NELZ, nn, FACTEMP, FACVEL, AX, BY,
     +                   Dz, TLITOS, scrust, GLit_ter, GLit, visco,
     +                   u, v, epuntzz, Tmoho, Qsurface, elevation, 
     +			 sediment, kpuntssd, D_Bodies, Points_t)

 9876     CONTINUE 

      WRITE(6,"(/5X,'END OF UHURU')")
      STOP
      END PROGRAM UHURU
 
C **********************************************************************
C               SUBROUTINE GRUIXOSINI 
C AQUESTA SUBRUTINA CONSTRUEIX UNA MALLA I ASSIGNA A CADA PUNT
C EL GRUIX CORTICAL I LITOSFERIC 
C GLit : gruix Litosfera. Profunditat a la que canvia la conductivitat.
C scrust : gruix cortical.
  
      SUBROUTINE GRUIXOSINI (m, n, nn, AX, BY, scrust, GLit, sediment)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*55 zname
      DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: XCORD,YCORD
      DIMENSION scrust(nn),GLit(nn),sediment(nn),x_pol(100),y_pol(100)
       
      ALLOCATE(XCORD(nn)) 
      ALLOCATE(YCORD(nn)) 
	XCORD=0.D0
	YCORD=0.D0
!	scrust=0.D0
!	GLit=0.D0
!	sediment=0.D0
      NDATA=0
      DO kxy=1,nn
            READ(11,*,IOSTAT=IOS) XCORD(kxy),YCORD(kxy),
     +             		  scrust(kxy),GLit(kxy)
!     +             		  scrust(kxy),GLit(kxy),sediment(kxy)
            IF (IOS < 0) THEN
               GO TO 120
            ELSE IF (IOS > 0) THEN
               WRITE (*, 90) NDATA
 90            FORMAT (/' ERROR: Encountered bad data after ',
     +                  I6,' successful READs.')
               STOP
            ENDIF
            NDATA=NDATA+1
      END DO
      
      IF (NDATA==0) PRINT*,' FITXER NO LLEGIT. NDATA:',NDATA
 120  WRITE(6,"(/5X,'Reading of data completed:',I7,' data points.')")
     +			NDATA
      IF(NDATA.NE.nn) THEN
            WRITE(6,"(' The number of read points:',I7,' is not equal',
     +		    ' to (n+1)*(m+1)=',I7//'STOP PROGRAMA')") NDATA, nn
            STOP
      ENDIF
      Xmin_file=VAL_MIN(XCORD,nn)
      Xmax_file=VAL_MAX(XCORD,nn)
      Ymin_file=VAL_MIN(YCORD,nn)
      Ymax_file=VAL_MAX(YCORD,nn)
      Dx_calc=XCORD(2)-XCORD(1)
      Dy_calc=YCORD(n+2)-YCORD(1)
      
      Xmin=XCORD(1)
      Xmax=XCORD(n+1)
      Ymin=YCORD(1)
      Ymax=YCORD(1+m*(n+1))
      AX=ABS(Xmax-Xmin)
      BY=ABS(Ymax-Ymin)
      Dx=AX/n
      Dy=BY/m
      
      IF(Xmin_file/=Xmin.OR.Xmax_file/=Xmax.OR.Ymin_file/=Ymin.OR.
     +           Ymax_file/=Ymax.OR.Dx_calc/=Dx.OR.Dy_calc/=Dy) THEN
          PRINT*,'The file (x,y,s,L) is not well arranged. STOP PROGRAM'
	  STOP
      ENDIF

       WRITE(6,"(5X,'Xmin:',F10.3,' km    Xmax:',F10.3,' km',/
     +		 5X,'Ymin:',F10.3,' km    Ymax:',F10.3,' km',/
     +           5X,'AX:  ',F10.3,' km    BY:  ',F10.3,' km',/
     +           5X,'Dx:  ',F10.3,' km    Dy:  ',F10.3,' km')")
     +			Xmin/1.D3,Xmax/1.D3,Ymin/1.D3,Ymax/1.D3,
     +			AX/1.D3,BY/1.D3,Dx/1.D3,Dy/1.D3

!!!!!!!  Change crust thickness inside a polygone   !!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      OPEN(54,FILE='fort.54',STATUS='OLD',ACTION='READ',IOSTAT=Ifile_P)
      IF(Ifile_P==0) THEN  
	 PRINT*,'Change crust and lithosphere thick inside the poligon'
	 PRINT*,'Add to the crustal thickness (metres) ?'
	 READ*,Add_crust
	 PRINT*,'Add to the lithosphere thickness (metres) ?'
	 READ*,Add_lithos     
	 !OPEN(54,STATUS='OLD',ACTION='READ')
	 npoints=0
	 DO i=1,100
	     READ(54,*,END=440) x,y
	     npoints=npoints+1
	     x_pol(i)=x
	     y_pol(i)=y
	 END DO
 440  	 DO iy=0,m
	    PY=iy*Dy
	    DO ix=0,n
               kxy=ix+1+iy*(n+1)
               PX=ix*Dx
               CALL outin (PX,PY,npoints,x_pol,y_pol,iadentro)
               IF(iadentro==1) THEN
	       		scrust(kxy)=scrust(kxy)+Add_crust
	       		GLit(kxy)=GLit(kxy)+Add_lithos
	       END IF	
	    END DO
	 END DO
	 zname='Input crustal thickness'
	 Dif=0.05
	 ds_allow=5.D-7
	 iwr_filt=-1
      	 CALL fit_gradients (m, n, nn, Dx, Dy, Dif, 
     +	    		              ds_allow, scrust, zname, iwr_filt)
	 OPEN(55,FILE='Grid_outUhuru.in')
	   WRITE(55,"('#  Lithosphere and crust modified')")
	   WRITE(55,"(2I5,2F8.1,' 3 (n m Dx Dy INPUT_DATA)')") n,m,Dx,Dy
	   DO kxy=1,nn
		WRITE(55,"(4F10.1)") XCORD(kxy),YCORD(kxy),
     +				scrust(kxy),GLit(kxy)
	   END DO
	 CLOSE (55)
      END IF		      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      scminim=VAL_MIN(scrust,nn)
      scmaxim=VAL_MAX(scrust,nn)
      GLitmin=VAL_MIN(GLit,nn)
      GLitmax=VAL_MAX(GLit,nn)
      sedmin=VAL_MIN(sediment,nn)
      sedmax=VAL_MAX(sediment,nn)
  
      WRITE(6,75) scminim/1.D3,scmaxim/1.D3,GLitmin/1.D3,GLitmax/1.D3,
     +			sedmin,sedmax  
 75   FORMAT(3X,'Crustal thickness (scrust),      minim:',F8.3,' km,',5X,
     +        'maxim:',F8.3,' km'/
     +       3X,'Lithospheric thickness (GLit),   minim:',F8.3,' km,',5X,
     +        'maxim:',F8.3,' km'/
     +       3X,'Sediment thickness (sediment),   minim:',F9.2,' m,',5X,
     +        'maxim:',F9.2,' m')

      DEALLOCATE(XCORD)       
      DEALLOCATE(YCORD)       
      RETURN 
      END SUBROUTINE GRUIXOSINI
 
!!  **********************************************************************
!!  Temperature at steady state, 1-D. Calculation with finite differences
!!  Boundary Conditions: to the surface (z=0) and asthenosphere (z=ZASTH)

      SUBROUTINE TEMP_STEADY_DF (m, n, nn, NELZ, Dz, DBANDAKZ, CONDUC,
     +                      HSURF, HEXP, PHEAT_m, TSURF, TBOTT, scrust,
     +                      GLit_ter, GLit, sediment, TEMPE, TMOHOlim,
     +                      TLITOS, TISOTER, Tmoho, Qsurface, AX, BY,
     +                      dGLit_allow, Dif_K_L, irel_ZL, DGL_ter)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Temp_z,T_plane
      DIMENSION CONDUC(0:3),sediment(nn),scrust(nn),GLit(nn),
     +          TEMPE(nn,0:NELZ),Tmoho(nn),GLit_ter(nn),Qsurface(nn)
      CHARACTER*55 zname
      PARAMETER (iswitch=0)		! iswitch = 1	=> keep in a file K(z) and H(z)
      DATA NMAXTM/0/,HFmin/500.D8/,HFmax/0.D0/

      WRITE(6,"(/'CALCULATION OF THE GEOTHERM AT STEADY STATE')")       
      ALLOCATE(Temp_z(0:NELZ)) 
      ALLOCATE(T_plane(nn)) 
	Temp_z=0.D0
	T_plane=0.D0
      Dx=AX/n
      Dy=BY/m
      Ndepth=NELZ+1
      NZLIT1=0
      NZLIT2=0
      BANDAKZ=Dz*DBANDAKZ
      TLITOS=0.D0

      itlit=0
      Point_y: DO iy=0,m
         Point_x: DO ix=0,n
                kxy=ix+1+iy*(n+1)
		Zsedim=sediment(kxy)
             	ZMOHO=scrust(kxy)+sediment(kxy)
             	ZLITOS=GLit(kxy)
                DO iz=0,NELZ
                    Temp_z(iz)=0.D0
                END DO
                CALL Temp_steadystate_1D (Ndepth, Dz, BANDAKZ, CONDUC,
     +                    	   HSURF, HEXP, PHEAT_m, TSURF, TBOTT,
     +				   Zsedim, ZMOHO, ZLITOS, 
     +                             Temp_z, NZLIT1, NZLIT2, iswitch)

               !IF(ix==1.AND.iy==1) WRITE(25,*) ' #  50  100 100 10'	!! File to plot with Temperature_1D.job
               DO iz=0,NELZ
                    TEMPE(kxy,iz)=Temp_z(iz)
                    !IF(ix==1.AND.iy==1) WRITE(25,*) iz*Dz,Temp_z(iz)	!! File to plot with Temperature_1D.job
               END DO

         END DO Point_x
      END DO Point_y
 
 ! Control of the lateral variations of the temperature:
      nfilters=0
      dT_allow=dGLit_allow
      Dif_T=Dif_K_L
      DO iz=2,NELZ-2
	 DO kxy=1,nn
	   T_plane(kxy)=TEMPE(kxy,iz)
	 END DO	
	 zname='Temerature steadystate z'
	 !CALL smooth_horizontal (m, n, nn, Dx, Dy, T_plane)
	 iwr_filt=0
	 CALL fit_gradients (m, n, nn, Dx, Dy, Dif_T, dT_allow,
     +				T_plane, zname, iwr_filt)
         IF(iwr_filt>0) nfilters=nfilters+1
	 DO kxy=1,nn
	   TEMPE(kxy,iz)=T_plane(kxy)
	 END DO	 
      END DO
      IF(nfilters/=0) THEN      
	   WRITE(6,"(8X,'Apply the diffusive filter to ',I4,
     +		'  Temperature layers with constant depth'/
     +		8X,'maximum gradient allowed:',1P,G10.3,
     +		3X,'Constant diffusive filter applied:',1P,G10.3)")
     +		nfilters, dT_allow, Dif_T
      END IF

      DO kxy=1,nn
          DO iz=0,NELZ
              Temp_z(iz)=TEMPE(kxy,iz)
          END DO
         Zsedim=sediment(kxy)
         ZMOHO=scrust(kxy)+sediment(kxy)
	 ZLITOS=GLit(kxy)
         GLit_ter(kxy)=Depth_isotherm (TISOTER,Temp_z,NELZ,Dz)
	 GLit_ter(kxy)=MAX(GLit_ter(kxy),ZMOHO)
	 Tmoho(kxy) = Temperature_DepthZ (ZMOHO,Temp_z,NELZ,Dz)		!! Moho Temperature
		IF(Tmoho(kxy).GT.TMOHOlim) THEN
			Tmoho(kxy)=TMOHOlim
			NMAXTM=NMAXTM+1
		ENDIF
	 T0=TEMPE(kxy,0)					!! Surface Heat Flow = - K*(dT/dz)			
	 T1=TEMPE(kxy,1)
	 Qsurface(kxy)=Surface_Heat_Flow (Zsedim,ZMOHO,CONDUC,T0,T1,Dz)
	 TLITOSP = Temperature_DepthZ (ZLITOS,Temp_z,NELZ,Dz)	!! Temperature at GLit (fixed)
	 TLITOS=TLITOS+TLITOSP
	 itlit=itlit+1
      END DO
      TLITOS=TLITOS/itlit
!       TISOTER=TLITOS	! ==> GLit_ter=GLit
      IF(NZLIT2.NE.0) PRINT*,'        ZLITOS2 > ZASTH  to ',NZLIT2,
     +                       ' points of the grid'     
      IF(NZLIT1.NE.0) PRINT*,'        ZLITOS1 < ZMOHO  tp ',NZLIT1,
     +                       ' points of the grid'          

      zname='Thermal lithosphere thickness (GLit_ter)'
      !CALL smooth_horizontal (m, n, nn, Dx, Dy, GLit_ter)
      iwr_filt=-1
      CALL fit_gradients (m, n, nn, Dx, Dy, Dif_K_L, dGLit_allow,
     +				GLit_ter, zname, iwr_filt)
      
      IF(irel_ZL==2) GLit=GLit_ter+DGL_ter
      IF(irel_ZL==1) THEN
	   DO kxy=1,nn
	      DO iz=0,NELZ
		  Temp_z(iz)=TEMPE(kxy,iz)
	      END DO
	      GLit(kxy)=Depth_isotherm(TLITOS,Temp_z,NELZ,Dz)
	   END DO
	   zname='Lithosphere thickness (GLit)'
	   !CALL smooth_horizontal (m, n, nn, Dx, Dy, GLit)
	   iwr_filt=-1
	   CALL fit_gradients (m, n, nn, Dx, Dy, Dif_K_L, dGLit_allow,
     +				GLit, zname, iwr_filt)
      ENDIF
      DO kxy=1,nn		! lithospheric thickness > crustal + sediment thickness
	 ZMOHO=scrust(kxy)+sediment(kxy)
	 GLit_ter(kxy)=MAX(GLit_ter(kxy),ZMOHO)
	 GLit(kxy)=MAX(GLit(kxy),ZMOHO)
      END DO

      GLtermin=VAL_MIN(GLit_ter,nn)
      GLtermax=VAL_MAX(GLit_ter,nn) 
      GLitmin= VAL_MIN(GLit,nn)    
      GLitmax= VAL_MAX(GLit,nn)    
      Tmmin= VAL_MIN(Tmoho,nn)    
      Tmmax= VAL_MAX(Tmoho,nn)
      Qsurface_min=VAL_MIN(Qsurface,nn)
      Qsurface_max=VAL_MAX(Qsurface,nn)
         
      IF(NMAXTM.NE.0) WRITE(6,182) TMOHOlim,NMAXTM
 182  FORMAT(10X,'The moho temperature has been limited to',F8.1,' K'/
     +   	10X,' in',I5,' points')

      WRITE(6,365) TISOTER,GLtermin/1.D3,GLtermax/1.D3,TLITOS,
     +     GLitmin/1.D3,GLitmax/1.D3,Tmmin,Tmmax,Qsurface_min*1.D3,
     +     Qsurface_max*1.D3
 365  FORMAT(3X,'GLit_ter (isotherm depth',F7.1,' K),  minima:',F9.4,
     +  		' km,    maxima:',F9.4,' km'/
     +        3X,'GLit (isotherm depth',F7.1,' K),      minima:',F9.4,
     +			' km,    maxima:',F9.4,' km'/
     +        3X,'Moho Temperature,                    minima:',F9.3,
     +			' K,     maxima:',F9.3,' K'/
     +	      3X 'Surface Heat Flow,                   minim:',F8.2,
     +			' mW/m2,   maxim:',F8.2,' mW/m2')
       
      DEALLOCATE(Temp_z)
      DEALLOCATE(T_plane)
      RETURN
      END SUBROUTINE TEMP_STEADY_DF


C **********************************************************************
      SUBROUTINE TEMPEDT (TEMPSMa, Dtsegons, m,n,nn, Dz, NELZ, DBANDAKZ,
     +                CONDUC, THDIFF, HSURF, HEXP, PHEAT_m, TLITOS,
     +                elevation, sediment, scrust, GLit, TEMPE1, TEMPE2,
     +                epuntzz, u, v, AX,BY, GLit_ter, Tmoho, TISOTER, 
     +		      TMOHOlim, Qsurface, dGLit_allow, Dif_K_L, irel_ZL,
     +		      DGL_ter, Iremoval, Time_removal, Zremoval,
     +		      TISO_rem, Zctall, Niremoved, A_removal, B_removal)
!!     New geotherm, after a Dt
!!  Niremoved : Number of geotherms that have been convective removed
     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TADVEC, Temp_z,
     +					T_k1,T_k2,T_k3,T_k4,T_plane
      DIMENSION CONDUC(0:3),epuntzz(nn),u(nn),v(nn),sediment(nn),
     +          elevation(nn),scrust(nn),Tmoho(nn),Qsurface(nn),
     +		TEMPE1(nn,0:NELZ),TEMPE2(nn,0:NELZ),
     +		GLit(nn),GLit_ter(nn),x_pol(100),y_pol(100),
     +		A_removal(nn), B_removal(nn) 
      CHARACTER*55 zname
      PARAMETER(FACTEMP=3.1536D7, TETA=0.5D0, Tmoho_min=300.D0, 
     +		zero=1.D-20)
      PARAMETER (iswitch=0)	! iswitch = 1  => keep in a file K(z) and H(z)

       WRITE(6,"(/'CALCULATION OF THE NEW GEOTHERM')") 	
	
      ALLOCATE(Temp_z(0:NELZ)) 
      ALLOCATE(T_k1(0:NELZ)) 
      ALLOCATE(T_k2(0:NELZ)) 
      ALLOCATE(T_k3(0:NELZ)) 
      ALLOCATE(T_k4(0:NELZ)) 
      ALLOCATE(TADVEC(NELZ)) 
      ALLOCATE(T_plane(nn)) 
	Temp_z=0.D0
	T_k1=0.D0
	T_k2=0.D0
	T_k3=0.D0
	T_k4=0.D0
	TADVEC=0.D0
	T_plane=0.D0
      NZLIT1=0
      NZLIT2=0

      DtMy=(Dtsegons/FACTEMP)*1.D-6
      

      ISwich_startR=0
      Time_int=ABS(TEMPSMa-Time_removal)		
      IF(Time_int<1D-3 .AND. Iremoval/=0) ISwich_startR=1

      ISwich_LongR=0
  !   t2=Time_removal+10.D0					!! Pluma mantelica durante 10My
      IF(Iremoval==3) t2=1.D4					!! Continuos Convective removal
      IF(TEMPSMa>Time_removal .AND. TEMPSMa<t2 .AND. Iremoval/=0) 
     +		ISwich_LongR=1


!      t1=Time_removal-1D-3
!      t2=Time_removal+1D-3
!  !   t2=Time_removal+10.D0					!! Pluma mantelica durante 10My
!      IF(Iremoval==3) t2=1.D4					!! Continuos Convective removal
!      IF(TEMPSMa>t1 .AND. TEMPSMa<t2 .AND. Iremoval/=0) ISwich_startR=1
      IF(ISwich_startR==1 .AND. Iremoval==22) THEN		!! read from a file the area of lithosphere removal
	 OPEN(44,STATUS='OLD',ACTION='READ')
	    npoints=0
	    DO i=1,100
		READ(44,*,END=440) x,y
		npoints=npoints+1
		x_pol(i)=x
		y_pol(i)=y
	    END DO
 440  	    WRITE(6,"(3X,'Read file where the lithosphere ',
     + 			'will be removed:',I4)") npoints
      END IF		

      Ndepth=NELZ+1
      Dx=AX/n
      Dy=BY/m
      ddx=2.D0*Dx
      ddy=2.D0*Dy
      ZASTH=Dz*NELZ
      BANDAKZ=Dz*DBANDAKZ
!      OPEN(25,FILE='T.tz')
!      OPEN(26,FILE='T2.tz')
      Niremoved=0
      Point_iy: DO iy=1,m-1
         Point_ix:  DO ix=1,n-1
                  kxy=ix+1+iy*(n+1)
	          Zsedim=sediment(kxy)
                  ZMOHO=scrust(kxy)+sediment(kxy)
                  ZLITOS=GLit(kxy)
                    k1=ix+1+(iy+1)*(n+1)
                    k2=ix+1+(iy-1)*(n+1)
                    k3=(ix+1)+1+iy*(n+1)
                    k4=(ix-1)+1+iy*(n+1)

                  DO iz=0,NELZ
                      T_k1(iz)=TEMPE1(k1,iz)
                      T_k2(iz)=TEMPE1(k2,iz)
                      T_k3(iz)=TEMPE1(k3,iz)
		      T_k4(iz)=TEMPE1(k4,iz)
                  END DO
		    
             Temp_z(0)=TEMPE1(kxy,0)
             Advective: DO iz=1,NELZ-1
		     Temp_z(iz)=TEMPE1(kxy,iz)
	      		zdepth=iz*Dz-elevation(kxy)

			ziz_k1=zdepth+elevation(k1)
			TEMPE1k1=Temperature_DepthZ(ziz_k1,T_k1,NELZ,Dz)
			IF(ziz_k1<=0) TEMPE1k1=T_k1(0)
			IF(ziz_k1>=ZASTH) TEMPE1k1=T_k1(NELZ)
     			
			ziz_k2=zdepth+elevation(k2)
                	TEMPE1k2=Temperature_DepthZ(ziz_k2,T_k2,NELZ,Dz)
			IF(ziz_k2<=0) TEMPE1k2=T_k2(0)
			IF(ziz_k2>=ZASTH) TEMPE1k2=T_k2(NELZ)     
     
			ziz_k3=zdepth+elevation(k3)
                	TEMPE1k3=Temperature_DepthZ(ziz_k3,T_k3,NELZ,Dz)
			IF(ziz_k3<=0) TEMPE1k3=T_k3(0)
			IF(ziz_k3>=ZASTH) TEMPE1k3=T_k3(NELZ)
     
			ziz_k4=zdepth+elevation(k4)
                	TEMPE1k4=Temperature_DepthZ(ziz_k4,T_k4,NELZ,Dz)
			IF(ziz_k4<=0) TEMPE1k4=T_k4(0)
			IF(ziz_k4>=ZASTH) TEMPE1k4=T_k4(NELZ)
			
			!dTdx=(TEMPE1k3-TEMPE1k4)/ddx		!! Derivada centrada
     			!dTdy=(TEMPE1k1-TEMPE1k2)/ddy		!! Derivada centrada


     			sign_x=SIGN(1.D0,u(kxy))			!! Derivada segons el fluxe del fluid		
     			sign_y=SIGN(1.D0,v(kxy))			!! Derivada segons el fluxe del fluid
     			dTdx_sup=(TEMPE1k3-Temp_z(iz))/Dx
     			dTdx_inf=(Temp_z(iz)-TEMPE1k4)/Dx
     			dTdy_sup=(TEMPE1k1-Temp_z(iz))/Dy
     			dTdy_inf=(Temp_z(iz)-TEMPE1k2)/Dy			
     			dTdx=dTdx_inf*MAX(sign_x,0D0)-dTdx_sup*MIN(sign_x,0D0)
     			dTdy=dTdy_inf*MAX(sign_y,0D0)-dTdy_sup*MIN(sign_y,0D0)			
			
                	TADVECX=u(kxy)*dTdx
                	TADVECY=v(kxy)*dTdy
                	TADVEC(iz)=-TADVECX-TADVECY
                !	TADVEC(iz)=0.D0					!! Nul advective term
             END DO Advective
	      
             Temp_z(NELZ)=TEMPE1(kxy,NELZ)
	       
!	      TADVECmaxP=VAL_MAX(TADVEC,nn)
!	      TADVECminP=VAL_MIN(TADVEC,nn)

	     HSURF_fe=HSURF
	     HEXP_fe=HEXP
	     IF(ZMOHO>40.D3) HEXP_fe=55.D3	!!! Tibetan Plateau: Heat Production function of ZMOHO - TECTONICS 2006
	     !!topo=elevation(kxy)
	     !!IF(topo>3500.0) HEXP_fe=HEXP*topo/3.5D3		!! Heat Production function of elevation
	     !!IF(topo>5000.0) HEXP_fe=HEXP*5.D3/3.5D3		!! Heat Production function of elevation
	     !!!IF(topo<0.0.AND.topo>-500.0) HSURF_fe=HSURF/1.5	!! Heat Production function of elevation
	     !!!IF(topo<=-500.0) HSURF_fe=HSURF/2.0		!! Heat Production function of elevation

             strain_rate_zz=epuntzz(kxy)
             CALL Temp_Dt1D (Dtsegons, Ndepth, Dz, BANDAKZ, CONDUC, 
     +                        HSURF_fe,HEXP_fe, PHEAT_m, THDIFF, Zsedim,
     +                        ZMOHO, ZLITOS, strain_rate_zz, TADVEC,
     +                        Temp_z, NZLIT1, NZLIT2, iswitch)

		A_removal_P=A_removal(kxy)
		B_removal_P=B_removal(kxy)
		IF(ISwich_startR==1 .OR. ISwich_LongR==1 .OR.
     +				ABS(B_removal_P)>zero) THEN
		   Iremoved=0
		   Iremoval_P=Iremoval
		   IF(Iremoval==22) THEN			!! read from file
		      Iremoval_P=2
		      IF(ISwich_startR==1) THEN
			 PX=ix*Dx
			 PY=iy*Dy
			 CALL outin (PX,PY,npoints,x_pol,y_pol,iadentro)
			 IF(iadentro/=1) Iremoval_P=0
		      END IF
		   END IF
		    
		    hlitos=Depth_isotherm (TISOTER,Temp_z,NELZ,Dz)
      		    CALL Convective_Removal_Type (TEMPSMa, Time_removal,
     +			  Iremoval_P, Zremoval, TISO_rem, Zctall, ZMOHO,
     +			  hlitos, Temp_z, TISOTER, DGL_ter, NELZ, Dz,
     +			  Iremoved, A_removal_P, B_removal_P)
     		    A_removal(kxy)=A_removal_P
     		    B_removal(kxy)=B_removal_P

!		    CALL Convective_Removal_Type (TEMPSMa, DtMy,
!     +				Time_removal, DTime_removal, Iremoval_P,
!     +			     	Zremoval, TISO_rem, Zctall,ZMOHO,hlitos,
!     +				Temp_z,TISOTER,DGL_ter,NELZ,Dz,Iremoved)
		    Niremoved=Niremoved+Iremoved
		END IF	

             DO iz=0,NELZ
                TEMPE2(kxy,iz)=Temp_z(iz)
!               	IF(ix.EQ.25.AND.iy.EQ.1) 
!     +    			WRITE(25,*) iz*Dz,Temp_z(iz)
!               	IF(ix.EQ.26.AND.iy.EQ.1) 
!     +    			WRITE(26,*) iz*Dz,Temp_z(iz)
             END DO
!	     TADVECmax=MAX(TADVECmax,TADVECmaxP)
!	     TADVECmin=MIN(TADVECmin,TADVECminP)
         END DO Point_ix
      END DO Point_iy

!! Temperature to the boundaries (ix=0,n; iy=0,m): dT/dx=0 and dT/dy=0
!!	(We need them to calculate the advective term) 
       TBoundy: DO iy=1,m-1
                   kxyw=1+iy*(n+1)
                   kxye=n+1+iy*(n+1)
                   DO iz=0,NELZ
                      TEMPE2(kxyw,iz)=TEMPE2(kxyw+1,iz)
                      TEMPE2(kxye,iz)=TEMPE2(kxye-1,iz)
		   END DO  
		   Tmoho(kxyw)=Tmoho(kxyw+1) 
		   Tmoho(kxye)=Tmoho(kxye-1) 
            END DO TBoundy 
       TBoundx: DO  ix=0,n
                   kxys=ix+1
                   kxyn=ix+1+m*(n+1)
                   DO iz=0,NELZ
                      TEMPE2(kxys,iz)=TEMPE2(kxys+n+1,iz)
                      TEMPE2(kxyn,iz)=TEMPE2(kxyn-n-1,iz)
		   END DO
		   Tmoho(kxys)=Tmoho(kxys+n+1)
		   Tmoho(kxyn)=Tmoho(kxyn-n-1)
            END DO TBoundx

      IF(NZLIT2.NE.0) PRINT*,'               ZLITOS2 > ZASTH  to',
     +                       NZLIT2,' points of the grid'     
      IF(NZLIT1.NE.0) PRINT*,'               ZLITOS1 < ZMOHO  to',
     +                       NZLIT1,' points of the grid'          
 55   FORMAT(1X,1F7.2,1P,1E12.2)
      CLOSE(25)
      CLOSE(26)
 ! Control of the lateral variations of the temperature:
      nfilters=0
      dT_allow=dGLit_allow
      Dif_T=Dif_K_L
      DO iz=2,NELZ-2
	 DO kxy=1,nn
	   T_plane(kxy)=TEMPE2(kxy,iz)
	 END DO	
	 zname='Temerature z'
	 !CALL smooth_horizontal (m, n, nn, Dx, Dy, T_plane) 
	 iwr_filt=0
	 CALL fit_gradients (m, n, nn, Dx, Dy, Dif_T, dT_allow,
     +				T_plane, zname, iwr_filt)
         IF(iwr_filt>0) nfilters=nfilters+1
	 DO kxy=1,nn
	   TEMPE2(kxy,iz)=T_plane(kxy)
	 END DO	 
      END DO
      IF(nfilters/=0) THEN      
	   WRITE(6,"(8X,'Apply the diffusive filter to ',I4,
     +		'  Temperature layers with constant depth'/
     +		8X,'maximum gradient allowed:',1P,G10.3,
     +		3X,'Constant diffusive filter applied:',1P,G10.3)")
     +		nfilters, dT_allow, Dif_T
      END IF

      NMAXTM=0
      DO kxy=1,nn		! Surface heat flow, Glit and GLit_ter
          Zsedim=sediment(kxy)
          ZMOHO=scrust(kxy)+sediment(kxy)
          T0=TEMPE2(kxy,0)
          T1=TEMPE2(kxy,1)
          Qsurface(kxy)=Surface_Heat_Flow (Zsedim,ZMOHO,CONDUC,
     +                                            T0,T1,Dz)		
          DO iz=0,NELZ
              Temp_z(iz)=TEMPE2(kxy,iz)
          END DO
	  Tmoho(kxy) = Temperature_DepthZ (ZMOHO,Temp_z,NELZ,Dz)	  
               IF(Tmoho(kxy)>TMOHOlim) THEN
                    Tmoho(kxy)=TMOHOlim
                    NMAXTM=NMAXTM+1
               ENDIF
     	       IF(Tmoho(kxy).LT.Tmoho_min)
     +     	      PRINT*,'kxy:',kxy,' Tmoho(kxy):',Tmoho(kxy)
          GLit_ter(kxy)=Depth_isotherm (TISOTER,Temp_z,NELZ,Dz)
 	  GLit_ter(kxy)=MAX(GLit_ter(kxy),ZMOHO)  	    
          IF(irel_ZL==1) GLit(kxy)=Depth_isotherm(TLITOS,Temp_z,NELZ,Dz)
	  GLit(kxy)=MAX(GLit(kxy),ZMOHO)	   	    
      END DO

      zname='Thermal lithosphere thickness (GLit_ter)'
      !CALL smooth_horizontal (m, n, nn, Dx, Dy, GLit_ter)
      iwr_filt=-1
      CALL fit_gradients (m, n, nn, Dx, Dy, Dif_K_L, dGLit_allow,
     +	    			GLit_ter, zname, iwr_filt)

      IF(irel_ZL==2) GLit=GLit_ter+DGL_ter
      IF(irel_ZL==1) THEN
           zname='Lithosphere thickness (GLit)'
	   !CALL smooth_horizontal (m, n, nn, Dx, Dy, GLit)
	   iwr_filt=-1
           CALL fit_gradients (m, n, nn, Dx, Dy, Dif_K_L, dGLit_allow,
     +			   GLit, zname, iwr_filt)
      ENDIF

      DO kxy=1,nn		! lithospheric thickness > crustal + sediment thickness
	 ZMOHO=scrust(kxy)+sediment(kxy)
	 GLit_ter(kxy)=MAX(GLit_ter(kxy),ZMOHO)
	 GLit(kxy)=MAX(GLit(kxy),ZMOHO)
      END DO

      GLtermin=VAL_MIN(GLit_ter,nn)
      GLtermax=VAL_MAX(GLit_ter,nn) 
      GLitmin= VAL_MIN(GLit,nn)    
      GLitmax= VAL_MAX(GLit,nn)    
      Tmmin= VAL_MIN(Tmoho,nn)    
      Tmmax= VAL_MAX(Tmoho,nn)
      Qsurface_min=VAL_MIN(Qsurface,nn)
      Qsurface_max=VAL_MAX(Qsurface,nn)

       IF(NMAXTM.NE.0) WRITE(6,182) TMOHOlim,NMAXTM
 182  FORMAT(/10X,'La Temperatura a la Moho ha estat limitada per ',/
     +   10X,'la Temeratura limit fixada a:',F8.1,' K, en',I4,' punts'/)
	      
      WRITE(6,"(3X,'GLit_ter (',F7.1,' K ),',2X,
     +       'minima:',F9.4,' km',6X,' maxima:',F9.4,' km'/
     +	        3X,'GLit    (',F7.1,'K ), ',3X,
     +       'minima:',F9.4,' km',6X,' maxima:',F9.4,' km'/
     +          3X,'Temperatura moho,       minima:',F8.2,' K,',8X
     +       'maxima:',F8.2,' K'/
     +		3X,'Surface Heat Flow,      minim:',F9.4,
     +       		' mW/m2',4X,' maxim:',F9.4,' mW/m2')") 
     +		TISOTER,GLtermin/1.D3,GLtermax/1.D3,
     +          TLITOS,GLitmin/1.D3,GLitmax/1.D3,
     +		Tmmin,Tmmax,Qsurface_min*1.D3,Qsurface_max*1.D3
      IF(GLitmin<=0.OR.GLitmax<=0) THEN
      	 WRITE(6,"(3X,'Negative lithospheric thickness (GLit)...',
     +		' Problem with the filter?'/' Stop the program')")
         STOP
      ENDIF
      IF(GLtermin<=0.OR.GLtermax<=0) THEN
      	 WRITE(6,"(3X,'Negative thermal lithospheric thickness ',
     +   '(GLit_ter)... Problem with the filter?'/' Stop the program')")
         STOP
      ENDIF
      
      IF(ALLOCATED(T_k1)) DEALLOCATE(T_k1)     
      IF(ALLOCATED(T_k2)) DEALLOCATE(T_k2)     
      IF(ALLOCATED(T_k3)) DEALLOCATE(T_k3)     
      IF(ALLOCATED(T_k4)) DEALLOCATE(T_k4)     
      IF(ALLOCATED(TADVEC)) DEALLOCATE(TADVEC)     
      IF(ALLOCATED(Temp_z)) DEALLOCATE(Temp_z)
      IF(ALLOCATED(T_plane)) DEALLOCATE(T_plane)
      RETURN
      END SUBROUTINE TEMPEDT

	  
CCC ********************************************************************
CCC       SUBROUTINE HEATFLOW_TOPO 
CCC        DADES D'ENTRADA:  x, y, elevation, surface Heat Flow. [S.I.]
CCC	   DADES DE SORTIDA: GRUIX CORTICAL I LITOSFERIC
C **********************************************************************
C  A partir de la topografia i el flux de calor superficial itera fins
C  trobar el gruix litosferic i cortical.

      SUBROUTINE HEATFLOW_TOPO (m, n, nn, Dz, NELZ, AX, BY,
     +                          roc, roalfa, RHOAST, RHOH2O, 
     +                          CONDUC, HSURF, HEXP, PHEAT_m, 
     +                          Topo, scrust, GLit, GLit_ter, TEMPE,
     +                          TMOHOlim, TSURF, TBOTT, TISOTER, TLITOS,
     +                          Qsurface, Tmoho)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XCORD, YCORD, 
     +				Topo_error
      DIMENSION CONDUC(0:3),Topo(nn),scrust(nn),GLit_ter(nn),GLit(nn),
     +           TEMPE(nn,0:NELZ),Qsurface(nn), Tmoho(nn)
      PARAMETER(tallTopo=0.001D0,Qlsea=40.D-3,scrust_ini=30.D3,
     +			crlimit=80.D3,crlmin=4.D3)
        
      ALLOCATE(XCORD(nn)) 
      ALLOCATE(YCORD(nn)) 
      ALLOCATE(Topo_error(nn)) 
	XCORD=0.D0
	YCORD=0.D0
	Topo_error=0.D0
      iteLQmax=15
            
      NDATA=0
      Lectura: DO kxy=1,nn
                 !scrust(kxy)=scrust_ini
                 READ(11,*,END=120) XCORD(kxy),YCORD(kxy),
     +                                 Topo(kxy),Qsurface(kxy)
                 NDATA=NDATA+1
                 s=scrust_ini+6.0*Topo(kxy)
		 scrust(kxy)=MIN(s,crlimit)
		 scrust(kxy)=MAX(s,crlmin)
		 !IF(Qsurface(kxy)<0.050) scrust(kxy)=20.D3
      END DO Lectura
 120  WRITE(6,"(/5X,'Reading of data completed:',I7,' data points.'/)")
     +			NDATA

      IF(NDATA>nn) THEN
	   WRITE(6,"('Number of data points read:',I7,
     +			'is higher than nn',I7,/
     +			'Increase n and m from file fort.11',/
     +			'STOP PROGRAM')") NDATA,nn
	   STOP
      ENDIF

C ---------  CONTROLO ELS VALORS MAXIMS I MINIMS  --------------
C --------  I QUE EL FITXER ESTIGUI BEN ORDENAT  ------------
      NQlimit=0
      DO kxy=1,NDATA       
             IF(Topo(kxy)>0.D0) THEN
                   !Qlimit=0.053D0 	     
                   Qlimit=0.05D0+0.01D0*(Topo(kxy)/7.D3) 
              ELSE
                   Qlimit=Qlsea
             ENDIF
             IF(Qsurface(kxy)<Qlimit) THEN 
                     Qsurface(kxy)=Qlimit
                     NQlimit=NQlimit+1  
             ENDIF 
      END DO  
      Qminim=VAL_MIN(Qsurface,NDATA)
      Qmaxim=VAL_MAX(Qsurface,NDATA)
      Topomin=VAL_MIN(Topo,NDATA)
      Topomax=VAL_MAX(Topo,NDATA)

      IF(NQlimit.NE.0) WRITE(6,"(/5X,'HEAT FLOW LIMITED TO ',F7.2,
     +      ' mW.m-2  IN ',I5,' POINTS')") Qlimit*1.D3,NQlimit
   
      WRITE(6,"(3X,'Surface Heat Flow, minimum:',F7.2,' mW/m2',3X,
     +          ' maximum:',F7.2,' mW/m2'/
     +          3X,'Elevation,         minima:',F8.1,' m',
     +          7X,' maxima:',F8.1,' m')") 
     +          Qminim*1.D3,Qmaxim*1.D3,Topomin,Topomax

      Xmin_file=VAL_MIN(XCORD,NDATA)
      Xmax_file=VAL_MAX(XCORD,NDATA)
      Ymin_file=VAL_MIN(YCORD,NDATA)
      Ymax_file=VAL_MAX(YCORD,NDATA)
      Dx_calc=XCORD(2)-XCORD(1)
      Dy_calc=YCORD(n+2)-YCORD(1)
      
      Xmin=XCORD(1)
      Xmax=XCORD(n+1)
      Ymin=YCORD(1)
      Ymax=YCORD(1+m*(n+1))
      AX=ABS(Xmax-Xmin)
      BY=ABS(Ymax-Ymin)
      Dx=AX/n
      Dy=BY/m
      
      WRITE(6,"(3X,'Xmin:',F10.3,' km    Xmax:',F10.3,' km',/
     +		3X,'Ymin:',F10.3,' km    Ymax:',F10.3,' km',/
     +		3X,'AX:',F9.3,' km       BY: ',F9.3,' km',/
     +		3X,'Dx:',F9.3,' km       Dy: ',F9.3,' km')") 
     +			Xmin/1.D3,Xmax/1.D3,Ymin/1.D3,Ymax/1.D3,
     +			AX/1.D3,BY/1.D3,Dx/1.D3,Dy/1.D3
        
      iteLQ=1
C     PRINT*,' Surface Heat production function of the elevation'
      WRITE(6,"(/3X,'ITERATION     Delta Topo')")
 10   CONTINUE

	  CALL Temp_TQsurf (n,m,NDATA,Dz,NELZ,CONDUC,HSURF,HEXP,PHEAT_m,
     +                      TSURF, TBOTT, TISOTER, TLITOS, Topo, scrust,
     +                      GLit,GLit_ter,Qsurface,TEMPE,Tmoho,TMOHOlim)

	  CALL crustal_thickness (iteLQ, NDATA, TISOTER, Tmoho, scrust,
     +					GLit_ter, roc, roalfa, RHOAST,
     +					RHOH2O, Topo, DTopo, Topo_error)                    
            
	  WRITE(6,"(I8,4X,F15.4)") iteLQ,DTopo
                iteLQ=iteLQ+1
               IF(DTopo.LE.tallTopo) THEN
		    WRITE(6,"(/' The elevation is converged')")
                    GOTO 99
               ENDIF     
               IF(iteLQ.GE.iteLQmax) THEN 
		    WRITE(6,"(' maxima iteration')")
                    GOTO 99
                 ELSE
                    GOTO 10 
               ENDIF
 99   CONTINUE

      crustmin=VAL_MIN(scrust,NDATA)
      crustmax=VAL_MAX(scrust,NDATA)
      zisomin=VAL_MIN(GLit_ter,NDATA)
      zisomax=VAL_MAX(GLit_ter,NDATA)
      WRITE(6,"(/3X,'Crustal thickness,       minim:',F9.1,' m',4X, 
     +       'maxim: ',F10.1,' m'/
     +           3X,'Lithospheric thickness,  minim:',F9.1,' m',4X,
     +       'maxim: ',F10.1,' m'/)") crustmin,crustmax,zisomin,zisomax

      IF(Xmin_file/=Xmin.OR.Xmax_file/=Xmax.OR.Ymin_file/=Ymin.OR.
     +   Ymax_file/=Ymax.OR.Dx_calc/=Dx.OR.Dy_calc/=Dy.OR.NDATA/=nn.OR.
     +		Xmin/=0.OR.Ymin/=0)THEN
!      IF(Xmin_file/=Xmin.OR.Xmax_file/=Xmax.OR.Ymin_file/=Ymin.OR.
!     +   Ymax_file/=Ymax.OR.Dx_calc/=Dx.OR.Dy_calc/=Dy.OR.NDATA/=nn)THEN
           PRINT*,'  The file (x,y,e,Q) is not well arranged or '
           PRINT*,' the number of data points read:',NDATA,
     +			'is different to (n+1)*(m+1)=',nn
           PRINT*,'  Saving in a file the Steady State and STOP PROGRAM'
           OPEN(40,FILE='Qe_sL_steady.res')
           WRITE(40,"('#  x, y, Crust.thick(km), Litos.thick(km), ',
     +		'Tmoho(K), Qlitos(mW/m2), Qsurface(mW/m2), topo(m), ',
     +		'error (data-calculated, m)')")
           DO k=1,NDATA          
	       QMOHO=Qsurface(k)-HSURF*HEXP*(1.D0-DEXP(-scrust(k)/HEXP))
	       QBASE=QMOHO-(GLit_ter(k)-scrust(k))*PHEAT_m	!! Heat flow at the base of lithosphere
	       WRITE(40,"(2F12.3,2F10.3,3F8.1,2F10.1)") XCORD(k),
     +		     YCORD(k),scrust(k)/1D3,GLit_ter(k)/1D3,Tmoho(k),
     +		     QBASE*1D3,Qsurface(k)*1D3,Topo(k),Topo_error(k)
           END DO
	   STOP
      ENDIF
                    
      IF(ALLOCATED(XCORD)) DEALLOCATE(XCORD)     
      IF(ALLOCATED(YCORD)) DEALLOCATE(YCORD)     
      IF(ALLOCATED(Topo_error)) DEALLOCATE(Topo_error)     
      RETURN 
      END SUBROUTINE HEATFLOW_TOPO
C **********************************************************************
C **********************************************************************
CC	CRITERI DE LACHENBRUCH AND MORGAN (1990)
CC Calcula el gruix cortical a partir de la topografia i el gruix 
CC   litosferic termic (GLit_ter) .
  
      SUBROUTINE crustal_thickness (iteLQ, nn, TISOTER, Tmoho,
     +                         scrust, GLit_ter,roc, roalfa, RHOAST,
     +                          RHOH2O, Topo, DTopo, Topo_error)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Topoini
      DIMENSION scrust(nn),Topo(nn),GLit_ter(nn),Tmoho(nn),
     +			Topo_error(nn)
      PARAMETER(crlimit=80.D3,crlmin=4.D3)

      ALLOCATE(Topoini(nn)) 
	Topoini=0.D0
       sedim=0.D3
       DO kxy=1,nn
          rom=RHOAST*(1.D0+(roalfa/2.D0)*(TISOTER-Tmoho(kxy)))
	  hcrust=scrust(kxy)
          hmantle=GLit_ter(kxy)-scrust(kxy)-sedim
	  Topoini(kxy)=elevation_isostasy(RHOH2O,rosed,roc,rom,RHOAST,
     +                                    sedim, hcrust, hmantle)          
       END DO            

       NLITCR=0
       NCRUSTL=0
       NCRMIN=0
       DO kxy=1,nn
             rom=RHOAST*(1.D0+(roalfa/2.D0)*(TISOTER-Tmoho(kxy)))
             hmantle=GLit_ter(kxy)-scrust(kxy)-sedim
	     elevation_P=Topo(kxy)
	     scrust(kxy)=crust_isostasy (RHOH2O, roc, rom, RHOAST,
     +                                      elevation_P, hmantle)	     
                IF(scrust(kxy).GT.crlimit) THEN
                      scrust(kxy)=crlimit
                      NCRUSTL=NCRUSTL+1
                ENDIF 
                IF(scrust(kxy).LT.crlmin) THEN
                      scrust(kxy)=crlmin
                      NCRMIN=NCRMIN+1
                ENDIF            
                IF(scrust(kxy).GT.GLit_ter(kxy)) THEN
                    GLit_ter(kxy)=scrust(kxy)
                    NLITCR= NLITCR+1
                ENDIF
       END DO
       IF(NCRUSTL.NE.0) WRITE(6,"(30X,'Crustal thickness limited to',
     +          F7.2,' km  at',I5,' points')") crlimit/1.D3,NCRUSTL
 
       IF(NCRMIN.NE.0) WRITE(6,"(30X,'Crustal thickness limited to',
     +                F7.2,' km  at',I5,' points')") crlmin/1.D3,NCRMIN

       IF(NLITCR.NE.0) WRITE(6,"(30X,'lithosphere limited by the crust',
     +                 ' at',I5,' points')") NLITCR
 
       NDATA=0
       DTopo=0.D0
       DO kxy=1,nn
C	    DTopo=DTopo+DABS((Topo(kxy)-Topoini(kxy))/Topo(kxy))
	    Topo_error(kxy)=Topo(kxy)-Topoini(kxy)
	    DTopo=DTopo+DABS(Topo_error(kxy))
	    NDATA=NDATA+1 
       END DO  
       DTopo=DTopo/NDATA

      IF(ALLOCATED(Topoini)) DEALLOCATE(Topoini)     
      RETURN 
      END SUBROUTINE crustal_thickness
C **********************************************************************
C **********************************************************************
C  CALCULA LA DISTRIBUCIO DE TEMPERATURES 
C   PROHEAT=HSURF*DEXP(-Z/HEXP)
C  !!!!!   No funcionara be per al llarg del temps.  !!!!!!!! 
C  !!!!!!!          No calcula T a l'astenosfera !!!!!!!!   

      SUBROUTINE Temp_TQsurf (n,m,nn, Dz, NELZ,CONDUC, HSURF, HEXP, 
     +                         PHEAT_m,TSURF, TBOTT, TISOTER, TLITOS,
     +                           Topo, scrust, GLit, GLit_ter, 
     +                           QSURF, TEMPE, Tmoho, TMOHOlim)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Temp_z
      DIMENSION CONDUC(0:3),QSURF(nn),Topo(nn),scrust(nn),GLit(nn),
     +		GLit_ter(nn),Tmoho(nn),TEMPE(nn,0:NELZ)
      PARAMETER(DepthH0=-500.D0, Tmoho_min=300.D0)
      
      iswitch=0			! iswitch = 1  => keep in a file K(z) and H(z)
      ALLOCATE(Temp_z(0:NELZ)) 
	Temp_z=0.D0
      ZASTH=Dz*NELZ
      NLIMZ=0
      NMAXTM=0
      OPEN(21,FILE='Heat_Production.res')
      WRITE(21,"(1X,'# ix  iy   HSURF   HEXP [H=HSURF*exp(-z/HEXP)] ')")
      DO ix=0,n
      !DO kxy=1,nn
        DO iy=0,m
           kxy=ix+1+iy*(n+1)
	   Qsurface=QSURF(kxy)
	   Zsedim=0.D3
	   ZMOHO=scrust(kxy)+Zsedim
	   ZLITOS=ZASTH
	   HSURF_fe=HSURF
	   HEXP_fe=HEXP
	   !IF(Topo(kxy)>4000.D0) HEXP_fe=25.D3	!HEXP*Topo(kxy)/3.5D3			!! Heat Production function of elevation
	   !IF(Topo(kxy)>6000.D0) HEXP_fe=40.D3	!HEXP*5.0D3/3.5D3			!! Heat Production function of elevation
	   !IF(Topo(kxy)<0.D0.AND.Topo(kxy)>DepthH0) HSURF_fe=HSURF/1.5D0	!! Heat Production function of elevation
	   !!IF(Topo(kxy)<=DepthH0) HSURF_fe=HSURF/2.D0				!! Heat Production function of elevation
           WRITE(21,"(2I5,2G15.6)") ix, iy, HSURF_fe, HEXP_fe
	   CALL Temp_steadystate_TQsurf_1D (NELZ, Dz, CONDUC,
     +                                 HSURF_fe, HEXP_fe, PHEAT_m,
     +			               Zsedim, ZMOHO, ZLITOS, Qsurface,
     +                                 TSURF, TLITOS, Temp_z, iswitch)
           IF(iswitch==999)  THEN
		WRITE(6,"(6X,'elevation: ',F9.2,' m',/
     +			/3X,'Error - STOP THE PROGRAM')") Topo(kxy)
		STOP
           END IF
	   GLit_ter(kxy)=Depth_isotherm (TISOTER,Temp_z,NELZ,Dz)
	   GLit(kxy)=GLit_ter(kxy)+20.D3				!! GLit	is 20 km deeper than GLit_ter
  !!	   GLit(kxy)=Depth_isotherm (TLITOS,Temp_z,NELZ,Dz)
  !!	   IF(GLit_ter(kxy)>GLit(kxy)) PRINT*,kxy,'GLit_ter > GLit'
               
	   Tmoho(kxy) = Temperature_DepthZ (ZMOHO,Temp_z,NELZ,Dz)
	   IF(Tmoho(kxy).GT.TMOHOlim) THEN
		Tmoho(kxy)=TMOHOlim
		NMAXTM=NMAXTM+1
	   ENDIF
	   IF(Tmoho(kxy).LT.Tmoho_min) THEN
		PRINT*,'kxy:',kxy,' Tmoho(kxy):',Tmoho(kxy)
!		DO 151 iz=0,NELZ
! 151			PRINT*,'TEMPE2(kxy,iz):',TEMPE2(kxy,iz)
	   ENDIF
	       
	   DO iz=0,NELZ
		TEMPE(kxy,iz)=Temp_z(iz)
	   END DO
        END DO
      END DO
      CLOSE(21)
      IF(NMAXTM.NE.0) WRITE(6,62) TMOHOlim,NMAXTM
 62   FORMAT(/10X,'La Temperatura a la Moho ha estat limitada per ',/
     +   10X,'la Temeratura limit fixada a:',F8.1,' K, en',I4,' punts'/)

      IF(NLIMZ.NE.0) WRITE(6,69) TISOTER,ZASTH/1.D3,NLIMZ
 69   FORMAT(32X,'LA PROFUNDITAT DE LA ISOTERMA',F8.1,' K',/32X,
     +       'HA ESTAT LIMITADA PER LA PROFUNDITAT',/32X,
     +       'MAXIMA',F8.2,' km EN',I5,' PUNTS.')

     
      IF (ALLOCATED(Temp_z)) DEALLOCATE(Temp_z)
      RETURN
      END SUBROUTINE Temp_TQsurf

CCC ********************************************************************
CCC              SUBROUTINE CRUST_TOPO 
CCC        DADES D'ENTRADA:  GRUIX CORTICAL I ELEVACIO
CCC	   DADES DE SORTIDA: GRUIX LITOSFERIC I FLUX DE CALOR
CCC ********************************************************************
CCC  A partir de la topografia i el gruix cortical itera fins trobar
CCC  la distribucio de temperatures i el gruix litosferic

      SUBROUTINE CRUST_TOPO (m, n, nn, Ndepth, Dz, NELZ,
     +                        AX, BY, roc, RHOAST, RHOH2O,
     +                        CONDUC, PROHEATS, HEXP, TSURF, 
     +                        TBOTT, TISOTER, scrust, 
     +                        GLit, GLit_ter, TEMPE, TLITOS, elevation)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(tallTopo=0.001D0,crust_limit=5.D3,roalfa=3.5D-5,
     +           QSUPINI=65.D-3,DELEV=25.D0,JCONVMAX=100)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: QSURF,Topo,
     +		XCORD,YCORD,TZ,CV
      DIMENSION CONDUC(0:3),PROHEAT(2),scrust(nn),GLit(nn),GLit_ter(nn),
     +           TEMPE(nn,0:NELZ),elevation(nn)
      CHARACTER*70 FILEQT
           
      ALLOCATE(QSURF(nn)) 
      ALLOCATE(Topo(nn)) 
      ALLOCATE(XCORD(nn)) 
      ALLOCATE(YCORD(nn)) 
      ALLOCATE(TZ(0:NELZ))
      ALLOCATE(CV(JCONVMAX)) 
	QSURF=0.D0
	Topo=0.D0
	XCORD=0.D0
	YCORD=0.D0
	TZ=0.D0
 	CV=0.D0
      iteLQmax=1 
        PRINT*,'   SUBROUTINE CRUST_TOPO  '
CC        OPEN(1,FILE='Delevacio.xy') 
C       FILEQT='/model/ivone/DADES/Golf/QesL/no_Mercator/es_input.xy'
CC       FILEQT='/model/ivone/DADES/Golf/QesL/Topo_crust_total.xy'
CC       FILEQT='/model/ivone/DADES/Golf/QesL/Topo_crust_2dominis.xy'
CC       FILEQT='/model/ivone/DADES/Golf/Topo_crust_Peninsula51x51.xy'
       WRITE (6,3) FILEQT
 3     FORMAT (/1X,'LLEGEIXO EL FITXER:  ',A55,/
     +    6X,'( x, y, Topografia, gruix cortical ) ',/
     +    4X,'on (x,y) es una malla regular. Tot esta en unitats S.I.')
       
C       OPEN(11,FILE=FILEQT)
      NDATA=0
      DO kxy=1,nn
            READ(11,*,IOSTAT=IOS) XCORD(kxy),YCORD(kxy),
     +                 Topo(kxy),scrust(kxy)
            IF (IOS < 0) THEN
               GO TO 120
            ELSE IF (IOS > 0) THEN
               WRITE (*, 90) NDATA
 90            FORMAT (/' ERROR: Encountered bad data after ',
     +                  I6,' successful READs.')
               STOP
            ENDIF
            NDATA=NDATA+1
      END DO

 100  CONTINUE
      IF (NDATA.EQ.0) PRINT*,' FITXER NO LLEGIT. NDATA:',NDATA
 120  CONTINUE
      WRITE(6,130)NDATA
 130  FORMAT(/5X,'Reading of data completed:',I7,' data points.'/)
      IF(NDATA.NE.nn) THEN
            WRITE(6,135)NDATA,nn
 135        FORMAT(' El numero de punts llegits:',I7,' no coincideix',
     &          'amb (n+1)*(m+1)=',I7 //'PROGRAMA ATURAT')
            STOP
      ENDIF
      CLOSE(15)
 66    FORMAT(1X,4F15.3)
CC ---------  CONTROLO ELS VALORS MAXIMS I MINIMS  --------------
        Ncrlimit=0
        scminim=500.D3
        scmaxim=0.D0
        DO 71 kxy=1,nn 
             IF(scrust(kxy).LT.crust_limit) THEN 
                     scrust(kxy)=crust_limit
                     Ncrlimit=Ncrlimit+1  
             ENDIF 
             scminim=MIN(scrust(kxy),scminim)
             scmaxim=MAX(scrust(kxy),scmaxim)
 71    CONTINUE
   
      IF(Ncrlimit.NE.0) WRITE(6,72) crust_limit,Ncrlimit
 72   FORMAT(/5X,'EL GRUIX CORTICAL HA ESTAT LIMITAT A',F9.2,
     +      ' m  EN',I4,' punts')
      WRITE(6,75) scminim, scmaxim
 75   FORMAT(3X,'GRUIX CORTICAL MINIM:',F9.2,' m',2X,
     +        'GRUIX CORTICAL MAXIM:',F9.2,' m'/)
         
        Xmin=XCORD(1)
        Xmax=XCORD(n+1)
        Ymin=YCORD(1)
        Ymax=YCORD(1+m*(n+1))
        AX=ABS(Xmax-Xmin)
        BY=ABS(Ymax-Ymin)
         Dx=AX/n
         Dy=BY/m
          WRITE(6,"(5X,'Xmin:',F10.3,' km    Xmax:',F10.3,' km',/
     +		5X,'Ymin:',F10.3,' km    Ymax:',F10.3,' km',/
     +		5X,'AX:',F9.3,' km       BY: ',F9.3,' km',/
     +		5X,'Dx:',F9.3,' km       Dy: ',F9.3,' km')") 
     +		Xmin/1.D3,Xmax/1.D3,Ymin/1.D3,Ymax/1.D3,
     +          AX/1.D3,BY/1.D3,Dx/1.D3,Dy/1.D3
        
 111          iteLQ=1
C              PRINT*,'   '
C              PRINT*,'   ITERATION     Delta Topo'
 10           CONTINUE
                                                            
CCC -------------------------------------------------------------------
CCC	PER CADA PUNT BUSCO LA DISTRIBUCIO DE TEMPERATURES 
CCC		I EL GRUIX LITOSFERIC		
CCC -------------------------------------------------------------------
CC   BUCLE PER CADA PUNT DE LA MALLA

       QSUPmin=1.D3
       QSUPmax=0.D0
       ERRImin=1.D20
       ERRImax=0.D0
       ERRImig=0.D0
       ITEmax=0
       PRINT*,' bucle 305'
       DO 305 iy=0,m
          DO 305 ix=0,n
              kxy=ix+1+iy*(n+1)
              QSUP=QSUPINI  
              JCRL=0
              NLITCR=0
      	      F=0.25D0
              ZMOHO=scrust(kxy)
              TopoP=Topo(kxy)
CC    BUCLE DE CONVERGENCIA

              DO 307 JCONV=1,JCONVMAX 
	      	    PRINT*,' ix,iy,JCONV:',ix,iy,JCONV
	      	    PRINT*,' CALL TEMP_CRUST_TOPO'
                    CALL TEMP_CRUST_TOPO (Dz, NELZ, CONDUC, PROHEATS, 
     +                                    HEXP, TSURF, TISOTER, TopoP,
     +                                    QSUP, ZMOHO, GLit_terP, TZ, 
     +                                    RHOAST, rom, roalfa, NLITCR)
                     
	      	    PRINT*,' CALL TOPO_LITOS'
                    CALL TOPO_LITOS (JCONVMAX, RHOAST, RHOH2O, roc,rom,
     +                              ZMOHO, GLit_terP, TopoP, DELEV,
     +                              QSUP, JCRL, CV, JCONV, F, ERRI, 
     +                              ERRISIG, elevatP)
	PRINT*,'RHOAST,RHOH2O,roc,rom:',RHOAST,RHOH2O,roc,rom
       PRINT*,'ZMOHO,GLit_terP,TopoP,DELEV:',ZMOHO,GLit_terP,TopoP,DELEV
	PRINT*,'QSUP,JCRL,CV(JCONV),JCONV:',QSUP,JCRL,CV(JCONV),JCONV
	PRINT*,'F,ERRI,ERRISIG,elevatP:',F,ERRI,ERRISIG,elevatP
                    IF(JCRL.EQ.1)  GOTO 308

 307          CONTINUE
               ITEmax=ITEmax+1
C              PRINT*,' ITERACIO MAXIMA '
 308          CONTINUE
             QSUPmin=MIN(QSUP,QSUPmin)
             QSUPmax=MAX(QSUP,QSUPmax)
             ERRImin=MIN(ERRI,ERRImin)
             ERRImax=MAX(ERRI,ERRImax)
             ERRImig=ERRImig+ERRI             
C              PRINT*,'Punt:(',ix,iy,')  iteracio:',JCONV,
C     +               'Diferencia elevacio:',ERRI
C              WRITE(1,*) ix*Dx,iy*Dy,ERRISIG
            IF(NLITCR.EQ.1) PRINT*,'            EL GRUIX LITOSFERIC ',
     +                             'HA LIMITAT PER EL GRUIX CORTICAL'
              GLit_ter(kxy)=GLit_terP
              GLit(kxy)=GLit_terP
              elevation(kxy)=elevatP
              DO 311 iz=0,NELZ
 311               TEMPE(kxy,iz)=TZ(iz)
                
 305   CONTINUE
         ERRImig=ERRImig/((n-1)*(m-1))

      WRITE(6,315) JCONVMAX,ITEmax,ERRImin,ERRImax,ERRImig
 315   FORMAT(/'  HA ARRIBAT AL NUMERO MAXIM DE ITERACCIONS (',I5,
     +        ')  A',I5,' PUNTS'//
     +        '  DIFERENCIA MINIMA EN LA ELEVACIO:',F9.2,' m'/
     +        '  DIFERENCIA MAXIM EN LA ELEVACIO:',F9.2,' m'/
     +        '  DIFERENCIA MITJA EN LA ELEVACIO:',F9.2,' m'/)
CCC		 ---------------------------

         crustmin=200.D3
         crustmax=0.D0
         zisomin=NELZ*Dz
         zisomax=0.D0
         DO 210 kxy=1,nn 
                crustmin=MIN(crustmin,scrust(kxy))
                crustmax=MAX(crustmax,scrust(kxy))
                zisomin=MIN(zisomin,GLit_ter(kxy))
                zisomax=MAX(zisomax,GLit_ter(kxy))
 210     CONTINUE      
      WRITE(6,65) crustmin,crustmax,zisomin,zisomax,
     +            QSUPmin*1.D3,QSUPmax*1.D3
 65   FORMAT(/3X,'GRUIX CORTICAL MINIM:  ',F9.1,' m',3X,
     +        'GRUIX CORTICAL MAXIM: ',F10.1,' m'/
     +        3X,'GRUIX LITOSFERIC MINIM:',F9.1,' m',3X,
     +        'GRUIX LITOSFERIC MAXIM:',F10.1,' m'/
     +        3X,'FLUX DE CALOR MINIM:',F7.2,' mW.m-2',5X,
     +        'FLUX DE CALOR MAXIM:',F7.2,' mW.m-2'/)

      CLOSE(1)

      IF (ALLOCATED(QSURF)) DEALLOCATE(QSURF)     
      IF (ALLOCATED(Topo)) DEALLOCATE(Topo)     
      IF (ALLOCATED(XCORD)) DEALLOCATE(XCORD)     
      IF (ALLOCATED(YCORD)) DEALLOCATE(YCORD)     
      IF (ALLOCATED(TZ)) DEALLOCATE(TZ)     
      IF (ALLOCATED(CV)) DEALLOCATE(CV)     

      RETURN 
      END SUBROUTINE CRUST_TOPO

CC ********************************************************************
CCC *******************************************************************
CCC             SUBROUTINE TEMP_CRUST_TOPO
CCC   CALCULA LA DISTRIBUCIO DE TEMPERATURES DE L'ESTAT ESTACIONARI 
CCC	EN FONDARIA. SOLUCIO ANALITICA PER UNA PRODUCCIO DE CALOR 
CCC	EXPONENCIAL A L'ESCORCA I CONSTANT AL MANTELL.
CCC   PROHEAT(a l'escorca)=HSURF*DEXP(-Z/HEXP)
CCC   COM A C.C. TINC LA TEMPERATURA I EL FLUX DE CALOR A LA SUPERFÖCIE
CCC   OUTPUT: GRUIX LITOSFERIC I DENSITAT MITJA DEL MANTELL

       SUBROUTINE TEMP_CRUST_TOPO (Dz, NELZ, CONDUC, PROHEATS, HEXP,
     +                              TSURF, TISOTER, TopoP, QSUP, ZMOHO, 
     +                              GLit_terP, TZ, RHOAST, rom, 
     +                              roalfa, NLITCR)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION CONDUC(0:3),PROHEAT(2),TZ(0:NELZ)
       PARAMETER (DepthH0=-500.D0)
       
       HSURF=PROHEATS
       IF(TopoP.LT.0.D0.AND.TopoP.GT.DepthH0) HSURF=HSURF/2
       IF(TopoP.LT.DepthH0) HSURF=HSURF/4.D0
       PROHEAT(2)=0.D0
C       Tm1=QSUP-HSURF*HEXP
C       Tm2=HEXP*HEXP*HSURF
C       TmohoP=TSURF+
C     +      (1.D0/CONDUC(1))*(Tm1*ZMOHO+Tm2*(1.D0-DEXP(-ZMOHO/HEXP)))
       TmohoP=TSURF+(1.D0/CONDUC(1))*((QSUP-HSURF*HEXP)*ZMOHO+
     +              HEXP*HEXP*HSURF*(1.D0-DEXP(-ZMOHO/HEXP)))
       QMOHO=QSUP-HSURF*HEXP*(1.D0-DEXP(-ZMOHO/HEXP))
        DO 15 iz=0,NELZ
            Z=iz*Dz
            IF(Z.LE.ZMOHO) TZ(iz)=TSURF+
     +                 (1.D0/CONDUC(1))*((QSUP-HSURF*HEXP)*Z+
     +                  HSURF*HEXP*HEXP*(1.D0-DEXP(-Z/HEXP)))
            IF(Z.GT.ZMOHO) THEN
                   Zlit=Z-ZMOHO
                   TZ(iz)=TmohoP+QMOHO*(Zlit/CONDUC(2))
     +                    -PROHEAT(2)*Zlit*Zlit/(2.D0*CONDUC(2))
                   IF(TZ(iz-1).LE.TISOTER.AND.TZ(iz).GT.TISOTER) THEN
                        izlitos=iz
                        diz=(TISOTER-TZ(iz-1))/(TZ(iz)-TZ(iz-1))
                        GLit_terP=(DBLE(iz-1)+diz)*Dz
                        IF(ZMOHO.GT.GLit_terP) THEN
                            PRINT*,'ZMOHO > GLit_ter -> ZMOHO=GLit_ter'
                            GLit_terP=ZMOHO
                            NLITCR=NLITCR+1
                        ENDIF
                        GOTO 16
                   ENDIF
            ENDIF
 15     CONTINUE
        GLit_terP=NELZ*Dz
 16    CONTINUE
       rom=RHOAST*(1.D0+(roalfa/2.D0)*(TISOTER-TmohoP))

       RETURN
       END SUBROUTINE TEMP_CRUST_TOPO
CC ********************************************************************
CC               SUBROUTINE  TOPO_LITOS
CC   COMPARA L'ELEVACIO OBSERVADA AMB LA QUE ES TINDRIA SEGONS EL 
CC	GRUIX LITOSFERIC TROBAT A PARTIR DE LA DISTRIBUCIO DE 
CC	TEMPERATURES EN FONDARIA.
  
      SUBROUTINE TOPO_LITOS (JCONVMAX, RHOAST, RHOH2O, roc, rom, 
     +                       ZMOHO, GLit_terP, TopoP, DELEV, QSUP,
     +                       JCRL,CV,JCONV,F,ERRI,ERRISIG,elevatP)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CVA
      DIMENSION CV(JCONVMAX)
      PARAMETER(gruixw=2.7D3,gruixc=5.D3,H0=2400.D0,QSUPmin_con=55.D-3,
     +           QSUPmin_lim=43.D-3,QSUPmax_lim=130.D-3)
      
      ALLOCATE(CVA(JCONVMAX))
	CVA=0.D0
      DELEV10=DELEV*30.D0
       
C       sc1=(RHOAST-RHOH2O)*gruixw
C       sc2=(RHOAST-ro1)*gruixc        
C       sc3=(rom-RHOAST)*GLit_terP
C       sc4=(rom-roc)*ZMOHO
C       FE=sc4-sc1-sc2-sc3ini
C         IF(FE.GE.0.D0) THEN
C              Topoini=FE/RHOAST
C           ELSE
C              Topoini=FE/(RHOAST-RHOH2O)
C         ENDIF

       RAW=RHOAST/(RHOAST-RHOH2O)
       EC=((RHOAST-roc)/RHOAST)*ZMOHO
       IF(TopoP.LT.0.)  TOP=TopoP/RAW
       IF(TopoP.GE.0.)  TOP=TopoP
       EMO=TOP-EC+H0
       EMC=((RHOAST-rom)/RHOAST)*(GLit_terP-ZMOHO)
       ERRISIG=EMO-EMC
       ERRI=ABS(EMO-EMC)
       CV(1)=EMO-EMC
CC  Elevacio segons els gruixos cortical i litosferic 
	rolitos=((roc*ZMOHO)+(rom*(GLit_terP-ZMOHO)))/GLit_terP
       elevatP=(((RHOAST-rolitos)/RHOAST)*GLit_terP)-H0
       IF(elevatP.LT.0.D0) elevatP=elevatP*RAW
      IF (ERRI.GT.DELEV)    THEN      
        IF (JCONV.GT.1) THEN
          CV(JCONV)=EMO-EMC
          CVA(JCONV)=CV(JCONV)-CV(JCONV-1)
          ERR1=ABS(CVA(JCONV))
          ERR2=ABS(CVA(JCONV-1))
      	  IF (CV(JCONV).GT.0.)   THEN
            IF(CVA(JCONV).GT.0.)  F=-F
            IF(CV(JCONV).GT.0. .AND. CV(JCONV-1).LT.0. .OR.
     &           CV(JCONV).LT.0. .AND. CV(JCONV-1).GT.0.) F=F/2.D0   
            IF(ERR1.GT.ERR2 .AND. ERRI.LT.DELEV10)  F=F/2.D0
            QSUP=QSUP+QSUP*F
          ELSE
            IF(CVA(JCONV).LT.0.)  F=-F 
            IF(CV(JCONV).GT.0. .AND. CV(JCONV-1).LT.0. .OR.
     &           CV(JCONV).LT.0. .AND. CV(JCONV-1).GT.0.) F=F/2.D0
            IF(ERR1.GT.ERR2 .AND. ERRI.LT.DELEV10)  F=F/2.D0
            QSUP=QSUP+QSUP*F
          ENDIF
        ENDIF
      ELSE
            JCRL=1
      ENDIF 

       IF(TopoP.GT.0.D0) THEN
             Qlimit=0.05D0+0.01D0*(TopoP/7.D3) 
         ELSE
             Qlimit=QSUPmin_lim
        ENDIF

      IF(TopoP.GE.500.D0.AND.QSUP.LE.QSUPmin_con) QSUP=QSUPmin_con
      IF(QSUP.LT.Qlimit) QSUP=Qlimit
      IF(QSUP.GT.QSUPmax_lim) QSUP=QSUPmax_lim      

      IF (ALLOCATED(CVA)) DEALLOCATE(CVA)     
      RETURN 
      END SUBROUTINE TOPO_LITOS

CCC *******************************************************************
CCC *******************************************************************
CCC                 SUBROUTINE ISOSTASY
CCCC	elevation calculation. Local isostasy

      SUBROUTINE ISOSTASY (m, n, nn, RHOH2O, rosed, roc, RHOAST,
     +			 roalfa, elevation, sediment, scrust, GLit_ter,
     +			 Tmoho, TISOTER, Dx, Dy, kpuntssd, D_Bodies)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x_pol1,y_pol1
      DIMENSION elevation(nn),sediment(nn),scrust(nn),GLit_ter(nn),
     +		  Tmoho(nn), D_Bodies(kpuntssd,2)

       ALLOCATE(x_pol1(kpuntssd),STAT=istat1) 
       ALLOCATE(y_pol1(kpuntssd),STAT=istat2) 
	x_pol1=0.D0
	y_pol1=0.D0
      WRITE(6,"(/'LOCAL ISOSTASY')") 

	!kpunts1=D_Bodies(1,1)			!! BODY 1		Lateral variations of density
	!qFlit1=D_Bodies(1,2)
	!i=1
  	!DO kp=2,kpunts1+1
        !    x_pol1(i)=D_Bodies(kp,1)
        !    y_pol1(i)=D_Bodies(kp,2)
	!    i=i+1
	!END DO

      DO iy=0,m
         DO ix=0,n
           kxy=ix+1+iy*(n+1)
	   sedim=sediment(kxy)
	   s=scrust(kxy)
	   hmantle=GLit_ter(kxy)-s-sedim
	   Tm=Tmoho(kxy)
	   !!!! Also on SUBROUTINE velvis  !!!!
	   rom=RHOAST*(1+((roalfa/2.D0)*(TISOTER-Tm)))		!! temperature depending
	!   rom=RHOAST						!! constant value
 	   ro_asth=RHOAST	! Isostasy level on the asthenosphere
	   roc_P=roc
	     !IF(hmantle<40.D3) roc_P=roc-100.D0	   
          	    !PY=iy*Dy
          	    !PX=ix*Dx
          	    !CALL outin (PX,PY,kpunts1,x_pol1,y_pol1,iadentro)	!! Lateral variations of density
               	    !IF(iadentro==1) roc_P=roc-100.D0	   
	   elevation(kxy)=elevation_isostasy(RHOH2O, rosed, roc_P, rom,
     +                                 ro_asth, sedim, s, hmantle)
         END DO 
      END DO 
         
      elevmin=VAL_MIN(elevation,nn)
      elevmax=VAL_MAX(elevation,nn)
      WRITE(6,22) elevmin,elevmax
 22   FORMAT(3X,'Elevation,   minima:',F9.2,
     + 			' m,   maxima:',F9.2,' m'/)

      IF (ALLOCATED(x_pol1)) DEALLOCATE(x_pol1,STAT=istat)       
      IF (ALLOCATED(y_pol1)) DEALLOCATE(y_pol1,STAT=istat)       
      RETURN 
      END SUBROUTINE ISOSTASY

!!! *******************************************************************
      SUBROUTINE DEFORMA_Bodies (Dtsegons, TEMPSMa, m, n, AX, BY, u, v,
     +				  nn, kpuntssd, D_Bodies, Bodyfile)

!!  Bodies (defined by points in counterclockwise direction) with
!!	    different strength and under time deformation
!!	    strength inside the body = qFlit * strength from the stress envelope
!!     All variables under S.I.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION u(nn), v(nn), D_Bodies(kpuntssd,2)
      CHARACTER*84 Bodyfile

        Dx=AX/n
        Dy=BY/m
	nbodies=0
      IF(TEMPSMa==0.D0) THEN
          D_Bodies=0.D0
	  OPEN(12,FILE=Bodyfile,STATUS='OLD',ACTION='READ')
	  READ(12,*,END=440) kpunts1,qFlit1	!! BODY 1
	  D_Bodies(1,1)=kpunts1
	  D_Bodies(1,2)=qFlit1
	  WRITE (6,"('Body 1: ',I4,' nodes.  qFlit1:',F11.2)")
     +             kpunts1,qFlit1
	  nbodies=nbodies+1
	  DO I=2,kpunts1+1
     	     READ(12,*) x,y
		D_Bodies(I,1)=x
		D_Bodies(I,2)=y
	  END DO
	  
	  READ(12,*,END=440) kpunts2,qFlit2	!! BODY 2
	  WRITE (6,"('Body 2: ',I4,' nodes.  qFlit2:',F11.2)")
     +             kpunts2,qFlit2
	  D_Bodies(kpunts1+2,1)=kpunts2
	  D_Bodies(kpunts1+2,2)=qFlit2
	  nbodies=nbodies+1
	  DO I=kpunts1+3,kpunts1+kpunts2+2
     	     READ(12,*) x,y
		D_Bodies(I,1)=x
		D_Bodies(I,2)=y
	  END DO	
	    
	  READ(12,*,END=440) kpunts3,qFlit3	!! BODY 3
	  WRITE (6,"('Body 3: ',I4,' nodes.  qFlit3:',F11.2)")
     +             kpunts3,qFlit3
	  D_Bodies(kpunts1+kpunts2+3,1)=kpunts3
	  D_Bodies(kpunts1+kpunts2+3,2)=qFlit3
	  nbodies=nbodies+1
	  DO I=kpunts1+kpunts2+4,kpunts1+kpunts2+kpunts3+3
     	     READ(12,*) x,y
		D_Bodies(I,1)=x
		D_Bodies(I,2)=y
	  END DO	  

 440	  WRITE(6,"(3X,'Number of bodies:',I3)") nbodies
 	  CLOSE(12)
      ELSE
      	  kpunts1=D_Bodies(1,1)
	  DO kp=2,kpunts1+1					!! BODY 1
                position_x=D_Bodies(kp,1)
		position_y=D_Bodies(kp,2)
		CALL new_position (position_x, position_y,
     +				    m, n, nn, Dx, Dy, u, v, Dtsegons)
                D_Bodies(kp,1)=position_x
                D_Bodies(kp,2)=position_y
          END DO
	  kpunts2=D_Bodies(kpunts1+2,1)
          DO kp=kpunts1+3,kpunts1+kpunts2+2			!! BODY 2
                position_x=D_Bodies(kp,1)
		position_y=D_Bodies(kp,2)
		CALL new_position (position_x, position_y,
     +					m, n, nn, Dx, Dy, u, v, Dtsegons)
                D_Bodies(kp,1)=position_x
                D_Bodies(kp,2)=position_y
          END DO
	  kpunts3=D_Bodies(kpunts1+kpunts2+3,1)
          DO kp=kpunts1+kpunts2+4,kpunts1+kpunts2+kpunts3+3	!! BODY 3
                position_x=D_Bodies(kp,1)
		position_y=D_Bodies(kp,2)
		CALL new_position (position_x, position_y,
     +					m, n, nn, Dx, Dy, u, v, Dtsegons)
                D_Bodies(kp,1)=position_x
                D_Bodies(kp,2)=position_y
          END DO
!! IF Body 2 is an indenter and I want to add nodes at the boundaries:
!		kf=0
!		ki=0
! 		IF(posicio2(kpunts2-1,2).GT.Dy) kf=1
! 		IF(posicio2(2,2).GT.Dy) ki=1
!		
!             	posicio2(kpunts2+kf+ki,1)=posicio2(kpunts2,1)
!             	posicio2(kpunts2+kf+ki,2)=posicio2(kpunts2,2)
! 		IF(kf.EQ.1) THEN
!             		posicio2(kpunts2+ki,1)=posicio2(kpunts2,1)
!             		posicio2(kpunts2+ki,2)=0.D0			
!		ENDIF	
! 		IF(ki.EQ.1) THEN
!	     		DO 49 kp=kpunts2-1,2,-1
!                		posicio2(kp+ki,1)=posicio2(kp,1)
! 49				posicio2(kp+ki,2)=posicio2(kp,2)
!             		posicio2(2,1)=posicio2(1,1)
!             		posicio2(2,2)=0.D0
!             		posicio2(1,1)=posicio2(1,1)
!             		posicio2(1,2)=posicio2(1,2)
!		ENDIF	
!	     	kpunts2=kpunts2+kf+ki
      ENDIF
      RETURN
      END SUBROUTINE DEFORMA_Bodies

!!! *******************************************************************
      SUBROUTINE DEFORMA_Points (Dtsegons, TEMPSMa, m, n, AX, BY, u, v,
     +				  nn, kpuntssd, Points_t, file_P)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION u(nn),v(nn),Points_t(0:kpuntssd,2)
      CHARACTER*84 file_P

      Dx=AX/n
      Dy=BY/m
      IF(TEMPSMa==0.D0) THEN
	  OPEN(10,FILE=file_P,STATUS='OLD',ACTION='READ')
	  npoints=0
	  DO I=1,kpuntssd
     	     READ(10,*,END=440) x,y
		npoints=npoints+1
		Points_t(I,1)=x
		Points_t(I,2)=y
	  END DO
 440  	  WRITE(6,"('DEFORMA_Points: Total Points:',I4)") npoints
	  Points_t(0,1)=npoints
      ELSE
          kppoints=Points_t(0,1)
	  kp_new=0
	  DO kp=1,kppoints
		position_x=Points_t(kp,1)
		position_y=Points_t(kp,2)
		CALL new_position (position_x, position_y,
     +				    m, n, nn, Dx, Dy, u, v, Dtsegons)
     
		IF(position_x>0 .AND. position_x<AX .AND. position_y>0
     +			.AND. position_y<BY) THEN
		    kp_new=kp_new+1
		    Points_t(kp_new,1)=position_x
		    Points_t(kp_new,2)=position_y
		ENDIF    
          END DO
	  Points_t(0,1)=kp_new
      ENDIF
      RETURN
      END SUBROUTINE DEFORMA_Points
!!*************************************************************************
      SUBROUTINE FINITE_ROTATION (Dtsegons, m, n, AX, BY, u, v, nn,
     +			   kpuntssd, Points_t, NumPoints, fin_rotat)
!!  Finite rotation: Vorticity Time Integration 

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION u(nn),v(nn),Points_t(0:kpuntssd,2),fin_rotat(NumPoints)
      PARAMETER (PI=3.1415926535897932D0)
      Dx=AX/n
      Dy=BY/m
      kppoints=Points_t(0,1)
      DO kp=1,kppoints
           zix=Points_t(kp,1)/Dx
           ziy=Points_t(kp,2)/Dy
           ix=aint(zix)
           iy=aint(ziy)
           kxy=ix+1+iy*(n+1)
           vx=(v(kxy+1)-v(kxy))/Dx
           uy=(u(kxy+n+1)-u(kxy))/Dy
           vorticity=vx-uy
           fin_rotat(kp)=fin_rotat(kp)+(Dtsegons*vorticity/2D0)*180D0/PI
      END DO
      rotmax=VAL_MAX(fin_rotat,NumPoints)
      rotmin=VAL_MIN(fin_rotat,NumPoints)
      WRITE(6,"(3X,'Finite Rotation (',I4,' points),  minima:',F12.4,
     +  ' degrees,   maxima:',F12.4,' degrees')") kppoints,rotmin,rotmax


      END SUBROUTINE FINITE_ROTATION

C **********************************************************************
C                 SUBROUTINE  CRUSTALTHICK_SYS
C          Troba el crustal thickness
C  Tinc en compte el nou i old crustal thickness -> 
C                                               SISTEMA D'EQUACIONS 
C  No faig les derivades en els punts de les vores, si no en els
C     punts del costat interiors.
C  Treballo amb VARIABLES DIMENSIONALS, S.I., metres i segons.
C  Dt(s) = Dtany(anys)*FACTEMP , FACTEMP=3.1536D7
C  scrust1 : old crustal thickness (t)
C  scrust2 : new crustal thickness (t+Dt)
C  u,v (m/s) : velocitats horizontals.

!       SUBROUTINE CRUSTALTHICK_SYS ( Dtany, m, n, nn, FACTEMP, AX, BY, 
!     +                               u, v, epuntzz, scrust1, scrust2)

!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!        DIMENSION u(nn),v(nn),epuntzz(nn),scrust1(nn),scrust2(nn),
!     +            ACRUST(nn,5),BCRUST(nn)
!        PARAMETER(TETA=0.0D0)

!         Dx=AX/n
!         Dy=BY/m
!         Dt=Dtany*FACTEMP
!         DO 3 kxy=1,nn
!             DO 3 lc=1,5
!                   ACRUST(kxy,lc)=0.D0
!                   BCRUST(kxy)=0.D0
! 3       CONTINUE
C  ------------- FIXO EL GRUIX DE L'ESCORA ALS CONTORN ---------------
C  NO TINC EN COMPTE TETA, PER TETA IGUAL A ZERO NO TINDRIA EQUACIà
!      DO 5 ix=0,n
C        iy=0 -> Dcrust/dy=0
!             kxy=ix+1
!             ACRUST(kxy,3)=1.D0
!             ACRUST(kxy,5)=-1.D0
!             BCRUST(kxy)=0.D0
C       iy=m -> Dcrust/dy=0
!             kxy=ix+1+m*(n+1)
!             ACRUST(kxy,3)=1.D0
!             ACRUST(kxy,2)=-1.D0
!             BCRUST(kxy)=0.D0
! 5    CONTINUE
!      DO 7 iy=1,m-1
C        ix=n -> Dcrust/dx=0
!             kxy=n+1+iy*(n+1)
!             ACRUST(kxy,3)=1.D0
!             ACRUST(kxy,1)=-1.D0
!             BCRUST(kxy)=0.D0
C       ix=0 -> Dcrust/dx=0
!             kxy=1+iy*(n+1)
!             ACRUST(kxy,3)=1.D0
!             ACRUST(kxy,4)=-1.D0
!             BCRUST(kxy)=0.D0
! 7    CONTINUE
C  -------------------------------------------------------------------
!      DO 13 iy=1,m-1
!         DO 13 ix=1,n-1
!             kxy=ix+1+iy*(n+1)
!             ACR4=(TETA*Dt*u(kxy))/(2.D0*Dx)
!             ACR5=(TETA*Dt*v(kxy))/(2.D0*Dy) 
!             ACR3=1.D0-TETA*Dt*epuntzz(kxy)
!             ACR1=-ACR4
!             ACR2=-ACR5
!             B1=epuntzz(kxy)*scrust1(kxy)
!             B2=(u(kxy)/(2.D0*Dx))*(scrust1(kxy+1)-scrust1(kxy-1))
!             B3=(v(kxy)/(2.D0*Dy))*(scrust1(kxy+n+1)-scrust1(kxy-n-1))
C   Ho normalitzo perque la diagonal sigui 1. 
!             ACRUST(kxy,1)=ACR1/ACR3
!             ACRUST(kxy,2)=ACR2/ACR3
!             ACRUST(kxy,3)=1.D0
!             ACRUST(kxy,4)=ACR4/ACR3
!             ACRUST(kxy,5)=ACR5/ACR3
!             BCRUST(kxy)=(scrust1(kxy)+(1.D0-TETA)*Dt*(B1-B2-B3))/ACR3
!13    CONTINUE
!
!       CALL sistbanda(ACRUST,BCRUST,nn,5,2,scrust2)

!      RETURN
!      END SUBROUTINE CRUSTALTHICK_SYS

!! **********************************************************************

      SUBROUTINE VISCOSITAT (TEMPSMa, npasmax, m, n, NELZ, nn, AX, BY,
     +                      Dz, u, v, epuntzz, elevation, sediment,
     +                      scrust,GLit_ter, TEMPE, visco, nitermax, 
     +                      strainrate, TISOTER, IRheology_type, Qarray,
     +                      enarray, Aarray, TIPPAR, kpuntssd, D_Bodies, 
     +                      vissup, visinf, ITSR, visco_cnst, Flitos, 
     +                      dvis_allow, Dif_K, iwrite, file_ep, ifixv)

!! ITSR= 0 - 3	
!!	ITSR=-1 -> Constant lithosphere strength to every where, no change in time.
!!	ITSR=0 -> Constant viscosity to every where, no change in time.
!!	ITSR=1 -> Calculate viscosity using a constant strain rate (strainrate) 
!! 	ITSR=2 -> Calculate viscosity using the efective strain rate at each point on the stress envelope 
!!			and on the constitutive equation
!!	ITSR=3 -> Calculate viscosity using the efective strain rate at each point on the stress envelope 
!!			and the constant strain rate (strainrate) on the constitutive equation.

! dvis_allow=[d(vis)/dx/]vis  [e.g. (((4-1)*1.D23)/30000m)/1.d23 = 1D-4 ]; dvis_allow=1.5D-5 (Poster EGS2002)
! Dif_K: Cnst. of diffusive filter,  Dif_K=1.D8  Dif_K<3.D8 (?)
!		major discretitzacio horizontal -> Dif_K ha de ser menor
! iswitch=1 => keep in a file the stress envelope and the strength

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*40 TIPPAR
      CHARACTER*55 zname
      CHARACTER*20 file_ep
      PARAMETER(FACVEL=3.1536D10,FACTEMP=3.1536D7)
      PARAMETER (epv_max=1.D-15, epv_min=1.D-18, iswitch=0)  
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x_pol1, y_pol1,
     +		x_pol2,y_pol2, x_pol3,y_pol3, T, ep_eff, ep_visc, ibody
      DIMENSION u(nn),v(nn),scrust(nn),GLit_ter(nn),TEMPE(nn,0:NELZ),
     +      elevation(nn),sediment(nn),visco(nn),Flitos(nn),epuntzz(nn),
     +      D_Bodies(kpuntssd,2),Qarray(3),enarray(3),Aarray(3)
!     +		Tmoho_tmp(nn), Qsurface_tmp(nn),vis_input(nn)
     
      WRITE(6,"('CALCULATION OF THE VISCOSITY')") 

       ALLOCATE(x_pol1(kpuntssd),STAT=istat1) 
       ALLOCATE(y_pol1(kpuntssd),STAT=istat2) 
       ALLOCATE(x_pol2(kpuntssd),STAT=istat3) 
       ALLOCATE(y_pol2(kpuntssd),STAT=istat4) 
       ALLOCATE(x_pol3(kpuntssd),STAT=istat5) 
       ALLOCATE(y_pol3(kpuntssd),STAT=istat6) 
       ALLOCATE(T(0:NELZ),STAT=istat7) 
       ALLOCATE(ep_visc(nn),STAT=istat8) 
       ALLOCATE(ep_eff(nn),STAT=istat9) 
       ALLOCATE(ibody(nn),STAT=istat10) 
	x_pol1=0.D0
	y_pol1=0.D0
	x_pol2=0.D0
	y_pol2=0.D0
	x_pol3=0.D0
	y_pol3=0.D0
	T=0.D0
	ep_visc=0.D0
	ep_eff=0.D0
	ibody=0.D0
       IF(istat1>0.OR.istat2>0.OR.istat3>0.OR.istat4>0.OR.istat5>0.OR.
     +       istat6>0.OR.istat7>0.OR.istat8>0.OR.istat9>0.OR.istat10>0)
     +       STOP " Error allocation on SUBR.VISCOSITAT" 

!      vis_input=visco
!      Tmoho_tmp=700.D0
!      Qsurface_tmp=70.D-3
      Flitmin=1.D35
      Flitmax=-1.D35
      
      ZASTH=NELZ*Dz
      Dx=AX/n
      Dy=BY/m
      ddx=2.D0*Dx
      ddy=2.D0*Dy
      IF(ifixv==1) OPEN(17,FILE='P_fixVelocity.res',STATUS='REPLACE')
      IF(TEMPSMa==0.D0.AND.ITSR/=0) OPEN(14,FILE='Flit_Zmec_xy.res')

	kpunts1=D_Bodies(1,1)			!! BODY 1
	qFlit1=D_Bodies(1,2)
	i=1
  	DO kp=2,kpunts1+1
            x_pol1(i)=D_Bodies(kp,1)
            y_pol1(i)=D_Bodies(kp,2)
	    i=i+1
	END DO

	kpunts2=D_Bodies(kpunts1+2,1)		!! BODY 2
	qFlit2=D_Bodies(kpunts1+2,2)
	i=1
	DO kp=kpunts1+3,kpunts1+kpunts2+2	
	    x_pol2(i)=D_Bodies(kp,1)
	    y_pol2(i)=D_Bodies(kp,2)
	    i=i+1
	END DO

	kpunts3=D_Bodies(kpunts1+kpunts2+3,1)	!! BODY 3
	qFlit3=D_Bodies(kpunts1+kpunts2+3,2)
	i=1
	DO kp=kpunts1+kpunts2+4,kpunts1+kpunts2+kpunts3+3	
	    x_pol3(i)=D_Bodies(kp,1)
	    y_pol3(i)=D_Bodies(kp,2)
	    i=i+1
	END DO

      IF(TEMPSMa.EQ.0.D0) THEN 
	   IF(kpunts1/=0) WRITE(6,"(3X,'BODY 1: kpunts1=',I3,
     +	     ', STRENGTH MULTIPLIED BY qFlit1=',F11.2)")kpunts1,qFlit1
	   IF(kpunts2/=0) WRITE(6,"(3X,'BODY 2: kpunts2=',I3,
     +	     ', STRENGTH MULTIPLIED BY qFlit2=',F11.2)")kpunts2,qFlit2
	   IF(kpunts3/=0) WRITE(6,"(3X,'BODY 3: kpunts3=',I3,
     +	     ', STRENGTH MULTIPLIED BY qFlit3=',F11.2)")kpunts3,qFlit3
      ENDIF
CCC -------------------------------------------------------------------
	
      IF(ITSR==0) THEN		!! Fixed viscosity
      	 DO iy=0,m
            PY=iy*Dy
            DO ix=0,n
	     TIPPAR='Constant viscosity:'
               kxy=ix+1+iy*(n+1)
               PX=ix*Dx
	       visco(kxy)=visco_cnst
		IF(kpunts1/=0) THEN
           	    CALL outin (PX,PY,kpunts1,x_pol1,y_pol1,iadentro)
               	    IF(iadentro==1) THEN
			visco(kxy)=visco_cnst*qFlit1
			IF(ifixv==1.AND.ix/=0.AND.ix/=n.AND.iy/=0.
     +				AND.iy/=m) WRITE(17,"(2I6,2F11.2)")
     +					ix, iy, PX/1.D3, PY/1.D3	!! Points where the velocity is fixed
		    ENDIF
		ENDIF
		IF(kpunts2/=0) THEN
           	    CALL outin (PX,PY,kpunts2,x_pol2,y_pol2,iadentro)
               	    IF(iadentro==1) visco(kxy)=visco_cnst*qFlit2
		ENDIF		
		IF(kpunts3/=0) THEN
           	    CALL outin (PX,PY,kpunts3,x_pol3,y_pol3,iadentro)
               	    IF(iadentro==1) visco(kxy)=visco_cnst*qFlit3
		ENDIF
               IF(visco(kxy).GT.vissup) visco(kxy)=vissup
               IF(visco(kxy).LT.visinf) visco(kxy)=visinf		
	    END DO
	 END DO
      ENDIF

      IF(ITSR/=0) THEN		!! stress envelope => viscosity
	 ep_eff=strainrate
	 IF(iwrite==1) THEN
      		OPEN(30,FILE=file_ep)
		WRITE(30,"(' # Time:',F8.3,'  x[km]  y[km]  ',
     +			'ep_eff[s-1]  ep_visc[s-1]')") TEMPSMa
	 ENDIF
	 IF(ITSR/=1.AND.TEMPSMa/=0.D0) THEN
	     DO iy=0,m
		DO ix=0,n
		   kxy=ix+1+iy*(n+1)
		   ep_eff(kxy)=effective_strainrate (ix, iy, m, n, nn,
     +							Dx, Dy, u, v)
                END DO  
             END DO  
	     !zname='effective strain rate, for stress envelope'
	     !CALL smooth_horizontal (m, n, nn, Dx, Dy, ep_eff)
	     !iwr_filt=-1
	     !CALL fit_gradients (m, n, nn, Dx, Dy, Dif_K, dvis_allow,
    !+					ep_eff, zname, iwr_filt)
	 ENDIF
	 ibody=99
         DO iy=0,m
            PY=iy*Dy
            DO ix=0,n
               kxy=ix+1+iy*(n+1)
               PX=ix*Dx
	       IF(ITSR==-1) THEN
		  Flitos(kxy)=visco_cnst
	       ELSE
		  ZMOHO=scrust(kxy)+sediment(kxy)
		  ZLITOS=GLit_ter(kxy)
		  DO iz=0,NELZ
		     T(iz)=TEMPE(kxy,iz)
		  END DO  
		  ep_eff_P=ep_eff(kxy)
		  IRheology_type_P=IRheology_type
		  !IF(ZMOHO>40) IRheology_type_P=110	!! OJO  change rheological parameters
		  CALL Ch_tree_lib (NELZ, Dz, T, ep_eff_P, ZMOHO, 
     +				   ZLITOS, Flitcomp, Flitext, ZMECANIC,
     +				   IRheology_type_P,Qarray, enarray, 
     +				   Aarray, TIPPAR, iswitch)
		  Flitos(kxy)=(ABS(Flitcomp)+Flitext)/2.D0
		  !IF(epuntzz(kxy)>0.0) THEN
		  !	Flitos(kxy)=ABS(Flitcomp)	!! epuntzz compressive -> Flitcomp
		  !ELSE
		  !	Flitos(kxy)=Flitext		!! epuntzz extensive -> Flitext
		  !ENDIF
	       ENDIF
		IF(kpunts1/=0) THEN
           	    CALL outin (PX,PY,kpunts1,x_pol1,y_pol1,iadentro)
               	    IF(iadentro==1) THEN
		    	ibody(kxy)=1
		    	Flitos(kxy)=Flitos(kxy)*qFlit1
			IF(ifixv==1.AND.ix/=0.AND.ix/=n.AND.iy/=0.
     +				AND.iy/=m) WRITE(17,"(2I6,2F11.2)")
     +					ix, iy, PX/1.D3, PY/1.D3	!! Points where the velocity is fixed
		    ENDIF
		ENDIF
		IF(kpunts2/=0) THEN
           	    CALL outin (PX,PY,kpunts2,x_pol2,y_pol2,iadentro)
               	    IF(iadentro==1) THEN
		    	ibody(kxy)=2
		    	Flitos(kxy)=Flitos(kxy)*qFlit2
		    ENDIF
		ENDIF		
		IF(kpunts3/=0) THEN
           	    CALL outin (PX,PY,kpunts3,x_pol3,y_pol3,iadentro)
               	    IF(iadentro==1) THEN
		    	ibody(kxy)=3
		    	Flitos(kxy)=Flitos(kxy)*qFlit3
		    ENDIF
		ENDIF		

	       IF(TEMPSMa==0.D0) WRITE(14,"(1X,2F13.1,1P,E13.4,I5)")
     +			        ix*Dx,iy*Dy,Flitos(kxy),ibody(kxy)
               Flitmin=MIN(Flitos(kxy),Flitmin)
               Flitmax=MAX(Flitos(kxy),Flitmax)
            END DO
         END DO
	 ep_effmin=VAL_MIN(ep_eff,nn)
	 ep_effmax=VAL_MAX(ep_eff,nn)
         WRITE(6,"(3X,'Strain rate used on the stress envelope  minim:',
     +			1P,G12.4,' s-1,     maxim:',1P,G12.4,' s-1'/
     +			3X,'Lithospheric strength,  minim:',1P,G12.4,
     +        		' N/m,     maxim:',1P,G12.4,' N/m')")
     +			ep_effmin, ep_effmax, Flitmin, Flitmax
!! Lithosphere Strength --> Viscosity
	 !zname='effective strain rate, for viscosity'
	 ep_visc=ep_eff 
	 !Dif_smooth=1.6*Dif_K
	 !dvis_al=dvis_allow/5.0
	 !iwr_filt=-1
	 !CALL fit_gradients (m, n, nn, Dx, Dy, Dif_smooth, dvis_al,	!! smoother strain rate for the Strength-Viscosity relation
    !+					ep_visc, zname, iwr_filt)
         DO iy=0,m
            DO ix=0,n
               kxy=ix+1+iy*(n+1)
	       !IF(ep_visc(kxy)>epv_max) ep_visc(kxy)=epv_max
	       !IF(ep_visc(kxy)<epv_min) ep_visc(kxy)=epv_min
	       !ep_visc(kxy)=MIN(ep_visc(kxy),epv_max)
	       !ep_visc(kxy)=MAX(ep_visc(kxy),epv_min)
	       IF(iwrite==1) WRITE(30,"(2F11.3,2G18.6)")
     +		       ix*Dx/1.D3, iy*Dy/1.D3, ep_eff(kxy), ep_visc(kxy)
	    END DO
	 END DO
	 ep_viscmin=VAL_MIN(ep_visc,nn)
	 ep_viscmax=VAL_MAX(ep_visc,nn)
         WRITE(6,"(3X,'Strain rate used on viscosity,       minim:',
     +			1P,G12.4,' s-1,     maxim:',1P,G12.4,' s-1')")
     +			ep_viscmin, ep_viscmax
         DO kxy=1,nn		
	     thick_layer=ZASTH	!+elevation(kxy)     thick_layer: GLit_ter o Zmecanica o ZASTH+elevation o thick_layer=ZLITOS ??
	     IF(ITSR==2. OR. ITSR==-1) THEN
		 visco(kxy)=Flitos(kxy)/(2.0*thick_layer*ep_visc(kxy))	!! Factor per incrementar la visco (*10.D0)
	      ELSE
		 visco(kxy)=Flitos(kxy)/(2.0*thick_layer*strainrate)
	     ENDIF
	     !IF(visco(kxy)>vissup .AND. ibody(kxy)/=1) visco(kxy)=vissup
	     visco(kxy)=MIN(visco(kxy),vissup)
	     visco(kxy)=MAX(visco(kxy),visinf)
	!!     IF(ibody(kxy)==1) visco(kxy)=1.D25		!! TIBETAN PLATEAU
         END DO
      ENDIF

      DO kxy=1,nn
	 ZMOHO=scrust(kxy)+sediment(kxy)
	 IF(ZMOHO>50.D3) visco(kxy)=visco(kxy)/1.D1	!! Tibet,  Low viscosity for thicker crust
      END DO
      zname='viscosity'
      iwr_filt=-1							!! Filter Viscosity 
      CALL fit_gradients (m, n, nn, Dx, Dy, Dif_K, dvis_allow,
     +				visco, zname, iwr_filt)			!! Filter Viscosity

      DO kxy=1,nn
	 IF(ibody(kxy)/=1) THEN				!! Tibet,  Body 1 = Indenter
	 	visco(kxy)=MIN(visco(kxy),vissup)
	 	visco(kxy)=MAX(visco(kxy),visinf)
	 ENDIF
      END DO
      vismax=VAL_MAX(visco,nn)
      vismin=VAL_MIN(visco,nn)

      IF(TEMPSMa==0.D0) THEN
	 IF(ITSR==0) THEN 
	     WRITE(6,"(3X,'Constant viscosity:',1P,G12.3,' Pa.s')")
     + 				visco_cnst
	  ELSE	 
      	     WRITE(6,"(3X,'Rheological parameters used: ',A40)") TIPPAR
      	     IF(ITSR==1) WRITE(6,
     +           "(4X,'cnst strain rate:',1P,G12.3,' s-1')") strainrate
      	     IF(ITSR==2)	 
     +	          PRINT*,'  strain rate used: effective strain rate'
      	     IF(ITSR==3) WRITE(6,"(4X,'effective strain rate ',
     +		     'on the stress envelope to calculate Flit'/
     +		     5X,'cnst strain rate to calculate the viscosity:',
     +		     1P,G12.3,' s-1')") strainrate
         ENDIF
      ENDIF
      WRITE(6,"(3X,'Viscosity,              minima:',1P,G12.4,
     +		' Pa.s,   maxima:',1P,G12.4,' Pa.s')") vismin, vismax     
      CLOSE(14)
      CLOSE(17)  
      CLOSE(30)  
      IF (ALLOCATED(x_pol1)) DEALLOCATE(x_pol1,STAT=istat)       
      IF (ALLOCATED(y_pol1)) DEALLOCATE(y_pol1,STAT=istat)       
      IF (ALLOCATED(x_pol1)) DEALLOCATE(x_pol2,STAT=istat)       
      IF (ALLOCATED(y_pol1)) DEALLOCATE(y_pol2,STAT=istat)       
      IF (ALLOCATED(x_pol1)) DEALLOCATE(x_pol3,STAT=istat)       
      IF (ALLOCATED(y_pol1)) DEALLOCATE(y_pol3,STAT=istat)       
      DEALLOCATE (ibody,STAT=istat) 
      DEALLOCATE (ep_eff,STAT=istat) 
      DEALLOCATE (ep_visc,STAT=istat) 
      DEALLOCATE (T,STAT=istat) 
      RETURN
      END SUBROUTINE VISCOSITAT

CC ********************************************************************
CC ********************************************************************
      SUBROUTINE velvis (TEMPSMa, Dt, AX, BY, n, m, nn, nitermax, 
     +			 tallmax, alfa, g, Flitos, visco, vissup, 
     +			 visinf, elevation, sediment, s, GL, u, v,
     +			 epuntzz, RHOH2O, rosed, roc, RHOAST, roalfa,
     +			 TISOTER, Tmoho, ZASTH, BCfile, iwrite, 
     +			 file_Pav, kpuntssd, D_Bodies, ifixv)

CC  Calculation of velocity field of a thin sheet, from a given 
CC      lithospheric structure and viscosity. 
CC	viscosity is velocity depending -> iteration
CC
CC   ix=0,...,n.  iy=0,...,m.
CC   'x' positive to the rigth, 'y' positive upward.
CC   vel [m/s]	Array of the velocity field: u_1, v_1, u_2, v_2, ... u_nn, v_nn.
CC			From west to east and south to north.
CC   vel(nincogn)=(u,v)(nincogn)
CC   velold : velocity field at the first iteration
CC   vel = alfa * vel + (1-alfa) * velold    (alfa=.5 is recommended)
CC   ep_limit [s-1]  Imposed limit for calculated vert.str.rate 'epuntzz'.
CC
!!   To high strain rate -> decrease of the Dt
!!  As small are Dx and Dy -> smaller should be Dt
!!     (strain maxim) = Dt*(strain rate maxim) < quotadef

!!  s [m]	Array of crustal thickness. Does not include sediments.
!!  GL [m]	Array of lithospheric thickness. Includes sedims. and crust.
!!  z_compens [m]	Depth of compensation (e.g., ZASTH, 300e3 m)
 

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*84 TITOL_BC, BCfile
      CHARACTER*20 file_Pav
      PARAMETER (FACTEMP=3.1536D7, DP_sup=10000.D6, DP_inf=-10000.D6)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: average_Szz,
     +						     x_pol1, y_pol1
      DIMENSION u(nn), v(nn), Flitos(nn), visco(nn), elevation(nn),
     +           sediment(nn), s(nn), GL(nn), epuntzz(nn), Tmoho(nn),
     +		 D_Bodies(kpuntssd,2)
      PARAMETER (ep_limit=5.D-15)

      WRITE(6,"(/'CALCULATION OF THE VELOCITY FIELD')") 

       ALLOCATE(x_pol1(kpuntssd),STAT=istat1) 
       ALLOCATE(y_pol1(kpuntssd),STAT=istat2) 
       ALLOCATE(average_Szz(nn),STAT=istat) 
	x_pol1=0.D0
	y_pol1=0.D0
	average_Szz=0.D0
      IF (istat>0) STOP " Error allocation average_Szz(:)" 
      
      rommin=500.D30
      rommax=0.D0
      Dx=AX/n
      Dy=BY/m
      z_compens=ZASTH
  !    z_compens=0.D3 
  !    DO kxy=1,nn			! z_compens: depth until where I calculate the GPE	OJO No funciona!!!
  !        z_lithos=GL(kxy)-elevation(kxy)							No extension a continental 
  !        z_compens=MAX(z_compens,z_lithos)							lithosphere root remove !!!!
  !    END DO
  !	z_compens=z_compens+5.D3	! 5 km deeper than the deepest lithosphere-asthenosphere boundary	

      IF(iwrite==1) THEN
      	     OPEN(30,FILE=file_Pav)
	     WRITE(30,"(' # Time:',F8.3,' My,  Z_comp:',F9.3,' km.  ',
     +			'x[km] y[km] Szz[Pa] Flit[N/m] rom[Kg/m3]')")
     +			TEMPSMa,z_compens/1.D3
      ENDIF

	!kpunts1=D_Bodies(1,1)			!! BODY 1   Lateral variations of density
	!qFlit1=D_Bodies(1,2)
	!i=1
  	!DO kp=2,kpunts1+1
        !    x_pol1(i)=D_Bodies(kp,1)
        !    y_pol1(i)=D_Bodies(kp,2)
	!    i=i+1
	!END DO

      DO iy=0,m
         DO ix=0,n
             kixy=ix+1+iy*(n+1)
	     elev=elevation(kixy)
	     sed=sediment(kixy)
	     crust=s(kixy)
	     hlitos=GL(kixy)
	   !!!! Also on SUBROUTINE ISOSTASY !!!    Controlar utilitzar el mateix criteri a graficth.f  ****
	     rom=RHOAST*(1+((roalfa/2.D0)*(TISOTER-Tmoho(kixy))))		!! temperature depending
 	  !   rom=RHOAST							!! constant value
 	     ro_asth=RHOAST							! Isostasy level on the asthenosphere	     
	     roc_P=roc
	     !hm=hlitos-crust-sed
	     !IF(hm<40.D3) roc_P=roc-100.D0
          	    !PY=iy*Dy
          	    !PX=ix*Dx
          	    !CALL outin (PX,PY,kpunts1,x_pol1,y_pol1,iadentro)		!! Lateral variations of density
               	    !IF(iadentro==1) roc_P=roc-100.D0
     	     average_Szz(kixy)=P_average(elev, sed, crust, hlitos,		!! RW_PARAMETRES.f
     +					     z_compens, RHOH2O, rosed,
     +					     roc_P, rom, ro_asth, g)
!              average_Szz(kixy)=1.5*average_Szz(kixy)
!     	     CALL P_averagelib (elev, sed, crust, hlitos,
!     +				z_compens, RHOH2O, rosed,
!     +				roc, rom, ro_asth, g, P_average)
!	     average_Szz(kixy)=P_average
	     IF(iwrite==1) WRITE(30,"(2F11.3,F16.2,E16.6,F12.3)")
     +	        ix*Dx/1.D3,iy*Dy/1.D3,average_Szz(kixy),Flitos(kixy),rom
	     rommin=MIN(rommin,rom)
	     rommax=MAX(rommax,rom)
         END DO
      END DO
      CLOSE(30)

      GLmin=VAL_MIN(GL,nn)
      GLmax=VAL_MAX(GL,nn)
      elevmin=VAL_MIN(elevation,nn)
      elevmax=VAL_MAX(elevation,nn)
      vismin=VAL_MIN(visco,nn)
      vismax=VAL_MAX(visco,nn)
      Pmin=VAL_MIN(average_Szz,nn)
      Pmax=VAL_MAX(average_Szz,nn)
      WRITE(6,"(3X,'Density of lithospheric mantle, minim:',F9.2,
     +	     ' kg/m3,  maxim:',F9.2,' kg/m3.  (asthenosphere:',F9.2,')'/
     +	     3X,'Depth of considered compensation: ',F8.2,' km'/
     +	     3X,'Depth average vertical stress = (Integrated vertical ',
     +	     'stress)/(layer thickness)'/
     +	     8X,'minim:',1P,E14.6,' N/m3,   maxim:',1P,E14.6,' N/m3')")
     +	     	rommin,rommax,ro_asth,z_compens/1.D3,Pmin,Pmax
     
      nbanda=4*(n+1)+7 
      nincogn=2*nn
      DO kxy=1,nn		! OJO utilitzar la mateixa que a la subroutine VISCOSITAT thick_layer
	  thick_layer=ZASTH	!+elevation(kxy) !thick_layer: GLit_ter o Zmecanica o ZASTH+elevation o thick_layer=ZLITOS ??
	  Flitos(kxy)=Flitos(kxy)/(2.0*thick_layer)
      END DO
      CALL VELOCITYFIELD (TEMPSMa, Dt, AX,BY, n,m, nn, vissup, visinf,
     +			 Flitos, visco, nincogn,tallmax, alfa, nitermax,
     +			 average_Szz, u, v,BCfile, nbanda, ifixv)
!      CALL vertical_strain_rate_2D (AX, BY, n, m, nn, u, v, epuntzz)

	epuntzz_min=1.D30
	epuntzz_max=-1.D30
	DO iy=0,m
	   DO ix=0,n
      		CALL vertical_strain_rate (ix, iy, m, n, nn, Dx, Dy,	!! Vertical Strain Rate lib/RW_PARAMETRES.f
     +					    u, v, ep_limit, epuntzz_P)
		kxy=ix+1+iy*(n+1)
     		epuntzz(kxy)=epuntzz_P
		!IF(iy==0.OR.iy==1.OR.iy==2) epuntzz(kxy)=0.0		!Fixo el strain rate a les bores
		!IF(iy==m.OR.iy==(m-1).OR.iy==(m-2)) epuntzz(kxy)=0.0
          	epuntzz_max=MAX(epuntzz(kxy),epuntzz_max)
              	epuntzz_min=MIN(epuntzz(kxy),epuntzz_min)
	   END DO 
	END DO
	WRITE(6,"(4X,'Vertical strain rate: minim:',1P,G12.4,
     +	  ' s-1,    maxim:',1P,G12.4,' s-1,  limit:',1P,G9.2,' s-1')")
     +		epuntzz_min, epuntzz_max, ep_limit
       
      IF (ALLOCATED(average_Szz)) 
     +			DEALLOCATE(average_Szz,STAT=istat)       
      IF (ALLOCATED(x_pol1)) DEALLOCATE(x_pol1,STAT=istat)       
      IF (ALLOCATED(y_pol1)) DEALLOCATE(y_pol1,STAT=istat)       
      RETURN
      END SUBROUTINE velvis
              
CC ********************************************************************
CC ********************************************************************
CC		SUBROUTINE call_uhuru_tisc
!!  Call 'call_surf_proc' SURFACE PROCESES.
!!  Input: elevation, sediment and crustal thickness
!!  Output: sediment and crustal thickness
!!  noSed er.:  erosio de roca mare (l'erosio de roques que no son sediments)
!!  sed.incr.:  increment de sediments
!!  outp.seds:  sediments que surten per les vores.

      SUBROUTINE call_uhuru_tisc (iw_erosion,iflexure, iwrite, Dtsegons,
     +		      m, n, nn, nn_tao, rosed, roc, roalfa, RHOAST, 
     +		      elevation, elevation0, sediment, scrust, scrust0,
     +		      erosion, nresults, Dx, Dy, hmantle0,
     +		      GLit_ter, Tmoho, TISOTER, ZASTH)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER*2 write_files
      LOGICAL exist_file
      CHARACTER*115 string
      CHARACTER*15 file_st,file_bas,file_xyw, file_tisc
      REAL topo_tisc(0:nn_tao), sedi_tisc(0:nn_tao), q_load(0:nn_tao),
     +	   topo_2eros(0:nn_tao),sed0(0:nn_tao), Dt_tisc, 
     +	   deflection(0:nn_tao)
      DIMENSION elevation(nn), elevation0(nn), sediment(nn), scrust(nn),
     +	    scrust0(nn),erosion(nn),Tmoho(nn),GLit_ter(nn),hmantle0(nn)
      PARAMETER (g=9.8D0)
	
      defl_min=1.D100
      defl_max=-1.D100
      topo_min=1.D100
      topo_max=-1.D100

      Nx_tao=n+1
!  Definition 'tao' arrays, topo_tisc, sedi_tisc:
      DO iy=0,m
	 DO ix=0,n
	     ktao=(m-iy)*Nx_tao+ix
	     kxy=ix+1+iy*(n+1)
	     rom=RHOAST*(1+((roalfa/2.D0)*(TISOTER-Tmoho(kxy))))
	     hmantle=GLit_ter(kxy)-scrust(kxy)-sediment(kxy)
	     topo_tisc(ktao)=REAL(elevation(kxy))
	     sedi_tisc(ktao)=REAL(sediment(kxy))
	     q=roc*g*(scrust(kxy)-scrust0(kxy))+(rosed*g*sediment(kxy))
     +			+(hmantle-hmantle0(kxy))*g*(rom-RHOAST)
	     q_load(ktao)=REAL(q)
	     topo_2eros(ktao)=topo_tisc(ktao)
	     sed0(ktao)=REAL(sediment(kxy))
         END DO
      END DO

      write_files=iwrite
      Dt_tisc=REAL(Dtsegons)
      IF(iflexure==1) THEN      
	  WRITE(6,"(/'FLEXURE CALCULATION ')") 
	  CALL CALL_FLEXURE (q_load, deflection, write_files)
	  DO iy=0,m
	     DO ix=0,n
	        ktao=(m-iy)*Nx_tao+ix
	        kxy=ix+1+iy*(n+1)
	        t=elevation0(kxy)+scrust(kxy)-scrust0(kxy)+
     +				sediment(kxy)-deflection(ktao)
	        topo_tisc(ktao)=REAL(t)
	        topo_2eros(ktao)=topo_tisc(ktao)
	        d=deflection(ktao)
                topo_min=MIN(t,topo_min)
                topo_max=MAX(t,topo_max)
                defl_min=MIN(d,defl_min)
                defl_max=MAX(d,defl_max)
	     END DO
	  END DO
	  WRITE(6,"(2X,'Deflection,',9X,'minima:',F9.2,
     +		' m,    maxima:',F9.2,' m'/2X,'Topo flexure:',
     +		7X,'minima:',F9.2,' m,    maxima:',F9.2,' m')")
     +		defl_min, defl_max, topo_min, topo_max
      ENDIF   

      IF(iw_erosion==1) THEN      
	  WRITE(6,"(/'CALCULATION OF THE SURFACE PROCESSES')") 
	  CALL call_surf_proc (Dt_tisc,topo_tisc,sedi_tisc,write_files)
      ENDIF
		
      DO iy=0,m
	 DO ix=0,n
	     ktao=(m-iy)*Nx_tao+ix
	     kxy=ix+1+iy*(n+1)
	     Dtopo=topo_tisc(ktao)-topo_2eros(ktao)
		 IF(Dtopo>0.0) THEN
		     sediment(kxy)=sed0(ktao)+Dtopo
		 ELSE
                     erosion(kxy)=erosion(kxy)-Dtopo
		     IF(ABS(Dtopo)>sed0(ktao)) THEN
			  sediment(kxy)=0.D0
			  scrust(kxy)=scrust(kxy)+Dtopo+sed0(ktao)
		     ELSE
			  sediment(kxy)=sed0(ktao)+Dtopo
		     ENDIF
		ENDIF
	     elevation(kxy)=elevation(kxy)+Dtopo	!elevation(kxy)=topo_tisc(ktao)	
         END DO
      END DO
      elevmin=VAL_MIN(elevation,nn)
      elevmax=VAL_MAX(elevation,nn)
      sedmin=VAL_MIN(sediment,nn)
      sedmax=VAL_MAX(sediment,nn)
      crustmin=VAL_MIN(scrust,nn)
      crustmax=VAL_MAX(scrust,nn)
      WRITE(6,"(/
     +   2X,'Elevation,',10X,'minima:',F9.2,' m,    maxima:',F9.2,' m'/
     +   2X,'Sediments,',10X,'minim:',F9.2,' m,     maxim:',F9.2,' m'/
     +   2X,'Crustal thickness,  minim:',F9.3,' km,    maxim:',
     +   F9.4,' km'/)") elevmin, elevmax, sedmin,sedmax,
     +       		crustmin/1.D3, crustmax/1.D3
     
      IF(iwrite==1) THEN
      
         nres_dec1=AINT((nresults-48.D0)/10.D0)
         nres_unit=nresults-10*nres_dec1
         nres_dec=nres_dec1+48
         IF(nres_dec==48) THEN
           file_tisc='results_tisc'//CHAR(nres_unit)
           file_st='res_tisc'//CHAR(nres_unit)//'.st'
           file_xyw='res_tisc'//CHAR(nres_unit)//'.xyw'
         ELSE
           file_tisc='results_tisc'//CHAR(nres_dec)//CHAR(nres_unit)
           file_st='res_tisc'//CHAR(nres_dec)//CHAR(nres_unit)//'.st'
           file_xyw='res_tisc'//CHAR(nres_dec)//CHAR(nres_unit)//'.xyw'
         ENDIF

         OPEN(30,FILE=file_tisc)
	    WRITE(30,"('#    x[km]     y[km]   deflection[m]   ', 
     +		  'e_flexure[m]   e_tisc[m]   e_uhuru[m]    load/roc*g',
     +		  '     crust[m]   init_crust[m]   sediments[m]')")
	    DO iy=0,m
		DO ix=0,n
		   ktao=(m-iy)*Nx_tao+ix
		   kxy=ix+1+iy*(n+1)
		   q=q_load(ktao)/(roc*g)
		   WRITE(30,"(2F10.3,8F14.2)")  ix*Dx/1D3, iy*Dy/1D3,
     +				deflection(ktao), topo_2eros(ktao),
     +				topo_tisc(ktao), elevation(kxy), q, 
     +				scrust(kxy), scrust0(kxy), sediment(kxy)
		END DO
	    END DO
	 CLOSE (30)    

!         OPEN(30,FILE='res_tisc.st',STATUS='OLD',ACTION='READ',
!     +			IOSTAT=iexist)
!	 IF(iexist==0) THEN
!	    OPEN(31,FILE=file_st)
!	    DO i=1,200000
!		READ(30,"(A60)",END=40) string
!		WRITE(31,"(A60)") string
!            END DO
!	 ENDIF	  
! 40  	 CLOSE(30)
!	 CLOSE(31)

!         OPEN(30,FILE='res_tisc.xyw',STATUS='OLD',ACTION='READ',
!     +			IOSTAT=iexist)
!	 IF(iexist==0) THEN
!	    OPEN(31,FILE=file_xyw)
!	    DO i=1,200000
!		READ(30,"(A115)",END=60) string
!		WRITE(31,"(A115)") string
!            END DO
!	 ENDIF	  
! 60  	 CLOSE(30)
!	 CLOSE(31)
      
          INQUIRE(FILE='res_tisc.st',EXIST=exist_file)
	    IF(exist_file) call system ('mv res_tisc.st '//file_st)
!!	    IF(exist_file) print *, system ('mv res_tisc.st '//file_st)
          INQUIRE(FILE='res_tisc.xyw',EXIST=exist_file)
	    IF(exist_file) call system ('mv res_tisc.xyw '//file_xyw)
!!	    IF(exist_file) print *, system ('mv res_tisc.xyw '//file_xyw)
      ENDIF

      RETURN
      END SUBROUTINE call_uhuru_tisc

       
!! ********************************************************************
      SUBROUTINE new_position (position_x, position_y,
     +				m, n, nn, Dx, Dy, u, v, Dtsegons)  
!  Calculate the new position of the point (position_x, position_y) 
!	under the velocity field (u,v)
!  Input: position_x, position_y, m, n, nn, Dx, Dy, u, v, Dtsegons	
!  Output: position_x, position_y

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION u(nn), v(nn)
	
      zix=position_x/Dx
      ziy=position_y/Dy
      ix=aint(zix)
      iy=aint(ziy)
      k=ix+1+iy*(n+1)
      IF(ix<0.or.ix>n.or.iy<0.or.iy>m) RETURN
      d1=dsqrt(Dx*((zix-ix)**2)+Dy*((ziy-iy)**2))
      d2=dsqrt(Dx*((zix-ix-1)**2)+Dy*((ziy-iy)**2))
      d3=dsqrt(Dx*((zix-ix)**2)+Dy*((ziy-iy-1)**2))
      d4=dsqrt(Dx*((zix-ix-1)**2)+Dy*((ziy-iy-1)**2))
      tnum1=d2*d3*d4*u(k)+d1*d3*d4*u(k+1)
      tnum2=d1*d2*d4*u(k+n+1)+d1*d2*d3*u(k+n+2)
      tden=d2*d3*d4+d1*d3*d4+d1*d2*d4+d1*d2*d3
      upunt=(tnum1+tnum2)/tden
      tnum1=d2*d3*d4*v(k)+d1*d3*d4*v(k+1)
      tnum2=d1*d2*d4*v(k+n+1)+d1*d2*d3*v(k+n+2)
      vpunt=(tnum1+tnum2)/tden
      despx=upunt*Dtsegons	!! Displacement (despx,despy) in meters
      despy=vpunt*Dtsegons
      position_x=zix*Dx+despx
      position_y=ziy*Dy+despy

      END SUBROUTINE new_position 

!!! *******************************************************************
      SUBROUTINE TIME_BCFILE (BCtime, nsBC, BCMy, BCfNew, NBC)
!!  Read from a file the time(Myear) and the new boundary conditions file
!!     BCMy[NBC] (array of the time (myear))	BCfNew[NBC] (array of character)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BCMy(nsBC),BCfNew(nsBC)
      CHARACTER*84 BCtime,BCfNew,b

      OPEN(13,FILE=BCtime,STATUS='OLD',ACTION='READ')
      NBC=0  
      DO i=1,nsBC
	  READ(13,*,END=440) a,b
	  NBC=NBC+1
	  BCMy(i)=a
	  BCfNew(i)=b
      END DO	  
 440  WRITE(6,"(3X,'Different Boundary conditions:',I3)") NBC
      CLOSE(13)

      RETURN
      END SUBROUTINE TIME_BCFILE
!! **********************************************************************
      SUBROUTINE fit_gradients (m, n, nn, Dx, Dy, Dif_K, dz_allow,
     +				zvalue, zname, iwr_filt)
! Control of the horizontal gradients of the function zvalue
! If iwr_filt=-1  => Write on the screen
! If iwr_filt=0	  => No write on the screen and return the number of filters
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION zvalue(nn)
      CHARACTER*55 zname
      PARAMETER (max_it=20)

      ddx=2.D0*Dx
      ddy=2.D0*Dy
      nfilters=0
      DO iter=1,max_it
        dz_max=0.D0
        DO iy=1,m-1
          DO ix=1,n-1
           kxy=ix+1+iy*(n+1)
	   dz_dx=ABS((zvalue(kxy+1)-zvalue(kxy-1))/ddx)
	   dz_dx_pc=dz_dx/zvalue(kxy)
	   dz_dy=ABS((zvalue(kxy+n+1)-zvalue(kxy-n-1))/ddy)
	   dz_dy_pc=dz_dy/zvalue(kxy)
           dz_max=MAX(dz_max,dz_dx_pc,dz_dy_pc)         
          END DO
        END DO
!	   WRITE(6,"(8X,'SUBR.fit_gradients, ',A55,':  dz_max=',G13.6,
!     +			'   dz_allow:',G13.6)") zname, dz_max, dz_allow
         IF(dz_max > dz_allow) THEN
     	       nfilters=nfilters+1   
               CALL diffusive_filter (m, n, nn, Dx, Dy, Dif_K, zvalue)	!! lib/lib_uhuru/mathematics.f
               CYCLE
          ELSE
               EXIT 
         ENDIF
      END DO
      IF(nfilters/=0 .AND. iwr_filt==-1) THEN
	   z_min=VAL_MIN(zvalue,nn)
	   z_max=VAL_MAX(zvalue,nn)
	   WRITE(6,"(8X,'Apply ',I4,' times (',I3,' max) the diffusive',
     +	     ' filter to ',A55/11X,'maximum gradient allowed:',1P,G10.3,
     +	      3X,'Constant diffusive filter applied:',1P,G10.3/
     +	      11X,'minim value:',1P,G12.4,5X,'maxim value: ',1P,G12.4)")
     +		nfilters, max_it, zname,dz_allow, Dif_K, z_min, z_max
       ELSE
	   iwr_filt=nfilters
      ENDIF

      RETURN
      END SUBROUTINE fit_gradients

!! **********************************************************************
!!      SUBROUTINE fit_gradients_3D (m, n, NELZ, nn, Dx, Dy, Dz, Dif_K,
!     +					dT_allow, Tvalue, zname)
!! Control of the gradients of the function Tvalue
!      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      DIMENSION Tvalue(nn,0:NELZ)
!      CHARACTER*55 zname
!      PARAMETER (max_it=5)

!      ddx=2.D0*Dx
!      ddy=2.D0*Dy
!      ddz=2.D0*Dz
!      nfilters=0
!      DO iter=1,max_it
!        dT_max=0.D0
!	DO iz=1,NELZ-1
!	   DO iy=1,m-1
!	      DO ix=1,n-1
!		kxy=ix+1+iy*(n+1)
!		dT_dx=ABS((Tvalue(kxy+1,iz)-Tvalue(kxy-1,iz))/ddx)
!		dT_dx_pc=dT_dx/Tvalue(kxy,iz)
!		dT_dy=ABS((Tvalue(kxy+n+1,iz)-Tvalue(kxy-n-1,iz))/ddy)
!		dT_dy_pc=dT_dy/Tvalue(kxy,iz)
!		dT_max=MAX(dT_max,dT_dx_pc,dT_dy_pc)         
!		!dT_dz=ABS((Tvalue(kxy,iz+1)-Tvalue(kxy,iz-1))/ddz)
!		!dT_dz_pc=dT_dz/Tvalue(kxy,iz)
!		!dT_max=MAX(dT_max,dT_dx_pc,dT_dy_pc,dT_dz_pc)         
!	      END DO
!	   END DO
!	END DO
!        IF(dT_max > dT_allow) THEN
!     	       nfilters=nfilters+1   
!               CALL diffusive_filter_3D (m, n, NELZ, nn, Dx, Dy, Dz,		!! lib/lib_uhuru/mathematics.f
!     +						Dif_K, Tvalue)
!               CYCLE
!          ELSE
!               EXIT 
!        ENDIF
!      END DO
!      IF(nfilters/=0) THEN
!	 Tmin=1.D200
!	 DO kxy=1,nn
!	    DO iz=0,NELZ
!		Tmin=MIN(Tvalue(kxy,iz),Tmin)
!	    END DO
!	 END DO
!	 Tmax=-1.D200
!	 DO kxy=1,nn
!	    DO iz=0,NELZ
!		Tmax=MAX(Tvalue(kxy,iz),Tmin)
!	    END DO
!	 END DO

!	   WRITE(6,"(8X,'Apply ',I4,' times (',I3,' max) the diffusive',
!     +	     ' filter to ',A55/11X,'maximum gradient allowed:',1P,G10.3,
!     +	      3X,'Constant diffusive filter applied:',1P,G10.3/
!     +	      11X,'minim value:',1P,G12.4,5X,'maxim value: ',1P,G12.4)")
!     +		nfilters, max_it, zname,dT_allow, Dif_K, Tmin, Tmax
!      ENDIF
!
!      RETURN
!      END SUBROUTINE fit_gradients_3D

!! **********************************************************************
      SUBROUTINE RE_MESH (m, n, nn, NELZ, kpuntssd, Dx, Dy, elevation,
     +			  scrust, sediment, GLit, GLit_ter, u, v,
     +			  epuntzz, visco, Tmoho, Qsurface,
     +			  erosion, TEMPE, D_Bodies, Points_t)

! Move the model to the south, deleting the firts row (iy=0) and adding a new row to the north

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION elevation(nn), scrust(nn), sediment(nn), GLit(nn),
     +		GLit_ter(nn), u(nn), v(nn), epuntzz(nn), visco(nn),
     +		erosion(nn), Tmoho(nn), Qsurface(nn), TEMPE(nn,0:NELZ), 
     +		D_Bodies(kpuntssd,2), Points_t(0:kpuntssd,2)

!      DO iy=0,m
!	   DO ix=0,n
!		kxy=ix+1+iy*(n+1)
!		WRITE(6,"(3X,2I3,'  elevation:',F9.1,'  visco:',G13.6)")
!     +				ix,iy,elevation(kxy), visco(kxy)
!	   END DO
!      END DO
!      kpunts1=D_Bodies(1,1)
!       DO kp=2,kpunts1+1					!! BODY 1
!      WRITE(6,"(3X,'D_Bodies 1:',2F12.1)") D_Bodies(kp,1),D_Bodies(kp,2)
!       END DO
!	kpunts2=D_Bodies(kpunts1+2,1)
!       DO kp=kpunts1+3,kpunts1+kpunts2+2			!! BODY 2
!      WRITE(6,"(3X,'D_Bodies 2:',2F12.1)") D_Bodies(kp,1),D_Bodies(kp,2)
!       END DO
!       kpunts3=D_Bodies(kpunts1+kpunts2+3,1)
!       DO kp=kpunts1+kpunts2+4,kpunts1+kpunts2+kpunts3+3	!! BODY 3
!      WRITE(6,"(3X,'D_Bodies 3:',2F12.1)") D_Bodies(kp,1),D_Bodies(kp,2)
!       END DO

      DO ix=0,n
	DO iy=0,m-1
	   kxy=ix+1+iy*(n+1)
	   kxy_sup=ix+1+(iy+1)*(n+1)
	   elevation(kxy)=elevation(kxy_sup)
	   scrust(kxy)=scrust(kxy_sup)
	   sediment(kxy)=sediment(kxy_sup)
	   GLit(kxy)=GLit(kxy_sup)
	   GLit_ter(kxy)=GLit_ter(kxy_sup)
	   u(kxy)=u(kxy_sup)
	   v(kxy)=v(kxy_sup)
	   epuntzz(kxy)=epuntzz(kxy_sup)
	   visco(kxy)=visco(kxy_sup)
	   Tmoho(kxy)=Tmoho(kxy_sup)
	   Qsurface(kxy)=Qsurface(kxy_sup)
	   erosion(kxy)=erosion(kxy_sup)
	   DO iz=0,NELZ
		TEMPE(kxy,iz)=TEMPE(kxy_sup,iz)
	   END DO
	END DO
      END DO
      DO ix=0,n	!! iy=m
	   kxy=ix+1+m*(n+1)
	   kxy_inf=ix+1+(m-1)*(n+1)
	   elevation(kxy)=elevation(kxy_inf)
	   scrust(kxy)=scrust(kxy_inf)
	   sediment(kxy)=sediment(kxy_inf)
	   GLit(kxy)=GLit(kxy_inf)
	   GLit_ter(kxy)=GLit_ter(kxy_inf)
	   u(kxy)=u(kxy_inf)
	   v(kxy)=v(kxy_inf)
	   epuntzz(kxy)=epuntzz(kxy_inf)
	   visco(kxy)=visco(kxy_inf)
	   Tmoho(kxy)=Tmoho(kxy_inf)
	   Qsurface(kxy)=Qsurface(kxy_inf)
	   erosion(kxy)=erosion(kxy_inf)
	   DO iz=0,NELZ
		TEMPE(kxy,iz)=TEMPE(kxy_inf,iz)
	   END DO
      END DO

	!!! New y position of the different points
      	  kpunts1=D_Bodies(1,1)
	  DO kp=2,kpunts1+1					!! BODY 1
		IF(D_Bodies(kp,2)>0.0 .AND. D_Bodies(kp,1)>0.0) 
     +		     D_Bodies(kp,2)=MAX(0.D0,D_Bodies(kp,2)-Dy)
          END DO
	  kpunts2=D_Bodies(kpunts1+2,1)
          DO kp=kpunts1+3,kpunts1+kpunts2+2			!! BODY 2
		IF(D_Bodies(kp,2)>0.0 .AND. D_Bodies(kp,1)>0.0) 
     +		     D_Bodies(kp,2)=MAX(0.D0,D_Bodies(kp,2)-Dy)
          END DO
	  kpunts3=D_Bodies(kpunts1+kpunts2+3,1)
          DO kp=kpunts1+kpunts2+4,kpunts1+kpunts2+kpunts3+3	!! BODY 3
		IF(D_Bodies(kp,2)>0.0 .AND. D_Bodies(kp,1)>0.0) 
     +		     D_Bodies(kp,2)=MAX(0.D0,D_Bodies(kp,2)-Dy)
          END DO

          kppoints=Points_t(0,1)				!! Points
	  DO kp=1,kppoints
		Points_t(kp,2)=Points_t(kp,2)-Dy
          END DO


!	PRINT*,'   -------------------------------------------'
!	PRINT*,'  After remeshing'
!		DO iy=0,m
!		DO ix=0,n
!			kxy=ix+1+iy*(n+1)
!	WRITE(6,"(3X,'after',2I3,'  elevation:',F9.1,'  visco:',G13.6)")
!     +				ix,iy,elevation(kxy), visco(kxy)
!		END DO
!		END DO
!	kpunts1=D_Bodies(1,1)
!	DO kp=2,kpunts1+1					!! BODY 1
!      WRITE(6,"(3X,'D_Bodies 1:',2F12.1)") D_Bodies(kp,1),D_Bodies(kp,2)
!          END DO
!	  kpunts2=D_Bodies(kpunts1+2,1)
!          DO kp=kpunts1+3,kpunts1+kpunts2+2			!! BODY 2
!      WRITE(6,"(3X,'D_Bodies 2:',2F12.1)") D_Bodies(kp,1),D_Bodies(kp,2)
!          END DO
!	  kpunts3=D_Bodies(kpunts1+kpunts2+3,1)
!          DO kp=kpunts1+kpunts2+4,kpunts1+kpunts2+kpunts3+3	!! BODY 3
!      WRITE(6,"(3X,'D_Bodies 3:',2F12.1)") D_Bodies(kp,1),D_Bodies(kp,2)
!          END DO


      RETURN
      END SUBROUTINE RE_MESH





