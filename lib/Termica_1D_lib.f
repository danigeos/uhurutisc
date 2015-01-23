
      SUBROUTINE Temp_steadystate_1D (Ndepth, Dz, BANDAKZ, CONDUC, 
     +                          HSURF, HEXP, PHEAT_m, TSURF, TBOTT, 
     +			        Zsedim, ZMOHO, ZLITOS,
     +                          Temp_z, NZLIT1, NZLIT2, iswitch)

!  Temperature at steady state, 1-D. Calculation with finite differences.
!  Boundary Conditions: Temperature at the surface and asthenosphere
!         T(z=0)=TSURF and T(z=ZASTH)=TBOTT

! Ndepth :      Number of vertical points 0, ..., NELZ.
! Dz [m] :      Vertical increment.
! BANDAKZ [m] : Conductivity transition zone, up and down of ZLITOS (2*BANDAKZ),
!                 between lithospheric mantle and asthenosphere.
! CONDUC(0:3)[W/m*K]: Conductivities sediments, crust, lith.mantle, asth.
! HSURF [W/m**3] :   Surface crustal heat production, H=HSURF*exp(-z/HEXP)
! HEXP [m] :         Exponent crustal heat production. 
! PHEAT_m [W/m**3] : Lithospheric mantle heat production.
! TSURF [K]  :       Surface Temperature (z=0)    
! TBOTT [K]  :       Bottom Temperature (z=ZASTH)
! Zsedim [m] :       Bottom of the sediments.
! ZMOHO [m]  :       Bottom of the crust (moho depth)
! ZLITOS [m]  :      Bottom of the lithospheric mantle.
! Temp_z(0:Ndepth)[K]: Geotherm (output)
! NZLIT1, NZLIT2:      Points out of the transition zone (output)

!  iswitch = 1	=> keep in a file K(z) and H(z)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ASTEADY(Ndepth,3),BSTEADY(Ndepth),CONDUC(0:3),
     +           xtemp(Ndepth),Temp_z(0:Ndepth)
      NELZ=Ndepth-1
      ZASTH=Dz*NELZ      

      CALL Transition_lith_asth (ZMOHO, ZLITOS, ZASTH, BANDAKZ,
     +                            ZLITOS1, ZLITOS2, NZLIT1, NZLIT2)

      IF(iswitch == 1) OPEN(20,FILE='H_K_steadystate.ter')
      BSTEADY=0.D0
      ASTEADY=0.D0
!!!------------- Construction of the matrix -----------------
!!	Surface BC,  iz=0 (Temperature(0)=TSURF) :
      ASTEADY(1,2)=1.D0
      BSTEADY(1)=TSURF

!!	Internal Points,  iz=1,NELZ-1 :
      Depth: DO iz=1,NELZ-1
                  Z=iz*Dz
                  Zp=(iz+1)*Dz
                  Zn=(iz-1)*Dz		  
                  CONDI=Conductivity (Z, CONDUC,
     +                          Zsedim, ZMOHO, ZLITOS, ZLITOS1, ZLITOS2)
                  CONDIp=Conductivity (Zp, CONDUC,
     +                          Zsedim, ZMOHO, ZLITOS, ZLITOS1, ZLITOS2)
                  CONDIn=Conductivity (Zn, CONDUC,
     +                          Zsedim, ZMOHO, ZLITOS, ZLITOS1, ZLITOS2)
                  PROHEATI=Heat_Production (Z, HSURF, HEXP, PHEAT_m,
     +                                        Zsedim, ZMOHO, ZLITOS)
                  IF(iswitch == 1) WRITE(20,"(F10.2,G14.5,F8.3)") 
     +						Z,PROHEATI,CONDI
                  FAC1=(CONDIp-CONDIn)/(4.D0*Dz*Dz)
                  FAC2=CONDI/(Dz*Dz)
                  ASTEADY(iz+1,1)=FAC2-FAC1
                  ASTEADY(iz+1,2)=-2.D0*FAC2
                  ASTEADY(iz+1,3)=FAC1+FAC2
                  BSTEADY(iz+1)=-PROHEATI
      END DO Depth

!!  Bottom BC,  iz=NELZ  (Temperature(NELZ)=TBOTT) :
      ASTEADY(NELZ+1,2)=1.D0
      BSTEADY(NELZ+1)=TBOTT
!!!----------------------------------------------------------------
      CALL sistbanda(ASTEADY,BSTEADY,Ndepth,3,1,xtemp)
      DO iz=0,NELZ
      	Temp_z(iz)=xtemp(iz+1)
      END DO
      IF(iswitch == 1) CLOSE(20)
      RETURN
      END SUBROUTINE Temp_steadystate_1D

!! **********************************************************************

      SUBROUTINE Temp_Dt1D (Dtemps, Ndepth, Dz, BANDAKZ, CONDUC, 
     +                     HSURF, HEXP, PHEAT_m, THDIFF, Zsedim,
     +                     ZMOHO, ZLITOS, strain_rate_zz, TADVEC,
     +                     Temp_z, Dsed_eros, NZLIT1, NZLIT2, iswitch)
      
!!	Temperature variation on Dt, 1-D.  
!! Boundary Conditions: Surface (iz=0) and bottom (iz=NELZ) temperature

! Dtemps [s] :  Time increment
! Ndepth :      Number of vertical points 0, ..., NELZ.
! Dz [m] :      Vertical increment.
! BANDAKZ [m] : Conductivity transition zone, up and down of ZLITOS (2*BANDAKZ),
!                 between lithospheric mantle and asthenosphere.
! CONDUC(0:3)[W/m*K]: Conductivities sediments, crust, lith.mantle, asth.
! HSURF [W/m**3] :   Surface crustal heat production, H=HSURF*exp(-z/HEXP)
! HEXP [m] :         Exponent crustal heat production. 
! PHEAT_m [W/m**3] : Lithospheric mantle heat production.
! THDIFF [m**2/s] :  Thermal Diffusivity
! Zsedim [m] :       Bottom of the sediments.
! ZMOHO [m]  :       Bottom of the crust (moho depth)
! ZLITOS [m]  :      Bottom of the lithospheric mantle.
! strain_rate_zz [s-1] : Vertical strain rate
! TADVEC : 	     Horizontal Advective term, -u*(dT/dx)-v(dT/dy)
! Temp_z(0:Ndepth)[K]: Geotherm (input and output)
! NZLIT1, NZLIT2:      Points out of the transition zone (input and output)
! Dsed_eros:	Change of sediment thickness(+) or erosion thickness(-) in this time step

!  iswitch = 1	=> keep in a file K(z) and H(z)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(TETA=0.5D0)
      DIMENSION ATERMICA(Ndepth,3),BTERMICA(Ndepth),TADVEC(Ndepth),
     +           CONDUC(0:3),xtemp(Ndepth),Temp_z(0:Ndepth), 
     +		 Temp_new(0:Ndepth)
     
      NELZ=Ndepth-1
      ZASTH=Dz*NELZ
      CALL Transition_lith_asth (ZMOHO, ZLITOS, ZASTH, BANDAKZ,
     +                            ZLITOS1, ZLITOS2, NZLIT1, NZLIT2)
	     
!!  New surface, depending of the Increment of sediment thickness(+) or erosion thickness(-) in this time step (Dsed_eros)
      
      IF(Dsed_eros>0) THEN		!! Sedimentation
	 DO iz=1, NELZ-1
	     Z=iz*Dz-Dsed_eros
	     IF(Z<0)THEN
		Temp_new(iz)=Temp_z(0)
	     ELSE
		Temp_new(iz)=Temperature_DepthZ(Z,Temp_z,NELZ,Dz)
	     ENDIF
	 END DO
	   	
      ELSE				!! Erosion
	 DO iz=1, NELZ-1
	     Z=ABS(Dsed_eros)+iz*Dz
	     Temp_new(iz)=Temperature_DepthZ(Z,Temp_z,NELZ,Dz)
	 END DO
      ENDIF
		
      DO iz=1, NELZ-1
	 Temp_z(iz)=Temp_new(iz)
      END DO
!!!!!!
      BTERMICA=0.D0
      ATERMICA=0.D0

      A2=THDIFF/(Dz*Dz)
!!  Surface Boundary Condition, iz=0 (Temperature(0,t)=Temperature(0,t+Dt) :
      ATERMICA(1,2)=1.D0
      BTERMICA(1)=Temp_z(0)
      
      IF(iswitch == 1) OPEN(20,FILE='H_K_Dt.ter')
!!  Internal Points,  iz=1,NELZ-1 :
      Depth: DO iz=1,NELZ-1
                Z=iz*Dz
                Zp=(iz+1)*Dz
                Zn=(iz-1)*Dz
                CONDI=Conductivity (Z, CONDUC,
     +                          Zsedim, ZMOHO, ZLITOS, ZLITOS1, ZLITOS2)
                CONDIp=Conductivity (Zp, CONDUC,
     +                          Zsedim, ZMOHO, ZLITOS, ZLITOS1, ZLITOS2)
                CONDIn=Conductivity (Zn, CONDUC,
     +                          Zsedim, ZMOHO, ZLITOS, ZLITOS1, ZLITOS2)
                PROHEATI=Heat_Production (Z, HSURF, HEXP, PHEAT_m,
     +                                        Zsedim, ZMOHO, ZLITOS)
                IF(iswitch == 1) WRITE(20,"(F10.2,G14.5,F8.3)") 
     +						Z,PROHEATI,CONDI
                A1=-(strain_rate_zz*Z)/(2.D0*Dz)
                A3=(A2/(4.D0*CONDI))*(CONDIp-CONDIn)
                DTPT=(1.D0-TETA)*(Temp_z(iz+1)-Temp_z(iz-1))
                BT1=A1*DTPT
                BT2=A2*(1.D0-TETA)*(Temp_z(iz+1)-
     +                      2.D0*Temp_z(iz)+Temp_z(iz-1))
                BT3=A3*DTPT
                BT4=(THDIFF/CONDI)*PROHEATI
                BT=BT1+BT2+BT3+BT4+TADVEC(iz)
                ATER2=1.D0+2.D0*A2*TETA*Dtemps
                ATERMICA(iz+1,1)=(TETA*(A1-A2+A3)*Dtemps)/ATER2
                ATERMICA(iz+1,2)=1.D0
                ATERMICA(iz+1,3)=(-TETA*(A1+A2+A3)*Dtemps)/ATER2
                BTERMICA(iz+1)=(Temp_z(iz)+BT*Dtemps)/ATER2
       END DO Depth 
!!  Bottom Boudary Conditions,  iz=NELZ 
!!		(Temperature(NELZ,t)=Temperature(NELZ,t+Dt)) :
       ATERMICA(NELZ+1,2)=1.D0
       BTERMICA(NELZ+1)=Temp_z(NELZ)

       CALL sistbanda(ATERMICA,BTERMICA,Ndepth,3,1,xtemp)
       !Tcont_max=xtemp(NELZ+1)-100.D0
       !ierror_T=0 
       Temp_z(0)=xtemp(1)
       DO iz=1,NELZ
           Temp_z(iz)=xtemp(iz+1)
	   Temp_z(iz)=MIN(Temp_z(iz),xtemp(NELZ+1))	   
	   Temp_z(iz)=MAX(Temp_z(iz),Temp_z(iz-1))
	   !IF(Temp_z(iz)<Tcont_max .AND. Temp_z(iz)<Temp_z(iz-1)) THEN
	   !IF(Temp_z(iz)<Temp_z(iz-1)) THEN
	   !	ierror_T=1
	   !	Temp_z(iz)=Temp_z(iz-1)
	   !END IF
       END DO
      IF(iswitch == 1) CLOSE(20)
      !IF(ierror_T==1) THEN
	!    WRITE(20,"(1X,2F10.2)") (iz*Dz,Temp_z(iz), iz=0,NELZ)
      !END IF
      RETURN
      END SUBROUTINE Temp_Dt1D

!! **********************************************************************
      SUBROUTINE Temp_steadystate_TQsurf_1D (NELZ, Dz, CONDUC,
     +                                 HSURF, HEXP, PHEAT_m, 
     +			               Zsedim, ZMOHO, ZLITOS, Qsurface,
     +                                 TSURF, TLITOS, Temp_z, iswitch)

!  Temperature at steady state, 1-D. Calculation with finite differences.
!  Boundary Conditions: Temperature and heat flow at the surface
!          T(z=0)=TSURF and Q(z=0)=Qsurface
!  No considero sediments

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CONDUC(0:3),Temp_z(0:NELZ)

      IF(iswitch == 1) OPEN(20,FILE='H_K_steadystate.ter')
      ZASTH=Dz*NELZ
      ZLITOS1=ZLITOS
      ZLITOS2=ZLITOS

!!  Boundary Conditions (TSURF,Qsurface),  iz=0,1 :
             Temp_z(0)=TSURF
             Temp_z(1)=Temp_z(0)+(Qsurface/CONDUC(1))*Dz 
!!  Internal Points,  iz=2,NELZ :
      Depth: DO iz=2,NELZ
                  Z=iz*Dz
                  Zn=(iz-1)*Dz		  
                  CONDI=Conductivity (Z, CONDUC,
     +                          Zsedim, ZMOHO, ZLITOS, ZLITOS1, ZLITOS2)
                  CONDIn=Conductivity (Zn, CONDUC,
     +                          Zsedim, ZMOHO, ZLITOS, ZLITOS1, ZLITOS2)
                  PROHEATI=Heat_Production (Z, HSURF, HEXP, PHEAT_m,
     +                                        Zsedim, ZMOHO, ZLITOS)
                  IF(iswitch == 1) WRITE(20,*) Z,PROHEATI,CONDI

                  IF(Temp_z(iz-1)>TLITOS) CONDI=CONDUC(3) 
                  IF(Temp_z(iz-2)>TLITOS) CONDIn=CONDUC(3) 
                  T1=(CONDI+CONDIn)*Temp_z(iz-1)
                  T2=CONDIn*Temp_z(iz-2)
                  T3=Dz*Dz*PROHEATI 
                  Temp_z(iz)=(T1-T2-T3)/CONDI
      END DO Depth
	IF(Temp_z(NELZ)<TLITOS) THEN
	      iswitch=999
     	    WRITE(6,"(/6X,'The deeper temperature calculated',F7.1,' K'/
     +		6X,'is lower than the thermal lithosphere (TLITOS):',
     +		F7.1,' K',/6X,'You have to give a higher value to ',
     +	      'ZASTH, actually',F8.2,' km'/6X,'Qsurface:',F7.2,' mW/m2'/
     +	    6X,'ZMOHO:',F7.2,' km.  ZLITOS:',F8.2,' km')") Temp_z(NELZ),
     +	       TLITOS,ZASTH/1D3,Qsurface*1D3,ZMOHO/1D3,ZLITOS/1D3
        ENDIF

      IF(iswitch == 1) CLOSE(20)
      RETURN
      END SUBROUTINE Temp_steadystate_TQsurf_1D

!! **********************************************************************
      SUBROUTINE Transition_lith_asth (ZMOHO, ZLITOS, ZASTH, BANDAKZ,
     +                                 ZLITOS1, ZLITOS2, NZLIT1, NZLIT2)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
        ZLITOS1=ZLITOS-BANDAKZ    ! Transition zone between mantle
        ZLITOS2=ZLITOS+BANDAKZ    !  lithosphere and asthenosphere.
        IF(ZLITOS1.LT.ZMOHO) THEN 
            ZLITOS1=ZMOHO
            NZLIT1=NZLIT1+1
        ENDIF
        IF(ZLITOS2.GT.ZASTH) THEN 
            ZLITOS2=ZASTH
            NZLIT2=NZLIT2+1
        ENDIF
      RETURN
      END SUBROUTINE Transition_lith_asth

!! ***************************************************************************
      SUBROUTINE Convective_Removal_Type (TEMPSMa, Time_removal,
     +			  Iremoval, Zremoval, TISO_rem, Zctall, ZMOHO, 
     +			  hlitos, Temp_z, TISOTER, DGL_ter, NELZ, Dz,
     +			  Iremoved, A_removal, B_removal)

  !    SUBROUTINE Convective_Removal_Type (TEMPSMa, DtMy, Time_removal,
  !   +			DTime_removal,Iremoval, Zremoval, TISO_rem, 
  !   +			Zctall, ZMOHO, hlitos, Temp_z, TISOTER,
  !   +				DGL_ter, NELZ, Dz, Iremoved)

!! Mode of Convective removal of the lithosphere:
!!  Iremoval=0 => No convective removal
!!  Iremoval=1 => At time=Time_removal, convective removal where lithosphere is deeper than Zctall at L=L(T=TISO_rem)
!!  Iremoval=2 => At time=Time_removal, convective removal where lithosphere is deeper than Zctall at Zremoval
!!  Iremoval=3 => Starting at time=Time_removal, continuos convective removal to Zremoval where lithosphere is deeper than Zctall until Time=t2
!!  Iremoval=4 => At time=Time_removal, convective removal where crust is deeper than Zctall at L=Zremoval
!!  Iremoval=5 => At time=Time_removal, convective removal where crust is deeper than Zctall at L=lineal transition 
!!  Iremoval=6 => At time=Time_removal, convective removal where crust is deeper than Zctall at L=L(T=TISO_rem)

!!  DTime_removal/=0  =>  Slow Lithosphere root removal, in DTime_removal [m.y.]   =0(instantaneous)
!!				Works with LITHOSPHERE = Zremoval, not with an isotherm. (Iremoval/=1,6)

!!  If return Iremoved=1 => the lithosphere has been convective removed

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Temp_z(0:NELZ)
      PARAMETER (zero=1.D-20)
      DATA DTime_removal /0.D0/

      IF(Iremoval==1.OR.Iremoval==6) DTime_removal=0.D0

      Time_removal_step=Time_removal
      Zremoval_step=Zremoval	    
      Zctall_P=Zctall

      ISwich_startR=0
      Time_int=ABS(TEMPSMa-Time_removal)		
      IF(Time_int<1D-3 .AND. Iremoval/=0) ISwich_startR=1
      IF(ISwich_startR==1 .AND. DTime_removal>0) THEN
	   B_removal=(Zremoval-hlitos)/DTime_removal
	   A_removal=hlitos-B_removal*Time_removal
	   Time_removal_step=TEMPSMa
	   Zremoval_step=A_removal+B_removal*TEMPSMa
      ENDIF

      Time_end=Time_removal+DTime_removal
      IF(ISwich_startR/=1 .AND. TEMPSMa>=Time_removal .AND.
     +		TEMPSMa<=Time_end .AND. ABS(B_removal)>zero) THEN
	   Time_removal_step=TEMPSMa
	   Zremoval_step=A_removal+B_removal*TEMPSMa
	   Zctall_P=9.D3
      ENDIF
      IF(ISwich_startR/=1.AND.DTime_removal>0.AND.ABS(B_removal)<zero)
     +		RETURN
      
!	Time_end=Time_removal+DTime_removal		
!        IF(TEMPSMa>=Time_removal .AND. TEMPSMa<Time_end) THEN
!		Time_removal_step=TEMPSMa
!		steps=(Time_end-TEMPSMa)/DtMy
!		IF(steps<=1) THEN
!		     Zremoval_step=Zremoval
!		 ELSE
!		     Delta_Litos=(Zremoval-hlitos)/steps
!		     Zremoval_step=hlitos+Delta_Litos
!		ENDIF
!	ENDIF

      t1=Time_removal_step-1D-3
      t2=Time_removal_step+1D-3
  !    t2=Time_removal_step+10.D0		!! Pluma mantelica durante 10My
      IF(Iremoval==3) t2=1.D4						!! Continuos Convective removal
      IF(TEMPSMa>t1 .AND. TEMPSMa<t2 .AND. Iremoval/=0) THEN
      
         ZISOTER=Depth_isotherm (TISOTER,Temp_z,NELZ,Dz)
         Zl=DGL_ter+ZISOTER						!! Zl=ZLITOS
	 Zrem=Zremoval_step	 
	 IF(Iremoval==1 .OR. Iremoval==6)  				!! New lithosphere depth where T=TISO_rem		
     +		Zrem=Depth_isotherm (TISO_rem,Temp_z,NELZ,Dz)

         IF(Iremoval>=1.AND.Iremoval<=3 .AND. ZISOTER>Zctall_P) THEN	!! Convective removal, for deeper lithospheres (Zremoval_step)
		CALL Convective_Removal(Dz,NELZ,Zrem,Zl,Temp_z)
!		WRITE(6,"(5X,' ConvectiveRemoval: Lithosphere to',F6.1,
!     +		    ' km at',F6.1,' My')") Zrem/1.D3,TEMPSMa
     		Iremoved=1
         END IF

         IF(Iremoval>=4 .AND. ZMOHO>Zctall_P) THEN			!! Convective removal, under thick crust (Zctall)
		IF(Iremoval==5) THEN 					!! lineal transition between 65.D3 and Zcm. All lineal => Zcm=Zctall
		     Zcm=55.D3					
		     A=(120.D3-Zremoval_step)/(65.D3-Zcm)
		     B=Zremoval_step-A*Zcm
		     IF(ZMOHO>Zcm) Zrem=A*ZMOHO+B
		     Zrem=ZMOHO+35.D3
		END IF	
		CALL Convective_Removal(Dz,NELZ,Zrem,Zl,Temp_z)	
		WRITE(6,"(5X,' ConvectiveRemoval: Lithosphere to',F6.1,
     +		    ' km, where moho is deeper than ',F6.1,' km at',
     +		    F6.1,' My')") Zrem/1.D3,Zctall_P/1.D3,TEMPSMa
     		Iremoved=1
         END IF

      END IF

      IF(Iremoved==0) THEN 
	   B_removal=0.D0
	   A_removal=0.D0
      END IF
      
      RETURN
      END SUBROUTINE Convective_Removal_Type

!! **********************************************************************
      SUBROUTINE Convective_Removal (Dz, NELZ, Zremoval, Zl, 
     +					Temp_z)

!! Convective removal of the lithosphere deper than Zremoval
!!	==> T(z) = TISOTER  if   Zremoval < z < Zl 
!!		Zl=ZLITOS, conductivity-depth-change from mantle to asthenosphere
!! Same process that Figure 3 on:
!!	 Platt and England, 1994, American Journal of Science, v 294. 

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Temp_z(0:NELZ)
      
      izremov1=1+(Zremoval/Dz)
      Tzremoval=Temperature_DepthZ (Zremoval,Temp_z,NELZ,Dz)
      izremov2=Zl/Dz
      Tzl=Temperature_DepthZ (Zl,Temp_z,NELZ,Dz)

      iztrans=izremov1+10
      iztrans=MIN(iztrans,izremov2)
      Ztrans=iztrans*Dz

      A=(Tzl-Tzremoval)/(Ztrans-Zremoval)
      B=Tzremoval-A*Zremoval
      DO iz=izremov1,iztrans-1
           Temp_z(iz)=A*iz*Dz+B
      END DO
      DO iz=iztrans,izremov2
           Temp_z(iz)=Tzl
      END DO

      END SUBROUTINE Convective_Removal

!! **********************************************************************
!!!!   FUNCTIONS !!!!!

      Double Precision FUNCTION crust_isostasy (RHOH2O, roc, rom,
     +                                      RHOAST, elevation, hmantle)
!  Calcul of the crustal thickness from the lithospheric mantle thickness
!    and densisties, assuming local isostasy.

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (H0=2400.D0)  

          TE=elevation
          IF(elevation.LT.0.D0) TE=elevation*((RHOAST-RHOH2O)/RHOAST)
          crust_isostasy=(RHOAST/(RHOAST-roc))*(TE-
     +             ((RHOAST-rom)/RHOAST)*hmantle+H0)
      END FUNCTION

!!  ---------------------------------------------------------------------
      Double Precision FUNCTION elevation_isostasy (RHOH2O, rosed, roc,
     +                                      rom, RHOAST, sediment, 
     +                                      hcrust, hmantle)     
!  Calcul of the elevation from the crustal and lithospheric mantle thickness
!    and densisties, assuming local isostasy.

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (H0=2400.D0)  	

	FACELEV=RHOAST/(RHOAST-RHOH2O)
	hlitosfere=hmantle+hcrust+sediment
        rol=(rosed*sediment+roc*hcrust+rom*hmantle)/hlitosfere
	e=(((RHOAST-rol)/RHOAST)*hlitosfere)-H0
	IF(e.LT.0.D0) e=FACELEV*e
        elevation_isostasy=e

      END FUNCTION

!!  ---------------------------------------------------------------------
      Double Precision FUNCTION Temperature_DepthZ(Z,Temp_z,NELZ,Dz)
!        Find the temperature at depth Z, 
!	   from the geotherm Temp_z(iz) iz=0,...,NELZ

      	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Temp_z(0:NELZ)
          diz=Z/Dz
          aiz=AINT(diz)
          iz=INT(diz)

	  IF(iz>(NELZ-1)) THEN
		Temperature_DepthZ=Temp_z(NELZ)
		RETURN
	  ENDIF

          Temperature_DepthZ=(diz-aiz)*Temp_z(iz+1)+
     +                    (aiz+1-diz)*Temp_z(iz)
      END FUNCTION

!!  ---------------------------------------------------------------------
      Double Precision FUNCTION Depth_isotherm (T_iso,Temp_z,NELZ,Dz)
!!        Find the depth of the isotherm T_iso
      	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Temp_z(0:NELZ)
        DO iz=0,NELZ-1
           IF(Temp_z(iz).LE.T_iso.AND.Temp_z(iz+1).GT.T_iso) THEN
	   	d1=((T_iso-Temp_z(iz))/(Temp_z(iz+1)-Temp_z(iz)))*Dz
     		Depth_isotherm=DBLE(iz)*Dz+d1
		EXIT
           ENDIF
        END DO
      END FUNCTION

!!  ---------------------------------------------------------------------
      Double Precision FUNCTION Conductivity (Z, CONDUC,
     +                          Zsedim, ZMOHO, ZLITOS, ZLITOS1, ZLITOS2)
!!       Conductivity at depth Z
!! Zsedim, ZMOHO, ZLITOS: Sediment, crustal and lithospheric thicknesses
!! ZLITOS1, ZLITOS2 : Transition zone between the litosphere and the asthenosphere
!!     ZLITOS1 <= ZLITOS <= ZLITOS2
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION CONDUC(0:3)
       PARAMETER(PI=3.1415926535897932D0)
        IF(ZLITOS1/=ZLITOS2)COBS=PI/(2.D0*(ZLITOS2-ZLITOS1))	
        cof=CONDUC(3)-CONDUC(2)
        IF(Z.LT.Zsedim) Conductivity=CONDUC(0)
        IF(Z.GE.Zsedim.AND.Z.LT.ZMOHO) Conductivity=CONDUC(1)
        IF(Z.GE.ZMOHO.AND.Z.LT.ZLITOS) Conductivity=CONDUC(2)
        IF(Z.GE.ZLITOS) Conductivity=CONDUC(3)
        IF(Z.GT.ZLITOS1.AND.Z.LT.ZLITOS2)
     +           Conductivity=CONDUC(2)+cof*(DSIN(COBS*(Z-ZLITOS1)))**2

      END FUNCTION


!!  ---------------------------------------------------------------------
       Double Precision FUNCTION Heat_Production (Z,HSURF,HEXP,PHEAT_m,
     +                                            Zsedim,ZMOHO,ZLITOS)
!!       Heat Production at depth Z
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
          IF(Z.LE.Zsedim) Heat_Production=HSURF*DEXP(-Z/HEXP)
          IF(Z.GT.Zsedim.AND.Z.LE.ZMOHO) 
     +             Heat_Production=HSURF*DEXP(-Z/HEXP)
          IF(Z.GT.ZMOHO.AND.Z.LE.ZLITOS) Heat_Production=PHEAT_m
          IF(Z.GT.ZLITOS) Heat_Production=0.D0
      END FUNCTION

!!  ---------------------------------------------------------------------
      Double Precision FUNCTION Surface_Heat_Flow (Zsedim,ZMOHO,CONDUC,
     +                                            T0,T1,Dz)
!!       Surface Heat Flow
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
        DIMENSION CONDUC(0:3)
          CONDUC_z0=CONDUC(2)
          IF(Zsedim.NE.0) CONDUC_z0=CONDUC(0)
          IF(Zsedim.EQ.0.AND.ZMOHO.NE.0)
     +	                        CONDUC_z0=CONDUC(1)
          Surface_Heat_Flow=-(-CONDUC_z0*(T1-T0)/Dz)

      END FUNCTION
!!  ---------------------------------------------------------------------
