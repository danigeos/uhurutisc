******************************************************************************************************************************
	UHURU PROGRAM, files and parameters
******************************************************************************************************************************

Input files: 	fort.11			Initial geometry: x, y, crust and lithosphere thickness (Grid.in)
		parametres.in		Input parameters - open(8)
		fort.4			Velocity Boundary conditions: ix,iy,type of BC,condition1,condition2 (BC.in)
		fort.12 (optional)	Bodies: Zones with diferent strength (bodies.in)
		fort.13 (optional)	Time and files with different velocity boundary conditions (Time_BCfiles.in)
		fort.10 (optional)	Points which uhuru follow with time. (Points.in)
		fort.44 (optional)	Contour of the area where the lithosphere will be removed
		fort.54 (optional) => Re-write the 'Grid.in' file to 'Grid_outUhuru.in'

Change the input model, Grid.in, defining a contour where you want to change the crust or lithosphere thicknes.
	1. copy the input grid to change (Grid.in) to 'fort.11'
	2. Create the contour inside where the crust or lithosphere will be changed and save in the file 'fort.54'. 
		ATENTION: File with just corners of the poligon, not need to know the number of corners
	3. Check you don't have any 'fort.4' in the directory
	4. Run UHURU
	5. The programm will ask the meters you want to add (+) or sustract(-) to the actual crust and/or lithosphere thicknesses.
	6. The output file 'Grid_outUhuru.in' has the same format than 'Grid.in'


________________________________  parametres.in file ________________________________________________________________________________
				      36 or 39 lines
------------------------------------------------------------------------------------------------------------------------------------
TFPARA 			- Title [max 80 characters]
RHOH2O rosed roc RHOAST - Densities: water, sediments, crust, asthenosphere, [kg/m**3]
roalfa			- Volumetric thermal expansion of lithospheric mantle, [1/K]
TSURF TBOTT ZASTH	- Surface (z=0) and bottom (z=ZASTH) temperature [K] and asthenosphere depth [m]
				ZASTH: Could be replaced on the program, in order to get Dz in meters with only 2 decimals.
TISOTER			- Isotherm that define the base of the lithosphere (ZISOTER), [K]
THDIFF			- Thermal diffusivity, [m**2/s]
CONDUC(0,1,2,3)		- Conductivities: sediments, crust, lithospheric mantel and asthenosphere, [W/m*K]
HSURF HEXP		- Surface [W/m**3] and exponent [m] of the crustal heat production, H=HSURF*exp(-z/HEXP)
PHEAT_m			- Lithospheric mantle heat production, [W/m**3]
IRheology_type		- Type of rheological parameters (forward explained). if IRheology_type=99 ==> 3 more lines
   Qarray(1,2,3)	- Power law activation energy: upper crust, lower crust, lithospheric mantle, [J/mol]
   enarray(1,2,3)	- Power law exponent: upper crust, lower crust, lithospheric mantle
   Aarray(1,2,3)	- Pre-exponencial constant: upper crust, lower crust, lithospheric mantle, [MPa-n s-1]
------------------------------------------------------------------------------------------------------------------------------------
n m NELZ		- Number of columns (x), rows (y) and vertical (z)
strainrate		- Reference strain rate (s-1)
nitermax		- Maximum number of iterations, for convergence between viscosity and velocity, within velocity solution in each time step
				If 0 =>	Lineal problem: viscosity doesn't depend on strain rate.
				If 1 =>	No lineal problem (but linear equation on each time step): 
					Viscosity depends on previous velocity. No iteration to find a consistent velocity and viscosity.
				If >1=>	No lineal problem: Iterating method recalculating the viscosity with the new velocity.
tallmax			- Criterium of convergence. If SUM(Dvel/vel) < tallmax => convergence. (p.e., =0.01).
alfa			- Relation between initial and new velocity when nitermax>1 => v[t+Dt]=alfa*v[t+Dt]+(1-alfa)*v[t]
					0 explicit;  1 implicit;  0.5 recommended				    
tinicial		- Beginning of calculation, [years]
Dtany			- Size of time step, [years]
npasmax			- Number of time steps
NPASINT			- Number of steps that will keep in files.
------------------------------------------------------------------------------------------------------------------------------------
ITSR visco_cnst		- Viscosity and strain rate dependence (ITSR = 0-3) :	
			     If -1 => Constant lithosphere strength to every where, no change in time => read constant strength [N/m]
			     If 0 => Constant viscosity to every where, no change in time => read visco_cnst [Pa s]
			     If 1 => Calculate viscosity using a constant strain rate (strainrate from the parameters file) 
			     If 2 => Calculate viscosity using the efective strain rate at each point on the stress envelope 
					and on the constitutive equation
			     If 3 => Calculate viscosity using the efective strain rate at each point on the stress envelope 
					and use the constant strain rate (strainrate) on the constitutive equation.
Iremoval Time_removal Zremoval/TISO_rem Zctall - At time=Time_removal, Lithosphere Root Removal (LR) 
						   where L/s deeper than Zctall at Zremoval/TISO_rem.  
			     If 0 => No convective removal/home/ivone/Projectes-Varios/Asia_Uhuru_tisc/ReMesh/185x125/M6
			     If 1 => LR at ISOTHERM L=L(T=TISO_rem) where LITHOSPHERE > Zctall 
			     If 2 => LR at LITHOSPHERE = Zremoval where LITHOSPHERE > Zctall
			     If 22 => = 2 but reading the convective point from a file (removal inside a band, uhuru.f) => you must
					assign a low value to Zctall (< GLit)
			     If 3 => = 2 but Continuous removal, start Time_removal and end t2 (SUBROUTINE Convective_Removal_Type)
			     If 4 => LR at LITHOSPHERE = Zremoval where CRUST > Zctall
			     If 5 => = 4 but LITHOSPHERE = Zremoval lineal transition 
			     If 6 => LR  ISOTHERM L=L(T=TISO_rem) where CRUST > Zctall
			 Then:	If Iremoval=2-5 => Read : Time_removal [My],  Zremoval [km], Zctall [km]
				If Iremoval=1,6 => Read : Time_removal [My],  TISO_rem [K],  Zctall [km]
     Slow Lithosphere root removal =>	at SUBROUTINE Convective_Removal_Type DTime_removal/=0 [m.y.] =0(instantaneous)
					Works with LITHOSPHERE = Zremoval, not with an isotherm. (Iremoval/=1,6)

iremesh			- If 1 => Remeshing, when the indenter go norther than Y0+Dy	
			  Of 0 => No remeshing			
dvis_allow Dif_K	- Control of the lateral variations of some variables (viscosity, lithosphere thickness):
				maximum gradient permited [dvis_allow=[d(vis)/dx/]vis] and diffusive filter (Dif_K)
-------------- ELASTIC THICKNESS AND SURFACE PROCESSES --------------------------------------------------------------------------------------
Te_flexure		- Elastic thickness [m],  if Te_flexure=0 => local isostasy
hydro_model		- hydro_model=1 for rain proportional to elevation; hydro_model=2 for orographic precipitation (wind-dependent); 
				hydro_model=3 for conservative orographic precipitation (wind-dependent);
Kerosdif 		- Diffusive transport erosion coefficient, [m2/a] (0 <= Kerosdif <=5000). if 0 => no diffusive erosion
rain			- (if hydro_model=1) Background precipitation (or runoff, water going to the drainage system), [l/m2/yr=mm/yr]  (rain~500)
			- (if hydro_model=2) Precipitation at T=0 C at plains or in absence of wind.
			- (if hydro_model=3) Precipitation in a saturated column (related to turbulence) P = rain * Water_col / Water_sat
Krain			- (if hydro_model=1) Proportionality of runoff with altitude, [l/m2/a/km]=[mm/a/km] (Krain~500)
			- (if hydro_model>=2: mean wind velocity [m/s] (Krain>0 and usually Krain<10).
windazimut relative_humidity	- windazimut (if hydro_model>=2): azimut of wind flow direction in degrees counted clockwise from north (positive y).
				- relative_humidity (if hydro_model==3): humidity factor of incomming wind (0<relative_humidity<1)
evaporation			- Evaporation rate at lakes [l/m2/yr]=[mm/yr]. evaporation>runoff can significantly slow down the hydrological calculations. 
					If hydro_model=1,2: laterally constant evaporation.
					If hydro_model=3: evaporation caused by dry air and 0 wind speed. E is then calculated as:
						E=evaporation * (1+beta*wind) * (Wmax-Wcol)/Wmax
erodability  erodability_sed	- erodability (20e3 - 200e3 m)  AND  erodability_sed (10e3 - 100e3 m)
				- (If erosed_model=6): erodability ~2e-6 and erodability_sed ~8e-6, in meters_rock/s / (kg/m3 * m/s2 * m)^a == meters_rock/s / (Pa)^a	
K_river_cap  l_fluv_sedim	- K_river_cap (60-80)  AND  l_fluv_sedim (10e3 - 100e3 m)
CXrain  CYrain		- (if hydro_model=1) Distances in meters over which rain duplicates in X and Y directions, positive mean increasing north and east. 
			- (if hydro_model>=2) CXrain: Distances in meters over which orographic precipitations smooths out. (5e3 - 20e3 m o 0).  CYrain: Not used.
________________________________END of parametres.in file_______________________________________________________________________________

	No rivers 		=>  rain=0 and Krain=0
	No surface porcesses	=>  Kerosdif=0, rain=0 and Krain=0

If we call a file with the initial lithosphere structure (Grid.in), this parameters are not taken into account:
	INPUT_DATA2, n, m, AX, BY

________________________________  Grid.in file  ________________________________________________________________________________
# TITLE
n m Dx Dy INPUT_DATA
x y DATA	(x=0,...,n*Dx) (y=0,...,m*Dy)
	   DATA two variables:	If INPUT_DATA=1 => elevation [m] and surface heat flow [W/m2]	=> Only for neotectonic studies
				If INPUT_DATA=2 => elevation [m] and crustal thickness [m]	=> Only for neotectonic studies
				If INPUT_DATA=3 => crustal [m] and lithospheric [m] thickness	=> Also temporal studies


________________________________  BC.in file  ________________________________________________________________________________
  ITBC=1    velocity fixed [Vx(m/s), Vy(m/s)]			->  vel_x=t1  and  vel_y=t2
  ITBC=12   stress xx and xy fixed (Est, West free)		-> tau_xx=t1  and  tau_xy=t2
  ITBC=13   stress xx and yy fixed   		   		-> tau_xx=t1  and  tau_yy=t2
  ITBC=23   stress xy and yy fixed (North, South free)	-> tau_xy=t1  and  tau_yy=t2
  ITBC=4    free slip (vel normal=0, tau_xy=0)		-> don't use t1 and t2
  ITBC=5    velocity fixed [modul(mm/yr), azimut(degree)]	-> modul=t1  and  azimut=t2


________________________________  Time_BCfiles.in file (optional) ________________________________________________________________



________________________________  Bodies.in file (optional) ________________________________________________________________________________
kpunts1 qFlit1		- number of vertex of the polygon and factor to multiply the strenght of the inside points (qFlit1=strenght_inside/strength)
x_pol1(1) y_pol1(1)	- Vertex in counterclockwise order
... 
x_pol1(kpunts1) y_pol1(kpunts1)
kpunts2 qFlit2		- number of vertex of the second polygon and strenght factor (qFlit2=strenght_inside/strength)
x_pol2(1) y_pol2(1)	- Vertex in counterclockwise order
...
x_pol2(kpunts2) y_pol2(kpunts2)


Maximum of 5 bodies.


________________________________  Points.in file (optional) ________________________________________________________________________________






*************  IRheology_type:  Different rheological parameters  ******************************************************
   IRheology_type=1 ==>  Lynch&Morgan, 1987
           Qarray=(138.D3, 251.D3, 523.D3), enarray=(3.0D0, 3.0D0, 3.0D0), Aarray=(2.5D-8, 3.2D-3, 1.D3)
   IRheology_type=2 ==>  Braun&Beaumont, 1989
           Qarray=(151.D3, 239.D3, 498.D3), enarray=(1.8D0, 3.2D0, 4.5D0), Aarray=(1.D-2, 3.D-2, 1.9D5)
   IRheology_type=3 ==>  Fadaie&Ranalli, 1990
           Qarray=(219.D3, 268.D3, 535.D3), enarray=(2.4D0, 3.3D0, 3.6D0), Aarray=(1.3D-3, 3.2D-3, 3.2D4)
   IRheology_type=4 ==>  Buck, 1991
           Qarray=(149.D3, 238.D3, 500.D3), enarray=(2.9D0, 3.2D0, 3.D0), Aarray=(1.7D-7, 8.9D-4, 1.D3)
   IRheology_type=5 ==>  Liu&Furlong, 1993
           Qarray=(123.D3, 260.D3, 420.D3), enarray=(3.D0, 3.4D0, 3.D0), Aarray=(1.6D-9, 2.D-4, 1.9D3)
   IRheology_type=6 ==>  Lowe&Ranalli, 1993
           Qarray=(144.D3, 238.D3, 535.D3), enarray=(3.2D0, 3.2D0, 3.5D0), Aarray=(1.3D-9, 3.3D-4, 1.4D5)
           sigmaD=15.D9
   IRheology_type=7 ==>  Mareschal, 1994
           Qarray=(141.D3, 445.D3, 527.D3), enarray=(1.9D0, 4.2D0, 3.D0), Aarray=(2.D-4, 1.4D-4, 4.3D2)
           sigmaD=12.D9
   IRheology_type=8 ==>  Bassi, 1995 (dry)
           Qarray=(185.D3, 235.D3, 535.D3), enarray=(2.8D0, 3.9D0, 3.6D0), Aarray=(3.4D-6, 2.3D-6, 2.9D4)
   IRheology_type=9 ==>  Bassi, 1991 (wet
           Qarray=(150.D3, 238.D3, 445.D3), enarray=(1.8D0, 3.2D0, 3.4D0), Aarray=(2.9D-3, 3.3D-4, 1.4D4)	   
   IRheology_type=20 ==>  Choosed parameters, L&M + n=10
           Qarray=(138.D3, 251.D3, 523.D3), enarray=(10.0D0, 10.0D0, 10.0D0), Aarray=(2.5D-8, 3.2D-3, 1.D3)
   IRheology_type=21 ==>  Choosed parameters, L&M + n=1
           Qarray=(138.D3, 251.D3, 523.D3), enarray=(1.0D0, 1.0D0, 1.0D0), Aarray=(2.5D-8, 3.2D-3, 1.D3)
   IRheology_type=99 ==>  Choosed parameters
           Parameters from the file (parametres.in)	   
   IRheology_type>=100 ==>  parameters used on the Mediterrani from Shells (Peter+Kirby)
	IRheology_type=100 ==> JGR Mediterrani (Peter+Kirby)	Qarray=(100.D3, 100.D3, 456.D3), enarray=(3.D0, 3.D0, 3.D0)
	IRheology_type=110 ==> Peter+Kirby, soft, Q minimes	Qarray=(100.D3, 100.D3, 420.D3)
	IRheology_type=120 ==> Peter+Kirby, hard, Q maximes	Qarray=(445.D3, 445.D3, 535.D3)
******************************************************************************************************************


**********************  Parameters fixed on the source UHURU.f ****************************************************
	irel_ZL=0 ==> nivell de relaxacio fixe, ZLITOS_rel => GLit(kxy)=ZLITOS_rel+elevation(kxy)
	irel_ZL=1 ==> GLit(kxy)= depth of isotherm (TLITOS)
	irel_ZL=2 ==> GLit(kxy)=GLit_ter(kxy)+DGL_ter

	Ldepth_var=0 -> constant Lithospheric depth (Depth_lit) and no temporal variations
	iTemp_var=0  -> No Temperature temporal variations. TEMPE2=TEMPE1

	iflexure=1 -> Flexure calculation

 	ifixv=1 => keep in a file the points inside body 1 (subrutine VISCOSITAT) => velocity fixed in these points
			in subrutine VELOCITY_FIELD

Boundary conditions of thickness variations (crust, sediments and erosion):
	IBC_thicken=0 -> No temporal thickening variations on the boundaries, thickness=thickness_old
	IBC_thicken!=0 -> No lateral variations of the thickening on the boundaries, d(thickness)/dx=d(thickness)/dy=0

***********************************************************************************************************************

********************   SOME COMMENTS **************************************************
Al hacer el ReMeshing se propagaran las condiciones de contorno del norte hacia el interior del modelo. 
En particular, cuidado con las del flexure. En el Modelo del Tibet, las rayas N-S en erosion son debidas a la condición de contorno de flexure 5555 (provoca baja topo en boundary) + remeshing (lo propaga al sur)

Mancha de mayor erosión en la parte distal del plateau es debida a la presencia de sedimentos de alta erodabilidad (mas erosión)



