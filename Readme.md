# Introduction #

Input files: 	fort.11			Initial geometry: x, y, crust and lithosphere thickness (Grid.in)


> parametres.in		Input parameters - open(8)

> fort.4			Velocity Boundary conditions: ix,iy,type of BC,condition1,condition2 (BC.in)

> fort.12 (optional)	Bodies: Zones with diferent strength (bodies.in)

> fort.13 (optional)	Time and files with different velocity boundary conditions (Time\_BCfiles.in)

> fort.10 (optional)	Points which uhuru follow with time. (Points.in)

> fort.44 (optional)	Contour of the area where the lithosphere will be removed

> fort.54 (optional) => Re-write the 'Grid.in' file to 'Grid\_outUhuru.in'



Change the input model, Grid.in, defining a contour where you want to change the crust or lithosphere thicknes.

  1. copy the input grid to change (Grid.in) to 'fort.11'

> 2. Create the contour inside where the crust or lithosphere will be changed and save in the file 'fort.54'.

> ATENTION: File with just corners of the poligon, not need to know the number of corners

> 3. Check you don't have any 'fort.4' in the directory

> 4. Run UHURU

> 5. The programm will ask the meters you want to add (+) or sustract(-) to the actual crust and/or lithosphere thicknesses.

> 6. The output file 'Grid\_outUhuru.in' has the same format than 'Grid.in'





_PARAMETERS file_
> Name:	 parametres.in   (obligatory)		36 or 39 lines


TFPARA 			- Title [80 characters](max.md)

RHOH2O rosed roc RHOAST - Densities: water, sediments, crust, asthenosphere, [kg/m3]

roalfa				- Volumetric thermal expansion of lithospheric mantle, [1/K]

TSURF TBOTT ZASTH	- Surface (z=0) and bottom (z=ZASTH) temperature [K](K.md) and asthenosphere depth [m](m.md)

> ZASTH: Could be replaced on the program, in order to get Dz in meters with only 2 decimals.

TISOTER			- Isotherm that define the base of the lithosphere (ZISOTER), [K](K.md)

THDIFF			- Thermal diffusivity, [m2/s]

CONDUC(0,1,2,3)		- Conductivities: sediments, crust, lithospheric mantel and asthenosphere, [W/m\*K]

HSURF HEXP		- Surface [W/m3] and exponent [m](m.md) of the crustal heat production, H=HSURF\*exp(-z/HEXP)

PHEAT\_m			- Lithospheric mantle heat production, [W/m3]

IRheology\_type		- Type of rheological parameters (forward explained). if IRheology\_type=99 ==> 3 more lines

> Qarray(1,2,3)		- Power law activation energy: upper crust, lower crust, lithospheric mantle, [J/mol]

> enarray(1,2,3)		- Power law exponent: upper crust, lower crust, lithospheric mantle

> Aarray(1,2,3)		- Pre-exponencial constant: upper crust, lower crust, lithospheric mantle, [MPa-n s-1]


n m NELZ			- Number of columns (x), rows (y) and vertical (z)

strainrate			- Reference strain rate (s-1)

nitermax			- Maximum number of iterations, for convergence between viscosity and velocity, within velocity solution in each 					time step

> If 0 =>	Lineal problem: viscosity doesn't depend on strain rate.

> If 1 =>	No lineal problem (but linear equation on each time step):

> Viscosity depends on previous velocity. No iteration to find a consistent velocity and viscosity.

> If >1=>	No lineal problem: Iterating method recalculating the viscosity with the new velocity.

tallmax			- Criterium of convergence. If SUM(Dvel/vel) < tallmax => convergence. (p.e., =0.01).

alfa				- Relation between initial and new velocity when nitermax>1 => v[t+Dt]=alfa\*v[t+Dt]+(1-alfa)**v[t](t.md)**

> 0 explicit;  1 implicit;  0.5 recommended

tinicial				- Beginning of calculation, [years](years.md)

Dtany				- Size of time step, [years](years.md)

npasmax			- Number of time steps

NPASINT			- Number of steps that will keep in files.




ITSR visco\_cnst		- Viscosity and strain rate dependence (ITSR = 0-3) :

> If -1 => Constant lithosphere strength to every where, no change in time => read constant strength [N/m]

> If 0 => Constant viscosity to every where, no change in time => read visco\_cnst [s](Pa.md)

> If 1 => Calculate viscosity using a constant strain rate (strainrate from the parameters file)

> If 2 => Calculate viscosity using the efective strain rate at each point on the stress envelope

> and on the constitutive equation

> If 3 => Calculate viscosity using the efective strain rate at each point on the stress envelope

> and use the constant strain rate (strainrate) on the constitutive equation.

Iremoval Time\_removal Zremoval/TISO\_rem Zctall - At time=Time\_removal, Lithosphere Root Removal (LR)

> where L/s deeper than Zctall at Zremoval/TISO\_rem.

> If 0 => No convective removal/home/ivone/Projectes-Varios/Asia\_Uhuru\_tisc/ReMesh/185x125/M6

> If 1 => LR at ISOTHERM L=L(T=TISO\_rem) where LITHOSPHERE > Zctall

> If 2 => LR at LITHOSPHERE = Zremoval where LITHOSPHERE > Zctall

> If 22 => = 2 but reading the convective point from a file (removal inside a band, uhuru.f) => you must

> assign a low value to Zctall (< GLit)

> If 3 => = 2 but Continuous removal, start Time\_removal and end t2 (SUBROUTINE Convective\_Removal\_Type)

> If 4 => LR at LITHOSPHERE = Zremoval where CRUST > Zctall

> If 5 => = 4 but LITHOSPHERE = Zremoval lineal transition

> If 6 => LR  ISOTHERM L=L(T=TISO\_rem) where CRUST > Zctall

> Then:	If Iremoval=2-5 => Read : Time\_removal [My](My.md),  Zremoval [km](km.md), Zctall [km](km.md)

> If Iremoval=1,6 => Read : Time\_removal [My](My.md),  TISO\_rem [K](K.md),  Zctall [km](km.md)

> Slow Lithosphere root removal =>	at SUBROUTINE Convective\_Removal\_Type DTime\_removal/=0 [m.y.] =0(instantaneous)

> Works with LITHOSPHERE = Zremoval, not with an isotherm. (Iremoval/=1,6)



iremesh			- If 1 => Remeshing, when the indenter go norther than Y0+Dy

> - If 0 => No remeshing

dvis\_allow Dif\_K		- Control of the lateral variations of some variables (viscosity, lithosphere thickness):

> maximum gradient permited [dvis\_allow=[d(vis)/dx/]vis] and diffusive filter (Dif\_K)

-------------- ELASTIC THICKNESS AND SURFACE PROCESSES --------------------------------------------------------------------------------------

Te\_flexure			- Elastic thickness [m](m.md),  if Te\_flexure=0 => local isostasy

hydro\_model			- hydro\_model=1 for rain proportional to elevation; hydro\_model=2 for orographic precipitation (wind-dependent);

> hydro\_model=3 for conservative orographic precipitation (wind-dependent);

Kerosdif 			- Diffusive transport erosion coefficient, [m2/a] (0 <= Kerosdif <=5000). if 0 => no diffusive erosion

rain				- (if hydro\_model=1) Background precipitation (or runoff, water going to the drainage system), [l/m2/yr=mm/yr]  (rain~500)

> - (if hydro\_model=2) Precipitation at T=0 C at plains or in absence of wind.

> - (if hydro\_model=3) Precipitation in a saturated column (related to turbulence) P = rain **Water\_col / Water\_sat**

Krain				- (if hydro\_model=1) Proportionality of runoff with altitude, [l/m2/a/km]=[mm/a/km] (Krain~500)

> - (if hydro\_model>=2: mean wind velocity [m/s] (Krain>0 and usually Krain<10).

windazimut relative\_humidity  - windazimut (if hydro\_model>=2): azimut of wind flow direction in degrees counted clockwise from north (+ y).

> - relative\_humidity (if hydro\_model==3): humidity factor of incomming wind (0<relative\_humidity<1)

evaporation			- Evaporation rate at lakes [l/m2/yr]=[mm/yr]. evaporation>runoff can significantly slow down the
> hydrological calculations.

> If hydro\_model=1,2: laterally constant evaporation.

> If hydro\_model=3: evaporation caused by dry air and 0 wind speed. E is then calculated as:

> E=evaporation **(1+beta\*wind)** (Wmax-Wcol)/Wmax

erodability  erodability\_sed	- erodability (20e3 - 200e3 m)  AND  erodability\_sed (10e3 - 100e3 m)

> - (If erosed\_model=6): erodability ~2e-6 and erodability\_sed ~8e-6, in meters\_rock/s /
> > (kg/m3 **m/s2** m)<sup>a == meters_rock/s / (Pa)</sup>a

K\_river\_cap  l\_fluv\_sedim	- K\_river\_cap (60-80)  AND  l\_fluv\_sedim (10e3 - 100e3 m)

CXrain  CYrain		- (if hydro\_model=1) Distances in meters over which rain duplicates in X and Y directions,

> positive mean increasing north and east.

> - (if hydro\_model>=2) CXrain: Distances in meters over which orographic precipitations smooths out.
> > (5e3 - 20e3 m o 0).  CYrain: Not used.

_END of parametres.in file_



> No rivers 		=>  rain=0 and Krain=0

> No surface porcesses	=>  Kerosdif=0, rain=0 and Krain=0



If we call a file with the initial lithosphere structure (Grid.in), this parameters are not taken into account:

> INPUT\_DATA2, n, m, AX, BY



_INITIAL LITHOSPHERE STRUCTURE_		Name:   Grid.in file


# TITLE

n m Dx Dy INPUT\_DATA

x y DATA	(x=0,...,n\*Dx) (y=0,...,m\*Dy)

> DATA two variables:	If INPUT\_DATA=1 => elevation [m](m.md) and surface heat flow [W/m2]	=>  Only for neotectonic studies

> If INPUT\_DATA=2 => elevation [m](m.md) and crustal thickness [m](m.md)		=>  Only for neotectonic studies

> If INPUT\_DATA=3 => crustal [m](m.md) and lithospheric [m](m.md) thickness		=>  Also temporal studies





_VELOCITY BOUNDARY CONDITIONS__> Name:  BC.in file_


> ITBC=1    velocity fixed [Vx(m/s), Vy(m/s)]			=>  vel\_x=t1  and  vel\_y=t2

> ITBC=12   stress xx and xy fixed (Est, West free)			=> tau\_xx=t1  and  tau\_xy=t2

> ITBC=13   stress xx and yy fixed   		   			=> tau\_xx=t1  and  tau\_yy=t2

> ITBC=23   stress xy and yy fixed (North, South free)		=> tau\_xy=t1  and  tau\_yy=t2

> ITBC=4    free slip (vel normal=0, tau\_xy=0)			=> don't use t1 and t2

> ITBC=5    velocity fixed [modul(mm/yr), azimut(degree)]	=> modul=t1  and  azimut=t2





_CHANGE WITH TIME OF THE VELOCITY BC_
> Name:  Time\_BCfiles.in file (optional)







    BODIES WITH DIFFERENT STRENGTH 

> Name: Bodies.in file (optional)

kpunts1 qFlit1		- number of vertex of the polygon and factor to multiply the strenght of the inside points (qFlit1=strenght\_inside/strength)

x\_pol1(1) y\_pol1(1)	- Vertex in counterclockwise order

...

x\_pol1(kpunts1) y\_pol1(kpunts1)

kpunts2 qFlit2		- number of vertex of the second polygon and strenght factor (qFlit2=strenght\_inside/strength)

x\_pol2(1) y\_pol2(1)	- Vertex in counterclockwise order

...

x\_pol2(kpunts2) y\_pol2(kpunts2)




Maximum of 5 bodies.


  POINTS TO BE FOLLOWED – MARKERS  
> Name:  Points.in file (optional)







RHEOLOGICAL PARAMETERS:


**IRheology\_type:  Different rheological parameters**

> IRheology\_type=1 ==>  Lynch&Morgan, 1987

> Qarray=(138.D3, 251.D3, 523.D3), enarray=(3.0D0, 3.0D0, 3.0D0), Aarray=(2.5D-8, 3.2D-3, 1.D3)

> IRheology\_type=2 ==>  Braun&Beaumont, 1989

> Qarray=(151.D3, 239.D3, 498.D3), enarray=(1.8D0, 3.2D0, 4.5D0), Aarray=(1.D-2, 3.D-2, 1.9D5)

> IRheology\_type=3 ==>  Fadaie&Ranalli, 1990

> Qarray=(219.D3, 268.D3, 535.D3), enarray=(2.4D0, 3.3D0, 3.6D0), Aarray=(1.3D-3, 3.2D-3, 3.2D4)

> IRheology\_type=4 ==>  Buck, 1991

> Qarray=(149.D3, 238.D3, 500.D3), enarray=(2.9D0, 3.2D0, 3.D0), Aarray=(1.7D-7, 8.9D-4, 1.D3)

> IRheology\_type=5 ==>  Liu&Furlong, 1993

> Qarray=(123.D3, 260.D3, 420.D3), enarray=(3.D0, 3.4D0, 3.D0), Aarray=(1.6D-9, 2.D-4, 1.9D3)

> IRheology\_type=6 ==>  Lowe&Ranalli, 1993

> Qarray=(144.D3, 238.D3, 535.D3), enarray=(3.2D0, 3.2D0, 3.5D0), Aarray=(1.3D-9, 3.3D-4, 1.4D5)
, sigmaD=15.D9

> IRheology\_type=7 ==>  Mareschal, 1994

> Qarray=(141.D3, 445.D3, 527.D3), enarray=(1.9D0, 4.2D0, 3.D0), Aarray=(2.D-4, 1.4D-4, 4.3D2)
, sigmaD=12.D9

> IRheology\_type=8 ==>  Bassi, 1995 (dry)

> Qarray=(185.D3, 235.D3, 535.D3), enarray=(2.8D0, 3.9D0, 3.6D0), Aarray=(3.4D-6, 2.3D-6, 2.9D4)

> IRheology\_type=9 ==>  Bassi, 1991 (wet
)
> > Qarray=(150.D3, 238.D3, 445.D3), enarray=(1.8D0, 3.2D0, 3.4D0), Aarray=(2.9D-3, 3.3D-4, 1.4D4)


> IRheology\_type=20 ==>  Choosed parameters, L&M + n=10

> Qarray=(138.D3, 251.D3, 523.D3), enarray=(10.0D0, 10.0D0, 10.0D0), Aarray=(2.5D-8, 3.2D-3, 1.D3)

> IRheology\_type=21 ==>  Choosed parameters, L&M + n=1

> Qarray=(138.D3, 251.D3, 523.D3), enarray=(1.0D0, 1.0D0, 1.0D0), Aarray=(2.5D-8, 3.2D-3, 1.D3)

> IRheology\_type=99 ==>  Choosed parameters

> Parameters from the file (parametres.in)

> IRheology\_type>=100 ==>  parameters used on the Mediterrani from Shells (Peter+Kirby)

> IRheology\_type=100 ==> JGR Mediterrani (Peter+Kirby)	Qarray=(100.D3, 100.D3, 456.D3), enarray=(3.D0, 3.D0, 3.D0)

> IRheology\_type=110 ==> Peter+Kirby, soft, Q minimes	Qarray=(100.D3, 100.D3, 420.D3)

> IRheology\_type=120 ==> Peter+Kirby, hard, Q maximes	Qarray=(445.D3, 445.D3, 535.D3)







  Parameters fixed on the source UHURU.f 

> irel\_ZL=0 ==> nivell de relaxacio fixe, ZLITOS\_rel => GLit(kxy)=ZLITOS\_rel+elevation(kxy)

> irel\_ZL=1 ==> GLit(kxy)= depth of isotherm (TLITOS)

> irel\_ZL=2 ==> GLit(kxy)=GLit\_ter(kxy)+DGL\_ter



> Ldepth\_var=0 -> constant Lithospheric depth (Depth\_lit) and no temporal variations

> iTemp\_var=0  -> No Temperature temporal variations. TEMPE2=TEMPE1



> iflexure=1 -> Flexure calculation



> ifixv=1 => keep in a file the points inside body 1 (subrutine VISCOSITAT) => velocity fixed in these points

> in subrutine VELOCITY\_FIELD



Boundary conditions of thickness variations (crust, sediments and erosion):

> IBC\_thicken=0 -> No temporal thickening variations on the boundaries, thickness=thickness\_old

> IBC\_thicken!=0 -> No lateral variations of the thickening on the boundaries, d(thickness)/dx=d(thickness)/dy=0







  Parameters fixed on the source call\_tisc.c 
see file tisc/doc/tisc.info.txt
Flexure Boundary Conditions


   SOME COMMENTS 

Al hacer el ReMeshing se propagaran las condiciones de contorno del norte hacia el interior del modelo.

En particular, cuidado con las del flexure. En el Modelo del Tibet, las rayas N-S en erosion son debidas a la condición de contorno de flexure 5555 (provoca baja topo en boundary) + remeshing (lo propaga al sur)
. Con la BC 6666 desaparece, mejora.


Mancha de mayor erosión en la parte distal del plateau es debida a la presencia de sedimentos de alta erodabilidad (mas erosión)