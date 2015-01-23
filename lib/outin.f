CCCCCCCCC*CCCCCCCCC*CCCCCCCCC*CCCCCCCCC*CCCCCCCCC*CCCCCCCCC*CCCCCCCCC*
C                SUBRUTINA  outin.f 
C     Determina si un punto se encuentra dentro de un poligono
C    (x_pol(i),y_pol(i)) : vertexs del poligon.  
C        i=1,np_pol.     np_pol=numero de vertexs del poligon.
C         Sentit antihorari i sense repetir l'ultim punt
C    (x,y) : punto que ens interessa saber si esta a dins o no.
C    iadentro : output. Si iadentro=1 -> el punt esta dins el poligon.

       SUBROUTINE outin (x,y,np_pol,x_pol,y_pol,iadentro)
       implicit integer (i-n)
       implicit double precision (a-h,o-z)
       PARAMETER (zero=1D-20)

       dimension x_pol(np_pol),y_pol(np_pol)

        circulacion=0
        pol_long=0
        parcial_circ=0
       do 10 i=1,np_pol
            x0=x_pol(i)
            y0=y_pol(i)

C	Si x0=xf y y0=yf:
            if (i.lt.np_pol) then
                 xf=x_pol(i+1)
                 yf=y_pol(i+1)
            else
                 xf=x_pol(1)
                 yf=y_pol(1)
            endif
            
C
            D2=(x0-x)**2+(y0-y)**2
	    E2=(xf-x)**2+(yf-y)**2
            rL2=(xf-x0)**2+(yf-y0)**2
            rL=sqrt(rL2)
	    D=sqrt(D2)
	    E=sqrt(E2)
            pol_long=pol_long+rL
            rK=((x0-x)*(yf-y0)-(y0-y)*(xf-x0))/rL
C		if (rK.eq.0) print*, 'rK=0 en outin.f'
             b=2/rL*((y0-y)*(yf-y0)+(x0-x)*(xf-x0))
C            if ((y0-y)*(yf-y0)+(x0-x)*(xf-x0).eq.0.) print*,'b=0'
            disc=b*b-4*D2
            if (disc.ge.0) then
                 sq=sqrt(disc)
            else
                 sq=sqrt(-disc)
            endif
C
C    Si (x,y) se encuentra en el segmento (xf,yf)(x0,y0) ==>
C    el pto (x,y) se considera dentro del poligono.
	    if ((D.lt.zero).OR.(E.lt.zero).OR.(ABS(D+E-rL).lt.rL/1e6))
     +		     then
		  iadentro=1
		  return
            endif
            if (disc.gt.0) 
     +               parcial_circ=rK/sq*(log((2*rL+b-sq)/(2*rL+b+sq))-
     +           log((b-sq)/(b+sq)))
            if (disc.lt.0)
     +             parcial_circ=rK/sq*(atan((2*rL+b)/sq)-atan(b/sq))
            if (disc==0) parcial_circ=1/D-1/E
	    !if (ABS(disc).lt.zero) parcial_circ=1/D-1/E

            if (rL.lt.zero) parcial_circ=0

            circulacion=circulacion+parcial_circ
C
 10    continue
C
       if (abs(circulacion).gt.1) then
         iadentro=1
       else
         iadentro=0
       endif
       
 100  continue  

C       if (iadentro.eq.1) then
C		print*,'El punto se encuentra DENTRO del poligono'
C	else
C		print*,'El punto se encuentra FUERA del poligono'
C	endif

       RETURN
       END

