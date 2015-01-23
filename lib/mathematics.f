C ****************************************************************
!!   SUBROUTINES : sistbanda, equation_system,  diffusive_filter, smooth_horizontal 
!!   FUNCTIONS: VAL_MIN, VAL_MAX
****************************************************************

      SUBROUTINE sistbanda (a,b,nfiles,nbanda,Li,x)
 
C   SISTEMA DE n INCOGNITES I n EQUACIONS LINEALMENT INDEPENDENTS
C                MATRIU a EN BANDA
C                  SOLUCIO UNICA

C                   a . x = b
C  a: Matriu del sistema d'equacions, a(nfiles,nbanda)
C  nfiles: numero d'equacions = numero d'incognites x(nfiles) = numero de files d'a
C  nbanda: numwero de columnes d'a   ( nbanda = Li+1+Ls )
C  b: terme independent, nfiles files
C  Li: banda inferior, termes per sota/esquerra de la diagonal
C  Ls: banda superior, termes per sobre/dreta de la diagonal
C      Els termes de la diagonal estan en la columna id

      implicit double precision (a-h,o-z)
      dimension a(nfiles,nbanda),b(nfiles),x(nfiles)
      x=0.D0
      id=Li+1
      Ls=nbanda-(Li+1)
      do 5 j=1,nfiles-1
          if(a(j,id).eq.0.D0) then
             DO l=j+1,j+Li
                 jj=id+j-l
                 if(a(l,jj).ne.0.D0) then
                     DO k1=1,nbanda
                        k2=k1+j-l
                        ap=a(l,k2)
                        a(l,k2)=a(j,k1)
                        a(j,k1)=ap
                     END DO
                     bp=b(l)
                     b(l)=b(j)
                     b(j)=bp
                     goto 54
                 endif
             END DO
             print*,'no he trobat com substituir el terme',j
             print*,'  SISTEMA INDETERMINAT '
             stop 
           endif
54        continue

          do 7 i=j+1,j+Li
                IF(i.GT.nfiles) GOTO 5
                kj=id+j-i
                if(a(i,kj).eq.0.D0) goto 7
                fac=a(i,kj)/a(j,id)
                do 9 k1=kj,nbanda
                     k2=k1+i-j
                     if(k2.gt.nbanda) then
                          aaa=0.d0
                       else
                          aaa=a(j,k2)
                     endif
                     a(i,k1)=aaa*fac-a(i,k1)
9               continue
                b(i)=b(j)*fac-b(i)
7         continue
5     continue

c ******** aillar les incognites ********************

       if(a(nfiles,id).eq.0.D0) then
           print*,'The last number of the diagonal is nul'
           stop
       endif
       x(nfiles)=(b(nfiles))/(a(nfiles,id))
       DO kk=1,nfiles-1
          k=nfiles-kk
          if(a(k,id)==0.D0) then
              print*,'The diagonal number of the row ',k,' is nul'
              stop
          endif
          c=0.d0
          DO l=id+1,nbanda
              ll=l+k-id
              c=c+a(k,l)*x(ll)
          END DO
         x(k)=(b(k)-c)/a(k,id)
       END DO

      RETURN
      END SUBROUTINE sistbanda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      SUBROUTINE equation_system (N, A, b, x)

!!  Linear equation solution by Gauss-Jordan elimination. 
!!  A * x =b
!!  A: Matrix of N rows and N columns. N: number of unknowns.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N), b(N), x(N)
      
      LOOP1: DO nf=1,N-1		
         IF(A(nf,nf)==0.D0) THEN
	    DO l=nf+1,N
	       IF(A(l,nf)/=0.D0) THEN
                  DO k=1,N
		     ap=A(l,k)
		     A(l,k)=A(nf,k)
		     A(nf,k)=ap
                  END DO
		  bp=b(l)    
                  b(l)=b(nf)
		  b(nf)=bp
		  goto 54
               ENDIF
	    END DO    
            PRINT*,'no he trobat com substituir el terme',nf
            PRINT*,'  SISTEMA INDETERMINAT '
            STOP 
         ENDIF
54       continue

         !!! Triangulacio de la matriu
         loop2: DO i=nf+1,N
	    if(A(i,nf)==0.D0) CYCLE loop2		
            fac=A(i,nf)/A(nf,nf)
            DO k1=nf,N
               A(i,k1)=A(nf,k1)*fac-A(i,k1)
            END DO
            b(i)=b(nf)*fac-b(i)
         END DO loop2

      END DO LOOP1		

!!! Aillar incognites
       IF(A(N,N)==0.D0) THEN
           PRINT*,'The last equation is nul'
           STOP
       ENDIF
       x(N)=(b(N))/(A(N,N))

       DO kk=1,N-1
          k=N-kk
          IF(A(k,k)==0.D0) THEN
              PRINT*,'The diagonal number of the row ',k,' is nul'
              STOP
          ENDIF
          c=0.d0
          DO l=k+1,N
              c=c+A(k,l)*x(l)
          END DO
	  x(k)=(b(k)-c)/A(k,k)
       END DO

      END SUBROUTINE equation_system

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE diffusive_filter (m, n, nn, Dx, Dy, Dif_K, zeta)

!   Filtre difusiu d'una magnitut (zeta) depenent de x e y => zeta(x,y) 
!	Dzeta = Dif_K_xy * (d2zeta/dx2 + d2zeta/dy2)
! 	Dif_K: constant of diffusivity
!	zeta: Input and output
!	nn=(n+1)*(m+1)


      implicit double precision (a-h,o-z)
      dimension zeta(nn),Delta_z(nn) 
      
      Dx2=Dx*Dx
      Dy2=Dy*Dy
      Dif_K_xy=Dif_K*Dx*Dy

      DO iy=1,m-1
         DO ix=1,n-1
            kxy=ix+1+iy*(n+1)
            kxy_u=ix+1+(iy-1)*(n+1)	! Point (ix, iy-1)
            kxy_d=ix+1+(iy+1)*(n+1)	! Point (ix, iy+1)
	    D2zD2x=(zeta(kxy+1)-(2.D0*zeta(kxy))+zeta(kxy-1))/Dx2   ! derivada segona de zeta en x
	    D2zD2y=(zeta(kxy_d)-(2.D0*zeta(kxy))+zeta(kxy_u))/Dy2   ! derivada segona de zeta en y
	    Delta_z_x=Dif_K_xy*D2zD2x
	    DZmean_x=((zeta(kxy+1)+zeta(kxy-1))/2.D0)-zeta(kxy)	    	    
	    IF(D2zD2x<0) Delta_z_x=MAX(DZmean_x,Delta_z_x)
	    IF(D2zD2x>0) Delta_z_x=MIN(DZmean_x,Delta_z_x)	    
	    Delta_z_y=Dif_K_xy*D2zD2y	    
	    DZmean_y=((zeta(kxy_d)+zeta(kxy_u))/2.D0)-zeta(kxy)	    	    
	    IF(D2zD2y<0) Delta_z_y=MAX(DZmean_y,Delta_z_y)
	    IF(D2zD2y>0) Delta_z_y=MIN(DZmean_y,Delta_z_y)	    	    
	    Delta_z(kxy)=Delta_z_x+Delta_z_y		!Delta_z(kxy)=Dif_K_xy*(D2zD2x+D2zD2y)
         END DO
      END DO
! to the boundaries     
      DO iy=1,m-1
	  kxyw=1+iy*(n+1)		! (0,iy)
	  kxye=n+1+iy*(n+1)		! (n,iy)
	  Delta_z(kxyw)=Delta_z(kxyw+1)
	  Delta_z(kxye)=Delta_z(kxye-1)
      END DO 
      DO ix=0,n
	  kxys=ix+1			! (ix,0)
	  kxys1=ix+1+n+1		! (ix,1)
	  kxyn=ix+1+m*(n+1)		! (ix,m)
	  kxyn1=ix+1+(m-1)*(n+1)	! (ix,m-1)
	  Delta_z(kxys)=Delta_z(kxys1)
	  Delta_z(kxyn)=Delta_z(kxyn1)
      END DO
      
      DO iy=0,m
         DO ix=0,n
            kxy=ix+1+iy*(n+1)
            zeta(kxy)=zeta(kxy)+Delta_z(kxy)
         END DO
      END DO
      
      RETURN
      END SUBROUTINE diffusive_filter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      SUBROUTINE diffusive_filter_3D (m,n,NELZ, nn, Dx,Dy,Dz, Dif_K, T)

!   Filtre difusiu d'una magnitut (T) depenent de x, y, z => T(x,y,z) 
!	DT = Dif_K_xy * (d2T/dx2 + d2T/dy2 + d2T/dz2)
! 	Dif_K: constant of diffusivity
!	T: Input and output
!	nn=(n+1)*(m+1)

!      implicit double precision (a-h,o-z)
!      dimension T(nn,0:NELZ),Delta_T(nn,0:NELZ) 
      
!      Delta_T=0.D0
!      Dx2=Dx*Dx
!      Dy2=Dy*Dy
!      Dz2=Dz*Dz
!      Dif_K_xy=Dif_K*Dx*Dy

!      DO iz=1,NELZ-1
!	 DO iy=1,m-1
!	    DO ix=1,n-1
!		kxy=ix+1+iy*(n+1)
!		kxy_u=ix+1+(iy-1)*(n+1)	! Point (ix, iy-1)
!		kxy_d=ix+1+(iy+1)*(n+1)	! Point (ix, iy+1)
!		D2TD2x=(T(kxy+1,iz)-(2.D0*T(kxy,iz))+T(kxy-1,iz))/Dx2   ! derivada segona de T en x
!		D2TD2y=(T(kxy_d,iz)-(2.D0*T(kxy,iz))+T(kxy_u,iz))/Dy2   ! derivada segona de T en y
!		D2TD2z=(T(kxy,iz+1)-(2.D0*T(kxy,iz))+T(kxy,iz-1))/Dz2   ! derivada segona de T en z
!		Delta_T(kxy,iz)=Dif_K_xy*(D2TD2x+D2TD2y+D2TD2z)    
!	    END DO
!	 END DO
!! to the boundaries     
!	 DO iy=1,m-1
!	    kxyw=1+iy*(n+1)		! (0,iy)
!	    kxye=n+1+iy*(n+1)		! (n,iy)
!	    Delta_T(kxyw,iz)=Delta_T(kxyw+1,iz)
!	    Delta_T(kxye,iz)=Delta_T(kxye-1,iz)
!	 END DO 
!	 DO ix=0,n
!	    kxys=ix+1			! (ix,0)
!	    kxys1=ix+1+n+1		! (ix,1)
!	    kxyn=ix+1+m*(n+1)		! (ix,m)
!	    kxyn1=ix+1+(m-1)*(n+1)	! (ix,m-1)
!	    Delta_T(kxys,iz)=Delta_T(kxys1,iz)
!	    Delta_T(kxyn,iz)=Delta_T(kxyn1,iz)
!	 END DO
!      END DO
      
!      DO iz=1,NELZ-1
!	 DO kxy=1,nn
!            T(kxy,iz)=T(kxy,iz)+Delta_T(kxy,iz)
!         END DO
!      END DO
      
!      RETURN
!      END SUBROUTINE diffusive_filter_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   suavitzo els valors segons els 8 valors dels voltants 
!! zeta: Input i output

      SUBROUTINE smooth_horizontal (m, n, nn, Dx, Dy, zeta)

      implicit double precision (a-h,o-z)
      dimension zeta_ini(nn),zeta(nn) 

      zeta_ini=zeta
      diagonal=dsqrt(Dx*Dx+Dy*Dy)
      denom=(4.D0/diagonal)+(2.D0/Dy)+(2.D0/Dx)
      DO iy=1,m-1
         DO ix=1,n-1
            kxy=ix+1+iy*(n+1)
            kxy_u=ix+1+(iy-1)*(n+1)
            kxy_d=ix+1+(iy+1)*(n+1)
            zD=(zeta_ini(kxy_d-1)+zeta_ini(kxy_d+1)+zeta_ini(kxy_u-1)+
     +		         zeta_ini(kxy_u+1))/diagonal
            zy=(zeta_ini(kxy_d)+zeta_ini(kxy_u))/Dy
            zx=(zeta_ini(kxy-1)+zeta_ini(kxy+1))/Dx
	    z_contorn=(zD+zy+zx)/denom
	    !zeta(kxy)=(z_contorn+zeta_ini(kxy))/2.D0	
	    zeta(kxy)=(z_contorn+29.D0*zeta_ini(kxy))/30.D0	
         END DO
      END DO
! to the boundaries     
      DO iy=1,m-1
	  kxyw=1+iy*(n+1)
	  kxye=n+1+iy*(n+1)
	  zeta(kxyw)=zeta(kxyw+1)
	  zeta(kxye)=zeta(kxye-1)
      END DO 
      DO ix=0,n
	  kxys=ix+1
	  kxyn=ix+1+m*(n+1)
	  zeta(kxys)=zeta(kxys+n+1)
	  zeta(kxyn)=zeta(kxyn-n-1)
      END DO
      
      RETURN
      END SUBROUTINE smooth_horizontal

!!------------------------------------------------------!!!
      Double Precision FUNCTION VAL_MAX (zeta, nn)

	!! It finds the maximum value from the array zeta(nn)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (z_sup=-1.D200)
      DIMENSION zeta(nn)
       zm=z_sup
      DO kxy=1,nn
	 zm=MAX(zeta(kxy),zm)
      END DO
      VAL_MAX=zm

      END FUNCTION 
!!------------------------------------------------------!!!
      Double Precision FUNCTION VAL_MIN (zeta, nn)

	!! It finds the minimum value from the array zeta(nn)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (z_inf=1.D200)
      DIMENSION zeta(nn)
      zm=z_inf
      DO kxy=1,nn
	 zm=MIN(zeta(kxy),zm)
      END DO
      VAL_MIN=zm

      END FUNCTION 

