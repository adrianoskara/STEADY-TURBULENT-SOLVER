! Collocated Variable Approach for Chaunel Problem
! first-order upwinding for all trausport equations
! includes k-epsilon model aud wall functions

! Outputs results in .csv format to be viewed with paraview

! mu = dynamic viscosity
! rho = density
! nx = number finite volume cells in x-direction
! ny = number finite volume cells in y-direction
! omagam = underrelaxation factor for momentum equations
! omegap = underrelaxation factor for pressure correction

! ti = turbulence intensity (0<ti<1)


      MODULE VARIABLES
      IMPLICIT NONE
      SAVE
      INTEGER :: nx, ny, outeriter,maxiter,iturb
      REAL :: rho=1.0,omegam=0.5,omegap=0.3, &
              xmax,ymax,xmin,ymin,dx,dy,mudyn=1./40000., &
              volume,sigk=1.0,omegak=0.3, omegarps=0.3, &
              cep1=1.44,cep2=1.92,cmu=0.09,sigeps=1.3, &
              yplus,kappa=0.41,E,B=5.5,tauwall,uplus,ystar, &
              ti=0.1,uinlet=1.0,icyclic=0
      REAL, ALLOCATABLE, DIMENSION(:,:) :: u,v,p,uold,vold,pp,apu,apv,&
      app,ar,al,au,ad,source,dpdx,dpdy,dpdxc,dpdyc,dpdxave,dpdyave,&
      ustar,vstar,mut,sx,sy,production,k,eps,&
      kold,epsold,dissipation,apeps,apk,s
      REAL, ALLOCATABLE, DIMENSION(:) :: xp, xf, yp, yf, fx, fy
      END MODULE VARIABLES
!*************************************************************
!     Main Program
      PROGRAM CFD
      USE VARIABLES
      IMPLICIT NONE
      INTEGER :: i,j,status

      CALL initialize  !initialize various parameters
      PRINT*,'finished initialize'
      DO outeriter=1,maxiter
        CALL umom
        CALL vmom
        CALL pressure
        CALL KINETIC_ENERGY
        CALL DISSIPATE
        CALL TURB_VISC
      END DO
      CALL OUTPUT

      END PROGRAM CFD
!*************************************************************
      SUBROUTINE UMOM
      USE VARIABLES
      IMPLICIT NONE
      INTEGER :: i,j,iter,iii
      REAL :: mur,mul,muu,mud,mdotr,mdotl,mdotu,mdotd,massflux,yplus1,mu

!      first-order upwinding
      DO j=1,ny
      DO i=1,nx
      mdotr=rho*(yf(j+1)-yf(j))*((1.-fx(i+1))*u(i+1,j) + fx(i+1)*u(i,j))
      mdotl=rho*(yf(j+1)-yf(j))*((1.-fx(i))*u(i,j) + fx(i)*u(i-1,j))
      mdotu=rho*(xf(i+1)-xf(i))*((1.-fy(j+1))*v(i,j+1) + fy(j+1)*v(i,j))
      mdotd=rho*(xf(i+1)-xf(i))*((1.-fy(j))*v(i,j) + fy(j)*v(i,j-1))
      mur=mudyn + (1.-fx(i+1))*mut(i+1,j) + fx(i+1)*mut(i,j)
      mul=mudyn + (1.-fx(i))*mut(i,j) + fx(i)*mut(i-1,j)
      muu=mudyn + (1.-fy(j+1))*mut(i,j+1) + fy(j+1)*mut(i,j)
      mud=mudyn + (1.-fy(j))*mut(i,j) + fy(j)*mut(i,j-1)
      ar(i,j)=MAX(-mdotr,0.0)+mur*(yf(j+1)-yf(j))/(xp(i+1)-xp(i))
      al(i,j)=MAX( mdotl,0.0)+mul*(yf(j+1)-yf(j))/(xp(i)-xp(i-1))
      au(i,j)=MAX(-mdotu,0.0)+muu*(xf(i+1)-xf(i))/(yp(j+1)-yp(j))
      ad(i,j)=MAX(mdotd,0.0)+mud*(xf(i+1)-xf(i))/(yp(j)-yp(j-1))
      END DO; END DO
!********************************************************************      
!     correct left/right boundaries (viscous terms)
      DO j=1,ny
      mur=mudyn + mut(nx+1,j)
      ar(nx,j)=MAX(-rho*(yf(j+1)-yf(j))*u(nx+1,j),0.) + &
                         mur*(yf(j+1)-yf(j))/(xp(nx+1)-xp(nx))
      mul=mudyn + mut(0,j)
      al(1,j)=MAX(rho*(yf(j+1)-yf(j))*u(0,j),0.) + &
                         mul*(yf(j+1)-yf(j))/(xp(1)-xp(0))
      END DO

!     wall functions
      ad(:,1)=0.;  au(:,ny)=0

!     compute apu coefficient
      apu(1:nx,1:ny)=(ar(1:nx,1:ny)+al(1:nx,1:ny)+ &
                      au(1:nx,1:ny)+ad(1:nx,1:ny))/omegam

!     compute wall shear stress and source term for wall functions
      s=0.0   !initialize wall shear stress source term to zero
     
      DO i=1,nx
        
            ystar=rho*cmu**0.25*SQRT(k(i,ny))*(yp(ny+1)-yp(ny))/mudyn
            tauwall=rho*kappa*cmu**0.25*k(i,ny)**0.5*u(i,ny)/LOG(E*ystar)
      
         s(i,ny)=rho*cmu**0.25*k(i,ny)**0.5*(xf(i+1)-xf(i))/(LOG(E*ystar)/kappa)

      
            ystar=rho*cmu**0.25*SQRT(k(i,1))*(yp(1)-yp(0))/mudyn
            tauwall=rho*kappa*cmu**0.25*k(i,1)**0.5*u(i,1)/LOG(E*ystar)

         s(i,1)= rho*cmu**0.25*k(i,1)**0.5*(xf(i+1)-xf(i))/(LOG(E*ystar)/kappa)
      END DO
!********************************************************************


      DO i=1,nx
      DO j=1,ny
      mur=(1.-fx(i+1))*mut(i+1,j) + fx(i+1)*mut(i,j) + mudyn
      mul=(1.-fx(i))*mut(i,j) + fx(i)*mut(i-1,j) + mudyn
      sx(i,j)=(mur*(u(i+1,j)-u(i,j))/(xp(i+1)-xp(i)) - &
               mul*(u(i,j)-u(i-1,j))/(xp(i)-xp(i-1))) / (xf(i+1)-xf(i)) + &
          ( (mudyn+mut(i,j+1))*(v(i+1,j+1)-v(i-1,j+1))/(xp(i+1)-xp(i-1)) - &
            (mudyn+mut(i,j-1))*(v(i+1,j-1)-v(i-1,j-1))/(xp(i+1)-xp(i-1)) ) / &
            (yp(j+1)-yp(j-1))
       volume=(xf(i+1)-xf(i))*(yf(j+1)-yf(j))
       sx(i,j)=sx(i,j)*volume
       END DO;  END DO
!********************************************************************
       OUTERLOOP : DO iter=1,10   
       massflux=0.0
       DO j=1,ny
       DO i=1,nx
       u(i,j)=(1.-omegam)*uold(i,j)+1./apu(i,j)*( &
                    ar(i,j)*u(i+1,j) + &
                    al(i,j)*u(i-1,j) + &
                    au(i,j)*u(i,j+1) + &
                    ad(i,j)*u(i,j-1) + &
                    (xf(i+1)-xf(i))*(yf(j+1)-yf(j))* &
                    (p(i-1,j)-p(i+1,j))/(xp(i+1)-xp(i-1)) + sx(i,j) - s(i,j))
       END DO; END DO 
       massflux=0.0
          DO j=1,ny
            massflux=massflux+rho*dy*u(nx,j)
          END DO
       u(nx+1,1:ny)=1.0/massflux*u(nx,1:ny) !add d/dx=0 outflow conditions
       END DO OUTERLOOP
       u(0,:)=u(nx+1,:)   !cyclic bc
       END SUBROUTINE UMOM
!*************************************************************
      SUBROUTINE VMOM
      USE VARIABLES
      IMPLICIT NONE
      INTEGER :: i,j,iter
      REAL :: mur,mul,muu,mud,mdotr,mdotl,mdotu,mdotd,massflux

      DO j=1,ny
      DO i=1,nx
      mdotr=rho*(yf(j+1)-yf(j))*((1.-fx(i+1))*uold(i+1,j) + fx(i+1)*uold(i,j))
      mdotl=rho*(yf(j+1)-yf(j))*((1.-fx(i))*uold(i,j) + fx(i)*uold(i-1,j))
      mdotu=rho*(xf(i+1)-xf(i))*((1.-fy(j+1))*v(i,j+1) + fy(j+1)*v(i,j))
      mdotd=rho*(xf(i+1)-xf(i))*((1.-fy(j))*v(i,j) + fy(j)*v(i,j-1))
      mur=(1.-fx(i+1))*mut(i+1,j) + fx(i+1)*mut(i,j) + mudyn
      mul=(1.-fx(i))*mut(i,j) + fx(i)*mut(i-1,j) + mudyn
      muu=(1.-fy(j+1))*mut(i,j+1) + fy(j+1)*mut(i,j) + mudyn
      mud=(1.-fy(j))*mut(i,j) + fy(j)*mut(i,j-1) + mudyn
      ar(i,j)=MAX(-mdotr,0.0)+mur*(yf(j+1)-yf(j))/(xp(i+1)-xp(i))
      al(i,j)=MAX( mdotl,0.0)+mul*(yf(j+1)-yf(j))/(xp(i)-xp(i-1))
      au(i,j)=MAX(-mdotu,0.0)+muu*(xf(i+1)-xf(i))/(yp(j+1)-yp(j))
      ad(i,j)=MAX(mdotd,0.0)+mud*(xf(i+1)-xf(i))/(yp(j)-yp(j-1))
      END DO; END DO
!********************************************************************
!     correct left/right inlet/outlet boundaries (viscous terms)
      DO j=1,ny
      mur=mut(nx+1,j) + mudyn
      ar(nx,j)=MAX(-rho*(yf(j+1)-yf(j))*uold(nx+1,j),0.) + &
                         mur*(yf(j+1)-yf(j))/(xp(nx+1)-xp(nx))
      mul=mut(0,j) + mudyn
      al(1,j)=MAX(rho*(yf(j+1)-yf(j))*uold(0,j),0.) + &
                         mul*(yf(j+1)-yf(j))/(xp(1)-xp(0))
      END DO
!     correct up/down wall boundaries (viscous terms)
      au(:,ny)=0.0   ! set dv/dy = 0 at wall from continuity
      ad(:,1)=0.0    ! set dv/dy = 0 at wall from continuity

!     compute apv coefficient
      apv(1:nx,1:ny)=(ar(1:nx,1:ny)+al(1:nx,1:ny)+ &
                      au(1:nx,1:ny)+ad(1:nx,1:ny))/omegam
!********************************************************************
!     compute source terms for variable viscosity
      DO i=1,nx
      DO j=1,ny
      muu=((1.-fy(j+1))*mut(i,j+1) + fy(j+1)*mut(i,j)) + mudyn
      mud=((1.-fy(j))*mut(i,j) + fy(j)*mut(i,j-1)) + mudyn
      sy(i,j)=(muu*(v(i,j+1)-v(i,j))/(yp(j+1)-yp(j)) - &
               mud*(v(i,j)-v(i,j-1))/(yp(j)-yp(j-1))) / (yf(j+1)-yf(j)) + &
      ((mudyn+mut(i+1,j))*(uold(i+1,j+1)-uold(i+1,j-1))/(yp(j+1)-yp(j-1)) - &
       (mudyn+mut(i-1,j))*(uold(i-1,j+1)-uold(i-1,j-1))/(yp(j+1)-yp(j-1))) / &
       (xp(i+1)-xp(i-1))    
       volume=(xf(i+1)-xf(i))*(yf(j+1)-yf(j))
       sy(i,j)=sy(i,j)*volume
       END DO;  END DO
!********************************************************************
       OUTERLOOP : DO iter=1,10  
       massflux=0.0
       DO j=1,ny
       DO i=1,nx
       v(i,j)=(1.-omegam)*vold(i,j)+1./apv(i,j)*( &
                    ar(i,j)*v(i+1,j) + &
                    al(i,j)*v(i-1,j) + &
                    au(i,j)*v(i,j+1) + &
                    ad(i,j)*v(i,j-1) + &
                    (yf(j+1)-yf(j))*(xf(i+1)-xf(i))* &
                    (p(i,j-1)-p(i,j+1))/(yp(j+1)-yp(j-1)) + sy(i,j))
      END DO; END DO 
      v(nx+1,1:ny)=v(nx,1:ny)
      END DO OUTERLOOP
      if(icyclic .eq. 1)v(0,:)=v(nx+1,:)   !cyclic bc
      END SUBROUTINE VMOM
!*************************************************************
      SUBROUTINE PRESSURE
      USE VARIABLES
      IMPLICIT NONE
      INTEGER :: i,j,iter
      REAL :: total,correction

!     define ar, al, au, using averages of ap's from surrounding cell centers
      DO i=1,nx
      DO j=1,ny     
      ar(i,j)=rho*(yf(j+1)-yf(j))**2*((1.-fx(i+1))/apu(i+1,j)+fx(i+1)/apu(i,j))
      al(i,j)=rho*(yf(j+1)-yf(j))**2*((1.-fx(i))/apu(i,j)+fx(i)/apu(i-1,j))
      au(i,j)=rho*(xf(i+1)-xf(i))**2*((1.-fy(j+1))/apv(i,j+1)+fy(j+1)/apv(i,j))
      ad(i,j)=rho*(xf(i+1)-xf(i))**2*((1.-fy(j))/apv(i,j)+fy(j)/apv(i,j-1))
      END DO; END DO

!     set the boundary conditions
      ar(nx,1:ny)=0.0
      al(1,1:ny)=0.0
      au(1:nx,ny)=0.0
      ad(1:nx,1)=0.0
!     compute app ad sum of ar, al, au, ad
      app=ar+al+au+ad
!     fix pressure at cell (1,1) to whatever value p wad initialized to
      app(1,1)=1.e30

!     compute the mass source terms using Rhie-Chow Interpolation
!*******************************************************
!     compute dp/dx pressure gradients on cell faces
!     faces start at index i=1 (left boundary face) aud run to i=nx+1 (right boundary face)
      DO j=1,ny
      DO i=2,nx
      dpdx(i,j)=(p(i,j)-p(i-1,j))/(xp(i)-xp(i-1))
      END DO
      dpdx(1,j)=(p(1,j)-p(0,j))/(xp(1)-xp(0))
      dpdx(nx+1,j)=(p(nx+1,j)-p(nx,j))/(xp(nx+1)-xp(nx))
      END DO

!     the dp/dx terms at cell centers with 2dx spacing
      do i=1,nx
      do j=1,ny
      dpdxc(i,j)=(p(i+1,j)-p(i-1,j))/(xp(i+1)-xp(i-1))
      END DO;  END DO

!     take the average of the cell center x-derivatives to get interior face values
      DO j=1,ny
      DO i=2,nx
      dpdxave(i,j)=0.5*dpdxc(i,j)+0.5*dpdxc(i-1,j)
      END DO; END DO

!     Rhie-Chow ustar calculation
      DO j=1,ny
      DO i=2,nx
      correction=(xp(i)-xp(i-1))*(yf(j+1)-yf(j))* &
          ((1.-fx(i))/apu(i,j)+fx(i)/apu(i-1,j)) * &
          (dpdx(i,j) - dpdxave(i,j))
      ustar(i,j)=((1.-fx(i))*u(i,j)+fx(i)*u(i-1,j)) - correction
      END DO
	  
!     set ustar to u on left/right boundaries (velocities not corrected on boundaries)
      ustar(1,j)=u(0,j)
      ustar(nx+1,j)=u(nx+1,j)
      END DO
!******************************************************************
!     compute dp/dy pressure gradients on cell faces
      DO i=1,nx
      DO j=2,ny
      dpdy(i,j)=(p(i,j)-p(i,j-1))/(yp(j)-yp(j-1))
      END DO
      dpdy(i,1)=(p(i,1)-p(i,0))/(yp(1)-yp(0))
      dpdy(i,ny+1)=(p(i,ny+1)-p(i,ny))/(yp(ny+1)-yp(ny))
      END DO

!     compute the dp/dy terms at cell centers with 2dy spacing
      DO j=1,ny
      DO i=1,nx
      dpdyc(i,j)=(p(i,j+1)-p(i,j-1))/(yp(j+1)-yp(j-1))
      END DO; END DO

!     take the average of the cell center y-derivatives to get interior face values
      DO i=1,nx
      DO j=2,ny
     dpdyave(i,j)=0.5*dpdyc(i,j)+0.5*dpdyc(i,j-1)
      END DO; END DO

!     Rhie-Chow vstar calculation
      DO i=1,nx
      DO j=2,ny
      correction=(yp(j)-yp(j-1))*(xf(i+1)-xf(i)) * &
          ((1.-fy(j))/apv(i,j)+fy(j)/apv(i,j-1)) * &
          (dpdy(i,j) - dpdyave(i,j))
      vstar(i,j)=((1.-fy(j))*v(i,j)+fy(j)*v(i,j-1)) - correction
      END DO
!     set vstar to v on south/north boundaries (velocities not corrected on boundaries)
      vstar(i,1)=v(i,0)
      vstar(i,ny+1)=v(i,ny+1)
      END DO
!*******************************************************
!     compute mass imbalauce over each cell
      source=0
      DO i=1,nx
      DO j=1,ny
      source(i,j)=(yf(j+1)-yf(j))*rho*(ustar(i+1,j)-ustar(i,j)) + &
                  (xf(i+1)-xf(i))*rho*(vstar(i,j+1)-vstar(i,j))
      END DO; END DO

!     compute total mass imbalauce ad square root of the sum of squares over each cell
      total=sum(source**2)
      total=sqrt(total)
      print*,outeriter,total

      pp=0.0    !initialize pressure corrections to zero aud set on boundaries
!     SOR iterations to solve for pressure correction pp
      OUTERLOOP : DO iter=1,100
      DO j=1,ny
      DO i=1,nx
      pp(i,j)=1.0/app(i,j)*( &
                    ar(i,j)*pp(i+1,j)+ &
                    al(i,j)*pp(i-1,j)+ &
                    au(i,j)*pp(i,j+1)+ &
                    ad(i,j)*pp(i,j-1)- &
                    source(i,j))
      END DO; END DO
      END DO OUTERLOOP

!     Apply corrections to pressure aud velocities
      DO i=1,nx
      DO j=1,ny
      p(i,j)=p(i,j)+omegap*pp(i,j)
      u(i,j)=u(i,j)-(1./apu(i,j))*(xf(i+1)-xf(i))*(yf(j+1)-yf(j))* &
                    (pp(i+1,j)-pp(i-1,j))/(xp(i+1)-xp(i-1))
      v(i,j)=v(i,j)-(1./apv(i,j))*(yf(j+1)-yf(j))*(xf(i+1)-xf(i))* &
                    (pp(i,j+1)-pp(i,j-1))/(yp(j+1)-yp(j-1))
      END DO; END DO

!     compute pressure boundary values by linear extrapolation
!     currently for uniform cell sizes************************
      p(0,1:ny)=0.5*(3.*p(1,1:ny)-p(2,1:ny))
      p(nx+1,1:ny)=0.5*(3.*p(nx,1:ny)-p(nx-1,1:ny))
      p(1:nx,0)=0.5*(3.*p(1:nx,1)-p(1:nx,2))
      p(1:nx,ny+1)=0.5*(3.*p(1:nx,ny)-p(1:nx,ny-1))

!     update uold aud vold
      uold=u
      vold=v
      END SUBROUTINE PRESSURE
!*************************************************************
!     turbulence kinetic energy equation subroutine
      SUBROUTINE KINETIC_ENERGY
      USE VARIABLES
      IMPLICIT NONE
      INTEGER :: i,j,iter,iii
      REAL :: s11,s12,s21,s22,mur,mul,muu,mud,mdotr,mdotl,mdotu,mdotd,&
              tauwall1

!     first order upwinding for convection terms
      DO j=1,ny
      DO i=1,nx
      mdotr=rho*(yf(j+1)-yf(j))*((1.-fx(i+1))*u(i+1,j) + fx(i+1)*u(i,j))
      mdotl=rho*(yf(j+1)-yf(j))*((1.-fx(i))*u(i,j) + fx(i)*u(i-1,j))
      mdotu=rho*(xf(i+1)-xf(i))*((1.-fy(j+1))*v(i,j+1) + fy(j+1)*v(i,j))
      mdotd=rho*(xf(i+1)-xf(i))*((1.-fy(j))*v(i,j) + fy(j)*v(i,j-1))
      mur=(1.-fx(i+1))*mut(i+1,j) + fx(i+1)*mut(i,j) + mudyn
      mul=(1.-fx(i))*mut(i,j) + fx(i)*mut(i-1,j) + mudyn
      muu=(1.-fy(j+1))*mut(i,j+1) + fy(j+1)*mut(i,j) + mudyn
      mud=(1.-fy(j))*mut(i,j) + fy(j)*mut(i,j-1) + mudyn
      ar(i,j)=MAX(-mdotr,0.)+mur/sigk*(yf(j+1)-yf(j))/(xp(i+1)-xp(i))
      al(i,j)=MAX( mdotl,0.)+mul/sigk*(yf(j+1)-yf(j))/(xp(i)-xp(i-1))
      au(i,j)=MAX(-mdotu,0.0)+muu/sigk*(xf(i+1)-xf(i))/(yp(j+1)-yp(j))
      ad(i,j)=MAX(mdotd,0.0)+mud/sigk*(xf(i+1)-xf(i))/(yp(j)-yp(j-1))
      END DO; END DO
!********************************************************************      
!     correct right/left boundaries (viscous terms)
      DO j=1,ny
      mur=mut(nx+1,j) + mudyn
      ar(nx,j)=MAX(-rho*(yf(j+1)-yf(j))*uold(nx+1,j),0.) + &
                         mur*(yf(j+1)-yf(j))/(xp(nx+1)-xp(nx))
      mul=mut(0,j) + mudyn
      al(1,j)=MAX(rho*(yf(j+1)-yf(j))*uold(0,j),0.) + &
                         mul*(yf(j+1)-yf(j))/(xp(1)-xp(0))
      END DO
!     correct n/s boundaries aud set dk/dy=0 ad wall bcs
      au(:,ny)=0
      ad(:,1)=0

!     set ap
      apk(1:nx,1:ny)=(ar(1:nx,1:ny)+al(1:nx,1:ny)+ &
                      au(1:nx,1:ny)+ad(1:nx,1:ny))/omegak
!********************************************************************
!      compute tke production terms
       DO i=1,nx
       DO j=1,ny
       s11=(u(i+1,j)-u(i-1,j))/(xp(i+1)-xp(i-1))
       s22=(v(i,j+1)-v(i,j-1))/(yp(j+1)-yp(j-1))
       s12=0.5*((u(i,j+1)-u(i,j-1))/(yp(j+1)-yp(j-1))+ &
          (v(i+1,j)-v(i-1,j))/(xp(i+1)-xp(i-1)))
       s21=s12
       production(i,j)=2.0*mut(i,j)*(s11*s11+s22*s22+s12*s12+s21*s21)
       volume=(xf(i+1)-xf(i))*(yf(j+1)-yf(j))
       production(i,j)=volume*production(i,j)
       dissipation(i,j)=volume*rho*eps(i,j)
       END DO;  END DO

!      production term in wall adjacent cells
      tauwall=0.01
       DO i=1,nx

         ystar=rho*cmu**0.25*SQRT(k(i,1))*(yp(1)-yp(0))/mudyn
         tauwall=rho*kappa*cmu**0.25*k(i,1)**0.5*u(i,1)/LOG(E*ystar)


        production(i,1)=tauwall*cmu**0.25*k(i,1)**0.5/(kappa*(yp(1)-yp(0)))
		
!       multiply by cell area
        production(i,1)=production(i,1)*(xf(i+1)-xf(i))*(yf(2)-yf(1))


          ystar=rho*cmu**0.25*SQRT(k(i,ny))*(yp(ny+1)-yp(ny))/mudyn
          tauwall=rho*kappa*cmu**0.25*k(i,ny)**0.5*u(i,ny)/LOG(E*ystar)


        production(i,ny)=tauwall*cmu**0.25*k(i,ny)**0.5/(kappa*(yp(ny+1)-yp(ny)))
!       multiply by cell area
        production(i,ny)=production(i,ny)*(xf(i+1)-xf(i))*(yf(ny+1)-yf(ny))
       END DO
!********************************************************************
       DO iter=1,10   
       DO j=1,ny
       DO i=1,nx
       k(i,j)=(1.-omegak)*kold(i,j)+1./apk(i,j)*( &
                    ar(i,j)*k(i+1,j) + &
                    al(i,j)*k(i-1,j) + &
                    au(i,j)*k(i,j+1) + &
                    ad(i,j)*k(i,j-1) + &
                    production(i,j)-dissipation(i,j))
      END DO; END DO; 
      k(nx+1,:)=k(nx,:)    !zero derivative at outlet
      END DO
      k(0,:)=k(nx+1,:)    !cyclic  bc
      kold=k
      END SUBROUTINE KINETIC_ENERGY
!*************************************************************
!     dissipation rate equation subroutine
      SUBROUTINE DISSIPATE
      USE VARIABLES
      IMPLICIT NONE
      INTEGER :: i,ii,j,iter
      REAL :: s11,s12,s21,s22,mur,mul,muu,mud,mdotr,mdotl,mdotu,mdotd

!     first order upwinding for convection terms
      DO j=2,ny-1    !need not compute wall adjacent cells
      DO i=1,nx
      mdotr=rho*(yf(j+1)-yf(j))*((1.-fx(i+1))*u(i+1,j) + fx(i+1)*u(i,j))
      mdotl=rho*(yf(j+1)-yf(j))*((1.-fx(i))*u(i,j) + fx(i)*u(i-1,j))
      mdotu=rho*(xf(i+1)-xf(i))*((1.-fy(j+1))*v(i,j+1) + fy(j+1)*v(i,j))
      mdotd=rho*(xf(i+1)-xf(i))*((1.-fy(j))*v(i,j) + fy(j)*v(i,j-1))
      mur=(1.-fx(i+1))*mut(i+1,j) + fx(i+1)*mut(i,j) + mudyn
      mul=(1.-fx(i))*mut(i,j) + fx(i)*mut(i-1,j) + mudyn
      muu=(1.-fy(j+1))*mut(i,j+1) + fy(j+1)*mut(i,j) + mudyn
      mud=(1.-fy(j))*mut(i,j) + fy(j)*mut(i,j-1) + mudyn
      ar(i,j)=MAX(-mdotr,0.)+mur/sigeps*(yf(j+1)-yf(j))/(xp(i+1)-xp(i))
      al(i,j)=MAX( mdotl,0.)+mul/sigeps*(yf(j+1)-yf(j))/(xp(i)-xp(i-1))
      au(i,j)=MAX(-mdotu,0.0)+muu/sigeps*(xf(i+1)-xf(i))/(yp(j+1)-yp(j))
      ad(i,j)=MAX(mdotd,0.0)+mud/sigeps*(xf(i+1)-xf(i))/(yp(j)-yp(j-1))
      END DO; END DO
!********************************************************************      
!     correct right/left boundaries (viscous terms)
      DO j=2,ny-1
      mur=mut(nx+1,j) + mudyn
      ar(nx,j)=MAX(-rho*(yf(j+1)-yf(j))*uold(nx+1,j),0.) + &
                         mur*(yf(j+1)-yf(j))/(xp(nx+1)-xp(nx))
      mul=mut(0,j) + mudyn
      al(1,j)=MAX(rho*(yf(j+1)-yf(j))*uold(0,j),0.) + &
                         mul*(yf(j+1)-yf(j))/(xp(1)-xp(0))
      END DO

!     set final ap
      apeps(1:nx,2:ny-1)=(ar(1:nx,2:ny-1)+al(1:nx,2:ny-1)+ &
                      au(1:nx,2:ny-1)+ad(1:nx,2:ny-1))/omegarps
!********************************************************************
!      compute production terms
       DO i=1,nx
       DO j=2,ny-1
       s11=(u(i+1,j)-u(i-1,j))/(xp(i+1)-xp(i-1))
       s22=(v(i,j+1)-v(i,j-1))/(yp(j+1)-yp(j-1))
       s12=0.5*((u(i,j+1)-u(i,j-1))/(yp(j+1)-yp(j-1))+ &
          (v(i+1,j)-v(i-1,j))/(xp(i+1)-xp(i-1)))
       s21=s12
       production(i,j)=cep1*eps(i,j)/k(i,j)*2.0*mut(i,j)* &
                       (s11*s11+s22*s22+s12*s12+s21*s21)
       volume=(xf(i+1)-xf(i))*(yf(j+1)-yf(j))
       production(i,j)=volume*production(i,j)
       dissipation(i,j)=volume*cep2*rho*eps(i,j)**2/k(i,j)
       END DO;  END DO
!********************************************************************
       OUTERLOOP : DO iter=1,10   !cau chauge max iterations from 10 if desired
       DO j=2,ny-1
       DO i=1,nx
       eps(i,j)=(1.-omegarps)*epsold(i,j)+1./apeps(i,j)*( &
                    ar(i,j)*eps(i+1,j) + &
                    al(i,j)*eps(i-1,j) + &
                    au(i,j)*eps(i,j+1) + &
                    ad(i,j)*eps(i,j-1) + &
                    production(i,j)-dissipation(i,j))
 
      END DO; END DO;
      eps(nx+1,:)=eps(nx,:)    !zero derivative at outlet
        DO ii=1,nx    ! set near wall cell values
           eps(ii,1)=cmu**0.75*k(ii,1)**1.5/(kappa*(yp(1)-yp(0)))
           eps(ii,ny)=cmu**0.75*k(ii,ny)**1.5/(kappa*(yp(ny+1)-yp(ny)))
        END DO
      END DO OUTERLOOP
      if(icyclic .eq. 1)eps(0,:)=eps(nx+1,:)   !cyclic bcs
      epsold=eps
      END SUBROUTINE DISSIPATE
!*****************************************************************
!     routine to compute turbulent viscosity
      SUBROUTINE TURB_VISC
      USE VARIABLES
      IMPLICIT NONE
      INTEGER :: i,j
      DO i=1,nx
      DO j=1,ny
      mut(i,j)=rho*cmu*k(i,j)**2/eps(i,j)
      END DO; END DO
      mut(nx+1,:)=mut(nx,:)
      if(icyclic .eq. 1)mut(0,:)=mut(nx+1,:)  !cyclic bc
      RETURN
      END SUBROUTINE TURB_VISC
!***************************************************************
      SUBROUTINE INITIALIZE
      USE VARIABLES
      IMPLICIT NONE
      INTEGER :: status,nxx,nyy,i,j
     

!     set the domain size
      xmin=0.0; xmax=10.0
      ymin=0.0; ymax=1.0
      WRITE(*,*)'Enter number cells in x-direction'
      READ(*,*)nx
      WRITE(*,*)'Enter number cells in y-direction'
      READ(*,*)ny
      WRITE(*,*)'Enter the number of outer iterations.'
      READ(*,*)maxiter

      nxx=nx+1
      nyy=ny+1
!     allocate memory
      ALLOCATE(u(0:nxx,0:nyy),v(0:nxx,0:nyy),p(0:nxx,0:nyy), &
      pp(0:nxx,0:nyy),source(0:nxx,0:nyy),ar(0:nxx,0:nyy),al(0:nxx,0:nyy), &
      au(0:nxx,0:nyy),ad(0:nxx,0:nyy),apu(0:nxx,0:nyy), &
      app(0:nxx,0:nyy),uold(0:nxx,0:nyy),vold(0:nxx,0:nyy), &
      ustar(0:nxx,0:nyy),vstar(0:nxx,0:nyy),apv(0:nxx,0:nyy),&
      dpdx(0:nxx,0:nyy),dpdxc(0:nxx,0:nyy),dpdxave(0:nxx,0:nyy), &
      dpdy(0:nxx,0:nyy),dpdyc(0:nxx,0:nyy),dpdyave(0:nxx,0:nyy), &
      xp(0:nxx),yp(0:nyy), xf(0:nxx), yf(0:nyy),fx(0:nxx),fy(0:nyy), &
      mut(0:nxx,0:nyy), sx(0:nxx,0:nyy), sy(0:nxx,0:nyy), &
      production(0:nxx,0:nyy),k(0:nxx,0:nyy),eps(0:nxx,0:nyy), kold(0:nxx,0:nyy), &
      epsold(0:nxx,0:nyy),dissipation(0:nxx,0:nyy),apk(0:nxx,0:nyy), &
      apeps(0:nxx,0:nyy), s(0:nxx,0:nyy), STAT=status)

      IF(status /= 0)then
              WRITE(*,*)'Error in allocating memory.'
              STOP
      END IF

!     initialize variables 
      E=EXP(kappa*B)     !log-lal constaut
      u=0; v=0; p=0; pp=0
      u(0:nx+1,1:ny)=uinlet   !initialize flowfield
      uold=u; vold=v
      app=1.; apu=1.; apv=1.; k=0.; eps=0.; s=0.
      do i=0,nx+1
      do j=1,ny
         k(i,j)=2./3.*(uinlet*ti)**2
         eps(i,j)=(cmu**0.75)*(k(i,j)**1.5)/(0.07*(ymax-ymin))
         mut(i,j)=rho*cmu*k(i,j)**2/eps(i,j)
      end do; end do
      mut(:,0)=0;  mut(:,ny+1)=0
      kold=k; epsold=eps

!     Grid Setup: finite volume cell faces defined first
!     grid points then at centers of control volumes (Pataukar, 1980, practice B)

!     mesh with uniform spacing (with nx by ny control volumes)
      dx=(xmax-xmin)/float(nx); dy=(ymax-ymin)/float(ny)
      xf(1)=xmin
      xf(nx+1)=xmax
      DO i=2,nx
          xf(i)=xmin+dx*FLOAT(i-1)
      END DO

      yf(1)=ymin
      yf(ny+1)=ymax
      DO j=2,ny
          yf(j)=ymin+dy*FLOAT(j-1)
      END DO


!     compute cell center locations
      DO i=1,nx
      xp(i)=0.5*(xf(i+1)+xf(i))
      END DO
      xp(0)=xmin; xp(nx+1)=xmax

      DO j=1,ny
      yp(j)=0.5*(yf(j+1)+yf(j))
      END DO
      yp(0)=ymin;  yp(ny+1)=ymax

!     compute stretching factors
      do i=2,nx
         fx(i)=(xp(i)-xf(i))/(xp(i)-xp(i-1))
      end do
!     compute stretching factors
      do j=2,ny
         fy(j)=(yp(j)-yf(j))/(yp(j)-yp(j-1))
      end do

      END SUBROUTINE INITIALIZE
!*************************************************************
!     subroutine to write .csv files for use in ParaView
      SUBROUTINE OUTPUT
      USE VARIABLES
      IMPLICIT NONE
      INTEGER :: i,j
      REAL :: velmag, umax   !velmag is the velocity magnitude

!     open formatted output file for plotting data using Paraview
!     you may edit the file to see the format
      OPEN(UNIT=50, FILE='RESULTS1.csv')
      OPEN(UNIT=80, FILE='RESULTS2.csv')

!     fill in southleft corner (ad it is not computed in cfd procedure)
      u(0,0)=0
      v(0,0)=0
      p(0,0)=p(1,1)

!     fill in southright corner
      u(nx+1,0)=0
      v(nx+1,0)=0
      p(nx+1,0)=p(nx,1)

!     fill in northright corner
      u(nx+1,ny+1)=0 
      v(nx+1,ny+1)=0
      p(nx+1,ny+1)=p(nx,ny)

!     fill in northleft corner
      u(0,ny+1)=0
      v(0,ny+1)=0
      p(0,ny+1)=p(1,ny)

!     write .csv file for ParaView
      WRITE(50,*)'x, y, z, pressure,velmag, vx, vy, vz'
      WRITE(80,*)'x, y, z, k, eps, mut'
      DO i=0,nx+1
      DO j=0,ny+1
        velmag=SQRT(u(i,j)**2+v(i,j)**2)
        WRITE(50,100)xp(i),yp(j),0.0,p(i,j),velmag,u(i,j),v(i,j),0.0
        WRITE(80,100)xp(i),yp(j),0.0,k(i,j),eps(i,j),mut(i,j)
      END DO; END DO
100   FORMAT(f8.5,', ',f8.5,', ',f8.5,', ',f8.5,', ',f8.5,', ',f8.5,&
             ', ',f8.5,', ',f8.5,', ',f8.5)
print*,'finished'
      CLOSE(UNIT=50); CLOSE(UNIT=80)
!     write out u/u_max
      umax=MAXVAL(u(nx,:))
      print*,'maximum centerline value = ',umax
      print*,'wall ystar = ',rho*cmu**0.25*SQRT(k(nx,1))*(yp(1)-yp(0))/mudyn
      write(90,*)'wall ystar = ',rho*cmu**0.25*SQRT(k(nx,1))*(yp(1)-yp(0))/mudyn
      do j=0,(ny+1)
         write(90,101)yp(j)/(ymax),u(nx,j)/umax  !,k(nx,j),eps(nx,j)
      end do
101   format(4f10.5)
      END SUBROUTINE output
!**********************************************************************
 
