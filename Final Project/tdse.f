*
      program prop1d
      implicit double precision (a-h,o-z)
*
*  This program uses the staggered leap-frog method to 
*  propagate the solution of the one-dimensional
*  Schroedinger Equation in time.
*
*  written by:    Klaus Bartschat
*                 Department of Physics and Astronomy
*                 Drake University
*                 Des Moines, Iowa 50311, USA
*
*  last update:   May 6, 2023
*
      parameter(nmax=80000)
      parameter(max=2400000)
      parameter(zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0)
      complex psi0,psi1,psi2
      complex czero
      dimension y(nmax),rpart(nmax),orbin(nmax),z(nmax)
      ! dimension envelope(max)
*     
      common/psi / psi0(nmax),psi1(nmax),psi2(nmax),psiin(nmax)
      common/pot / rval(nmax),vpot(nmax),field  !(max)
      common/rest/ tmax,dt,h,npts
*
 1002 format(1p10e14.6)
 1003 format(1p,5e14.5)
 1004 format('#    time         norm',
     >       '         surv')
*
      pi = acos(-one)
      czero = dcmplx(zero,zero)
*
* set up the radial mesh
*
      read(5,*) h,npts,iprinx
      do 10 n=1,npts+1
* 2n+1 points
       rval(n) = h*(n)
       vpot(n) = -1.0d0/(sqrt(1+(rval(n)-(npts*h)/2)**2))
 10   continue
*
* set up the potential
*
* HAVE TO SET ELECTRIC FIELD AS PER PAPER!!
      read(5,*) tmax, dt, iprint

      T = NINT(tmax/dt) !total number of timesteps

      omega = 0.148d0 !Trevin saw this on Barty's screen
      T_1 = 2*pi/omega !period
      nperiods = tmax/T_1
      ! do n=1, T
      !       envelope(n) = (sin(pi*(dt*(n-1))/tmax))**2
      !       field(n) = 0.1d0 * envelope(n) * sin(omega*(dt*(n-1))) !
      ! end do
*
* set up the initial state
*
*This will be changed. We can just read in bound state energy and other stuff
* wfn.out for him contained initial x_vals, dummy, orbin
******
      read(5,*) nin,ebound
!initialize psi0(n) and psi1(n) as czero for all points
      open(1,file='wfn.out')
      ! read(1,*) rin, dummy, orbin
      do n=1,nin+1
            read(1,*) rin, dummy, orbin(n)
      end do
      nstart = nin/2+1
      do n=1, nin+1
            psi0(n) = cmplx(orbin(n),0.0d0)
            psiin(n) = psi0(n)
            psi1(n) = psi0(n)*cmplx(cos(ebound*dt),-sin(ebount*dt))
      end do
******

*
* check the normalization
*
      itout  = 99
      ymax = 0.0d0
      do 50 n=1,npts
       y(n) = cabs(psi0(n)**2)
       if (y(n).ge.ymax) ymax = y(n)
 50   continue
      call integr(rval,y,npts,result,error)
*
* print the initial wavefunction and the potential
*
      open(99,file='potential.dat')
      if (ehart.gt.height*1.002) then 
       pfac = 0.8d0*ymax/abs(height)
      else 
       if (ehart.lt.height) then 
        pfac = 1.2d0*ymax/abs(height)
       else 
        pfac = ymax/abs(height)
       endif
      endif
*
      itout  = 100
      do 60 n=1,npts,iprinx
       write(100,1003) rval(n),cabs(psi0(n))**2
 60   continue
      do 61 n=1,npts
       write(itout-1,1003) rval(n),pfac*vpot(n)
 61   continue
      call integr(rval,y,npts,result,error)
*
* now propagate the solution
*

      !       envelope(n) = (sin(pi*(dt*(n-1))/tmax))**2
      !       field(n) = 0.1d0 * envelope(n) * sin(omega*(dt*(n-1))) !
      E0 = 0.10d0
      open(7,file='results.dat')
      write(7,1004)
      ntimes = nint(tmax/dt)
      do 110 nt=1,ntimes-1   !he has nt=1,times-1, we had up to ntimes originally
       ttt=nt*dt !new
       envelope = (sin(pi*ttt/tmax))**2 !new
       field = E0*envelope*sin(omega*ttt) !new
       call propagate !!!!!! field
       if (nt/iprint*iprint.eq.nt) then
        itout = itout+1
!         do 70 n=1,npts,iprinx
!          write(itout,1003) rval(n),cabs(psi2(n))**2
!  70     continue
*
* extract the physics (reflection and transmission coefficients)
*
**** NOW WE WANT CHECK NORM AND EXTRACT SURVIVAL PROBABILITY: CHANGE THIS SECTION
      ! nlast = 100 !DUMMY SO THAT WE DONT SEG FAULT
        do 80 n=1,npts
         y(n) = cabs(psi2(n))**2
 80     continue
        call integr(rval,y,npts,xnorm,error)
        ttt = nt*dt
        do n=1, npts
        y(n) = real(psi2(n)*psiin(n))
        z(n) = aimag(psi2(n)*psiin(n))

        end do
        call integr(rval,y,npts,surv_real,error)
        call integr(rval,z,npts,surv_imag,error)
        surv = surv_real**2 + surv_imag**2
        write(7,1002) ttt,xnorm,surv
       endif
*
* move the wavefunctions
*
       do 100 n=1,npts
        psi0(n) = psi1(n)
        psi1(n) = psi2(n)
        psi2(n) = czero
100    continue
110   continue
*
      stop
      end
************************************************************************
      subroutine propagate !ttt is current timestep
      implicit double precision (a-h,o-z)
*
      parameter(nmax=80000)
      parameter(zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0)
*
      complex psi0,psi1,psi2,Hpsi1
      complex czero,ci
      dimension Hpsi1(nmax)
*
      common/psi / psi0(nmax),psi1(nmax),psi2(nmax),psiin(nmax)
      common/pot / rval(nmax),vpot(nmax),field
      common/rest/ tmax,dt,h,npts
*
      czero = dcmplx(0.0d0,0.0d0)
      ci    = dcmplx(0.0d0,1.0d0)
*
* apply the hamiltonian to psi1; begin with the edges
*
      Hpsi1(1) = czero
      Hpsi1(npts) = czero
*
* use 3-point formula for the rest 
*    
      do 10 n = 2,npts-1
         Hpsi1(n) = -half/(h*h)*(psi1(n-1)-two*psi1(n)+psi1(n+1))
 10   continue
*
* now leapfrog to get psi2; set to zero at the edges to enforce
* the boundary conditions
*
      psi2(1) = czero
      psi2(npts) = czero
!  
      do n = 2,npts-1
       psi2(n) = 
     > psi0(n)-two*ci*dt*(Hpsi1(n)+(vpot(n)-rval(n)*field)*psi1(n))
      end do
*
 2000 format(i5,1p5d14.4)
*
      end
************************************************************************
      SUBROUTINE integr(X,Y,N,RESULT,ERROR)
************************************************************************
*                                                                      *
*  THIS SUBROUTINE INTEGRATES ARBITRARILY SPACED DATA                  *
*                                                                      *
*                      INPUT                                           *
*                      -----                                           *
*                                                                      *
*  X .........     VECTOR CONTAINING THE MESHPOINTS                    *
*                  (EITHER IN ASCENDING OR DESCENDING ORDER)           *
*  Y .........     VECTOR CONTAINING THE FUNCTION VALUES               *
*  N .........     NUMBER OF MESHPOINTS   (AT LEAST 4)                 *
*                                                                      *
*                      OUTPUT                                          *
*                      ------                                          *
*                                                                      *
*  RESULT ....     RESULT OF THE INTEGRATION                           *
*  ERROR .....     ESTIMATED ERROR                                     *
*                                                                      *
************************************************************************
*                                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,                    
     +           FIVE=5.0D0,SIX=6.0D0,TEN=10.0D0,TWELVE=12.0D0,                 
     +           SIXTY=60.0D0,HUNTWE=120.0D0)
      PARAMETER (IREAD=5,IWRITE=6)
      DIMENSION X(N),Y(N)

      RESULT = ZERO
      ERROR = ZERO
*
*  CHECK THAT WE HAVE ENOUGH POINTS
*
      IF (N.LT.4) THEN
       WRITE(IWRITE,1000) N
       STOP
      ENDIF
*
*  CHECK THAT THE MESHPOINTS AR IN EITHER ASCENDING OR DESCENDING ORDER
*
      H2 = X(2)-X(1)
      DO 10 I=3,N
       H3 = X(I)-X(I-1)
       IF(H2*H3.LE.ZERO) THEN
        WRITE(IWRITE,1001)
        STOP
       ENDIF
10    CONTINUE
*
*  START THE INTEGRATION
*
      D3 = (Y(2)-Y(1))/H2
      H3 = X(3)-X(2)
      D1 = (Y(3)-Y(2))/H3
      H1 = H2+H3
      D2 = (D1-D3)/H1
      H4 = X(4)-X(3)
      R1 = (Y(4)-Y(3))/H4
      R2 = (R1-D1)/(H4+H3)
      H1 = H1+H4
      R3 = (R2-D2)/H1
      RESULT = H2*(Y(1)+H2*(D3/TWO-H2*(D2/SIX-(H2+TWO*H3)*R3/TWELVE)))
      S = -(H2**3)*(H2*(THREE*H2+FIVE*H4)+TEN*H3*H1)/SIXTY
      R4 = ZERO
      NN = N-1
*
*  LOOP OVER POINTS 2 TO N-1
*
      DO 20 I=3,NN
       RESULT = RESULT+H3*((Y(I)+Y(I-1))/TWO-H3*H3*(D2+R2+(H2-H4)*R3) 
     +         /TWELVE)
       C = H3**3*(TWO*H3*H3+FIVE*(H3*(H4+H2)+TWO*H4*H2))/HUNTWE
       ERROR = ERROR+(C+S)*R4
       IF (I.NE.3) THEN
        S = C
       ELSE
        S = S+TWO*C
       ENDIF
       IF (I.EQ.(N-1)) GOTO 30
       H1 = H2
       H2 = H3
       H3 = H4
       D1 = R1
       D2 = R2
       D3 = R3
       H4 = X(I+2)-X(I+1)
       R1 = (Y(I+2)-Y(I+1))/H4
       R4 = H4+H3
       R2 = (R1-D1)/R4
       R4 = R4+H2
       R3 = (R2-D2)/R4
       R4 = R4+H1
       R4 = (R3-D3)/R4
20    CONTINUE
30    CONTINUE
*
*  FINISH INTEGRATION
*
      RESULT = RESULT+H4*(Y(N)-H4*(R1/TWO+H4*(R2/SIX+(TWO*H3+H4)*R3             
     +        /TWELVE)))
      ERROR = ERROR-H4**3*R4*(H4*(THREE*H4+FIVE*H2)+TEN*H3*(H2+H3+H4)) 
     +        /SIXTY+S*R4
*
1000  FORMAT(//,' ERROR IN SUBROUTINE INTEGR . N = ',I2,3X,                     
     +      'PROGRAM STOPS.',//)                                                
1001  FORMAT(//,' ERROR IN SUBROUTINE INTEGR . MESHPOINTS ARE OUT ',            
     +      'OF ORDER. PROGRAM STOPS.',//)                                      
      RETURN
      END
****************  input data follow
*      0.2   9001   40                           h,npts,iprinx
*      2     900.0   1.0   30.0    2.500         nbar,r1,ra,rd,height
*      1.88  100.0 500.0  600.0    0.002   500   xk0,xw,x0,tmax,dt,iprint
