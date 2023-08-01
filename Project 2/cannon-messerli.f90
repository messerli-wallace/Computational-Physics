!================================================
Program cannon 
!================================================
! a cannon ball trajectory using the Euler method
!------------------------------------------------
implicit none
integer, parameter :: nmax = 500000
double precision x(nmax),y(nmax),vx(nmax),vy(nmax),dt,vinit,B2m,theta,h,current_angle, &
                                                  top_angle, top_distance, temp, previous_dist
integer n, test, test_index, air_drag_type, gravity_adjustment, optimization_input
common / constants / dt, B2m, h, gravity_adjustment, air_drag_type, test

read (*,*) vinit          ! initial speed
read (*,*) dt             ! time step
read (*,*) B2m            ! drag/m, B2/m (see book)
read (*,*) theta          ! firing angle, ignored if test=1
read (*,*) h              ! height difference between launching and landing (m)
!note that when h>0, we are shooting uphill 
read (*,*) test           ! 1 to test for optimal angle
read (*,*) air_drag_type  ! 0 for Eq. 2.23, 1 for Eq. 2.24
read (*,*) gravity_adjustment ! 0 for constant gravity, else for varying gravity
read (*,*) optimization_input ! 0 for maximizing distance, else for acquiring 19.5 km

IF (test .NE. 1) THEN
  
  !Function call to begin trajectory calculation.
  Call calculate(x,y,vx,vy,dt,vinit,theta,B2m,n,h,air_drag_type,gravity_adjustment) 

  !Output function call.
  Call output(x,y,vx,vy,n,theta,B2m,dt,vinit,h,air_drag_type)
ELSE
  IF (optimization_input .EQ. 0) THEN


    top_distance = 0.0d0 !starting this variable for the first comparison
    DO test_index=1, 7000 !Test from 10-80 degrees
      current_angle = 10.0d0 + (test_index-1)*0.01d0 ! current angle

      !call calculate for current angle
      Call calculate(x,y,vx,vy,dt,vinit,current_angle,B2m,n,h,air_drag_type,gravity_adjustment)

      !check if current last distance is greater than the stored current max
      IF (x(n) .gt. top_distance) THEN
        top_distance = x(n)
        top_angle = current_angle
      END IF

    END DO
    print*, top_angle,top_distance !print max


  ELSE !calculate to get angle of 19.5 km distance
    top_distance = 19.5d0 !stopping point

    DO test_index=3000, 9000 !Test from 30-90 degrees
      current_angle = 10.0d0 + (test_index-1)*0.01d0 ! current angle

      !call calculate for current angle
      Call calculate(x,y,vx,vy,dt,vinit,current_angle,B2m,n,h,air_drag_type,gravity_adjustment)
      write(1,*) current_angle,x(n)
      !check if current last distance is greater than our desired amount
      !IF (x(n) .lt. top_distance) THEN
        !print*, current_angle, x(n) !print the straddling distances and their angle
        !print*, current_angle-0.1d0, previous_dist
        !exit
      !END IF
      !previous_dist = x(n) !store distance in case the next angle is the last one
    END DO

  END IF
END IF

End program cannon

!=====================================================
Subroutine calculate (x,y,vx,vy,dt,vinit,theta,B2m,n,h,air_drag_type,gravity_adjustment)
!=====================================================
implicit none
integer, parameter :: nmax = 500000
integer i,n, air_drag_type, gravity_adjustment
double precision x(nmax),y(nmax),dt,vinit,theta,B2m,vx(nmax),vy(nmax), &
                                        dvx,dvy,dx,dy,pi,h
!set initial position to the origin
x(1) = 0.0d0
y(1) = 0.0d0
!establish pi for calculations
pi = acos(-1.0d0)
!initial velocity component calculation
vx(1) = vinit*cos(theta/180.0d0*pi)
vy(1) = vinit*sin(theta/180.0d0*pi)
!print*, vx(1)
!Computation of trajectory
Do i = 2,nmax
  !Call deriv to calculate change in position and velocity
  Call deriv(x(i-1),y(i-1),vx(i-1),vy(i-1),B2m,dx,dy,dvx,dvy,h,air_drag_type,gravity_adjustment)
  
  !update position and velocity components
  x(i) = x(i-1) + dt*dx
  y(i) = y(i-1) + dt*dy
  vx(i) = vx(i-1) + dt*dvx
  vy(i) = vy(i-1) + dt*dvy
  n = i	   
  !print*, vx(i),vy(i)

  !stop when the cannonball lands at the desired height
  if (ABS(vy(i)) .NE. vy(i)) THEN
    if (y(i).le.h) exit
  End if
  
End do

End Subroutine calculate 

!======================================================
Subroutine deriv (x_in,y_in,vx_in,vy_in,B2m,dx,dy,dvx,dvy,h,air_drag_type,gravity_adjustment)
!======================================================
implicit none
integer, parameter :: nmax = 500000
double precision x_in,y_in,vx_in,vy_in,B2m,dx,dy,dvx,dvy,f_drag,h,rho,rho_0,temp,g
integer air_drag_type, gravity_adjustment

!check air_drag_type to determine air drag
IF (air_drag_type .EQ. 0) THEN  !Eq. 2.23
  f_drag = B2m*EXP(-y_in / (10000.0d0))*sqrt(vx_in**2+vy_in**2)
ELSE IF (air_drag_type .EQ. 1) THEN !Eq. 2.24
  temp = 290.0d0
  f_drag = B2m*((1.0d0-((6.5d0*(10d0**(-3.0d0)))*y_in/temp))**2.50d0)*sqrt(vx_in**2+vy_in**2)
ELSE
  f_drag = B2m*sqrt(vx_in**2+vy_in**2) !default
END IF

!Change gravity if gravity_adjustment is not 0
IF (gravity_adjustment .NE. 0) THEN
  g = (3.986019* (10.0d0**14))/(((6.3781d0*10.0d0**6.0d0)+y_in)**2.0d0)
ELSE
  g = 9.80665d0 !constant gravity
END IF

!calculate change in speed
dx = vx_in
dy = vy_in
dvx = -f_drag*vx_in
dvy = -f_drag*vy_in-g
End Subroutine deriv

!=============================================
Subroutine output (x,y,vx,vy,n,theta,B2m,dt,vinit,h,air_drag_type)
!=============================================
implicit none
integer, parameter :: nmax = 500000
double precision x(nmax),y(nmax),vx(nmax),vy(nmax),theta,B2m,dt,vinit,r,xl,h
integer n,i,air_drag_type
Character(80) :: AF
write(6,1000) theta
1000 format('theta = ',f6.3) 
write(AF,'(a,f5.2,I0)') 'cannon_',theta,air_drag_type
open(1,file=AF)

write(1,'(a)') '# Cannonball in air: Euler'
write(1,'(a)') '# ' 
write(1,'(a,f10.6)') '# initial speed    -> ', vinit
write(1,'(a,f10.6)') '# time step        -> ', dt
write(1,'(a,f10.6)') '# drag/m           -> ', B2m
write(1,'(a,f10.6)') '# firing angle     -> ', theta
write(1,'(a,f10.6)') '# height           -> ', h
write(1,'(a,I5)') '# air drag type    -> ', air_drag_type
write(1,'(a)') '# ' 
write(1,'(a)') '# Height vs Distance (in km)'
write(1,'(a)') '# ' 
write(1,'(a)') '# x position, y position, x velocity, y velocity' 
!
! Let's make sure the last point is zero rather than negative
!
r = -y(n-1)/y(n)
! print *,y(n-1),y(n),r
xl = (x(n-1)+r*x(n))/(r+1.0d0)
! print *,'range = ',xl
!
x(n) = xl
y(n) = h
write(1,'(4f14.6)') (x(i)/1000.0d0,y(i)/1000.0d0,vx(i)/1.0d0,vy(i)/1.0d0,i=1,n)
close(1)

End Subroutine output
