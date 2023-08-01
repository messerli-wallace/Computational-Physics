!================================================
Program baseball 
!================================================
! a baseball trajectory using the Euler method
!------------------------------------------------
implicit none
integer, parameter :: nmax = 500000
double precision :: field_altitude = 1609.344d0 ![m]=5280 ft above sea level
double precision x(nmax),y(nmax),z(nmax),vx(nmax),vy(nmax),vz(nmax),dt,vinit,B2m, &
theta, current_angle, omega, top_angle, h_0, final_x, final_y, final_z, &
center_0, headwind_speed, headwind_angle, theta_answer, omega_answer, &
previous_dist
integer n,i,test,ball_spin_index, theta_index
common / constants / dt, h_0, center_0, vinit, final_x, final_y, final_z, &
headwind_speed, headwind_angle, field_altitude

! Define constants given in the assignment
vinit = 31.2928d0         ! Initial speed of 70 miles/hour
h_0 = 1.58496d0           ! Initial height (y) the ball is thrown at [m]
center_0 = 0.33528d0      ! Initial distance (z) from the centerline [m]
final_x = 18.4404d0       ! 60.5ft, length between pitcher and home plate [m]
headwind_speed = 2.2352d0 ! Headwind speed [m/s]
headwind_angle = 28.0d0   ! Angle of headwind
dt = 0.01d0                ! Time step [s]
!read input for calculations
read(*,*) final_y        ! Final height [m]
read(*,*) final_z         ! Final distance from centerline [m]
read(*,*) test            ! Test argument: 0 for theta, 1 for omega, 3 to plot with following values
read(*,*) theta_answer    ! Theta that gets answer, after testing for theta
read(*,*) omega_answer    ! Omega that gets answer, after testing for omega 

!print*,vinit, h_0, center_0,final_x
!print*, headwind_speed,headwind_angle,dt
!print*,final_y, final_z, test, theta_answer

!Call calculate(x,y,z,vx,vy,vz,10.0d0,n,0.0d0)
!print*, (z(i),i=1,n) !works for x, y, and z
!Call output(x,y,z,vx,vy,vz,n,10.0d0,0.0d0)

IF (test.eq.0) THEN
  omega = 0.0d0
  DO theta_index = 1110, 10000 !Test from 1.110-10 degrees
    current_angle = (theta_index)*0.001d0 ! current angle

    !call calculate for current angle
    Call calculate(x,y,z,vx,vy,vz,current_angle,n,omega)
    write(1,*) current_angle,x(n),y(n) !write the angle and final distance to track in file
  END DO
END IF

IF (test.eq.1) THEN !testing for ball spin
  DO ball_spin_index = 0, 3000 !Test ball spin 0-300 rad/s with interval of 0.1
    omega = ball_spin_index * (0.1d0)
    !call calculate for current ball spin rate
    Call calculate(x,y,z,vx,vy,vz,theta_answer,n,omega)
    write(2,*) omega,x(n),y(n),z(n) !write omega and position
  END DO
END IF

IF (test.eq.2) THEN !getting plot for final answer
  Call calculate(x,y,z,vx,vy,vz,theta_answer,n,omega_answer)
  Call output(x,y,z,vx,vy,vz,n,theta_answer,omega_answer)
END IF
End program baseball

!=====================================================
Subroutine calculate (x,y,z,vx,vy,vz,theta,n,omega)
!=====================================================
implicit none
integer, parameter :: nmax = 500000
integer i,n
double precision x(nmax),y(nmax),z(nmax),dt,vinit,theta,vx(nmax),vy(nmax),vz(nmax), &
dvx,dvy,dvz,dx,dy,dz,pi,B2m,m,y_final, center_0, final_x, final_y, final_z, h_0, &
headwind_angle, headwind_speed, omega
common / constants / dt, h_0, center_0, vinit, final_x, final_y, final_z, &
headwind_speed, headwind_angle
!print*,final_x !works

!set initial position
x(1) = 0.0d0
y(1) = h_0
z(1) = center_0

!establish pi for calculations
pi = acos(-1.0d0)

!initial velocity component calculation
vx(1) = vinit*cos(theta/180.0d0*pi)
vy(1) = vinit*sin(theta/180.0d0*pi)
vz(1) = 0.0d0
!print*, vx(1)

!Computation of trajectory
Do i = 2,nmax
  !Call deriv to calculate change in position and velocity
  Call deriv(x(i-1),y(i-1),z(i-1),vx(i-1),vy(i-1),vz(i-1),dx,dy,dz,dvx,dvy,dvz,omega)
  
  !update position and velocity components
  x(i) = x(i-1) + dt*dx
  y(i) = y(i-1) + dt*dy
  z(i) = z(i-1) + dt*dz
  vx(i) = vx(i-1) + dt*dvx
  vy(i) = vy(i-1) + dt*dvy
  vz(i) = vz(i-1) + dt*dvz
  n = i
  !print*, x(i), y(i), z(i) !working

  !stop if the ball has reached the ground.
  !if (y(i).lt.0.0d0) exit

  !stop when the ball is at the desired x (at the mount)
  if (x(i).ge.final_x) exit
  
End do

! Interpolate final x to acquire the proper distance, only if the ball made it
IF (y(n).ge.0.0d0) THEN
  m = (y(n)-y(n-1))/(x(n)-x(n-1))
  y_final = m*(final_x-x(n-1))+y(n-1)
  x(n) = final_x
  y(n) = y_final
END IF
!print*,x(n),y(n),z(n) !working

End Subroutine calculate 

!======================================================
Subroutine deriv (x_in,y_in,z_in,vx_in,vy_in,vz_in,dx,dy,dz,dvx,dvy,dvz,omega)
!======================================================
implicit none
integer, parameter :: nmax = 500000
double precision x_in,y_in,z_in,vx_in,vy_in,vz_in,dx,dy,dz,dvx,dvy,dvz, &
f_drag,rho,rho_0,g,B2m, v_d, delta, magnitude_velocity, magnus_force, &
center_0, dt, field_altitude, final_x, final_y, final_z, h_0, &
headwind_angle, headwind_speed, omega, vinit
common / constants / dt, h_0, center_0, vinit, final_x, final_y, final_z, &
headwind_speed, headwind_angle, field_altitude

!determine air drag
v_d = 35 !m/s, from text
delta = 5 !m/s, from text
magnitude_velocity = sqrt(vx_in**2+vy_in**2+vz_in**2)
B2m = 0.0039d0 + 0.0058d0/(1+EXP((vy_in - v_d)/delta))
!print*,B2m
f_drag = B2m*EXP(-(y_in+field_altitude)/(10000.0d0))*magnitude_velocity
!print*,f_drag

!determine magnus force
magnus_force = (4.1d0*10.0d0**(-4.0d0))*vx_in*omega
!print*,magnus_force

g = 9.80665d0 !constant gravity
!calculate change in speed
dx = vx_in
dy = vy_in
dz = vz_in
dvx = -f_drag*vx_in
dvy = -f_drag*vy_in - g
dvz = -f_drag*vz_in - magnus_force
!print*,dvx,dvy,dvz !seems reasonable

End Subroutine deriv

!=============================================
Subroutine output (x,y,z,vx,vy,vz,n,theta,omega)
!=============================================
implicit none
integer, parameter :: nmax = 500000
double precision x(nmax),y(nmax),z(nmax),vx(nmax),vy(nmax),vz(nmax),&
dt,vinit,theta,omega
integer n,i
Character(80) :: AF
common / constants / dt, vinit

write(6,1000) theta
1000 format('theta = ',f6.3) 
write(AF,'(a,f5.3)') 'baseball_',theta
open(1,file=AF)

write(1,'(a)') '# Baseball in air: Euler'
write(1,'(a)') '# ' 
write(1,'(a,f10.6)') '# initial speed    -> ', vinit
write(1,'(a,f10.6)') '# time step        -> ', dt
write(1,'(a,f10.6)') '# firing angle     -> ', theta
write(1,'(a,f10.6)') '# spin speed[rad/s]-> ', omega
write(1,'(a,f10.6)') '# final height     -> ', y(n)
write(1,'(a)') '# ' 
write(1,'(a)') '#'
write(1,'(a)') '# ' 
write(1,'(a)') '# x position, y position, z position, x velocity, y velocity, z velocity' 


write(1,'(6f14.6)') (x(i),y(i),z(i),vx(i),vy(i),vz(i),i=1,n)
close(1)

End Subroutine output