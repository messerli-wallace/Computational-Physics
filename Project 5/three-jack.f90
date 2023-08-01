!================================================
   program threebody2 
!================================================
!  The Three-body Problem: 
!  The Effect of Jupiter on Earth 
!------------------------------------------------
   implicit double precision (a-h,o-z)
   Real :: m1,m2,ms
   integer i
   pi = acos(-1.0)

!  set up initial conditions (in AU units): 
   open(1,file='threebody-final.inp')
   read(1,*) x1   
   read(1,*) x2    
   read(1,*) m1
   read(1,*) m2
   read(1,*) ms
   read(1,*) dt
   read(1,*) nt
   close(1)


   y1 = 0.0; vx1 = 0.0
   vy1 = 2.0 * pi / sqrt(x1)
   y2 = 0.0; vx2 = 0.0    
   vy2 = 2.0 * pi / sqrt(x2)
   y3 = 0.0d0; vxs = 0.0d0
   print*,vy1,vy2
   print*,x1,x2
   x3 = -(m1*x1 + m2*x2)/ms !sun position based off CM=0
   print*,x3
   !vy3 = -2.0d0 * pi / sqrt(abs(x3))
   !vy3 = -2.10805e-4 !
   vy3 = -(m1*vy1+m2*vy2)/(ms)
   print*,vy3

   open(2,file='threebody.out')
   write(2,'(a,8(5x,a,4x))') '#', ' t','x_earth','y_earth',&
            'x_jupiter','y_jupiter','x_sun','y_sun'
   Call Calculate(x1,vx1,y1,vy1,x2,vx2,y2,vy2,x3,vx3,y3,vy3,ms,m1,m2,dt,nt)

   call system ('gnuplot -load three-f.plt') !gnuplot call
   end program threebody2

!===============================================================
  subroutine calculate (x_earth, vx_earth, y_earth, vy_earth, &
     x_jupiter,vx_jupiter,y_jupiter,vy_jupiter,&
     x_sun,vx_sun,y_sun,vy_sun,ms,me,mj,dt,nt) 
!===============================================================
! use the Euler-Cromer method; update the velocities first
!---------------------------------------------------------------
  implicit double precision (a-l,o-z)
  Real :: me,mj,ms
  integer i
  pi = acos(-1.0); pp = 4*pi**2;  pj = pp*(mj/ms); pe = pp*(me/ms) 
  t = 0.0
  do  i=1,nt
   res   = sqrt(((x_earth-x_sun))**2 + ((y_earth-y_sun))**2) 
   rjs  = sqrt((abs(x_jupiter-x_sun))**2 + (abs(y_jupiter-y_sun))**2) 
   rej = sqrt(((x_earth-x_jupiter))**2 + ((y_earth-y_jupiter))**2) 
   vx_earth = vx_earth - pp * (x_earth-x_sun) * dt / res**3   &
                       - pj * (x_earth-x_jupiter) * dt / rej**3 
   vy_earth = vy_earth - pp * (y_earth-y_sun) * dt / res**3   &
                       - pj * (y_earth-y_jupiter) * dt / rej**3 
   vx_jupiter = vx_jupiter - pp * (x_jupiter-x_sun) * dt / rjs**3  &
                           - pe * (x_jupiter-x_earth) * dt / rej**3 
   vy_jupiter = vy_jupiter - pp * (y_jupiter-y_sun) * dt / rjs**3  &
                           - pe * (y_jupiter-y_earth) * dt / rej**3 
   vx_sun = vx_sun - pe * (x_sun-x_earth) * dt / res**3 &
                   - pj * (x_sun-x_jupiter) * dt / rjs**3
   vy_sun = vy_sun - pe * (y_sun-y_earth) * dt / res**3 &
                   - pj * (y_sun-y_jupiter) * dt / rjs**3
   x_earth = x_earth + vx_earth * dt 
   y_earth = y_earth + vy_earth * dt 
   x_jupiter = x_jupiter + vx_jupiter * dt 
   y_jupiter = y_jupiter + vy_jupiter * dt 
   x_sun = x_sun + vx_sun * dt
   y_sun = y_sun + vy_sun * dt
   t = t + dt

   write(2,'(7f14.7)')  t,x_earth,y_earth, x_jupiter,y_jupiter,x_sun,y_sun 
  end do
  end subroutine calculate

