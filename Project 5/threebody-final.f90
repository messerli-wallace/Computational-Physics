!================================================
   program threebody2 
!================================================
!  The Three-body Problem: 
!  The Effect of Jupiter on Earth 
!------------------------------------------------
   Real :: m1,m2,ms

   pi = acos(-1.0)

!  set up initial conditions (in AU units): 
   open(1,file='threebody-final.inp')
   read(1,*) x1
   y1 = 0.0; vx1 = 0.0
   vy1 = 2.0 * pi / sqrt(x1)   
   read(1,*) x2
   y2 = 0.0; vx2 = 0.0    
   vy2 = 2.0 * pi / sqrt(x2)    
   read(1,*) m1
   read(1,*) m2
   read(1,*) ms
   read(1,*) dt
   read(1,*) nt
   close(1)

   open(2,file='threebody-final.out')
   write(2,'(a,6(5x,a,4x))') '#', ' t','x_earth','y_earth',&
                               'x_jupiter','y_jupiter' 
   Call Calculate(x1,vx1,y1,vy1,x2,vx2,y2,vy2,ms,m1,m2,dt,nt)

   end program threebody2

!===============================================================
  subroutine calculate (x_earth, vx_earth, y_earth, vy_earth, &
     x_jupiter,vx_jupiter,y_jupiter,vy_jupiter,ms,me,mj,dt,nt) 
!===============================================================
! use the Euler-Cromer method; update the velocities first
!---------------------------------------------------------------
  Real :: me,mj,ms
  pi = acos(-1.0); pp = 4*pi**2;  pj = pp*(mj/ms); pe = pp*(me/ms) 
  t = 0.0
  do  i=1,nt
   r   = sqrt(x_earth**2 + y_earth**2) 
   rj  = sqrt(x_jupiter**2 + y_jupiter**2) 
   rej = sqrt((x_earth-x_jupiter)**2 + (y_earth-y_jupiter)**2) 
   vx_earth = vx_earth - pp * x_earth * dt / r**3   &
                       - pj * (x_earth-x_jupiter) * dt / rej**3 
   vy_earth = vy_earth - pp * y_earth * dt / r**3   &
                       - pj * (y_earth-y_jupiter) * dt / rej**3 
   vx_jupiter = vx_jupiter - pp * x_jupiter * dt / rj**3  &
                           - pe * (x_jupiter-x_earth) * dt / rej**3 
   vy_jupiter = vy_jupiter - pp * y_jupiter * dt / rj**3  &
                           - pe * (y_jupiter-y_earth) * dt / rej**3 
   x_earth = x_earth + vx_earth * dt 
   y_earth = y_earth + vy_earth * dt 
   x_jupiter = x_jupiter + vx_jupiter * dt 
   y_jupiter = y_jupiter + vy_jupiter * dt 
   t = t +dt
   write(2,'(6f14.7)')  t,x_earth,y_earth, x_jupiter,y_jupiter 
  end do
  end subroutine calculate

