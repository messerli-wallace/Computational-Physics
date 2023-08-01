!================================================
   program kepler3
!================================================
!  Planetary motion (simplified circle) 
!------------------------------------------------
   implicit none
   double precision pi,twopi,pp, beta, ecc, v_0, euler_switch
   double precision a,x,y,x_s,y_s,vx,vy,vx_s,vy_s,r
   double precision dt,t,phi,phinew,phi_sum
   double precision ms,me
   integer n

   pi = acos(-1.0d0)
   twopi = 2.0d0*pi
   ms = 2.0e30 !mass of the sun
   me = 6.0e24 !mass of the earth
! Read initial conditions.  We always start on x-axis with non-zero vy.
   open(1,file='2-body.in')
   a = 1.0d0                 !           This is for the Earth in A.U.
   x = a                     !           This is just for cosmetics
   y = 0.0d0; vx = 0.0d0     !           We always start on the x-axisi, so vx = 0.0d0   
   !vy = twopi/sqrt(a)        !           This is v_y in A.U.  
   read(1,*) dt,n            !           time steps, number of orbits
   read(1,*) beta, ecc     !           beta and eccentricity
   read(1,*) v_0             !           initial v_y
   read(1,*) euler_switch    !    switch for Euler (0) or Euler-Cromer (1)
   
   !assign circular orbit if we pass zero velocity
   if (v_0.eq.0.0d0) then
      vy=twopi/sqrt(a)
   else
      vy=v_0
   end if

   print*, vy

   x_s = 0.0d0               ! Initial conditions for the sun
   y_s = 0.0d0
   vx_s = 0.0d0
   vy_s = 0.0d0

! Calculations and output: 

   pp = 4*pi*pi  
   phi = 0.0d0
   phi_sum = 0.0

   open(2,file = 'earth.out',status='unknown')
   write(2,1000)
   t = 0.d0
   
   if (euler_switch.eq.1) then
      !Euler-Cromer
      do 
      write(2,1001) t,x,y,x_s,y_s,vx,vy,vx_s,vy_s,phi_sum
      !r = sqrt(x**2 + y**2)
      r = sqrt((x-x_s)**2 + (y-y_s)**2)
      vx = vx - (pp*(x-x_s)*dt)/r**(beta+1)
      vy = vy - (pp*(y-y_s)*dt)/r**(beta+1)
      vx_s = vx_s - ((pp*me/ms)*(x_s-x))/r**(beta+1)
      vy_s = vy_s - ((pp*me/ms)*(y_s-y))/r**(beta+1)
      x = x + vx * dt
      y = y + vy * dt
      x_s = x_s + vx_s*dt
      y_s = y_s + vy_s*dt
      phinew = atan2(abs(y),abs(x))
      phi_sum = phi_sum + abs(phinew-phi)
      if(phi_sum.ge.2*pi*n) stop
      phi = phinew 
      t = t + dt
      enddo
   end if
   if (euler_switch.eq.0) then
      !Euler method
      do 
         write(2,1001) t,x,y,x_s,y_s,vx,vy,vx_s,vy_s,phi_sum
         !r = sqrt(x**2 + y**2)
         r = sqrt((x-x_s)**2 + (y-y_s)**2)
         x = x + vx * dt
         y = y + vy * dt
         x_s = x_s + vx_s*dt
         y_s = y_s + vy_s*dt
         vx = vx - (pp*(x-x_s)*dt)/r**(beta+1)
         vy = vy - (pp*(y-y_s)*dt)/r**(beta+1)
         vx_s = vx_s - (pp*(me/ms)*(x_s-x))/r**(beta+1)
         vy_s = vy_s - (pp*(me/ms)*(y_s-y))/r**(beta+1)
         phinew = atan2(abs(y),abs(x))
         phi_sum = phi_sum + abs(phinew-phi)
         if(phi_sum.ge.2*pi*n) stop
         phi = phinew 
         t = t + dt
         enddo
   end if
1000  format('#      t             x             y            x_s            &
           y_s            vx            vy            vx_s            vy_s         phi-total')
1001  format(1p10e14.6)

   end program kepler3

