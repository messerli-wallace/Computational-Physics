!================================================
   program kepler3
!================================================
!  Planetary motion (simplified circle) 
!------------------------------------------------
   implicit none
   double precision pi,twopi,pp, beta, ecc, v_0
   double precision a,x,y,vx,vy,r
   double precision dt,t,phi,phinew,phi_sum
   double precision ms,me, mp, mm,mpluto
   integer n, planet_switch, euler_switch

   pi = acos(-1.0d0)
   twopi = 2.0d0*pi
   ms = 2.0e30 !mass of the sun
   me = 6.0e24 !mass of the earth
   mm = 2.4e23 !mass of mercury
   mpluto = 1.3e22 !mass of pluto

! Read initial conditions.  We always start on x-axis with non-zero vy.
   open(1,file='2-body.in') 
   read(1,*) dt,n            !           time steps, number of orbits
   read(1,*) beta, ecc     !           beta and eccentricity
   read(1,*) v_0             !           initial vy
   read(1,*) euler_switch    !    switch for Euler (0) or Euler-Cromer (1)
   read(1,*) planet_switch   ! switch for planet 0=earth, 1=mercury, 2=pluto

   a=1.0d0
   mp = me
   vy=v_0
   x = 0.5d0 
   y = 0.0d0; vx = 0.0d0     !           We always start on the x-axisi, so vx = 0.0d0   

   !use pre-set planet values if we wish
   if (planet_switch.eq.0) then
      a=1.0d0 !earth
      ecc = 0.017d0
      mp = me
      x = a*(1+ecc)                     !     intitial x based off eccentricity
      vy = -2*pi*sqrt((1+ecc)/abs(x)/(1-ecc)*(1+(mp/ms))) !initial velocity
   else if (planet_switch.eq.1) then
      a=0.39d0 !mercury
      ecc = 0.206d0
      mp = mm
      x = a*(1+ecc)                     !     intitial x based off eccentricity
      vy = -2*pi*sqrt((1+ecc)/abs(x)/(1-ecc)*(1+(mp/ms))) !initial velocity
   else if (planet_switch.eq.2) then
      a=39.53d0 !pluto, perihelion
      ecc = 0.248d0
      mp = mpluto
      x = a*(1+ecc)                     !     intitial x based off eccentricity
      vy = -2*pi*sqrt((1+ecc)/abs(x)/(1-ecc)*(1+(mp/ms))) !initial velocity
   end if

   
   print*, vy

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
      write(2,1001) t,x,y,vx,vy,phi_sum
      r = sqrt(x**2 + y**2)
      vx = vx - (pp*x*dt)/r**(beta+1)
      vy = vy - (pp*y*dt)/r**(beta+1)
      x = x + vx * dt
      y = y + vy * dt
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
         write(2,1001) t,x,y,vx,vy,phi_sum
         r = sqrt(x**2 + y**2)
         x = x + vx * dt
         y = y + vy * dt
         vx = vx - (pp*x*dt)/r**(beta+1)
         vy = vy - (pp*y*dt)/r**(beta+1)
         phinew = atan2(abs(y),abs(x))
         phi_sum = phi_sum + abs(phinew-phi)
         if(phi_sum.ge.2*pi*n) stop
         phi = phinew 
         t = t + dt
         enddo
   end if
1000  format('#      t             x             y            &
           vx            vy            phi-total')
1001  format(1p10e14.6)

   end program kepler3

