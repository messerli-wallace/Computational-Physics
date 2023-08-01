program pendulum
implicit none
double precision th(4),om(4),t(4),Omega,k,q,pi, &
dt
integer nperiod, npp
common / constants / pi

print*,'initial angle, velocity, periods, points per period?'
read(5,*) th(1),om(1),nperiod,npp
t(1) = 0.0d0
print *,'k (damping const), q (force amp), Omega (force frequency)?'
read(5,*) k,q,Omega
pi = acos(-1.0d0)
dt = pi/FLOAT(npp)
print *, k,q,omega


open(1,file='data.out',status='unknown')
!1001 format('# damping, drive ampl., drive freq., dt :',1p4e14.6,/,&
!'# periods,points per period, total steps:',3i10)

write(1,1002)

1000    format(1p10ES16.6E3)
1002    format('#    time            phi         omega')

call calculate(th,om,t,dt,k,q,Omega,nperiod,npp)

end program pendulum


subroutine calculate(th,om,t,dt,k,q,Omega,nperiod,npp)
implicit none
double precision th(*),om(*),t(*),dt,k,q,Omega, &
dom(4), dth(4), pi
integer nperiod, npp, i, j
common / constants / pi

!4th-order Runge-Kutta
do i = 2,nperiod
    do j=1,npp
        call dv(om(1),th(1),t(1),k,q,Omega,dom(1),dth(1))
        om(2) = om(1)+0.5d0*dt*dom(1)
        th(2) = th(1)+0.5d0*dt*dth(1)
        t(2) = t(1)+0.5d0*dt
        call dv(om(2),th(2),t(2),k,q,Omega,dom(2),dth(2))
        t(3) = t(2)
        om(3) = om(1)+0.5d0*dt*dom(2)
        th(3) = th(1)+0.5d0*dt*dth(2)
        call dv(om(3),th(3),t(3),k,q,Omega,dom(3),dth(3))
        t(4) = t(1)+dt
        om(4) = om(1)+dt*dom(3)
        th(4) = th(1)+dt*dth(3)
        call dv(om(4),th(4),t(4),k,q,Omega,dom(4),dth(4))
        om(1)=om(1)+dt*(dom(1)+2.0d0*dom(2)+2.0d0*dom(3)+dom(4))/6.0d0
        th(1)=th(1)+dt*(dth(1)+2.0d0*dth(2)+2.0d0*dth(3)+dth(4))/6.0d0
        t(1) = t(1)+dt

        if (th(1).gt.pi) th(1) = th(1)-2.0d0*pi
        if (th(1).lt.-pi) th(1) = th(1)+2.0d0*pi
        
        if(j.eq.npp .and. i.gt.nperiod-64) write(1,1000) t(1),th(1),om(1) ! writing the values we need
    enddo
enddo
1000    format(1p10ES16.6E3)
!
End Subroutine calculate


Subroutine dv(om0,th0,t0,k,q,Omega,dom,dth)
implicit none
double precision om0,dth,k,Omega,t0,q,th0,dom,pi
common / constants / pi

dth = om0
dom = -k*om0 - (1.0d0+2.0d0*q*cos(Omega*t0))*sin(th0)

End Subroutine dv
