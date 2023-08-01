pi 

read x1

y1 = 0.0; vx1 = 0.0
vy1=2*pi/sqrt(x1)
read x2
y2=0.0;vx2=0.0
vy2=2*pi/sqrt(x2)

read m1
read m2
read m3
read dtread nt
close(1)

open(2,file-'threebody-final.out')
write(2,6(...)) 'stuff'

call calculate


calculate

!use Euler-Chromer
do i=1,nt
    r=sqrt()
    rj =sqrt()
    rej = sqrt()
    vx_earth = vx_earth - pp*x_earth

    !do for other two

    t=t+dt
    write 


