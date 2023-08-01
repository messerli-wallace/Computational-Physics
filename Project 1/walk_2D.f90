program randomwalk2D
    IMPLICIT NONE
    integer :: nsteps
    integer,parameter :: nwalkers = 1000
    integer k,T(8),seed,x_pos(nwalkers),y_pos(nwalkers),i,iout,inp
    real*8 x,y
    real*8 rsum,r(nwalkers)
    !
    ! Initiate the random-number generator
    !
    CALL DATE_AND_TIME(VALUES = T)
    SEED = T(1)+T(2)+T(3)+T(4)+T(5)+T(6)+T(7)+T(8)
    call srand(seed)
    !
    ! This is just to show how to set an output unit.  Never use 5, only use 6 if you want to come to the screen.
    !
    iout = 100
    !
    10     print *,'Enter 0 to stop the program.'
    read(*,*) inp
    if (inp.eq.0) stop 'done'
    x_pos=0
    y_pos=0
    !
    ! Files for step counts
    open(unit=11,file='step.100',status='unknown')
    open(unit=12,file='step.200',status='unknown')
    open(unit=13,file='step.500',status='unknown')
    open(unit=14,file='step.1000',status='unknown')
    open(unit=15,file='step.5000',status='unknown')
    open(unit=16,file='step.10000',status='unknown')

    !set nsteps to run the experiment
    nsteps = 10000

    1001 format(10i8)
    !
    Do i=1,nwalkers
      Do k=1,nsteps
        x=rand()
        if(x>0.5d0) x_pos(i) = x_pos(i) + 1
        if(x<0.5d0) x_pos(i) = x_pos(i) - 1
        y=rand()
        if(y>0.5d0) y_pos(i) = y_pos(i) + 1
        if(y<0.5d0) y_pos(i) = y_pos(i) - 1
  
      ! Write the positions to the files at various checkpoints
        if(k.eq.100) write(11,1001) i,x_pos(i),y_pos(i)
        if(k.eq.200) write(12,1001) i,x_pos(i),y_pos(i)
        if(k.eq.500) write(13,1001) i,x_pos(i),y_pos(i)
        if(k.eq.1000) write(14,1001) i,x_pos(i),y_pos(i)
        if(k.eq.5000) write(15,1001) i,x_pos(i),y_pos(i)
        if(k.eq.10000) write(16,1001) i,x_pos(i),y_pos(i)
      Enddo
      r(i) = sqrt((x_pos(i)**2)*1.0 + (y_pos(i)**2)*1.0)
      !
      ! fort.101 will contain the last position of our walkers
      ! Not needed because we have separate files for each checkpoint.
      !write(iout,*) i,x_pos(i),y_pos(i)
    Enddo
    !
    ! Now do some analysis of the results
    ! We sum up the absolute values of the last point for each walker and then divide by the number of walkers.
    !
    rsum = 0.0d0
    do i=1,nwalkers
      rsum = rsum+r(i)
    enddo
    rsum = rsum/dble(nwalkers)
    write(6,1000) rsum
    1000 format('average distance = ',f10.5)
    !
    ! Now go back and do it again
    !
    goto 10
      
  end program randomwalk2D