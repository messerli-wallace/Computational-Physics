      program potential
      implicit double precision(a-h,o-z)
      dimension pot0(0:100),pot1(0:100)
*
*  initialize
*
      do 10 i=1,99
        pot0(i) = 0.5d0
10    continue
      pot0(0) = 1.0d0
      pot0(100) = 0.0d0
40    continue
*
*   iterate
*
      h = 0.1d0
      niter = 0
      nitmax = 10000
70    niter = niter+1
100   do 110 i=1,99
        pot1(i) = 0.5d0*(pot0(i-1)+pot0(i+1))
110   continue
*
*    get the change and transfer values
*
      error = 0.0d0
      do 200 i=1,99
        error = error+abs(pot0(i)-pot1(i))
        pot0(i) = pot1(i)
200   continue
      if (error.gt.1.0d-5.and.niter.le.nitmax) goto 70
*
*   print final result to file
*
      write(6,1000) niter,error
      do 310 i=1,99
        write(1,1001) h*i,pot1(i)
310   continue
*
1000  format('niter = ',i5,'   error = ',1pe14.6)
1001  format(1p10e14.6)
*
      stop
      end
