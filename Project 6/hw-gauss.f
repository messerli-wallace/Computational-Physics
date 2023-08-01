      program potential
      implicit double precision(a-h,o-z)
      dimension pot0(0:400,0:400)
      dimension dx(0:400,0:400)
      dimension dy(0:400,0:400)
*
*  initialize
*
      h = 0.1d0
      niter = 0
      nitmax = 100000
      error_old=0.0d0
      error=0.0d0
      do 13 i=0,400
            do 12 j=0,400
                  pot0(i,j) = 0.0d0
12    continue
13    continue

*
* start
*
70    niter = niter+1
      error_old = error
      do 10 i=160,240
            do 11 j=180,220
*                 pot0(i,j) = 0.0d0
                  if (j.eq.180) pot0(i,j) = -1.0d0
                  if (j.eq.220) pot0(i,j) = 1.0d0
11    continue
10    continue


      error = 0.0d0
100   do 110 i=1,399
       do 120 j=1,399
         sum = pot0(i-1,j)+pot0(i+1,j)
     >        +pot0(i,j-1)+pot0(i,j+1)

         old = pot0(i,j)
         pot0(i,j) = sum/4.0d0
         error = error+abs(pot0(i,j)-old)
120    continue
110   continue
*
*      write(6,1000) niter, error
      write(2,*) niter,error
      if (abs(error_old-error).gt.1.0d-5.and.niter.le.nitmax) goto 70
*
*   print final result to file
*
      do 15 i=160,240
            do 14 j=180,220
*                 pot0(i,j) = 0.0d0
                  if (j.eq.180) pot0(i,j) = -1.0d0
                  if (j.eq.220) pot0(i,j) = 1.0d0
14    continue
15    continue
      do 310 i=0,400
       do 320 j=0,400
         write(1,1001) h*i,h*j,pot0(i,j)
       write(1,*)
320    continue
310   continue
      do 400 i=0,400
       write(3,1001) h*i,pot0(i,0)
400   continue
*
1000  format('niter = ',i5,'   error = ',1pe14.6)
1001  format(1p10e14.6)
*
      stop
      end
