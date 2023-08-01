! **********************************************************
Program Diagonalization
  implicit double precision (a-h,o-z)
  integer  npts,ntotstates, counter, q
  Real*8  dx,norm,value,lang,Zn
  real*8, allocatable :: r1(:),v1(:),r(:),v(:),hamiltonianD(:),hamiltonianU(:)
  real*8, allocatable :: W(:),WORK(:),wf(:,:,:),wf1(:),Z(:,:)
  INTEGER LWORK,INFO,lmax,nmax,i,LDZ
  character*1024, allocatable :: state(:)
  
  !Variables:
  !---------
  ! npts  = number of points
  ! dx    = stepsize
  ! nmax  = maximum principal number
  ! lmax  = maximum angular momentum
  
  print *,'input number of points and stepsize'
  read(5,*) npts,dx
  print *,'input maximum n and l'
  read(5,*) nmax,lmax
  ntotstates = (lmax+1)*nmax - lmax*(lmax+1)/2
  
  ! The next two statements are for the LAPACK routine
  LWORK = 2*npts-2
  LDZ = npts
  
  allocate(r1(npts1+1),v1(npts1+1),r(npts+1),v(npts+1),Z(npts,npts),hamiltonianD(npts),hamiltonianU(npts-1))
  allocate(W(npts+1), WORK(LWORK),wf(npts+1,nmax,0:lmax),wf1(npts+1))
  allocate(state(ntotstates))

  write(6,1000) npts,dx,npts*dx
  write(6,1001) nmax,lmax
  1000 format('npts = ',i6,'     dx = ',f10.6,'     rmax = ',f12.4)
  1001 format('nmax = ',i6,'     lmax = ',2i5)
  write(6,*)
  
  !Set up the grid and set potential
  do i=1,npts+1 !revisit, maybe 2npts+1
    ! counter = i - (npts/2)
    ! r(i)=(i-1 - (npts/2))*dx ! center it at npts/2
    r(i) = (i-1)*dx
    ! v(i) = 1/(r(i)-(npts/2)*dx)        !coulomb
    v(i) = -1.0d0/(sqrt(1+(r(i)-((npts*dx)/2))**2)) ! project 9
    ! write(1,*) r(i), v(i)
  enddo
  
  ! Start diagonalization of hamiltonian
  print*,'****************************************'
  print*,'Setting up and Diagonalizing Hamiltonian'
  print*,'****************************************'
  print*,'---------------'
  print*,'Energies (a.u.)'
  print*,'---------------'
  Do il = 0,lmax
   Z = 0.d0
   Do i=1,npts
    Z(i,i) = 1.0D0
   enddo
   lang = dble(il)
   hamiltonianU = 0.d0
   hamiltonianD = 0.d0
   Do i=1,npts-1
     hamiltonianD(i) =  1.d0/dx**2
     hamiltonianD(i) = hamiltonianD(i) + lang * (lang + 1.d0)/2.d0/r(i+1)**2 + v(i+1)
     if (i.lt.(npts-1)) hamiltonianU(i) = -0.5d0/dx**2
   Enddo
   CALL DSTEQR('V',npts-1,hamiltonianD,hamiltonianU,Z,LDZ,WORK,INFO)
   write(*,1002) il
   write(*,1003) (hamiltonianD(j),j=1,nmax-il)
   write(1,1002) il
   write(1,1003) (hamiltonianD(j),j=1,nmax-il)
  1002 format(/,'energies for l = ',i3,':',/)
  1003 format(1p10e14.6)
   write(*,*)
  ! Renormalize with Simpson rule 
   Do j=1,nmax-il
    wf1(:) = abs(Z(:,j))**2
    wf(:,j,il) = Z(:,j)
    call arsimd(npts,dx,wf1,norm)
    wf(:,j,il) = wf(:,j,il)/dsqrt(norm)
   Enddo
  Enddo
  
  ! Write output file with the states in order [s-state from n= 1,...,nmax, then p-state from n= 1,...,nmax-1,...]
  open(9,file='wfn.out')
  write(9,999) 0.d0,0.0d0,(0.d0,j=1,ntotstates)
  999 format(1p600e14.6)
  ! Ensure that all wavefunctions start with positive values

  do j=1,nmax-il+1
   do il=0,lmax
    if (wf(2,j,il).le.0.0d0) then
     do i=1,npts-1
      wf(i,j,il) = -wf(i,j,il)
     enddo
     wf(npts,j,il) = 0.0d0
    endif
   enddo
  enddo
  !


  !trevin stuff
  ! Do q=1,ntotstates
  !   write(state(q),"(i3.3,a)") q+2
  ! Enddo

  ! write(9,*) "#       1               2              ",(trim(state(q)),"             ",q=1,ntotstates)

  ! q = 0
  ! Do il=0,lmax
  !   Do j=il+1,nmax
  !     q = q + 1
  !     write(state(q),"(a,i2.2,a,i2.2)") " n = ",j," l = ",il
  !   Enddo
  ! Enddo

  ! write(9,*) "#    position       potential    ",(trim(state(q)),"  ",q=1,ntotstates)
  !End trevin stuff


  Do i=1,npts
  write(9,999) r(i+1),v(i+1),((wf(i,j,il),j=1,nmax-il),il=0,lmax)
  Enddo
  



  End Program
  
  SUBROUTINE ARSIMD(N,DEL,A,R)
  IMPLICIT real*8 (A-H,O-Z)
  DIMENSION A(N)
  L=N
  SUM=A(1)-A(L)
  DO I=2,L,2
    SUM=SUM+4.D0*A(I)+2.D0*A(I+1)
  ENDDO
  R=(DEL*SUM)/3.D0
  RETURN
  END
  