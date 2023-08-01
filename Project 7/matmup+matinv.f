      program matdemo
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)
      PARAMETER (IDMKTM = 410)
      
      dimension a(idmktm,idmktm),b(idmktm,idmktm),c(idmktm,idmktm)
      dimension d(idmktm,idmktm),e(idmktm,idmktm)
      
      print *,'input matrix dimension'
      read(5,*) ndim
      do i=1,ndim
       do j=i,ndim
        a(i,j) = sqrt(dble(i))*log(dble(j))
        a(j,i) = a(i,j)
        b(i,j) = sin(dble(i))*cos(dble(j))
        b(j,i) = b(i,j)
       enddo
      enddo
      
      print *
      print *,'Matrix A:'
      do i=1,ndim
       write(6,1000) (a(i,j),j=1,ndim)
      enddo
      print *
      print *,'Matrix B:'
      do i=1,ndim
       write(6,1000) (b(i,j),j=1,ndim)
      enddo
      
      call matmup(a,b,c,ndim,ndim,ndim)
      print *
      print *,'Matrix C = A x B:'
      do i=1,ndim
       write(6,1000) (c(i,j),j=1,ndim)
      enddo
        
      do i=1,ndim
       do j=1,ndim
        d(i,j) = c(i,j)
       enddo
      enddo
      call matinv(d,ndim,det)
      print *
      print *,'Inverse of Matrix C:'
      do i=1,ndim
       write(6,1000) (d(i,j),j=1,ndim)
      enddo
      
      call matmup(c,d,e,ndim,ndim,ndim)
      print *
      print *,'Check C * C^{-1}:'
      do i=1,ndim
       write(6,1000) (e(i,j),j=1,ndim)
      enddo
      
1000  format(1p10e13.5)
      stop
      end
      
************************************************************************
*                                                                      *
*   THE FOLLOWING SUBROUTINE CALCULATES THE PRODUCT MATRIX AMULT       *
*   OF THE MATRICES AMATL (LEFT) AND AMATR (RIGHT).                    *
*                                                                      *
*   IROWL  : NUMBER OF ROWS OF AMATL                                   *
*   ICOROW : NUMBER OF COLUMNS OF AMATL (OF ROWS OF AMATR)             *
*   AMULT IS A MATRIX WITH IROWL ROWS AND ICOLR COLUMNS.               *
*                                                                      *
************************************************************************
      SUBROUTINE MATMUP (AMATL,AMATR,AMULT,IROWL,ICOROW,ICOLR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER (IDMKTM = 410)
*
      DIMENSION AMATL(IDMKTM,IDMKTM),AMATR(IDMKTM,IDMKTM),
     +          AMULT(IDMKTM,IDMKTM)
*
      DO 4 I=1,IROWL
       DO 1 J=1,ICOLR
        AMULT(I,J) = 0.0D0
1      CONTINUE
       DO 3 K=1,ICOROW
        DO 2 J=1,ICOLR
         AMULT(I,J) = AMULT(I,J) + AMATL(I,K)*AMATR(K,J)
2       CONTINUE
3      CONTINUE
4     CONTINUE
*
      RETURN
      END
           
************************************************************************
*                                                                      *
*   THE FOLLOWING SUBROUTINE INVERTS A MATRIX USING GAUSS-JORDAN       *
*   ELIMINATION WITH COLUMN SHIFTING TO MAXIMIZE PIVOT-ELEMENTS.       *
*                                                                      *
*                     INPUT                                            *
*                     =====                                            *
*                                                                      *
*     A               N X N MATRIX TO BE INVERTED                      *
*                     (A IS OVERWRITTEN ON OUTPUT)                     *
*     N               DIMENSION OF ROW OR COLUMN OF A                  *
*                                                                      *
*                     OUTPUT                                           *
*                     ======                                           *
*                                                                      *
*     A               MATRIX INVERSE                                   *
*     DETA            THE ABSOLUTE VALUE OF THE DETERMINANT OF         *
*                     THE INPUT MATRIX A                               *
*                                                                      *
************************************************************************
      SUBROUTINE MATINV(A,N,DETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)
      PARAMETER (IDMKTM = 410)
*
      DIMENSION A(IDMKTM,IDMKTM),ICOL(IDMKTM)
*
*   ICOL IS USED TO RECORD COLUMN SHIFTS TO ENABLE
*   UNSCRAMBLING OF THE SOLUTION AT THE END.
*
*   INITIALIZE DETERMINANT AND ICOL.
*
      DETA = ONE
      DO 1 I=1,N
       ICOL(I) = I
1     CONTINUE
      DO 8 I=1,N
*
*   FIND LARGEST ELEMENT IN ROW I IN A NON-PIVOTED COLUMN
*
       AA = ZERO
       AB = ZERO
       J=I
       DO 2 K=I,N
        IF ((AB-ABS(A(I,K))).GE.ZERO) GOTO 2
        J = K
        AA = A(I,K)
        AB = ABS(AA)
2      CONTINUE
       IF (I.EQ.J) GOTO 4
*
*   SHIFT COLUMNS
*
       M = ICOL(J)
       ICOL(J) = ICOL(I)
       ICOL(I) = M
       DO 3 M=1,N
        B = A(M,I)
        A(M,I) = A(M,J)
        A(M,J) = B
3      CONTINUE
4      A(I,I) = ONE
       DETA = DETA*AA
       DO 5 J=1,N
        A(I,J) = A(I,J)/AA
5      CONTINUE
*
*   PERFORM ELIMINATION
*
       DO 7 J=1,N
        IF (I.EQ.J) GOTO 7
        AD = A(J,I)
        IF (A(J,I).EQ.ZERO) GOTO 7
        A(J,I) = ZERO
        DO 6 K=1,N
         AE = A(J,K)
         AF = A(I,K)
         AE = AE - AD*AF
         A(J,K) = AE
6      CONTINUE
7      CONTINUE
8     CONTINUE
*
*   UNSCRAMBLE THE SOLUTION.
*
      DO 12 I=1,N
       IF (ICOL(I).EQ.I) GOTO 12
       J = I
9      J = J+1
       IF (ICOL(J).EQ.I) GOTO 10
       IF (N.GT.J) GOTO 9
10     ICOL(J) = ICOL(I)
       DO 11 K=1,N
        B=A(I,K)
        A(I,K) = A(J,K)
        A(J,K) = B
11     CONTINUE
       ICOL(I) = I
12    CONTINUE
      DETA = ABS(DETA)
*
      RETURN
      END
