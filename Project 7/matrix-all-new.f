      program matdemo
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)
      PARAMETER (IDMKTM = 2000)
      
      dimension a(idmktm,idmktm),b(idmktm,idmktm),c(idmktm,idmktm)
      dimension d(idmktm,idmktm),e(idmktm,idmktm)
      
      print *,'input matrix dimension'
      read(5,*) ndim
      if (ndim.gt.2000) then
       print *,'ndim = ',ndim
       print *,'increase ndim'
       stop
      endif

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
      if (ndim.le.10) then
      do i=1,ndim
       write(6,1000) (a(i,j),j=1,ndim)
      enddo
      endif
      print *
      print *,'Matrix B:'
      if (ndim.le.10) then
      do i=1,ndim
       write(6,1000) (b(i,j),j=1,ndim)
      enddo
      endif
      
      call matmup(a,b,c,ndim,ndim,ndim)
      print *
      print *,'Matrix C = A x B:'
      if (ndim.le.10) then
      do i=1,ndim
       write(6,1000) (c(i,j),j=1,ndim)
      enddo
      endif
        
      do i=1,ndim
       do j=1,ndim
        d(i,j) = c(i,j)
       enddo
      enddo
      call matinv(d,ndim,det)
      print *
      print *,'Inverse of Matrix C:'
      if (ndim.le.10) then
      do i=1,ndim
       write(6,1000) (d(i,j),j=1,ndim)
      enddo
      endif
      
      call matmup(c,d,e,ndim,ndim,ndim)
      print *
      print *,'Check C * C^{-1}:'
      if (ndim.le.10) then
      do i=1,ndim
       write(6,1000) (e(i,j),j=1,ndim)
      enddo
      endif
      
      call cpu_time(start)
      CALL DIAG(A,NDIM,0,D,NR)     
      WRITE (6,9003)         
      WRITE (6,9002) (A(I,I),I=1,NDIM)                       
      WRITE (6,9004)         
      DO 5 I=max(5,NDIM-4),NDIM   
         WRITE (6,9002) (D(I,J),J=1,NDIM)               
    5 CONTINUE            
      call cpu_time(finish)
      print '("Time = ",f6.3," seconds.")',finish-start       
*           
 9002 FORMAT (5F14.7)        
 9003 FORMAT (/' *****     EIGENVALUES     *****',/)      
 9004 FORMAT (/' *****     EIGENVECTORS    *****',/)      
      
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
      PARAMETER (IDMKTM = 2000)
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
      PARAMETER (IDMKTM = 2000)
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
      
      SUBROUTINE DIAG(H,N,IEIG,U,NR)                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                 
*      
*  DIAG ... DIAGONALIZATION OF A REAL SYMMETRIC MATRIX BY 
*           THE JACOBI METHOD. 
*  N        ORDER OF THE SUBMATRIX TO BE DIAGONALISED  
*  IEIG     SET TO ZERO FOR EIGENVALUES AND EIGENVECTORS   
*  IEIG     NOT ZERO IF ONLY EIGENVALUES ARE NEEDED                                
*  U IS THE UNITARY MATRIX USED FOR FORMATION OF THE EIGENVECTORS.              
*  NR IS USED AS AN INPUT PARAMETER. DIAG DIAGONALIZES THE NR-TH SUBMATRIX.     
*  ON OUTPUT, NR IS THE NUMBER OF ROTATIONS.   
*  THE SUBROUTINE OPERATES ONLY ON THE ELEMENTS OF H THAT ARE TO THE            
*  RIGHT OF THE MAIN DIAGONAL.  THUS, ONLY A TRIANGULAR                  
*  SECTION NEED BE STORED IN THE ARRAY H.          
*      
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0)                     
      PARAMETER (RAPIN=7.450580596D-9,HDIN=1.0D38)        
      PARAMETER (NDIM=2000)    
*      
      DIMENSION H(NDIM,NDIM),U(NDIM,NDIM),X(NDIM),IQ(NDIM)                      
*      
      IF (IEIG.NE.0) GOTO 4  
      DO 3 I=1,N             
         DO 2 J=1,N          
            IF (I.EQ.J) GOTO 1      
            U(I,J)=ZERO      
            U(J,I)=ZERO      
            GOTO 2           
1           U(I,J)=ONE       
2        CONTINUE            
3     CONTINUE               
*      
4     NR=0                   
      IF (N.LE.1) GOTO 30    
*      
*  SCAN FOR LARGEST OFF DIAGONAL ELEMENT IN EACH ROW      
*  X(I) CONTAINS LARGEST ELEMENT IN ITH ROW               
*  IQ(I) HOLDS SECOND SUBSCRIPT DEFINING POSITION OF ELEMENT                    
*      
      DO 5 I=1,N-1           
         X(I)=ZERO           
         IPL1=I+1            
         DO 5 J=IPL1,N       
            IF (X(I).GT.ABS(H(I,J))) GOTO 5               
            X(I)=ABS(H(I,J)) 
            IQ(I)=J          
5     CONTINUE               
*      
*  SET INDICATOR FOR SHUT-OFF.RAP=2**-27,NR=NO. OF ROTATIONS                    
*      
      RAP=RAPIN              
      HDTEST=HDIN            
*      
*  FIND MAXIMUM OF X(I)'S FOR PIVOT ELEMENT AND           
*  TEST FOR END OF PROBLEM   
*      
6     DO 8 I=1,N-1           
         IF (I.LE.1) GOTO 7  
         IF (XMAX.GE.X(I)) GOTO 8   
7        XMAX=X(I)           
         IPIV=I              
         JPIV=IQ(I)          
8     CONTINUE               
*      
*  IS MAX(X(I)) EQUAL TO ZERO, IF LESS THAN HDTEST, REVISE HDTEST               
*      
      IF (XMAX.LE.ZERO) GOTO 30     
      IF (HDTEST.LE.ZERO) GOTO 9    
      IF (XMAX.GT.HDTEST) GOTO 11   
9     HDIMIN=ABS(H(1,1))     
      DO 10 I=2,N            
         IF (HDIMIN.LE.ABS(H(I,I))) GOTO 10               
         HDIMIN=ABS(H(I,I))  
10    CONTINUE               
*      
      HDTEST=HDIMIN*RAP      
*      
*  RETURN IF MAX.H(I,J)LESS THAN(2**-27)ABSF(H(K,K)-MIN)  
*      
      IF (HDTEST.GE.XMAX) GOTO 30   
11    NR=NR+1                
*      
*  COMPUTE TANGENT, SINE AND COSINE,H(I,I),H(J,J)         
*      
12    TANG=SIGN(TWO,(H(IPIV,IPIV)-H(JPIV,JPIV)))*H(IPIV,JPIV)/ (ABS(H(IP        
     -IV,IPIV)-H(JPIV,JPIV))+SQRT((H(IPIV,IPIV) -H(JPIV,JPIV))**2               
     -+FOUR*H(IPIV,JPIV)**2)) 
      COSINE=ONE/SQRT(ONE+TANG**2)  
      SINE=TANG*COSINE       
      HII=H(IPIV,IPIV)       
      H(IPIV,IPIV)=COSINE**2*(HII+TANG*(TWO*H(IPIV,JPIV)+TANG* H(JPIV,          
     -JPIV)))                
      H(JPIV,JPIV)=COSINE**2*(H(JPIV,JPIV)-TANG*(TWO*H(IPIV,JPIV)               
     -- TANG*HII))           
      H(IPIV,JPIV)=ZERO      
*      
*   PSEUDO RANK THE EIGENVALUES     
*   ADJUST SINE AND COS FOR COMPUTATION OF H(IK) AND U(IK)                      
*      
      IF (H(IPIV,IPIV).GE.H(JPIV,JPIV)) GOTO 13           
      HTEMP=H(IPIV,IPIV)     
      H(IPIV,IPIV)=H(JPIV,JPIV)     
      H(JPIV,JPIV)=HTEMP     
*    RECOMPUTE SINE AND COS  
      HTEMP=SIGN(ONE,-SINE)*COSINE  
      COSINE=ABS(SINE)       
      SINE=HTEMP             
13    CONTINUE               
*      
*  INSPECT THE IQS BETWEEN I+1 AND N-1 TO DETERMINE       
*  WHETHER A NEW MAXIMUM VALUE SHOULD BE COMPUTED SINCE   
*  THE PRESENT MAXIMUM IS IN THE I OR J ROW.              
*      
      DO 19 I=1,N-1          
         IF (I-IPIV) 15,19,14 
14       IF (I.EQ.JPIV) GOTO 19     
15       IF (IPIV.EQ.IQ(I)) GOTO 16                       
         IF (JPIV.NE.IQ(I)) GOTO 19                       
16       K=IQ(I)             
17       HTEMP=H(I,K)        
         H(I,K)=ZERO         
         IPL1=I+1            
         X(I)=ZERO           
*      
*  SEARCH IN DEPLETED ROW FOR NEW MAXIMUM                 
*      
         DO 18 J=IPL1,N      
            IF (X(I).GT.ABS(H(I,J))) GOTO 18              
            X(I)=ABS(H(I,J)) 
            IQ(I)=J          
18       CONTINUE            
         H(I,K)=HTEMP        
19    CONTINUE               
*      
      X(IPIV)=ZERO           
      X(JPIV)=ZERO           
*      
*  CHANGE THE OTHER ELEMENTS OF H   
*      
      DO 28 I=1,N            
         IF (I-IPIV) 20,28,23 
20       HTEMP=H(I,IPIV)     
         H(I,IPIV)=COSINE*HTEMP+SINE*H(I,JPIV)            
         IF (X(I).GE.ABS(H(I,IPIV))) GOTO 21              
         X(I)=ABS (H(I,IPIV)) 
         IQ(I)=IPIV          
21       H(I,JPIV)=-SINE*HTEMP+COSINE*H(I,JPIV)           
         IF (X(I).GE.ABS(H(I,JPIV))) GOTO 28              
22       X(I)=ABS (H(I,JPIV)) 
         IQ(I)=JPIV          
         GOTO 28             
*      
23       IF (I-JPIV) 24,28,26 
24       HTEMP=H(IPIV,I)     
         H(IPIV,I)=COSINE*HTEMP+SINE*H(I,JPIV)            
         IF (X(IPIV).GE.ABS(H(IPIV,I))) GOTO 25           
         X(IPIV)=ABS(H(IPIV,I))     
         IQ(IPIV)=I          
25       H(I,JPIV)=-SINE*HTEMP+COSINE*H(I,JPIV)           
         IF(X(I)-ABS(H(I,JPIV))) 22,28,28                 
*      
26       HTEMP=H(IPIV,I)     
         H(IPIV,I)=COSINE*HTEMP+SINE*H(JPIV,I)            
         IF (X(IPIV).GE.ABS(H(IPIV,I))) GOTO 27           
         X(IPIV)=ABS(H(IPIV,I))     
         IQ(IPIV)=I          
27       H(JPIV,I)=-SINE*HTEMP+COSINE*H(JPIV,I)           
         IF (X(JPIV).GE.ABS(H(JPIV,I))) GOTO 28           
         X(JPIV)=ABS(H(JPIV,I))     
         IQ(JPIV)=I          
28    CONTINUE               
*      
*  TEST FOR COMPUTATION OF EIGENVECTORS                   
*      
      IF (IEIG.NE.0) GOTO 6  
      DO 29 I=1,N            
         HTEMP=U(I,IPIV)     
         U(I,IPIV)=COSINE*HTEMP+SINE*U(I,JPIV)            
         U(I,JPIV)=-SINE*HTEMP+COSINE*U(I,JPIV)           
29    CONTINUE               
      GOTO 6                 
30    CONTINUE               
      RETURN                 
      END                    

