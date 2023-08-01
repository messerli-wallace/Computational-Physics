*                            
*  PROGRAM TO TEST DIAGONALISATION ROUTINE   
*      
      PROGRAM DITEST         
*      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                 
      PARAMETER (NDIM=20)    
      DIMENSION H(NDIM,NDIM),VECTOR(NDIM,NDIM)            
*      
      READ (5,9001) N        
C      
      DO 1 I=1,N             
         READ (5,*) (H(I,J),J=I,N)                     
    1 CONTINUE               
*      
      DO 2 J=1,N             
         DO 3 I=1,J          
            H(J,I)=H(I,J)    
    3    CONTINUE            
    2 CONTINUE               
      WRITE (6,9005)         
      DO  4 I=1,N            
         WRITE (6,9002) (H(I,J),J=1,N)                    
    4 CONTINUE               
*      
      CALL DIAG(H,N,0,VECTOR,NR)    
*      
      WRITE (6,9003)         
      WRITE (6,9002) (H(I,I),I=1,N)                       
      WRITE (6,9004)         
      DO 5 I=1,N             
         WRITE (6,9002) (VECTOR(I,J),J=1,N)               
    5 CONTINUE               
*      
 9001 FORMAT (12I5)          
 9002 FORMAT (5F14.7)        
 9003 FORMAT (/' *****     EIGENVALUES     *****',/)      
 9004 FORMAT (/' *****     EIGENVECTORS    *****',/)      
 9005 FORMAT (/' *****   ORIGINAL MATRIX   *****',/)      
*      
      STOP                   
      END                    
      SUBROUTINE DIAG(H,N,IEIG,U,NR)                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                 
*      
*  DIAG ... DIAGONALIZATION OF A REAL SYMMETRIC MATRIX BY 
*         THE JACOBI METHOD. 
*  N IS THE ORDER OF THE SUBMATRIX TO BE DIAGONALISED     
*  IEIG MUST BE SET UNEQUAL TO ZERO IF ONLY EIGENVALUES ARE                     
*         TO BE COMPUTED.    
*  IEIG MUST BE SET EQUAL TO ZERO IF EIGENVALUES AND EIGENVECTORS               
*         ARE TO BE COMPUTED. 
*  U IS THE UNITARY MATRIX USED FOR FORMATION OF THE EIGENVECTORS.              
*  NR IS USED AS AN INPUT PARAMETER.   DIAG DIAGONALIZES THE                    
*  NR-TH SUBMATRIX.     WHEN  NR  IS USED AS AN OUTPUT PARAMTER,                
*  NR IS THE NUMBER OF ROTATIONS.   
*  THE SUBROUTINE OPERATES ONLY ON THE ELEMENTS OF H THAT ARE TO THE            
*         RIGHT OF THE MAIN DIAGONAL.  THUS, ONLY A TRIANGULAR                  
*         SECTION NEED BE STORED IN THE ARRAY H.          
*      
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0)                     
      PARAMETER (RAPIN=7.450580596D-9,HDIN=1.0D38)        
      PARAMETER (NDIM=20)    
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
