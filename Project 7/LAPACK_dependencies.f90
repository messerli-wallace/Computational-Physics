program matrix

end program

!==================================================================================================
    SUBROUTINE MATMUP (MAT_L,MAT_R,MULT,ROW_L,COLROW,COL_R,NDIM)
!==================================================================================================
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                      !
        !   THE FOLLOWING SUBROUTINE CALCULATES THE PRODUCT MATRIX MULT        !
        !   OF THE MATRICES MAT_L (LEFT) AND MAT_R (RIGHT).                    !
        !                                                                      !
        !   ROW_L  : NUMBER OF ROWS OF MAT_L                                   !
        !   COLROW : NUMBER OF COLUMNS OF MAT_L (OF ROWS OF MAT_R)             !
        !   MULT IS A MATRIX WITH ROW_L ROWS AND COL_R COLUMNS.                !
        !                                                                      !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IMPLICIT NONE
    !        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        real*8,allocatable,dimension(:,:) :: MAT_L,MAT_R,MULT
        integer I,J,K,ROW_L,COLROW,COL_R,NDIM
    !
        MULT = 0.0d0
        DO I=1,ROW_L
            DO K=1,COLROW
                DO J=1,COL_R
                    MULT(I,J) = MULT(I,J) + MAT_L(I,K)*MAT_R(K,J)
                ENDDO
            ENDDO
        ENDDO
    !
        RETURN
    END SUBROUTINE


!==================================================================================================
    SUBROUTINE dgetri(N,A,LDA,IPIV,WORK,LWORK,INFO)
!==================================================================================================
    !
    !  -- LAPACK computational routine --
    !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    !
    !     .. Scalar Arguments ..
        INTEGER            INFO,LDA,LWORK,N
    !     ..
    !     .. Array Arguments ..
        INTEGER            IPIV(*)
        DOUBLE PRECISION   A(LDA,*), WORK(*)
    !     ..
    !
    !============================================================================
    !
    !     .. Parameters ..
        DOUBLE PRECISION    ZERO,ONE
        parameter(ZERO = 0.0d0,ONE = 1.0d0 )
    !     ..
    !     .. Local Scalars ..
        LOGICAL            LQUERY
        INTEGER            I,IWS,J,JB,JJ,JP,LDWORK,LWKOPT,NB,NBMIN,NN
    !     ..
    !     .. External Functions ..
        INTEGER            ILAENV
        EXTERNAL           ilaenv
    !     ..
    !     .. External Subroutines ..
        EXTERNAL           dgemm,dgemv,dswap,dtrsm,dtrtri,xerbla
    !     ..
    !     .. Intrinsic Functions ..
        INTRINSIC          max,min
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
        info    = 0
        nb      = ilaenv(1,'DGETRI',' ',n,-1,-1,-1)
        lwkopt  = n*nb
        work(1) = lwkopt
        lquery  = (lwork == -1)
        IF(n < 0) THEN
            info = -1
        ELSE IF(lda < max(1,n)) THEN
            info = -3
        ELSE IF(lwork < max(1,n) .AND. .NOT.lquery) THEN
            info = -6
        END IF
        IF(info /= 0) THEN
            CALL xerbla('DGETRI',-info)
            RETURN
        ELSE IF(lquery) THEN
            RETURN
        END IF
    !
    !     Quick return if possible
    !
        IF(n == 0) RETURN
    !
    !     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
    !     and the inverse is not computed.
    !
        CALL dtrtri('Upper','Non-unit',n,a,lda,info)
        IF(info > 0) RETURN
    !
        nbmin  = 2
        ldwork = n
        IF(nb > 1 .AND. nb < n) THEN
            iws = max(ldwork*nb,1)
            IF(lwork < iws) THEN
                nb    = lwork/ldwork
                nbmin = max(2,ilaenv(2,'DGETRI',' ',n,-1,-1,-1))
            END IF
        ELSE
            iws = n
        END IF
    !
    !     Solve the equation inv(A)*L = inv(U) for inv(A).
    !
        IF(nb < nbmin .OR. nb >= n) THEN
    !
    !        Use unblocked code.
    !
            DO j=n,1,-1
    !
    !           Copy current column of L to WORK and replace with zeros.
    !
                DO i=j+1,n
                    work(i) = a(i,j)
                    a(i,j)  = ZERO
                ENDDO
    !
    !           Compute current column of inv(A).
    !
                IF(j < n) CALL dgemv('No transpose',n,n-j,-ONE,a(1,j+1),lda,work(j+1),1,ONE,a(1,j),1)
            ENDDO
        ELSE
    !
    !        Use blocked code.
    !
            nn = ((n-1)/nb)*nb + 1
            DO j=nn,1,-nb
                jb = min(nb,n-j+1)
    !
    !           Copy current block column of L to WORK and replace with
    !           zeros.
    !
                DO jj=j,j+jb-1
                    DO i=jj+1,n
                        work(i+(jj-j)*ldwork) = a(i,jj)
                        a(i,jj) = ZERO
                    ENDDO
                ENDDO
    !
    !           Compute current block column of inv(A).
    !
                IF(j+jb <= n) CALL dgemm( 'No transpose','No transpose',n,jb,n-j-jb+1,-ONE,a(1,j+jb),lda,work(j+jb),&
                                        ldwork,ONE,a(1,j),lda)
                CALL dtrsm('Right','Lower','No transpose','Unit',n,jb,ONE,work(j),ldwork,a(1,j),lda)
            ENDDO
        END IF
    !
    !     Apply column interchanges.
    !
        DO j=n-1,1,-1
            jp = IPIV(j)
            IF(jp /= j) CALL dswap(n,a(1,j),1,a(1,jp),1)
        ENDDO
    !
        work(1) = iws
        RETURN
    !
    !     End of DGETRI
    !
    END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE dsyev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
    !==================================================================================================
        !
        !  -- LAPACK driver routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          JOBZ,UPLO
            INTEGER            INFO,LDA,LWORK,N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A(LDA,*),W(*),WORK(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ZERO,ONE
            parameter(ZERO = 0.0d0,ONE = 1.0d0)
        !     ..
        !     .. Local Scalars ..
            LOGICAL            LOWER,LQUERY,WANTZ
            INTEGER            IINFO,IMAX,INDE,INDTAU,INDWRK,ISCALE,LLWORK,LWKOPT,NB
            DOUBLE PRECISION   ANRM,BIGNUM,EPS,RMAX,RMIN,SAFMIN,SIGMA,SMLNUM
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            INTEGER            ILAENV
            DOUBLE PRECISION   DLAMCH,DLANSY
            EXTERNAL           lsame,ilaenv,dlamch,dlansy
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           dlascl,dorgtr,dscal,dsteqr,dsterf,dsytrd,xerbla
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          max,sqrt
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters.
        !
            wantz  = lsame(jobz,'V')
            lower  = lsame(uplo,'L')
            lquery = (lwork == -1)
        !
            info = 0
            IF(.NOT.(wantz .OR. lsame(jobz,'N'))) THEN
                info = -1
            ELSE IF(.NOT.(lower .OR. lsame(uplo,'U'))) THEN
                info = -2
            ELSE IF(n < 0) THEN
                info = -3
            ELSE IF(lda < max(1,n)) THEN
                info = -5
            END IF
        !
            IF(info == 0) THEN
                nb      = ilaenv(1,'DSYTRD',uplo,n,-1,-1,-1)
                lwkopt  = max(1,(nb+2)*n)
                work(1) = lwkopt
        !
                IF(lwork < max(1,3*n-1) .AND. .NOT.lquery) info = -8
            END IF
        !
            IF(info /= 0) THEN
                CALL xerbla('DSYEV ',-info)
                RETURN
            ELSE IF(lquery) THEN
                RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF(n == 0) RETURN
        !
            IF(n == 1) THEN
                w(1)    = a(1,1)
                work(1) = 2
                IF(wantz) a(1,1) = ONE
                RETURN
            END IF
        !
        !     Get machine constants.
        !
            safmin = dlamch('Safe minimum')
            eps    = dlamch('Precision')
            smlnum = safmin/eps
            bignum = ONE/smlnum
            rmin   = sqrt(smlnum)
            rmax   = sqrt(bignum)
        !
        !     Scale matrix to allowable range, if necessary.
        !
            anrm   = dlansy('M',uplo,n,a,lda,work)
            iscale = 0
            IF(anrm > ZERO .AND. anrm < rmin) THEN
                iscale = 1
                sigma  = rmin/anrm
            ELSE IF(anrm > rmax) THEN
                iscale = 1
                sigma  = rmax/anrm
            END IF
            IF(iscale == 1) CALL dlascl(uplo,0,0,ONE,sigma,n,n,a,lda,info)
        !
        !     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
        !
            inde   = 1
            indtau = inde + n
            indwrk = indtau + n
            llwork = lwork - indwrk + 1
            CALL dsytrd(uplo,n,a,lda,w,work(inde),work(indtau),work(indwrk),llwork,iinfo)
        !
        !     For eigenvalues only, call DSTERF.  For eigenvectors, first call
        !     DORGTR to generate the orthogonal matrix, then call DSTEQR.
        !
            IF(.NOT.wantz) THEN
                CALL dsterf(n,w,work(inde),info)
            ELSE
                CALL dorgtr(uplo,n,a,lda,work(indtau),work(indwrk),llwork,iinfo)
                CALL dsteqr(jobz,n,w,work(inde),a,lda,work(indtau),info)
            END IF
        !
        !     If matrix was scaled, then rescale eigenvalues appropriately.
        !
            IF(iscale == 1) THEN
                IF(info == 0) THEN
                    imax = n
                ELSE
                    imax = info - 1
                END IF
                CALL dscal(imax,ONE/sigma,w,1)
            END IF
        !
        !     Set WORK(1) to optimal workspace size.
        !
            work(1) = lwkopt
        !
            RETURN
        !
        !     End of DSYEV
        !
        END SUBROUTINE
    
    !==================================================================================================
        LOGICAL FUNCTION LSAME(CA,CB)
    !==================================================================================================
        !  -- Reference BLAS level1 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER CA,CB
        !     ..
        !
        ! =====================================================================
        !
        !     .. Intrinsic Functions ..
            INTRINSIC ICHAR
        !     ..
        !     .. Local Scalars ..
            INTEGER INTA,INTB,ZCODE
        !     ..
        !
        !     Test if the characters are equal
        !
            LSAME = CA == CB
            IF (LSAME) RETURN
        !
        !     Now test for equivalence if both characters are alphabetic.
        !
            ZCODE = ICHAR('Z')
        !
        !     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
        !     machines, on which ICHAR returns a value with bit 8 set.
        !     ICHAR('A') on Prime machines returns 193 which is the same as
        !     ICHAR('A') on an EBCDIC machine.
        !
            INTA = ICHAR(CA)
            INTB = ICHAR(CB)
        !
            IF (ZCODE == 90 .OR. ZCODE == 122) THEN
        !
        !        ASCII is assumed - ZCODE is the ASCII code of either lower or
        !        upper case 'Z'.
        !
                IF (INTA >= 97 .AND. INTA <= 122) INTA = INTA - 32
                IF (INTB >= 97 .AND. INTB <= 122) INTB = INTB - 32
        !
            ELSE IF (ZCODE == 233 .OR. ZCODE == 169) THEN
        !
        !        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
        !        upper case 'Z'.
        !
                IF (INTA >= 129 .AND. INTA <= 137 .OR. &
                    INTA >= 145 .AND. INTA <= 153 .OR. &
                    INTA >= 162 .AND. INTA <= 169) INTA = INTA + 64
                IF (INTB >= 129 .AND. INTB <= 137 .OR. &
                    INTB >= 145 .AND. INTB <= 153 .OR. &
                    INTB >= 162 .AND. INTB <= 169) INTB = INTB + 64
        !
            ELSE IF (ZCODE == 218 .OR. ZCODE == 250) THEN
        !
        !        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
        !        plus 128 of either lower or upper case 'Z'.
        !
                IF (INTA >= 225 .AND. INTA <= 250) INTA = INTA - 32
                IF (INTB >= 225 .AND. INTB <= 250) INTB = INTB - 32
            END IF
            LSAME = INTA  ==  INTB
        !
        !     RETURN
        !
        !     End of LSAME
        !
        END FUNCTION
    
    !==================================================================================================
        INTEGER FUNCTION ILAENV(ISPEC,NAME,OPTS,N1,N2,N3,N4)
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER*( * )    NAME, OPTS
            INTEGER            ISPEC, N1, N2, N3, N4
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
            INTEGER            I, IC, IZ, NB, NBMIN, NX
            LOGICAL            CNAME, SNAME, TWOSTAGE
            CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*16
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
        !     ..
        !     .. External Functions ..
            INTEGER            IEEECK, IPARMQ, IPARAM2STAGE
            EXTERNAL           IEEECK, IPARMQ, IPARAM2STAGE
        !     ..
        !     .. Executable Statements ..
        !
            GO TO (10,10,10,80,90,100,110,120,&
                    130,140,150,160,160,160,160,160,160)ISPEC
        !
        !     Invalid value for ISPEC
        !
            ILAENV = -1
            RETURN
        !
        10 CONTINUE
        !
        !     Convert NAME to upper case if the first character is lower case.
        !
            ILAENV = 1
            SUBNAM = NAME
            IC = ICHAR( SUBNAM( 1: 1 ) )
            IZ = ICHAR( 'Z' )
            IF( IZ == 90 .OR. IZ == 122 ) THEN
        !
        !        ASCII character set
        !
                IF( IC >= 97 .AND. IC <= 122 ) THEN
                    SUBNAM( 1: 1 ) = CHAR( IC-32 )
                    DO 20 I = 2, 6
                    IC = ICHAR( SUBNAM( I: I ) )
                    IF( IC >= 97 .AND. IC <= 122 ) &
                        SUBNAM( I: I ) = CHAR( IC-32 )
        20       CONTINUE
                END IF
        !
            ELSE IF( IZ == 233 .OR. IZ == 169 ) THEN
        !
        !        EBCDIC character set
        !
                IF( ( IC >= 129 .AND. IC <= 137 ) .OR. &
                    ( IC >= 145 .AND. IC <= 153 ) .OR. &
                    ( IC >= 162 .AND. IC <= 169 ) ) THEN
                    SUBNAM( 1: 1 ) = CHAR( IC+64 )
                    DO 30 I = 2, 6
                    IC = ICHAR( SUBNAM( I: I ) )
                    IF( ( IC >= 129 .AND. IC <= 137 ) .OR. &
                        ( IC >= 145 .AND. IC <= 153 ) .OR. &
                        ( IC >= 162 .AND. IC <= 169 ) )SUBNAM( I:&
                        I ) = CHAR( IC+64 )
        30       CONTINUE
                END IF
        !
            ELSE IF( IZ == 218 .OR. IZ == 250 ) THEN
        !
        !        Prime machines:  ASCII+128
        !
                IF( IC >= 225 .AND. IC <= 250 ) THEN
                    SUBNAM( 1: 1 ) = CHAR( IC-32 )
                    DO 40 I = 2, 6
                    IC = ICHAR( SUBNAM( I: I ) )
                    IF( IC >= 225 .AND. IC <= 250 ) &
                        SUBNAM( I: I ) = CHAR( IC-32 )
        40       CONTINUE
                END IF
            END IF
        !
            C1 = SUBNAM( 1: 1 )
            SNAME = C1 == 'S' .OR. C1 == 'D'
            CNAME = C1 == 'C' .OR. C1 == 'Z'
            IF( .NOT.( CNAME .OR. SNAME ) ) &
                RETURN
            C2 = SUBNAM( 2: 3 )
            C3 = SUBNAM( 4: 6 )
            C4 = C3( 2: 3 )
            TWOSTAGE = LEN( SUBNAM ) >= 11 &
                        .AND. SUBNAM( 11: 11 ) == '2'
        !
            GO TO ( 50, 60, 70 )ISPEC
        !
        50 CONTINUE
        !
        !     ISPEC = 1:  block size
        !
        !     In these examples, separate code is provided for setting NB for
        !     real and complex.  We assume that NB will take the same value in
        !     single or double precision.
        !
            NB = 1
        !
            IF( SUBNAM(2:6) == 'LAORH' ) THEN
        !
        !        This is for *LAORHR_GETRFNP routine
        !
                IF( SNAME ) THEN
                    NB = 32
                ELSE
                    NB = 32
                END IF
            ELSE IF( C2 == 'GE' ) THEN
                IF( C3 == 'TRF' ) THEN
                    IF( SNAME ) THEN
                    NB = 64
                    ELSE
                    NB = 64
                    END IF
                ELSE IF( C3 == 'QRF' .OR. C3 == 'RQF' .OR. C3 == 'LQF' .OR. &
                        C3 == 'QLF' ) THEN
                    IF( SNAME ) THEN
                    NB = 32
                    ELSE
                    NB = 32
                    END IF
                ELSE IF( C3 == 'QR ') THEN
                    IF( N3  ==  1) THEN
                    IF( SNAME ) THEN
        !     M*N
                        IF ((N1*N2 <= 131072).OR.(N1 <= 8192)) THEN
                            NB = N1
                        ELSE
                            NB = 32768/N2
                        END IF
                    ELSE
                        IF ((N1*N2 <= 131072).OR.(N1 <= 8192)) THEN
                            NB = N1
                        ELSE
                            NB = 32768/N2
                        END IF
                    END IF
                    ELSE
                    IF( SNAME ) THEN
                        NB = 1
                    ELSE
                        NB = 1
                    END IF
                    END IF
                ELSE IF( C3 == 'LQ ') THEN
                    IF( N3  ==  2) THEN
                    IF( SNAME ) THEN
        !     M*N
                        IF ((N1*N2 <= 131072).OR.(N1 <= 8192)) THEN
                            NB = N1
                        ELSE
                            NB = 32768/N2
                        END IF
                    ELSE
                        IF ((N1*N2 <= 131072).OR.(N1 <= 8192)) THEN
                            NB = N1
                        ELSE
                            NB = 32768/N2
                        END IF
                    END IF
                    ELSE
                    IF( SNAME ) THEN
                        NB = 1
                    ELSE
                        NB = 1
                    END IF
                    END IF
                ELSE IF( C3 == 'HRD' ) THEN
                    IF( SNAME ) THEN
                    NB = 32
                    ELSE
                    NB = 32
                    END IF
                ELSE IF( C3 == 'BRD' ) THEN
                    IF( SNAME ) THEN
                    NB = 32
                    ELSE
                    NB = 32
                    END IF
                ELSE IF( C3 == 'TRI' ) THEN
                    IF( SNAME ) THEN
                    NB = 64
                    ELSE
                    NB = 64
                    END IF
                END IF
            ELSE IF( C2 == 'PO' ) THEN
                IF( C3 == 'TRF' ) THEN
                    IF( SNAME ) THEN
                    NB = 64
                    ELSE
                    NB = 64
                    END IF
                END IF
            ELSE IF( C2 == 'SY' ) THEN
                IF( C3 == 'TRF' ) THEN
                    IF( SNAME ) THEN
                    IF( TWOSTAGE ) THEN
                        NB = 192
                    ELSE
                        NB = 64
                    END IF
                    ELSE
                    IF( TWOSTAGE ) THEN
                        NB = 192
                    ELSE
                        NB = 64
                    END IF
                    END IF
                ELSE IF( SNAME .AND. C3 == 'TRD' ) THEN
                    NB = 32
                ELSE IF( SNAME .AND. C3 == 'GST' ) THEN
                    NB = 64
                END IF
            ELSE IF( CNAME .AND. C2 == 'HE' ) THEN
                IF( C3 == 'TRF' ) THEN
                    IF( TWOSTAGE ) THEN
                    NB = 192
                    ELSE
                    NB = 64
                    END IF
                ELSE IF( C3 == 'TRD' ) THEN
                    NB = 32
                ELSE IF( C3 == 'GST' ) THEN
                    NB = 64
                END IF
            ELSE IF( SNAME .AND. C2 == 'OR' ) THEN
                IF( C3( 1: 1 ) == 'G' ) THEN
                    IF( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == &
                        'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) &
                        THEN
                    NB = 32
                    END IF
                ELSE IF( C3( 1: 1 ) == 'M' ) THEN
                    IF( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == &
                        'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) &
                        THEN
                    NB = 32
                    END IF
                END IF
            ELSE IF( CNAME .AND. C2 == 'UN' ) THEN
                IF( C3( 1: 1 ) == 'G' ) THEN
                    IF( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == &
                        'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) &
                        THEN
                    NB = 32
                    END IF
                ELSE IF( C3( 1: 1 ) == 'M' ) THEN
                    IF( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == &
                        'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) &
                        THEN
                    NB = 32
                    END IF
                END IF
            ELSE IF( C2 == 'GB' ) THEN
                IF( C3 == 'TRF' ) THEN
                    IF( SNAME ) THEN
                    IF( N4 <= 64 ) THEN
                        NB = 1
                    ELSE
                        NB = 32
                    END IF
                    ELSE
                    IF( N4 <= 64 ) THEN
                        NB = 1
                    ELSE
                        NB = 32
                    END IF
                    END IF
                END IF
            ELSE IF( C2 == 'PB' ) THEN
                IF( C3 == 'TRF' ) THEN
                    IF( SNAME ) THEN
                    IF( N2 <= 64 ) THEN
                        NB = 1
                    ELSE
                        NB = 32
                    END IF
                    ELSE
                    IF( N2 <= 64 ) THEN
                        NB = 1
                    ELSE
                        NB = 32
                    END IF
                    END IF
                END IF
            ELSE IF( C2 == 'TR' ) THEN
                IF( C3 == 'TRI' ) THEN
                    IF( SNAME ) THEN
                    NB = 64
                    ELSE
                    NB = 64
                    END IF
                ELSE IF ( C3 == 'EVC' ) THEN
                    IF( SNAME ) THEN
                    NB = 64
                    ELSE
                    NB = 64
                    END IF
                ELSE IF( C3 == 'SYL' ) THEN
        !           The upper bound is to prevent overly aggressive scaling.
                    IF( SNAME ) THEN
                    NB = MIN( MAX( 48, INT( ( MIN( N1, N2 ) * 16 ) / 100) ), &
                            240 )
                    ELSE
                    NB = MIN( MAX( 24, INT( ( MIN( N1, N2 ) * 8 ) / 100) ), &
                            80 )
                    END IF
                END IF
            ELSE IF( C2 == 'LA' ) THEN
                IF( C3 == 'UUM' ) THEN
                    IF( SNAME ) THEN
                    NB = 64
                    ELSE
                    NB = 64
                    END IF
                ELSE IF( C3 == 'TRS' ) THEN
                    IF( SNAME ) THEN
                    NB = 32
                    ELSE
                    NB = 32
                    END IF
                END IF
            ELSE IF( SNAME .AND. C2 == 'ST' ) THEN
                IF( C3 == 'EBZ' ) THEN
                    NB = 1
                END IF
            ELSE IF( C2 == 'GG' ) THEN
                NB = 32
                IF( C3 == 'HD3' ) THEN
                    IF( SNAME ) THEN
                    NB = 32
                    ELSE
                    NB = 32
                    END IF
                END IF
            END IF
            ILAENV = NB
            RETURN
        !
        60 CONTINUE
        !
        !     ISPEC = 2:  minimum block size
        !
            NBMIN = 2
            IF( C2 == 'GE' ) THEN
                IF( C3 == 'QRF' .OR. C3 == 'RQF' .OR. C3 == 'LQF' .OR. C3 == &
                    'QLF' ) THEN
                    IF( SNAME ) THEN
                    NBMIN = 2
                    ELSE
                    NBMIN = 2
                    END IF
                ELSE IF( C3 == 'HRD' ) THEN
                    IF( SNAME ) THEN
                    NBMIN = 2
                    ELSE
                    NBMIN = 2
                    END IF
                ELSE IF( C3 == 'BRD' ) THEN
                    IF( SNAME ) THEN
                    NBMIN = 2
                    ELSE
                    NBMIN = 2
                    END IF
                ELSE IF( C3 == 'TRI' ) THEN
                    IF( SNAME ) THEN
                    NBMIN = 2
                    ELSE
                    NBMIN = 2
                    END IF
                END IF
            ELSE IF( C2 == 'SY' ) THEN
                IF( C3 == 'TRF' ) THEN
                    IF( SNAME ) THEN
                    NBMIN = 8
                    ELSE
                    NBMIN = 8
                    END IF
                ELSE IF( SNAME .AND. C3 == 'TRD' ) THEN
                    NBMIN = 2
                END IF
            ELSE IF( CNAME .AND. C2 == 'HE' ) THEN
                IF( C3 == 'TRD' ) THEN
                    NBMIN = 2
                END IF
            ELSE IF( SNAME .AND. C2 == 'OR' ) THEN
                IF( C3( 1: 1 ) == 'G' ) THEN
                    IF( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == &
                        'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) &
                        THEN
                    NBMIN = 2
                    END IF
                ELSE IF( C3( 1: 1 ) == 'M' ) THEN
                    IF( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == &
                        'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) &
                        THEN
                    NBMIN = 2
                    END IF
                END IF
            ELSE IF( CNAME .AND. C2 == 'UN' ) THEN
                IF( C3( 1: 1 ) == 'G' ) THEN
                    IF( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == &
                        'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) &
                        THEN
                    NBMIN = 2
                    END IF
                ELSE IF( C3( 1: 1 ) == 'M' ) THEN
                    IF( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == &
                        'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) &
                        THEN
                    NBMIN = 2
                    END IF
                END IF
            ELSE IF( C2 == 'GG' ) THEN
                NBMIN = 2
                IF( C3 == 'HD3' ) THEN
                    NBMIN = 2
                END IF
            END IF
            ILAENV = NBMIN
            RETURN
        !
        70 CONTINUE
        !
        !     ISPEC = 3:  crossover point
        !
            NX = 0
            IF( C2 == 'GE' ) THEN
                IF( C3 == 'QRF' .OR. C3 == 'RQF' .OR. C3 == 'LQF' .OR. C3 == &
                    'QLF' ) THEN
                    IF( SNAME ) THEN
                    NX = 128
                    ELSE
                    NX = 128
                    END IF
                ELSE IF( C3 == 'HRD' ) THEN
                    IF( SNAME ) THEN
                    NX = 128
                    ELSE
                    NX = 128
                    END IF
                ELSE IF( C3 == 'BRD' ) THEN
                    IF( SNAME ) THEN
                    NX = 128
                    ELSE
                    NX = 128
                    END IF
                END IF
            ELSE IF( C2 == 'SY' ) THEN
                IF( SNAME .AND. C3 == 'TRD' ) THEN
                    NX = 32
                END IF
            ELSE IF( CNAME .AND. C2 == 'HE' ) THEN
                IF( C3 == 'TRD' ) THEN
                    NX = 32
                END IF
            ELSE IF( SNAME .AND. C2 == 'OR' ) THEN
                IF( C3( 1: 1 ) == 'G' ) THEN
                    IF( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == &
                        'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) &
                        THEN
                    NX = 128
                    END IF
                END IF
            ELSE IF( CNAME .AND. C2 == 'UN' ) THEN
                IF( C3( 1: 1 ) == 'G' ) THEN
                    IF( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == &
                        'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) &
                        THEN
                    NX = 128
                    END IF
                END IF
            ELSE IF( C2 == 'GG' ) THEN
                NX = 128
                IF( C3 == 'HD3' ) THEN
                    NX = 128
                END IF
            END IF
            ILAENV = NX
            RETURN
        !
        80 CONTINUE
        !
        !     ISPEC = 4:  number of shifts (used by xHSEQR)
        !
            ILAENV = 6
            RETURN
        !
        90 CONTINUE
        !
        !     ISPEC = 5:  minimum column dimension (not used)
        !
            ILAENV = 2
            RETURN
        !
        100 CONTINUE
        !
        !     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
        !
            ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
            RETURN
        !
        110 CONTINUE
        !
        !     ISPEC = 7:  number of processors (not used)
        !
            ILAENV = 1
            RETURN
        !
        120 CONTINUE
        !
        !     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
        !
            ILAENV = 50
            RETURN
        !
        130 CONTINUE
        !
        !     ISPEC = 9:  maximum size of the subproblems at the bottom of the
        !                 computation tree in the divide-and-conquer algorithm
        !                 (used by xGELSD and xGESDD)
        !
            ILAENV = 25
            RETURN
        !
        140 CONTINUE
        !
        !     ISPEC = 10: ieee and infinity NaN arithmetic can be trusted not to trap
        !
        !     ILAENV = 0
            ILAENV = 1
            IF( ILAENV == 1 ) THEN
                ILAENV = IEEECK( 1, 0.0, 1.0 )
            END IF
            RETURN
        !
        150 CONTINUE
        !
        !     ISPEC = 11: ieee infinity arithmetic can be trusted not to trap
        !
        !     ILAENV = 0
            ILAENV = 1
            IF( ILAENV == 1 ) THEN
                ILAENV = IEEECK( 0, 0.0, 1.0 )
            END IF
            RETURN
        !
        160 CONTINUE
        !
        !     12 <= ISPEC <= 17: xHSEQR or related subroutines.
        !
            ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
            RETURN
        !
        !     End of ILAENV
        !
        END FUNCTION
    
    !==================================================================================================
          SUBROUTINE XERBLA( SRNAME, INFO )
    !==================================================================================================
        !
        !  -- Reference BLAS level1 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER*(*)      SRNAME
            INTEGER            INFO
        !     ..
        !
        ! =====================================================================
        !
        !     .. Intrinsic Functions ..
            INTRINSIC          LEN_TRIM
        !     ..
        !     .. Executable Statements ..
        !
            WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
        !
            STOP
        !
        9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ', &
                    'an illegal value' )
        !
        !     End of XERBLA
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DSCAL(N,DA,DX,INCX)
    !==================================================================================================
        !
        !  -- Reference BLAS level1 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION DA
            INTEGER INCX,N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION DX(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
            INTEGER I,M,MP1,NINCX
        !     .. Parameters ..
            DOUBLE PRECISION ONE
            PARAMETER (ONE=1.0D+0)
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MOD
        !     ..
            IF (N <= 0 .OR. INCX <= 0 .OR. DA == ONE) RETURN
            IF (INCX == 1) THEN
        !
        !        code for increment equal to 1
        !
        !
        !        clean-up loop
        !
                M = MOD(N,5)
                IF (M /= 0) THEN
                    DO I = 1,M
                    DX(I) = DA*DX(I)
                    END DO
                    IF (N < 5) RETURN
                END IF
                MP1 = M + 1
                DO I = MP1,N,5
                    DX(I) = DA*DX(I)
                    DX(I+1) = DA*DX(I+1)
                    DX(I+2) = DA*DX(I+2)
                    DX(I+3) = DA*DX(I+3)
                    DX(I+4) = DA*DX(I+4)
                END DO
            ELSE
        !
        !        code for increment not equal to 1
        !
                NINCX = N*INCX
                DO I = 1,NINCX,INCX
                    DX(I) = DA*DX(I)
                END DO
            END IF
            RETURN
        !
        !     End of DSCAL
        !
            END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
    !==================================================================================================
        !
        !  -- Reference BLAS level1 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER INCX,INCY,N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION DX(*),DY(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
            DOUBLE PRECISION DTEMP
            INTEGER I,IX,IY,M,MP1
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MOD
        !     ..
            IF (N <= 0) RETURN
            IF (INCX == 1 .AND. INCY == 1) THEN
        !
        !       code for both increments equal to 1
        !
        !
        !       clean-up loop
        !
                M = MOD(N,3)
                IF (M /= 0) THEN
                    DO I = 1,M
                    DTEMP = DX(I)
                    DX(I) = DY(I)
                    DY(I) = DTEMP
                    END DO
                    IF (N < 3) RETURN
                END IF
                MP1 = M + 1
                DO I = MP1,N,3
                    DTEMP = DX(I)
                    DX(I) = DY(I)
                    DY(I) = DTEMP
                    DTEMP = DX(I+1)
                    DX(I+1) = DY(I+1)
                    DY(I+1) = DTEMP
                    DTEMP = DX(I+2)
                    DX(I+2) = DY(I+2)
                    DY(I+2) = DTEMP
                END DO
            ELSE
        !
        !       code for unequal increments or equal increments not equal
        !         to 1
        !
                IX = 1
                IY = 1
                IF (INCX < 0) IX = (-N+1)*INCX + 1
                IF (INCY < 0) IY = (-N+1)*INCY + 1
                DO I = 1,N
                    DTEMP = DX(IX)
                    DX(IX) = DY(IY)
                    DY(IY) = DTEMP
                    IX = IX + INCX
                    IY = IY + INCY
                END DO
            END IF
            RETURN
        !
        !     End of DSWAP
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    !==================================================================================================
        !
        !  -- Reference BLAS level3 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION ALPHA,BETA
            INTEGER K,LDA,LDB,LDC,M,N
            CHARACTER TRANSA,TRANSB
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. External Functions ..
            LOGICAL LSAME
            EXTERNAL LSAME
        !     ..
        !     .. External Subroutines ..
            EXTERNAL XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MAX
        !     ..
        !     .. Local Scalars ..
            DOUBLE PRECISION TEMP
            INTEGER I,INFO,J,L,NROWA,NROWB
            LOGICAL NOTA,NOTB
        !     ..
        !     .. Parameters ..
            DOUBLE PRECISION ONE,ZERO
            PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
        !     ..
        !
        !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
        !     transposed and set  NROWA and NROWB  as the number of rows of  A
        !     and  B  respectively.
        !
            NOTA = LSAME(TRANSA,'N')
            NOTB = LSAME(TRANSB,'N')
            IF (NOTA) THEN
                NROWA = M
            ELSE
                NROWA = K
            END IF
            IF (NOTB) THEN
                NROWB = K
            ELSE
                NROWB = N
            END IF
        !
        !     Test the input parameters.
        !
            INFO = 0
            IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. &
                (.NOT.LSAME(TRANSA,'T'))) THEN
                INFO = 1
            ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. &
                    (.NOT.LSAME(TRANSB,'T'))) THEN
                INFO = 2
            ELSE IF (M < 0) THEN
                INFO = 3
            ELSE IF (N < 0) THEN
                INFO = 4
            ELSE IF (K < 0) THEN
                INFO = 5
            ELSE IF (LDA < MAX(1,NROWA)) THEN
                INFO = 8
            ELSE IF (LDB < MAX(1,NROWB)) THEN
                INFO = 10
            ELSE IF (LDC < MAX(1,M)) THEN
                INFO = 13
            END IF
            IF (INFO /= 0) THEN
                CALL XERBLA('DGEMM ',INFO)
                RETURN
            END IF
        !
        !     Quick return if possible.
        !
            IF ((M == 0) .OR. (N == 0) .OR. &
                (((ALPHA == ZERO).OR. (K == 0)).AND. (BETA == ONE))) RETURN
        !
        !     And if  alpha == zero.
        !
            IF (ALPHA == ZERO) THEN
                IF (BETA == ZERO) THEN
                    DO 20 J = 1,N
                        DO 10 I = 1,M
                            C(I,J) = ZERO
        10             CONTINUE
        20         CONTINUE
                ELSE
                    DO 40 J = 1,N
                        DO 30 I = 1,M
                            C(I,J) = BETA*C(I,J)
        30             CONTINUE
        40         CONTINUE
                END IF
                RETURN
            END IF
        !
        !     Start the operations.
        !
            IF (NOTB) THEN
                IF (NOTA) THEN
        !
        !           Form  C := alpha*A*B + beta*C.
        !
                    DO 90 J = 1,N
                        IF (BETA == ZERO) THEN
                            DO 50 I = 1,M
                                C(I,J) = ZERO
        50                 CONTINUE
                        ELSE IF (BETA /= ONE) THEN
                            DO 60 I = 1,M
                                C(I,J) = BETA*C(I,J)
        60                 CONTINUE
                        END IF
                        DO 80 L = 1,K
                            TEMP = ALPHA*B(L,J)
                            DO 70 I = 1,M
                                C(I,J) = C(I,J) + TEMP*A(I,L)
        70                 CONTINUE
        80             CONTINUE
        90         CONTINUE
                ELSE
        !
        !           Form  C := alpha*A**T*B + beta*C
        !
                    DO 120 J = 1,N
                        DO 110 I = 1,M
                            TEMP = ZERO
                            DO 100 L = 1,K
                                TEMP = TEMP + A(L,I)*B(L,J)
        100                 CONTINUE
                            IF (BETA == ZERO) THEN
                                C(I,J) = ALPHA*TEMP
                            ELSE
                                C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                            END IF
        110             CONTINUE
        120         CONTINUE
                END IF
            ELSE
                IF (NOTA) THEN
        !
        !           Form  C := alpha*A*B**T + beta*C
        !
                    DO 170 J = 1,N
                        IF (BETA == ZERO) THEN
                            DO 130 I = 1,M
                                C(I,J) = ZERO
        130                 CONTINUE
                        ELSE IF (BETA /= ONE) THEN
                            DO 140 I = 1,M
                                C(I,J) = BETA*C(I,J)
        140                 CONTINUE
                        END IF
                        DO 160 L = 1,K
                            TEMP = ALPHA*B(J,L)
                            DO 150 I = 1,M
                                C(I,J) = C(I,J) + TEMP*A(I,L)
        150                 CONTINUE
        160             CONTINUE
        170         CONTINUE
                ELSE
        !
        !           Form  C := alpha*A**T*B**T + beta*C
        !
                    DO 200 J = 1,N
                        DO 190 I = 1,M
                            TEMP = ZERO
                            DO 180 L = 1,K
                                TEMP = TEMP + A(L,I)*B(J,L)
        180                 CONTINUE
                            IF (BETA == ZERO) THEN
                                C(I,J) = ALPHA*TEMP
                            ELSE
                                C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                            END IF
        190             CONTINUE
        200         CONTINUE
                END IF
            END IF
        !
            RETURN
        !
        !     End of DGEMM
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    !==================================================================================================
        !
        !  -- Reference BLAS level2 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION ALPHA,BETA
            INTEGER INCX,INCY,LDA,M,N
            CHARACTER TRANS
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION A(LDA,*),X(*),Y(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION ONE,ZERO
            PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
        !     ..
        !     .. Local Scalars ..
            DOUBLE PRECISION TEMP
            INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
        !     ..
        !     .. External Functions ..
            LOGICAL LSAME
            EXTERNAL LSAME
        !     ..
        !     .. External Subroutines ..
            EXTERNAL XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MAX
        !     ..
        !
        !     Test the input parameters.
        !
            INFO = 0
            IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
                .NOT.LSAME(TRANS,'C')) THEN
                INFO = 1
            ELSE IF (M < 0) THEN
                INFO = 2
            ELSE IF (N < 0) THEN
                INFO = 3
            ELSE IF (LDA < MAX(1,M)) THEN
                INFO = 6
            ELSE IF (INCX == 0) THEN
                INFO = 8
            ELSE IF (INCY == 0) THEN
                INFO = 11
            END IF
            IF (INFO /= 0) THEN
                CALL XERBLA('DGEMV ',INFO)
                RETURN
            END IF
        !
        !     Quick return if possible.
        !
            IF ((M == 0) .OR. (N == 0) .OR. &
                ((ALPHA == ZERO).AND. (BETA == ONE))) RETURN
        !
        !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
        !     up the start points in  X  and  Y.
        !
            IF (LSAME(TRANS,'N')) THEN
                LENX = N
                LENY = M
            ELSE
                LENX = M
                LENY = N
            END IF
            IF (INCX > 0) THEN
                KX = 1
            ELSE
                KX = 1 - (LENX-1)*INCX
            END IF
            IF (INCY > 0) THEN
                KY = 1
            ELSE
                KY = 1 - (LENY-1)*INCY
            END IF
        !
        !     Start the operations. In this version the elements of A are
        !     accessed sequentially with one pass through A.
        !
        !     First form  y := beta*y.
        !
            IF (BETA /= ONE) THEN
                IF (INCY == 1) THEN
                    IF (BETA == ZERO) THEN
                        DO 10 I = 1,LENY
                            Y(I) = ZERO
        10             CONTINUE
                    ELSE
                        DO 20 I = 1,LENY
                            Y(I) = BETA*Y(I)
        20             CONTINUE
                    END IF
                ELSE
                    IY = KY
                    IF (BETA == ZERO) THEN
                        DO 30 I = 1,LENY
                            Y(IY) = ZERO
                            IY = IY + INCY
        30             CONTINUE
                    ELSE
                        DO 40 I = 1,LENY
                            Y(IY) = BETA*Y(IY)
                            IY = IY + INCY
        40             CONTINUE
                    END IF
                END IF
            END IF
            IF (ALPHA == ZERO) RETURN
            IF (LSAME(TRANS,'N')) THEN
        !
        !        Form  y := alpha*A*x + y.
        !
                JX = KX
                IF (INCY == 1) THEN
                    DO 60 J = 1,N
                        TEMP = ALPHA*X(JX)
                        DO 50 I = 1,M
                            Y(I) = Y(I) + TEMP*A(I,J)
        50             CONTINUE
                        JX = JX + INCX
        60         CONTINUE
                ELSE
                    DO 80 J = 1,N
                        TEMP = ALPHA*X(JX)
                        IY = KY
                        DO 70 I = 1,M
                            Y(IY) = Y(IY) + TEMP*A(I,J)
                            IY = IY + INCY
        70             CONTINUE
                        JX = JX + INCX
        80         CONTINUE
                END IF
            ELSE
        !
        !        Form  y := alpha*A**T*x + y.
        !
                JY = KY
                IF (INCX == 1) THEN
                    DO 100 J = 1,N
                        TEMP = ZERO
                        DO 90 I = 1,M
                            TEMP = TEMP + A(I,J)*X(I)
        90             CONTINUE
                        Y(JY) = Y(JY) + ALPHA*TEMP
                        JY = JY + INCY
        100         CONTINUE
                ELSE
                    DO 120 J = 1,N
                        TEMP = ZERO
                        IX = KX
                        DO 110 I = 1,M
                            TEMP = TEMP + A(I,J)*X(IX)
                            IX = IX + INCX
        110             CONTINUE
                        Y(JY) = Y(JY) + ALPHA*TEMP
                        JY = JY + INCY
        120         CONTINUE
                END IF
            END IF
        !
            RETURN
        !
        !     End of DGEMV
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    !==================================================================================================
        !
        !  -- Reference BLAS level3 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION ALPHA
            INTEGER LDA,LDB,M,N
            CHARACTER DIAG,SIDE,TRANSA,UPLO
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION A(LDA,*),B(LDB,*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. External Functions ..
            LOGICAL LSAME
            EXTERNAL LSAME
        !     ..
        !     .. External Subroutines ..
            EXTERNAL XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MAX
        !     ..
        !     .. Local Scalars ..
            DOUBLE PRECISION TEMP
            INTEGER I,INFO,J,K,NROWA
            LOGICAL LSIDE,NOUNIT,UPPER
        !     ..
        !     .. Parameters ..
            DOUBLE PRECISION ONE,ZERO
            PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
        !     ..
        !
        !     Test the input parameters.
        !
            LSIDE = LSAME(SIDE,'L')
            IF (LSIDE) THEN
                NROWA = M
            ELSE
                NROWA = N
            END IF
            NOUNIT = LSAME(DIAG,'N')
            UPPER = LSAME(UPLO,'U')
        !
            INFO = 0
            IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
                INFO = 1
            ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
                INFO = 2
            ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
                    (.NOT.LSAME(TRANSA,'T')) .AND. &
                    (.NOT.LSAME(TRANSA,'C'))) THEN
                INFO = 3
            ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
                INFO = 4
            ELSE IF (M < 0) THEN
                INFO = 5
            ELSE IF (N < 0) THEN
                INFO = 6
            ELSE IF (LDA < MAX(1,NROWA)) THEN
                INFO = 9
            ELSE IF (LDB < MAX(1,M)) THEN
                INFO = 11
            END IF
            IF (INFO /= 0) THEN
                CALL XERBLA('DTRSM ',INFO)
                RETURN
            END IF
        !
        !     Quick return if possible.
        !
            IF (M == 0 .OR. N == 0) RETURN
        !
        !     And when  alpha == zero.
        !
            IF (ALPHA == ZERO) THEN
                DO 20 J = 1,N
                    DO 10 I = 1,M
                        B(I,J) = ZERO
        10         CONTINUE
        20     CONTINUE
                RETURN
            END IF
        !
        !     Start the operations.
        !
            IF (LSIDE) THEN
                IF (LSAME(TRANSA,'N')) THEN
        !
        !           Form  B := alpha*inv( A )*B.
        !
                    IF (UPPER) THEN
                        DO 60 J = 1,N
                            IF (ALPHA /= ONE) THEN
                                DO 30 I = 1,M
                                    B(I,J) = ALPHA*B(I,J)
        30                     CONTINUE
                            END IF
                            DO 50 K = M,1,-1
                                IF (B(K,J) /= ZERO) THEN
                                    IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                                    DO 40 I = 1,K - 1
                                        B(I,J) = B(I,J) - B(K,J)*A(I,K)
        40                         CONTINUE
                                END IF
        50                 CONTINUE
        60             CONTINUE
                    ELSE
                        DO 100 J = 1,N
                            IF (ALPHA /= ONE) THEN
                                DO 70 I = 1,M
                                    B(I,J) = ALPHA*B(I,J)
        70                     CONTINUE
                            END IF
                            DO 90 K = 1,M
                                IF (B(K,J) /= ZERO) THEN
                                    IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                                    DO 80 I = K + 1,M
                                        B(I,J) = B(I,J) - B(K,J)*A(I,K)
        80                         CONTINUE
                                END IF
        90                 CONTINUE
        100             CONTINUE
                    END IF
                ELSE
        !
        !           Form  B := alpha*inv( A**T )*B.
        !
                    IF (UPPER) THEN
                        DO 130 J = 1,N
                            DO 120 I = 1,M
                                TEMP = ALPHA*B(I,J)
                                DO 110 K = 1,I - 1
                                    TEMP = TEMP - A(K,I)*B(K,J)
        110                     CONTINUE
                                IF (NOUNIT) TEMP = TEMP/A(I,I)
                                B(I,J) = TEMP
        120                 CONTINUE
        130             CONTINUE
                    ELSE
                        DO 160 J = 1,N
                            DO 150 I = M,1,-1
                                TEMP = ALPHA*B(I,J)
                                DO 140 K = I + 1,M
                                    TEMP = TEMP - A(K,I)*B(K,J)
        140                     CONTINUE
                                IF (NOUNIT) TEMP = TEMP/A(I,I)
                                B(I,J) = TEMP
        150                 CONTINUE
        160             CONTINUE
                    END IF
                END IF
            ELSE
                IF (LSAME(TRANSA,'N')) THEN
        !
        !           Form  B := alpha*B*inv( A ).
        !
                    IF (UPPER) THEN
                        DO 210 J = 1,N
                            IF (ALPHA /= ONE) THEN
                                DO 170 I = 1,M
                                    B(I,J) = ALPHA*B(I,J)
        170                     CONTINUE
                            END IF
                            DO 190 K = 1,J - 1
                                IF (A(K,J) /= ZERO) THEN
                                    DO 180 I = 1,M
                                        B(I,J) = B(I,J) - A(K,J)*B(I,K)
        180                         CONTINUE
                                END IF
        190                 CONTINUE
                            IF (NOUNIT) THEN
                                TEMP = ONE/A(J,J)
                                DO 200 I = 1,M
                                    B(I,J) = TEMP*B(I,J)
        200                     CONTINUE
                            END IF
        210             CONTINUE
                    ELSE
                        DO 260 J = N,1,-1
                            IF (ALPHA /= ONE) THEN
                                DO 220 I = 1,M
                                    B(I,J) = ALPHA*B(I,J)
        220                     CONTINUE
                            END IF
                            DO 240 K = J + 1,N
                                IF (A(K,J) /= ZERO) THEN
                                    DO 230 I = 1,M
                                        B(I,J) = B(I,J) - A(K,J)*B(I,K)
        230                         CONTINUE
                                END IF
        240                 CONTINUE
                            IF (NOUNIT) THEN
                                TEMP = ONE/A(J,J)
                                DO 250 I = 1,M
                                    B(I,J) = TEMP*B(I,J)
        250                     CONTINUE
                            END IF
        260             CONTINUE
                    END IF
                ELSE
        !
        !           Form  B := alpha*B*inv( A**T ).
        !
                    IF (UPPER) THEN
                        DO 310 K = N,1,-1
                            IF (NOUNIT) THEN
                                TEMP = ONE/A(K,K)
                                DO 270 I = 1,M
                                    B(I,K) = TEMP*B(I,K)
        270                     CONTINUE
                            END IF
                            DO 290 J = 1,K - 1
                                IF (A(J,K) /= ZERO) THEN
                                    TEMP = A(J,K)
                                    DO 280 I = 1,M
                                        B(I,J) = B(I,J) - TEMP*B(I,K)
        280                         CONTINUE
                                END IF
        290                 CONTINUE
                            IF (ALPHA /= ONE) THEN
                                DO 300 I = 1,M
                                    B(I,K) = ALPHA*B(I,K)
        300                     CONTINUE
                            END IF
        310             CONTINUE
                    ELSE
                        DO 360 K = 1,N
                            IF (NOUNIT) THEN
                                TEMP = ONE/A(K,K)
                                DO 320 I = 1,M
                                    B(I,K) = TEMP*B(I,K)
        320                     CONTINUE
                            END IF
                            DO 340 J = K + 1,N
                                IF (A(J,K) /= ZERO) THEN
                                    TEMP = A(J,K)
                                    DO 330 I = 1,M
                                        B(I,J) = B(I,J) - TEMP*B(I,K)
        330                         CONTINUE
                                END IF
        340                 CONTINUE
                            IF (ALPHA /= ONE) THEN
                                DO 350 I = 1,M
                                    B(I,K) = ALPHA*B(I,K)
        350                     CONTINUE
                            END IF
        360             CONTINUE
                    END IF
                END IF
            END IF
        !
            RETURN
        !
        !     End of DTRSM
        !
        END SUBROUTINE
    
    !==================================================================================================
        DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine (version 3.3.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     Based on LAPACK DLAMCH but with Fortran 95 query functions
        !     See: http://www.cs.utk.edu/~luszczek/lapack/lamch.html
        !     and  http://www.netlib.org/lapack-dev/lapack-coding/program-style.html#id2537289
        !     July 2010
        !
        !     .. Scalar Arguments ..
            CHARACTER          CMACH
        !     ..
        !
        !  Purpose
        !  =======
        !
        !  DLAMCH determines double precision machine parameters.
        !
        !  Arguments
        !  =========
        !
        !  CMACH   (input) CHARACTER*1
        !          Specifies the value to be returned by DLAMCH:
        !          = 'E' or 'e',   DLAMCH := eps
        !          = 'S' or 's ,   DLAMCH := sfmin
        !          = 'B' or 'b',   DLAMCH := base
        !          = 'P' or 'p',   DLAMCH := eps*base
        !          = 'N' or 'n',   DLAMCH := t
        !          = 'R' or 'r',   DLAMCH := rnd
        !          = 'M' or 'm',   DLAMCH := emin
        !          = 'U' or 'u',   DLAMCH := rmin
        !          = 'L' or 'l',   DLAMCH := emax
        !          = 'O' or 'o',   DLAMCH := rmax
        !
        !          where
        !
        !          eps   = relative machine precision
        !          sfmin = safe minimum, such that 1/sfmin does not overflow
        !          base  = base of the machine
        !          prec  = eps*base
        !          t     = number of (base) digits in the mantissa
        !          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
        !          emin  = minimum exponent before (gradual) underflow
        !          rmin  = underflow threshold - base**(emin-1)
        !          emax  = largest exponent before overflow
        !          rmax  = overflow threshold  - (base**emax)*(1-eps)
        !
        ! =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE, ZERO
            PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, RMACH
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            EXTERNAL           LSAME
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT, &
                               MINEXPONENT, RADIX, TINY
        !     ..
        !     .. Executable Statements ..
        !
        !
        !     Assume rounding, not chopping. Always.
        !
            RND = ONE
        !
            IF( ONE == RND ) THEN
                EPS = EPSILON(ZERO) * 0.5
            ELSE
                EPS = EPSILON(ZERO)
            END IF
        !
            IF( LSAME( CMACH, 'E' ) ) THEN
                RMACH = EPS
            ELSE IF( LSAME( CMACH, 'S' ) ) THEN
                SFMIN = TINY(ZERO)
                SMALL = ONE / HUGE(ZERO)
                IF( SMALL >= SFMIN ) THEN
        !
        !           Use SMALL plus a bit, to avoid the possibility of rounding
        !           causing overflow when computing  1/sfmin.
        !
                    SFMIN = SMALL*( ONE+EPS )
                END IF
                RMACH = SFMIN
            ELSE IF( LSAME( CMACH, 'B' ) ) THEN
                RMACH = RADIX(ZERO)
            ELSE IF( LSAME( CMACH, 'P' ) ) THEN
                RMACH = EPS * RADIX(ZERO)
            ELSE IF( LSAME( CMACH, 'N' ) ) THEN
                RMACH = DIGITS(ZERO)
            ELSE IF( LSAME( CMACH, 'R' ) ) THEN
                RMACH = RND
            ELSE IF( LSAME( CMACH, 'M' ) ) THEN
                RMACH = MINEXPONENT(ZERO)
            ELSE IF( LSAME( CMACH, 'U' ) ) THEN
                RMACH = tiny(zero)
            ELSE IF( LSAME( CMACH, 'L' ) ) THEN
                RMACH = MAXEXPONENT(ZERO)
            ELSE IF( LSAME( CMACH, 'O' ) ) THEN
                RMACH = HUGE(ZERO)
            ELSE
                RMACH = ZERO
            END IF
        !
            DLAMCH = RMACH
            RETURN
        !
        !     End of DLAMCH
        !
        END FUNCTION
    
    !==================================================================================================
        DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          NORM, UPLO
            INTEGER            LDA, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * ), WORK( * )
        !     ..
        !
        ! =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE, ZERO
            PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            I, J
            DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DLASSQ
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME, DISNAN
            EXTERNAL           LSAME, DISNAN
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          ABS, SQRT
        !     ..
        !     .. Executable Statements ..
        !
            IF( N == 0 ) THEN
            VALUE = ZERO
            ELSE IF( LSAME( NORM, 'M' ) ) THEN
        !
        !        Find max(abs(A(i,j))).
        !
            VALUE = ZERO
            IF( LSAME( UPLO, 'U' ) ) THEN
                DO 20 J = 1, N
                    DO 10 I = 1, J
                        SUM = ABS( A( I, J ) )
                        IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
        10          CONTINUE
        20       CONTINUE
            ELSE
                DO 40 J = 1, N
                    DO 30 I = J, N
                        SUM = ABS( A( I, J ) )
                        IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
        30          CONTINUE
        40       CONTINUE
            END IF
            ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. &
                    ( NORM == '1' ) ) THEN
        !
        !        Find normI(A) ( = norm1(A), since A is symmetric).
        !
            VALUE = ZERO
            IF( LSAME( UPLO, 'U' ) ) THEN
                DO 60 J = 1, N
                    SUM = ZERO
                    DO 50 I = 1, J - 1
                        ABSA = ABS( A( I, J ) )
                        SUM = SUM + ABSA
                        WORK( I ) = WORK( I ) + ABSA
        50          CONTINUE
                    WORK( J ) = SUM + ABS( A( J, J ) )
        60       CONTINUE
                DO 70 I = 1, N
                    SUM = WORK( I )
                    IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
        70       CONTINUE
            ELSE
                DO 80 I = 1, N
                    WORK( I ) = ZERO
        80       CONTINUE
                DO 100 J = 1, N
                    SUM = WORK( J ) + ABS( A( J, J ) )
                    DO 90 I = J + 1, N
                        ABSA = ABS( A( I, J ) )
                        SUM = SUM + ABSA
                        WORK( I ) = WORK( I ) + ABSA
        90          CONTINUE
                    IF( VALUE  <  SUM .OR. DISNAN( SUM ) ) VALUE = SUM
        100       CONTINUE
            END IF
            ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
        !
        !        Find normF(A).
        !
            SCALE = ZERO
            SUM = ONE
            IF( LSAME( UPLO, 'U' ) ) THEN
                DO 110 J = 2, N
                    CALL DLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
        110       CONTINUE
            ELSE
                DO 120 J = 1, N - 1
                    CALL DLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
        120       CONTINUE
            END IF
            SUM = 2*SUM
            CALL DLASSQ( N, A, LDA+1, SCALE, SUM )
            VALUE = SCALE*SQRT( SUM )
            END IF
        !
            DLANSY = VALUE
            RETURN
        !
        !     End of DLANSY
        !
        END FUNCTION
    
    !==================================================================================================
        SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          TYPE
            INTEGER            INFO, KL, KU, LDA, M, N
            DOUBLE PRECISION   CFROM, CTO
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ZERO, ONE
            PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
        !     ..
        !     .. Local Scalars ..
            LOGICAL            DONE
            INTEGER            I, ITYPE, J, K1, K2, K3, K4
            DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME, DISNAN
            DOUBLE PRECISION   DLAMCH
            EXTERNAL           LSAME, DLAMCH, DISNAN
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          ABS, MAX, MIN
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           XERBLA
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
            INFO = 0
        !
            IF( LSAME( TYPE, 'G' ) ) THEN
            ITYPE = 0
            ELSE IF( LSAME( TYPE, 'L' ) ) THEN
            ITYPE = 1
            ELSE IF( LSAME( TYPE, 'U' ) ) THEN
            ITYPE = 2
            ELSE IF( LSAME( TYPE, 'H' ) ) THEN
            ITYPE = 3
            ELSE IF( LSAME( TYPE, 'B' ) ) THEN
            ITYPE = 4
            ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
            ITYPE = 5
            ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
            ITYPE = 6
            ELSE
            ITYPE = -1
            END IF
        !
            IF( ITYPE == -1 ) THEN
            INFO = -1
            ELSE IF( CFROM == ZERO .OR. DISNAN(CFROM) ) THEN
            INFO = -4
            ELSE IF( DISNAN(CTO) ) THEN
            INFO = -5
            ELSE IF( M < 0 ) THEN
            INFO = -6
            ELSE IF( N < 0 .OR. ( ITYPE == 4 .AND. N /= M ) .OR. &
                    ( ITYPE == 5 .AND. N /= M ) ) THEN
            INFO = -7
            ELSE IF( ITYPE <= 3 .AND. LDA < MAX( 1, M ) ) THEN
            INFO = -9
            ELSE IF( ITYPE >= 4 ) THEN
            IF( KL < 0 .OR. KL >  MAX( M-1, 0 ) ) THEN
                INFO = -2
            ELSE IF( KU < 0 .OR. KU >  MAX( N-1, 0 ) .OR. &
                    ( ( ITYPE == 4 .OR. ITYPE == 5 ) .AND. KL /= KU ) ) &
                    THEN
                INFO = -3
            ELSE IF( ( ITYPE == 4 .AND. LDA < KL+1 ) .OR. &
                    ( ITYPE == 5 .AND. LDA < KU+1 ) .OR. &
                    ( ITYPE == 6 .AND. LDA < 2*KL+KU+1 ) ) THEN
                INFO = -9
            END IF
            END IF
        !
            IF( INFO /= 0 ) THEN
            CALL XERBLA( 'DLASCL', -INFO )
            RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( N == 0 .OR. M == 0 ) &
                RETURN
        !
        !     Get machine parameters
        !
            SMLNUM = DLAMCH( 'S' )
            BIGNUM = ONE / SMLNUM
        !
            CFROMC = CFROM
            CTOC = CTO
        !
        10 CONTINUE
            CFROM1 = CFROMC*SMLNUM
            IF( CFROM1 == CFROMC ) THEN
        !        CFROMC is an inf.  Multiply by a correctly signed zero for
        !        finite CTOC, or a NaN if CTOC is infinite.
            MUL = CTOC / CFROMC
            DONE = .TRUE.
            CTO1 = CTOC
            ELSE
            CTO1 = CTOC / BIGNUM
            IF( CTO1 == CTOC ) THEN
        !           CTOC is either 0 or an inf.  In both cases, CTOC itself
        !           serves as the correct multiplication factor.
                MUL = CTOC
                DONE = .TRUE.
                CFROMC = ONE
            ELSE IF( ABS( CFROM1 ) >  ABS( CTOC ) .AND. CTOC /= ZERO ) THEN
                MUL = SMLNUM
                DONE = .FALSE.
                CFROMC = CFROM1
            ELSE IF( ABS( CTO1 ) >  ABS( CFROMC ) ) THEN
                MUL = BIGNUM
                DONE = .FALSE.
                CTOC = CTO1
            ELSE
                MUL = CTOC / CFROMC
                DONE = .TRUE.
                IF (MUL  ==  ONE) &
                    RETURN
            END IF
            END IF
        !
            IF( ITYPE == 0 ) THEN
        !
        !        Full matrix
        !
            DO 30 J = 1, N
                DO 20 I = 1, M
                    A( I, J ) = A( I, J )*MUL
        20       CONTINUE
        30    CONTINUE
        !
            ELSE IF( ITYPE == 1 ) THEN
        !
        !        Lower triangular matrix
        !
            DO 50 J = 1, N
                DO 40 I = J, M
                    A( I, J ) = A( I, J )*MUL
        40       CONTINUE
        50    CONTINUE
        !
            ELSE IF( ITYPE == 2 ) THEN
        !
        !        Upper triangular matrix
        !
            DO 70 J = 1, N
                DO 60 I = 1, MIN( J, M )
                    A( I, J ) = A( I, J )*MUL
        60       CONTINUE
        70    CONTINUE
        !
            ELSE IF( ITYPE == 3 ) THEN
        !
        !        Upper Hessenberg matrix
        !
            DO 90 J = 1, N
                DO 80 I = 1, MIN( J+1, M )
                    A( I, J ) = A( I, J )*MUL
        80       CONTINUE
        90    CONTINUE
        !
            ELSE IF( ITYPE == 4 ) THEN
        !
        !        Lower half of a symmetric band matrix
        !
            K3 = KL + 1
            K4 = N + 1
            DO 110 J = 1, N
                DO 100 I = 1, MIN( K3, K4-J )
                    A( I, J ) = A( I, J )*MUL
        100       CONTINUE
        110    CONTINUE
        !
            ELSE IF( ITYPE == 5 ) THEN
        !
        !        Upper half of a symmetric band matrix
        !
            K1 = KU + 2
            K3 = KU + 1
            DO 130 J = 1, N
                DO 120 I = MAX( K1-J, 1 ), K3
                    A( I, J ) = A( I, J )*MUL
        120       CONTINUE
        130    CONTINUE
        !
            ELSE IF( ITYPE == 6 ) THEN
        !
        !        Band matrix
        !
            K1 = KL + KU + 2
            K2 = KL + 1
            K3 = 2*KL + KU + 1
            K4 = KL + KU + 1 + M
            DO 150 J = 1, N
                DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
                    A( I, J ) = A( I, J )*MUL
        140       CONTINUE
        150    CONTINUE
        !
            END IF
        !
            IF( .NOT.DONE ) &
                GO TO 10
        !
            RETURN
        !
        !     End of DLASCL
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
    !==================================================================================================
        !
        !  -- LAPACK computational routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          UPLO
            INTEGER            INFO, LDA, LWORK, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ), &
                               WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE
            PARAMETER          ( ONE = 1.0D+0 )
        !     ..
        !     .. Local Scalars ..
            LOGICAL            LQUERY, UPPER
            INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB, &
                               NBMIN, NX
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DLATRD, DSYR2K, DSYTD2, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MAX
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            INTEGER            ILAENV
            EXTERNAL           LSAME, ILAENV
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters
        !
            INFO = 0
            UPPER = LSAME( UPLO, 'U' )
            LQUERY = ( LWORK == -1 )
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
                INFO = -1
            ELSE IF( N < 0 ) THEN
                INFO = -2
            ELSE IF( LDA < MAX( 1, N ) ) THEN
                INFO = -4
            ELSE IF( LWORK < 1 .AND. .NOT.LQUERY ) THEN
                INFO = -9
            END IF
        !
            IF( INFO == 0 ) THEN
        !
        !        Determine the block size.
        !
                NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
                LWKOPT = N*NB
                WORK( 1 ) = LWKOPT
            END IF
        !
            IF( INFO /= 0 ) THEN
                CALL XERBLA( 'DSYTRD', -INFO )
                RETURN
            ELSE IF( LQUERY ) THEN
                RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( N == 0 ) THEN
                WORK( 1 ) = 1
                RETURN
            END IF
        !
            NX = N
            IWS = 1
            IF( NB > 1 .AND. NB < N ) THEN
        !
        !        Determine when to cross over from blocked to unblocked code
        !        (last block is always handled by unblocked code).
        !
                NX = MAX( NB, ILAENV( 3, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
                IF( NX < N ) THEN
        !
        !           Determine if workspace is large enough for blocked code.
        !
                    LDWORK = N
                    IWS = LDWORK*NB
                    IF( LWORK < IWS ) THEN
        !
        !              Not enough workspace to use optimal NB:  determine the
        !              minimum value of NB, and reduce NB or force use of
        !              unblocked code by setting NX = N.
        !
                    NB = MAX( LWORK / LDWORK, 1 )
                    NBMIN = ILAENV( 2, 'DSYTRD', UPLO, N, -1, -1, -1 )
                    IF( NB < NBMIN ) &
                        NX = N
                    END IF
                ELSE
                    NX = N
                END IF
            ELSE
                NB = 1
            END IF
        !
            IF( UPPER ) THEN
        !
        !        Reduce the upper triangle of A.
        !        Columns 1:kk are handled by the unblocked method.
        !
                KK = N - ( ( N-NX+NB-1 ) / NB )*NB
                DO 20 I = N - NB + 1, KK + 1, -NB
        !
        !           Reduce columns i:i+nb-1 to tridiagonal form and form the
        !           matrix W which is needed to update the unreduced part of
        !           the matrix
        !
                    CALL DLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK, &
                                LDWORK )
        !
        !           Update the unreduced submatrix A(1:i-1,1:i-1), using an
        !           update of the form:  A := A - V*W**T - W*V**T
        !
                    CALL DSYR2K( UPLO, 'No transpose', I-1, NB, -ONE, A( 1, I ), &
                                LDA, WORK, LDWORK, ONE, A, LDA )
        !
        !           Copy superdiagonal elements back into A, and diagonal
        !           elements into D
        !
                    DO 10 J = I, I + NB - 1
                    A( J-1, J ) = E( J-1 )
                    D( J ) = A( J, J )
        10       CONTINUE
        20    CONTINUE
        !
        !        Use unblocked code to reduce the last or only block
        !
                CALL DSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
            ELSE
        !
        !        Reduce the lower triangle of A
        !
                DO 40 I = 1, N - NX, NB
        !
        !           Reduce columns i:i+nb-1 to tridiagonal form and form the
        !           matrix W which is needed to update the unreduced part of
        !           the matrix
        !
                    CALL DLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ), &
                                TAU( I ), WORK, LDWORK )
        !
        !           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
        !           an update of the form:  A := A - V*W**T - W*V**T
        !
                    CALL DSYR2K( UPLO, 'No transpose', N-I-NB+1, NB, -ONE, &
                                A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE, &
                                A( I+NB, I+NB ), LDA )
        !
        !           Copy subdiagonal elements back into A, and diagonal
        !           elements into D
        !
                    DO 30 J = I, I + NB - 1
                    A( J+1, J ) = E( J )
                    D( J ) = A( J, J )
        30       CONTINUE
        40    CONTINUE
        !
        !        Use unblocked code to reduce the last or only block
        !
                CALL DSYTD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ), &
                            TAU( I ), IINFO )
            END IF
        !
            WORK( 1 ) = LWKOPT
            RETURN
        !
        !     End of DSYTRD
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DSTERF( N, D, E, INFO )
    !==================================================================================================
        !
        !  -- LAPACK computational routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER            INFO, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   D( * ), E( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ZERO, ONE, TWO, THREE
            PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                               THREE = 3.0D0 )
            INTEGER            MAXIT
            PARAMETER          ( MAXIT = 30 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M, &
                               NMAXIT
            DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC, &
                               OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN, &
                               SIGMA, SSFMAX, SSFMIN, RMAX
        !     ..
        !     .. External Functions ..
            DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
            EXTERNAL           DLAMCH, DLANST, DLAPY2
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          ABS, SIGN, SQRT
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters.
        !
            INFO = 0
        !
        !     Quick return if possible
        !
            IF( N < 0 ) THEN
                INFO = -1
                CALL XERBLA( 'DSTERF', -INFO )
                RETURN
            END IF
            IF( N <= 1 ) &
                RETURN
        !
        !     Determine the unit roundoff for this environment.
        !
            EPS = DLAMCH( 'E' )
            EPS2 = EPS**2
            SAFMIN = DLAMCH( 'S' )
            SAFMAX = ONE / SAFMIN
            SSFMAX = SQRT( SAFMAX ) / THREE
            SSFMIN = SQRT( SAFMIN ) / EPS2
            RMAX = DLAMCH( 'O' )
        !
        !     Compute the eigenvalues of the tridiagonal matrix.
        !
            NMAXIT = N*MAXIT
            SIGMA = ZERO
            JTOT = 0
        !
        !     Determine where the matrix splits and choose QL or QR iteration
        !     for each block, according to whether top or bottom diagonal
        !     element is smaller.
        !
            L1 = 1
        !
        10 CONTINUE
            IF( L1 > N ) &
                GO TO 170
            IF( L1 > 1 ) &
                E( L1-1 ) = ZERO
            DO 20 M = L1, N - 1
                IF( ABS( E( M ) ) <= ( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ &
                    1 ) ) ) )*EPS ) THEN
                    E( M ) = ZERO
                    GO TO 30
                END IF
        20 CONTINUE
            M = N
        !
        30 CONTINUE
            L = L1
            LSV = L
            LEND = M
            LENDSV = LEND
            L1 = M + 1
            IF( LEND == L ) &
                GO TO 10
        !
        !     Scale submatrix in rows and columns L to LEND
        !
            ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) )
            ISCALE = 0
            IF( ANORM == ZERO ) &
                GO TO 10
            IF( (ANORM > SSFMAX) ) THEN
                ISCALE = 1
                CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, &
                            INFO )
                CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, &
                            INFO )
            ELSE IF( ANORM < SSFMIN ) THEN
                ISCALE = 2
                CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, &
                            INFO )
                CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, &
                            INFO )
            END IF
        !
            DO 40 I = L, LEND - 1
                E( I ) = E( I )**2
        40 CONTINUE
        !
        !     Choose between QL and QR iteration
        !
            IF( ABS( D( LEND ) ) < ABS( D( L ) ) ) THEN
                LEND = LSV
                L = LENDSV
            END IF
        !
            IF( LEND >= L ) THEN
        !
        !        QL Iteration
        !
        !        Look for small subdiagonal element.
        !
        50    CONTINUE
                IF( L /= LEND ) THEN
                    DO 60 M = L, LEND - 1
                    IF( ABS( E( M ) ) <= EPS2*ABS( D( M )*D( M+1 ) ) ) &
                        GO TO 70
        60       CONTINUE
                END IF
                M = LEND
        !
        70    CONTINUE
                IF( M < LEND ) &
                    E( M ) = ZERO
                P = D( L )
                IF( M == L ) &
                    GO TO 90
        !
        !        If remaining matrix is 2 by 2, use DLAE2 to compute its
        !        eigenvalues.
        !
                IF( M == L+1 ) THEN
                    RTE = SQRT( E( L ) )
                    CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
                    D( L ) = RT1
                    D( L+1 ) = RT2
                    E( L ) = ZERO
                    L = L + 2
                    IF( L <= LEND ) &
                        GO TO 50
                    GO TO 150
                END IF
        !
                IF( JTOT == NMAXIT ) &
                    GO TO 150
                JTOT = JTOT + 1
        !
        !        Form shift.
        !
                RTE = SQRT( E( L ) )
                SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
                R = DLAPY2( SIGMA, ONE )
                SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
        !
                C = ONE
                S = ZERO
                GAMMA = D( M ) - SIGMA
                P = GAMMA*GAMMA
        !
        !        Inner loop
        !
                DO 80 I = M - 1, L, -1
                    BB = E( I )
                    R = P + BB
                    IF( I /= M-1 ) &
                        E( I+1 ) = S*R
                    OLDC = C
                    C = P / R
                    S = BB / R
                    OLDGAM = GAMMA
                    ALPHA = D( I )
                    GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
                    D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
                    IF( C /= ZERO ) THEN
                    P = ( GAMMA*GAMMA ) / C
                    ELSE
                    P = OLDC*BB
                    END IF
        80    CONTINUE
        !
                E( L ) = S*P
                D( L ) = SIGMA + GAMMA
                GO TO 50
        !
        !        Eigenvalue found.
        !
        90    CONTINUE
                D( L ) = P
        !
                L = L + 1
                IF( L <= LEND ) &
                    GO TO 50
                GO TO 150
        !
            ELSE
        !
        !        QR Iteration
        !
        !        Look for small superdiagonal element.
        !
        100    CONTINUE
                DO 110 M = L, LEND + 1, -1
                    IF( ABS( E( M-1 ) ) <= EPS2*ABS( D( M )*D( M-1 ) ) ) &
                        GO TO 120
        110    CONTINUE
                M = LEND
        !
        120    CONTINUE
                IF( M > LEND ) &
                    E( M-1 ) = ZERO
                P = D( L )
                IF( M == L ) &
                    GO TO 140
        !
        !        If remaining matrix is 2 by 2, use DLAE2 to compute its
        !        eigenvalues.
        !
                IF( M == L-1 ) THEN
                    RTE = SQRT( E( L-1 ) )
                    CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
                    D( L ) = RT1
                    D( L-1 ) = RT2
                    E( L-1 ) = ZERO
                    L = L - 2
                    IF( L >= LEND ) &
                        GO TO 100
                    GO TO 150
                END IF
        !
                IF( JTOT == NMAXIT ) &
                    GO TO 150
                JTOT = JTOT + 1
        !
        !        Form shift.
        !
                RTE = SQRT( E( L-1 ) )
                SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
                R = DLAPY2( SIGMA, ONE )
                SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
        !
                C = ONE
                S = ZERO
                GAMMA = D( M ) - SIGMA
                P = GAMMA*GAMMA
        !
        !        Inner loop
        !
                DO 130 I = M, L - 1
                    BB = E( I )
                    R = P + BB
                    IF( I /= M ) &
                        E( I-1 ) = S*R
                    OLDC = C
                    C = P / R
                    S = BB / R
                    OLDGAM = GAMMA
                    ALPHA = D( I+1 )
                    GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
                    D( I ) = OLDGAM + ( ALPHA-GAMMA )
                    IF( C /= ZERO ) THEN
                    P = ( GAMMA*GAMMA ) / C
                    ELSE
                    P = OLDC*BB
                    END IF
        130    CONTINUE
        !
                E( L-1 ) = S*P
                D( L ) = SIGMA + GAMMA
                GO TO 100
        !
        !        Eigenvalue found.
        !
        140    CONTINUE
                D( L ) = P
        !
                L = L - 1
                IF( L >= LEND ) &
                    GO TO 100
                GO TO 150
        !
            END IF
        !
        !     Undo scaling if necessary
        !
        150 CONTINUE
            IF( ISCALE == 1 ) &
                CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, &
                            D( LSV ), N, INFO )
            IF( ISCALE == 2 ) &
                CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, &
                            D( LSV ), N, INFO )
        !
        !     Check for no convergence to an eigenvalue after a total
        !     of N*MAXIT iterations.
        !
            IF( JTOT < NMAXIT ) &
                GO TO 10
            DO 160 I = 1, N - 1
                IF( E( I ) /= ZERO ) &
                    INFO = INFO + 1
        160 CONTINUE
            GO TO 180
        !
        !     Sort eigenvalues in increasing order.
        !
        170 CONTINUE
            CALL DLASRT( 'I', N, D, INFO )
        !
        180 CONTINUE
            RETURN
        !
        !     End of DSTERF
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
    !==================================================================================================
        !
        !  -- LAPACK computational routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          UPLO
            INTEGER            INFO, LDA, LWORK, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ZERO, ONE
            PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
        !     ..
        !     .. Local Scalars ..
            LOGICAL            LQUERY, UPPER
            INTEGER            I, IINFO, J, LWKOPT, NB
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            INTEGER            ILAENV
            EXTERNAL           LSAME, ILAENV
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DORGQL, DORGQR, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MAX
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
            INFO = 0
            LQUERY = ( LWORK == -1 )
            UPPER = LSAME( UPLO, 'U' )
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
                INFO = -1
            ELSE IF( N < 0 ) THEN
                INFO = -2
            ELSE IF( LDA < MAX( 1, N ) ) THEN
                INFO = -4
            ELSE IF( LWORK < MAX( 1, N-1 ) .AND. .NOT.LQUERY ) THEN
                INFO = -7
            END IF
        !
            IF( INFO == 0 ) THEN
                IF( UPPER ) THEN
                    NB = ILAENV( 1, 'DORGQL', ' ', N-1, N-1, N-1, -1 )
                ELSE
                    NB = ILAENV( 1, 'DORGQR', ' ', N-1, N-1, N-1, -1 )
                END IF
                LWKOPT = MAX( 1, N-1 )*NB
                WORK( 1 ) = LWKOPT
            END IF
        !
            IF( INFO /= 0 ) THEN
                CALL XERBLA( 'DORGTR', -INFO )
                RETURN
            ELSE IF( LQUERY ) THEN
                RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( N == 0 ) THEN
                WORK( 1 ) = 1
                RETURN
            END IF
        !
            IF( UPPER ) THEN
        !
        !        Q was determined by a call to DSYTRD with UPLO = 'U'
        !
        !        Shift the vectors which define the elementary reflectors one
        !        column to the left, and set the last row and column of Q to
        !        those of the unit matrix
        !
                DO 20 J = 1, N - 1
                    DO 10 I = 1, J - 1
                    A( I, J ) = A( I, J+1 )
        10       CONTINUE
                    A( N, J ) = ZERO
        20    CONTINUE
                DO 30 I = 1, N - 1
                    A( I, N ) = ZERO
        30    CONTINUE
                A( N, N ) = ONE
        !
        !        Generate Q(1:n-1,1:n-1)
        !
                CALL DORGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
        !
            ELSE
        !
        !        Q was determined by a call to DSYTRD with UPLO = 'L'.
        !
        !        Shift the vectors which define the elementary reflectors one
        !        column to the right, and set the first row and column of Q to
        !        those of the unit matrix
        !
                DO 50 J = N, 2, -1
                    A( 1, J ) = ZERO
                    DO 40 I = J + 1, N
                    A( I, J ) = A( I, J-1 )
        40       CONTINUE
        50    CONTINUE
                A( 1, 1 ) = ONE
                DO 60 I = 2, N
                    A( I, 1 ) = ZERO
        60    CONTINUE
                IF( N > 1 ) THEN
        !
        !           Generate Q(2:n,2:n)
        !
                    CALL DORGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, &
                                LWORK, IINFO )
                END IF
            END IF
            WORK( 1 ) = LWKOPT
            RETURN
        !
        !     End of DORGTR
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
    !==================================================================================================
        !
        !  -- LAPACK computational routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          COMPZ
            INTEGER            INFO, LDZ, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ZERO, ONE, TWO, THREE
            PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                               THREE = 3.0D0 )
            INTEGER            MAXIT
            PARAMETER          ( MAXIT = 30 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, &
                               LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1, &
                               NM1, NMAXIT
            DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2, &
                               S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
            EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASET, DLASR, &
                               DLASRT, DSWAP, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          ABS, MAX, SIGN, SQRT
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters.
        !
            INFO = 0
        !
            IF( LSAME( COMPZ, 'N' ) ) THEN
                ICOMPZ = 0
            ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
                ICOMPZ = 1
            ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
                ICOMPZ = 2
            ELSE
                ICOMPZ = -1
            END IF
            IF( ICOMPZ < 0 ) THEN
                INFO = -1
            ELSE IF( N < 0 ) THEN
                INFO = -2
            ELSE IF( ( LDZ < 1 ) .OR. ( ICOMPZ > 0 .AND. LDZ < MAX( 1, &
                    N ) ) ) THEN
                INFO = -6
            END IF
            IF( INFO /= 0 ) THEN
                CALL XERBLA( 'DSTEQR', -INFO )
                RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( N == 0 ) &
                RETURN
        !
            IF( N == 1 ) THEN
                IF( ICOMPZ == 2 ) &
                    Z( 1, 1 ) = ONE
                RETURN
            END IF
        !
        !     Determine the unit roundoff and over/underflow thresholds.
        !
            EPS = DLAMCH( 'E' )
            EPS2 = EPS**2
            SAFMIN = DLAMCH( 'S' )
            SAFMAX = ONE / SAFMIN
            SSFMAX = SQRT( SAFMAX ) / THREE
            SSFMIN = SQRT( SAFMIN ) / EPS2
        !
        !     Compute the eigenvalues and eigenvectors of the tridiagonal
        !     matrix.
        !
            IF( ICOMPZ == 2 ) &
                CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
        !
            NMAXIT = N*MAXIT
            JTOT = 0
        !
        !     Determine where the matrix splits and choose QL or QR iteration
        !     for each block, according to whether top or bottom diagonal
        !     element is smaller.
        !
            L1 = 1
            NM1 = N - 1
        !
        10 CONTINUE
            IF( L1 > N ) &
                GO TO 160
            IF( L1 > 1 ) &
                E( L1-1 ) = ZERO
            IF( L1 <= NM1 ) THEN
                DO 20 M = L1, NM1
                    TST = ABS( E( M ) )
                    IF( TST == ZERO ) &
                    GO TO 30
                    IF( TST <= ( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ &
                        1 ) ) ) )*EPS ) THEN
                    E( M ) = ZERO
                    GO TO 30
                    END IF
        20    CONTINUE
            END IF
            M = N
        !
        30 CONTINUE
            L = L1
            LSV = L
            LEND = M
            LENDSV = LEND
            L1 = M + 1
            IF( LEND == L ) &
                GO TO 10
        !
        !     Scale submatrix in rows and columns L to LEND
        !
            ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) )
            ISCALE = 0
            IF( ANORM == ZERO ) &
                GO TO 10
            IF( ANORM > SSFMAX ) THEN
                ISCALE = 1
                CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, &
                            INFO )
                CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, &
                            INFO )
            ELSE IF( ANORM < SSFMIN ) THEN
                ISCALE = 2
                CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, &
                            INFO )
                CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, &
                            INFO )
            END IF
        !
        !     Choose between QL and QR iteration
        !
            IF( ABS( D( LEND ) ) < ABS( D( L ) ) ) THEN
                LEND = LSV
                L = LENDSV
            END IF
        !
            IF( LEND > L ) THEN
        !
        !        QL Iteration
        !
        !        Look for small subdiagonal element.
        !
        40    CONTINUE
                IF( L /= LEND ) THEN
                    LENDM1 = LEND - 1
                    DO 50 M = L, LENDM1
                    TST = ABS( E( M ) )**2
                    IF( TST <= ( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+ &
                        SAFMIN )GO TO 60
        50       CONTINUE
                END IF
        !
                M = LEND
        !
        60    CONTINUE
                IF( M < LEND ) &
                    E( M ) = ZERO
                P = D( L )
                IF( M == L ) &
                    GO TO 80
        !
        !        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
        !        to compute its eigensystem.
        !
                IF( M == L+1 ) THEN
                    IF( ICOMPZ > 0 ) THEN
                    CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
                    WORK( L ) = C
                    WORK( N-1+L ) = S
                    CALL DLASR( 'R', 'V', 'B', N, 2, WORK( L ), &
                                WORK( N-1+L ), Z( 1, L ), LDZ )
                    ELSE
                    CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
                    END IF
                    D( L ) = RT1
                    D( L+1 ) = RT2
                    E( L ) = ZERO
                    L = L + 2
                    IF( L <= LEND ) &
                        GO TO 40
                    GO TO 140
                END IF
        !
                IF( JTOT == NMAXIT ) &
                    GO TO 140
                JTOT = JTOT + 1
        !
        !        Form shift.
        !
                G = ( D( L+1 )-P ) / ( TWO*E( L ) )
                R = DLAPY2( G, ONE )
                G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
        !
                S = ONE
                C = ONE
                P = ZERO
        !
        !        Inner loop
        !
                MM1 = M - 1
                DO 70 I = MM1, L, -1
                    F = S*E( I )
                    B = C*E( I )
                    CALL DLARTG( G, F, C, S, R )
                    IF( I /= M-1 ) &
                        E( I+1 ) = R
                    G = D( I+1 ) - P
                    R = ( D( I )-G )*S + TWO*C*B
                    P = S*R
                    D( I+1 ) = G + P
                    G = C*R - B
        !
        !           If eigenvectors are desired, then save rotations.
        !
                    IF( ICOMPZ > 0 ) THEN
                    WORK( I ) = C
                    WORK( N-1+I ) = -S
                    END IF
        !
        70    CONTINUE
        !
        !        If eigenvectors are desired, then apply saved rotations.
        !
                IF( ICOMPZ > 0 ) THEN
                    MM = M - L + 1
                    CALL DLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ), &
                                Z( 1, L ), LDZ )
                END IF
        !
                D( L ) = D( L ) - P
                E( L ) = G
                GO TO 40
        !
        !        Eigenvalue found.
        !
        80    CONTINUE
                D( L ) = P
        !
                L = L + 1
                IF( L <= LEND ) &
                    GO TO 40
                GO TO 140
        !
            ELSE
        !
        !        QR Iteration
        !
        !        Look for small superdiagonal element.
        !
        90    CONTINUE
                IF( L /= LEND ) THEN
                    LENDP1 = LEND + 1
                    DO 100 M = L, LENDP1, -1
                    TST = ABS( E( M-1 ) )**2
                    IF( TST <= ( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+ &
                        SAFMIN )GO TO 110
        100       CONTINUE
                END IF
        !
                M = LEND
        !
        110    CONTINUE
                IF( M > LEND ) &
                    E( M-1 ) = ZERO
                P = D( L )
                IF( M == L ) &
                    GO TO 130
        !
        !        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
        !        to compute its eigensystem.
        !
                IF( M == L-1 ) THEN
                    IF( ICOMPZ > 0 ) THEN
                    CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
                    WORK( M ) = C
                    WORK( N-1+M ) = S
                    CALL DLASR( 'R', 'V', 'F', N, 2, WORK( M ), &
                                WORK( N-1+M ), Z( 1, L-1 ), LDZ )
                    ELSE
                    CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
                    END IF
                    D( L-1 ) = RT1
                    D( L ) = RT2
                    E( L-1 ) = ZERO
                    L = L - 2
                    IF( L >= LEND ) &
                        GO TO 90
                    GO TO 140
                END IF
        !
                IF( JTOT == NMAXIT ) &
                    GO TO 140
                JTOT = JTOT + 1
        !
        !        Form shift.
        !
                G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
                R = DLAPY2( G, ONE )
                G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
        !
                S = ONE
                C = ONE
                P = ZERO
        !
        !        Inner loop
        !
                LM1 = L - 1
                DO 120 I = M, LM1
                    F = S*E( I )
                    B = C*E( I )
                    CALL DLARTG( G, F, C, S, R )
                    IF( I /= M ) &
                        E( I-1 ) = R
                    G = D( I ) - P
                    R = ( D( I+1 )-G )*S + TWO*C*B
                    P = S*R
                    D( I ) = G + P
                    G = C*R - B
        !
        !           If eigenvectors are desired, then save rotations.
        !
                    IF( ICOMPZ > 0 ) THEN
                    WORK( I ) = C
                    WORK( N-1+I ) = S
                    END IF
        !
        120    CONTINUE
        !
        !        If eigenvectors are desired, then apply saved rotations.
        !
                IF( ICOMPZ > 0 ) THEN
                    MM = L - M + 1
                    CALL DLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ), &
                                Z( 1, M ), LDZ )
                END IF
        !
                D( L ) = D( L ) - P
                E( LM1 ) = G
                GO TO 90
        !
        !        Eigenvalue found.
        !
        130    CONTINUE
                D( L ) = P
        !
                L = L - 1
                IF( L >= LEND ) &
                    GO TO 90
                GO TO 140
        !
            END IF
        !
        !     Undo scaling if necessary
        !
        140 CONTINUE
            IF( ISCALE == 1 ) THEN
                CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, &
                            D( LSV ), N, INFO )
                CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ), &
                            N, INFO )
            ELSE IF( ISCALE == 2 ) THEN
                CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, &
                            D( LSV ), N, INFO )
                CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ), &
                            N, INFO )
            END IF
        !
        !     Check for no convergence to an eigenvalue after a total
        !     of N*MAXIT iterations.
        !
            IF( JTOT < NMAXIT ) &
                GO TO 10
            DO 150 I = 1, N - 1
                IF( E( I ) /= ZERO ) &
                    INFO = INFO + 1
        150 CONTINUE
            GO TO 190
        !
        !     Order eigenvalues and eigenvectors.
        !
        160 CONTINUE
            IF( ICOMPZ == 0 ) THEN
        !
        !        Use Quick Sort
        !
                CALL DLASRT( 'I', N, D, INFO )
        !
            ELSE
        !
        !        Use Selection Sort to minimize swaps of eigenvectors
        !
                DO 180 II = 2, N
                    I = II - 1
                    K = I
                    P = D( I )
                    DO 170 J = II, N
                    IF( D( J ) < P ) THEN
                        K = J
                        P = D( J )
                    END IF
        170       CONTINUE
                    IF( K /= I ) THEN
                    D( K ) = D( I )
                    D( I ) = P
                    CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
                    END IF
        180    CONTINUE
            END IF
        !
        190 CONTINUE
            RETURN
        !
        !     End of DSTEQR
        !
        END SUBROUTINE
    
    !==================================================================================================
        INTEGER FUNCTION IEEECK( ISPEC, ZERO, ONE )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER            ISPEC
            REAL               ONE, ZERO
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
            REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
                               NEGZRO, NEWZRO, POSINF
        !     ..
        !     .. Executable Statements ..
            IEEECK = 1
        !
            POSINF = ONE / ZERO
            IF( POSINF <= ONE ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            NEGINF = -ONE / ZERO
            IF( NEGINF >= ZERO ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            NEGZRO = ONE / ( NEGINF+ONE )
            IF( NEGZRO /= ZERO ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            NEGINF = ONE / NEGZRO
            IF( NEGINF >= ZERO ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            NEWZRO = NEGZRO + ZERO
            IF( NEWZRO /= ZERO ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            POSINF = ONE / NEWZRO
            IF( POSINF <= ONE ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            NEGINF = NEGINF*POSINF
            IF( NEGINF >= ZERO ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            POSINF = POSINF*POSINF
            IF( POSINF <= ONE ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
        !
        !
        !
        !     Return if we were only asked to check infinity arithmetic
        !
            IF( ISPEC == 0 ) &
                RETURN
        !
            NAN1 = POSINF + NEGINF
        !
            NAN2 = POSINF / NEGINF
        !
            NAN3 = POSINF / POSINF
        !
            NAN4 = POSINF*ZERO
        !
            NAN5 = NEGINF*NEGZRO
        !
            NAN6 = NAN5*ZERO
        !
            IF( NAN1 == NAN1 ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            IF( NAN2 == NAN2 ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            IF( NAN3 == NAN3 ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            IF( NAN4 == NAN4 ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            IF( NAN5 == NAN5 ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            IF( NAN6 == NAN6 ) THEN
                IEEECK = 0
                RETURN
            END IF
        !
            RETURN
        END FUNCTION
    
    !==================================================================================================
        INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER            IHI, ILO, ISPEC, LWORK, N
            CHARACTER          NAME*( * ), OPTS*( * )
        !
        !  ================================================================
        !     .. Parameters ..
            INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22, ICOST
            PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14, &
                               ISHFTS = 15, IACC22 = 16, ICOST = 17 )
            INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP, RCOST
            PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14, &
                               NIBBLE = 14, KNWSWP = 500, RCOST = 10 )
            REAL               TWO
            PARAMETER          ( TWO = 2.0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            NH, NS
            INTEGER            I, IC, IZ
            CHARACTER          SUBNAM*6
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          LOG, MAX, MOD, NINT, REAL
        !     ..
        !     .. Executable Statements ..
            IF( ( ISPEC == ISHFTS ) .OR. ( ISPEC == INWIN ) .OR. &
                ( ISPEC == IACC22 ) ) THEN
        !
        !        ==== Set the number simultaneous shifts ====
        !
                NH = IHI - ILO + 1
                NS = 2
                IF( NH >= 30 ) &
                    NS = 4
                IF( NH >= 60 ) &
                    NS = 10
                IF( NH >= 150 ) &
                    NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
                IF( NH >= 590 ) &
                    NS = 64
                IF( NH >= 3000 ) &
                    NS = 128
                IF( NH >= 6000 ) &
                    NS = 256
                NS = MAX( 2, NS-MOD( NS, 2 ) )
            END IF
        !
            IF( ISPEC == INMIN ) THEN
        !
        !
        !        ===== Matrices of order smaller than NMIN get sent
        !        .     to xLAHQR, the classic double shift algorithm.
        !        .     This must be at least 11. ====
        !
                IPARMQ = NMIN
        !
            ELSE IF( ISPEC == INIBL ) THEN
        !
        !        ==== INIBL: skip a multi-shift qr iteration and
        !        .    whenever aggressive early deflation finds
        !        .    at least (NIBBLE*(window size)/100) deflations. ====
        !
                IPARMQ = NIBBLE
        !
            ELSE IF( ISPEC == ISHFTS ) THEN
        !
        !        ==== NSHFTS: The number of simultaneous shifts =====
        !
                IPARMQ = NS
        !
            ELSE IF( ISPEC == INWIN ) THEN
        !
        !        ==== NW: deflation window size.  ====
        !
                IF( NH <= KNWSWP ) THEN
                    IPARMQ = NS
                ELSE
                    IPARMQ = 3*NS / 2
                END IF
        !
            ELSE IF( ISPEC == IACC22 ) THEN
        !
        !        ==== IACC22: Whether to accumulate reflections
        !        .     before updating the far-from-diagonal elements
        !        .     and whether to use 2-by-2 block structure while
        !        .     doing it.  A small amount of work could be saved
        !        .     by making this choice dependent also upon the
        !        .     NH=IHI-ILO+1.
        !
        !
        !        Convert NAME to upper case if the first character is lower case.
        !
                IPARMQ = 0
                SUBNAM = NAME
                IC = ICHAR( SUBNAM( 1: 1 ) )
                IZ = ICHAR( 'Z' )
                IF( IZ == 90 .OR. IZ == 122 ) THEN
        !
        !           ASCII character set
        !
                    IF( IC >= 97 .AND. IC <= 122 ) THEN
                    SUBNAM( 1: 1 ) = CHAR( IC-32 )
                    DO I = 2, 6
                        IC = ICHAR( SUBNAM( I: I ) )
                        IF( IC >= 97 .AND. IC <= 122 ) &
                            SUBNAM( I: I ) = CHAR( IC-32 )
                    END DO
                    END IF
        !
                ELSE IF( IZ == 233 .OR. IZ == 169 ) THEN
        !
        !           EBCDIC character set
        !
                    IF( ( IC >= 129 .AND. IC <= 137 ) .OR. &
                        ( IC >= 145 .AND. IC <= 153 ) .OR. &
                        ( IC >= 162 .AND. IC <= 169 ) ) THEN
                    SUBNAM( 1: 1 ) = CHAR( IC+64 )
                    DO I = 2, 6
                        IC = ICHAR( SUBNAM( I: I ) )
                        IF( ( IC >= 129 .AND. IC <= 137 ) .OR. &
                            ( IC >= 145 .AND. IC <= 153 ) .OR. &
                            ( IC >= 162 .AND. IC <= 169 ) )SUBNAM( I:&
                            I ) = CHAR( IC+64 )
                    END DO
                    END IF
        !
                ELSE IF( IZ == 218 .OR. IZ == 250 ) THEN
        !
        !           Prime machines:  ASCII+128
        !
                    IF( IC >= 225 .AND. IC <= 250 ) THEN
                    SUBNAM( 1: 1 ) = CHAR( IC-32 )
                    DO I = 2, 6
                        IC = ICHAR( SUBNAM( I: I ) )
                        IF( IC >= 225 .AND. IC <= 250 ) &
                            SUBNAM( I: I ) = CHAR( IC-32 )
                    END DO
                    END IF
                END IF
        !
                IF( SUBNAM( 2:6 ) == 'GGHRD' .OR. &
                    SUBNAM( 2:6 ) == 'GGHD3' ) THEN
                    IPARMQ = 1
                    IF( NH >= K22MIN ) &
                        IPARMQ = 2
                ELSE IF ( SUBNAM( 4:6 ) == 'EXC' ) THEN
                    IF( NH >= KACMIN ) &
                        IPARMQ = 1
                    IF( NH >= K22MIN ) &
                        IPARMQ = 2
                ELSE IF ( SUBNAM( 2:6 ) == 'HSEQR' .OR. &
                        SUBNAM( 2:5 ) == 'LAQR' ) THEN
                    IF( NS >= KACMIN ) &
                        IPARMQ = 1
                    IF( NS >= K22MIN ) &
                        IPARMQ = 2
                END IF
        !
            ELSE IF( ISPEC == ICOST ) THEN
        !
        !        === Relative cost of near-the-diagonal chase vs
        !            BLAS updates ===
        !
                IPARMQ = RCOST
            ELSE
        !        ===== invalid value of ispec =====
                IPARMQ = -1
        !
            END IF
        !
        !     ==== End of IPARMQ ====
        !
        END FUNCTION
    
    !==================================================================================================
        SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
    !==================================================================================================
        !
        !  -- LAPACK computational routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          DIAG, UPLO
            INTEGER            INFO, LDA, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE, ZERO
            PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            LOGICAL            NOUNIT, UPPER
            INTEGER            J, JB, NB, NN
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            INTEGER            ILAENV
            EXTERNAL           LSAME, ILAENV
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MAX, MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters.
        !
            INFO = 0
            UPPER = LSAME( UPLO, 'U' )
            NOUNIT = LSAME( DIAG, 'N' )
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
                INFO = -1
            ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
                INFO = -2
            ELSE IF( N < 0 ) THEN
                INFO = -3
            ELSE IF( LDA < MAX( 1, N ) ) THEN
                INFO = -5
            END IF
            IF( INFO /= 0 ) THEN
                CALL XERBLA( 'DTRTRI', -INFO )
                RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( N == 0 ) &
                RETURN
        !
        !     Check for singularity if non-unit.
        !
            IF( NOUNIT ) THEN
                DO 10 INFO = 1, N
                    IF( A( INFO, INFO ) == ZERO ) &
                       RETURN
        10    CONTINUE
                INFO = 0
            END IF
        !
        !     Determine the block size for this environment.
        !
            NB = ILAENV( 1, 'DTRTRI', UPLO // DIAG, N, -1, -1, -1 )
            IF( NB <= 1 .OR. NB >= N ) THEN
        !
        !        Use unblocked code
        !
                CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
            ELSE
        !
        !        Use blocked code
        !
                IF( UPPER ) THEN
        !
        !           Compute inverse of upper triangular matrix
        !
                    DO 20 J = 1, N, NB
                    JB = MIN( NB, N-J+1 )
        !
        !              Compute rows 1:j-1 of current block column
        !
                    CALL DTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1, &
                                JB, ONE, A, LDA, A( 1, J ), LDA )
                    CALL DTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1, &
                                JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
        !
        !              Compute inverse of current diagonal block
        !
                    CALL DTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
        20       CONTINUE
                ELSE
        !
        !           Compute inverse of lower triangular matrix
        !
                    NN = ( ( N-1 ) / NB )*NB + 1
                    DO 30 J = NN, 1, -NB
                    JB = MIN( NB, N-J+1 )
                    IF( J+JB <= N ) THEN
        !
        !                 Compute rows j+jb:n of current block column
        !
                        CALL DTRMM( 'Left', 'Lower', 'No transpose', DIAG, &
                                    N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA, &
                                    A( J+JB, J ), LDA )
                        CALL DTRSM( 'Right', 'Lower', 'No transpose', DIAG, &
                                    N-J-JB+1, JB, -ONE, A( J, J ), LDA, &
                                    A( J+JB, J ), LDA )
                    END IF
        !
        !              Compute inverse of current diagonal block
        !
                    CALL DTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
        30       CONTINUE
                END IF
            END IF
        !
            RETURN
        !
        !     End of DTRTRI
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          UPLO
            INTEGER            LDA, M, N
            DOUBLE PRECISION   ALPHA, BETA
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * )
        !     ..
        !
        ! =====================================================================
        !
        !     .. Local Scalars ..
            INTEGER            I, J
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            EXTERNAL           LSAME
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MIN
        !     ..
        !     .. Executable Statements ..
        !
            IF( LSAME( UPLO, 'U' ) ) THEN
        !
        !        Set the strictly upper triangular or trapezoidal part of the
        !        array to ALPHA.
        !
                DO 20 J = 2, N
                    DO 10 I = 1, MIN( J-1, M )
                    A( I, J ) = ALPHA
        10       CONTINUE
        20    CONTINUE
        !
            ELSE IF( LSAME( UPLO, 'L' ) ) THEN
        !
        !        Set the strictly lower triangular or trapezoidal part of the
        !        array to ALPHA.
        !
                DO 40 J = 1, MIN( M, N )
                    DO 30 I = J + 1, M
                    A( I, J ) = ALPHA
        30       CONTINUE
        40    CONTINUE
        !
            ELSE
        !
        !        Set the leading m-by-n submatrix to ALPHA.
        !
                DO 60 J = 1, N
                    DO 50 I = 1, M
                    A( I, J ) = ALPHA
        50       CONTINUE
        60    CONTINUE
            END IF
        !
        !     Set the first min(M,N) diagonal elements to BETA.
        !
            DO 70 I = 1, MIN( M, N )
                A( I, I ) = BETA
        70 CONTINUE
        !
            RETURN
        !
        !     End of DLASET
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          DIRECT, PIVOT, SIDE
            INTEGER            LDA, M, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE, ZERO
            PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            I, INFO, J
            DOUBLE PRECISION   CTEMP, STEMP, TEMP
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            EXTERNAL           LSAME
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MAX
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters
        !
            INFO = 0
            IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
                INFO = 1
            ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT, &
                    'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
                INFO = 2
            ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) ) &
                    THEN
                INFO = 3
            ELSE IF( M < 0 ) THEN
                INFO = 4
            ELSE IF( N < 0 ) THEN
                INFO = 5
            ELSE IF( LDA < MAX( 1, M ) ) THEN
                INFO = 9
            END IF
            IF( INFO /= 0 ) THEN
                CALL XERBLA( 'DLASR ', INFO )
                RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( ( M == 0 ) .OR. ( N == 0 ) ) &
                RETURN
            IF( LSAME( SIDE, 'L' ) ) THEN
        !
        !        Form  P * A
        !
                IF( LSAME( PIVOT, 'V' ) ) THEN
                    IF( LSAME( DIRECT, 'F' ) ) THEN
                    DO 20 J = 1, M - 1
                        CTEMP = C( J )
                        STEMP = S( J )
                        IF( ( CTEMP /= ONE ) .OR. ( STEMP /= ZERO ) ) THEN
                            DO 10 I = 1, N
                                TEMP = A( J+1, I )
                                A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                                A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
        10                CONTINUE
                        END IF
        20          CONTINUE
                    ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
                    DO 40 J = M - 1, 1, -1
                        CTEMP = C( J )
                        STEMP = S( J )
                        IF( ( CTEMP /= ONE ) .OR. ( STEMP /= ZERO ) ) THEN
                            DO 30 I = 1, N
                                TEMP = A( J+1, I )
                                A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                                A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
        30                CONTINUE
                        END IF
        40          CONTINUE
                    END IF
                ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
                    IF( LSAME( DIRECT, 'F' ) ) THEN
                    DO 60 J = 2, M
                        CTEMP = C( J-1 )
                        STEMP = S( J-1 )
                        IF( ( CTEMP /= ONE ) .OR. ( STEMP /= ZERO ) ) THEN
                            DO 50 I = 1, N
                                TEMP = A( J, I )
                                A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                                A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
        50                CONTINUE
                        END IF
        60          CONTINUE
                    ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
                    DO 80 J = M, 2, -1
                        CTEMP = C( J-1 )
                        STEMP = S( J-1 )
                        IF( ( CTEMP /= ONE ) .OR. ( STEMP /= ZERO ) ) THEN
                            DO 70 I = 1, N
                                TEMP = A( J, I )
                                A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                                A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
        70                CONTINUE
                        END IF
        80          CONTINUE
                    END IF
                ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
                    IF( LSAME( DIRECT, 'F' ) ) THEN
                    DO 100 J = 1, M - 1
                        CTEMP = C( J )
                        STEMP = S( J )
                        IF( ( CTEMP /= ONE ) .OR. ( STEMP /= ZERO ) ) THEN
                            DO 90 I = 1, N
                                TEMP = A( J, I )
                                A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                                A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
        90                CONTINUE
                        END IF
        100          CONTINUE
                    ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
                    DO 120 J = M - 1, 1, -1
                        CTEMP = C( J )
                        STEMP = S( J )
                        IF( ( CTEMP /= ONE ) .OR. ( STEMP /= ZERO ) ) THEN
                            DO 110 I = 1, N
                                TEMP = A( J, I )
                                A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                                A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
        110                CONTINUE
                        END IF
        120          CONTINUE
                    END IF
                END IF
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
        !
        !        Form A * P**T
        !
                IF( LSAME( PIVOT, 'V' ) ) THEN
                    IF( LSAME( DIRECT, 'F' ) ) THEN
                    DO 140 J = 1, N - 1
                        CTEMP = C( J )
                        STEMP = S( J )
                        IF( ( CTEMP /= ONE ) .OR. ( STEMP /= ZERO ) ) THEN
                            DO 130 I = 1, M
                                TEMP = A( I, J+1 )
                                A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                                A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
        130                CONTINUE
                        END IF
        140          CONTINUE
                    ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
                    DO 160 J = N - 1, 1, -1
                        CTEMP = C( J )
                        STEMP = S( J )
                        IF( ( CTEMP /= ONE ) .OR. ( STEMP /= ZERO ) ) THEN
                            DO 150 I = 1, M
                                TEMP = A( I, J+1 )
                                A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                                A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
        150                CONTINUE
                        END IF
        160          CONTINUE
                    END IF
                ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
                    IF( LSAME( DIRECT, 'F' ) ) THEN
                    DO 180 J = 2, N
                        CTEMP = C( J-1 )
                        STEMP = S( J-1 )
                        IF( ( CTEMP /= ONE ) .OR. ( STEMP /= ZERO ) ) THEN
                            DO 170 I = 1, M
                                TEMP = A( I, J )
                                A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                                A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
        170                CONTINUE
                        END IF
        180          CONTINUE
                    ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
                    DO 200 J = N, 2, -1
                        CTEMP = C( J-1 )
                        STEMP = S( J-1 )
                        IF( ( CTEMP /= ONE ) .OR. ( STEMP /= ZERO ) ) THEN
                            DO 190 I = 1, M
                                TEMP = A( I, J )
                                A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                                A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
        190                CONTINUE
                        END IF
        200          CONTINUE
                    END IF
                ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
                    IF( LSAME( DIRECT, 'F' ) ) THEN
                    DO 220 J = 1, N - 1
                        CTEMP = C( J )
                        STEMP = S( J )
                        IF( ( CTEMP /= ONE ) .OR. ( STEMP /= ZERO ) ) THEN
                            DO 210 I = 1, M
                                TEMP = A( I, J )
                                A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                                A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
        210                CONTINUE
                        END IF
        220          CONTINUE
                    ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
                    DO 240 J = N - 1, 1, -1
                        CTEMP = C( J )
                        STEMP = S( J )
                        IF( ( CTEMP /= ONE ) .OR. ( STEMP /= ZERO ) ) THEN
                            DO 230 I = 1, M
                                TEMP = A( I, J )
                                A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                                A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
        230                CONTINUE
                        END IF
        240          CONTINUE
                    END IF
                END IF
            END IF
        !
            RETURN
        !
        !     End of DLASR
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DLASRT( ID, N, D, INFO )
    !==================================================================================================
        !
        !  -- LAPACK computational routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          ID
            INTEGER            INFO, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   D( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            INTEGER            SELECT
            PARAMETER          ( SELECT = 20 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            DIR, ENDD, I, J, START, STKPNT
            DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
        !     ..
        !     .. Local Arrays ..
            INTEGER            STACK( 2, 32 )
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            EXTERNAL           LSAME
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           XERBLA
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters.
        !
            INFO = 0
            DIR = -1
            IF( LSAME( ID, 'D' ) ) THEN
                DIR = 0
            ELSE IF( LSAME( ID, 'I' ) ) THEN
                DIR = 1
            END IF
            IF( DIR == -1 ) THEN
                INFO = -1
            ELSE IF( N < 0 ) THEN
                INFO = -2
            END IF
            IF( INFO /= 0 ) THEN
                CALL XERBLA( 'DLASRT', -INFO )
                RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( N <= 1 ) &
                RETURN
        !
            STKPNT = 1
            STACK( 1, 1 ) = 1
            STACK( 2, 1 ) = N
        10 CONTINUE
            START = STACK( 1, STKPNT )
            ENDD = STACK( 2, STKPNT )
            STKPNT = STKPNT - 1
            IF( ENDD-START <= SELECT .AND. ENDD-START > 0 ) THEN
        !
        !        Do Insertion sort on D( START:ENDD )
        !
                IF( DIR == 0 ) THEN
        !
        !           Sort into decreasing order
        !
                    DO 30 I = START + 1, ENDD
                    DO 20 J = I, START + 1, -1
                        IF( D( J ) > D( J-1 ) ) THEN
                            DMNMX = D( J )
                            D( J ) = D( J-1 )
                            D( J-1 ) = DMNMX
                        ELSE
                            GO TO 30
                        END IF
        20          CONTINUE
        30       CONTINUE
        !
                ELSE
        !
        !           Sort into increasing order
        !
                    DO 50 I = START + 1, ENDD
                    DO 40 J = I, START + 1, -1
                        IF( D( J ) < D( J-1 ) ) THEN
                            DMNMX = D( J )
                            D( J ) = D( J-1 )
                            D( J-1 ) = DMNMX
                        ELSE
                            GO TO 50
                        END IF
        40          CONTINUE
        50       CONTINUE
        !
                END IF
        !
            ELSE IF( ENDD-START > SELECT ) THEN
        !
        !        Partition D( START:ENDD ) and stack parts, largest one first
        !
        !        Choose partition entry as median of 3
        !
                D1 = D( START )
                D2 = D( ENDD )
                I = ( START+ENDD ) / 2
                D3 = D( I )
                IF( D1 < D2 ) THEN
                    IF( D3 < D1 ) THEN
                    DMNMX = D1
                    ELSE IF( D3 < D2 ) THEN
                    DMNMX = D3
                    ELSE
                    DMNMX = D2
                    END IF
                ELSE
                    IF( D3 < D2 ) THEN
                    DMNMX = D2
                    ELSE IF( D3 < D1 ) THEN
                    DMNMX = D3
                    ELSE
                    DMNMX = D1
                    END IF
                END IF
        !
                IF( DIR == 0 ) THEN
        !
        !           Sort into decreasing order
        !
                    I = START - 1
                    J = ENDD + 1
        60       CONTINUE
        70       CONTINUE
                    J = J - 1
                    IF( D( J ) < DMNMX ) &
                        GO TO 70
        80       CONTINUE
                    I = I + 1
                    IF( D( I ) > DMNMX ) &
                        GO TO 80
                    IF( I < J ) THEN
                    TMP = D( I )
                    D( I ) = D( J )
                    D( J ) = TMP
                    GO TO 60
                    END IF
                    IF( J-START > ENDD-J-1 ) THEN
                    STKPNT = STKPNT + 1
                    STACK( 1, STKPNT ) = START
                    STACK( 2, STKPNT ) = J
                    STKPNT = STKPNT + 1
                    STACK( 1, STKPNT ) = J + 1
                    STACK( 2, STKPNT ) = ENDD
                    ELSE
                    STKPNT = STKPNT + 1
                    STACK( 1, STKPNT ) = J + 1
                    STACK( 2, STKPNT ) = ENDD
                    STKPNT = STKPNT + 1
                    STACK( 1, STKPNT ) = START
                    STACK( 2, STKPNT ) = J
                    END IF
                ELSE
        !
        !           Sort into increasing order
        !
                    I = START - 1
                    J = ENDD + 1
        90       CONTINUE
        100       CONTINUE
                    J = J - 1
                    IF( D( J ) > DMNMX ) &
                        GO TO 100
        110       CONTINUE
                    I = I + 1
                    IF( D( I ) < DMNMX ) &
                        GO TO 110
                    IF( I < J ) THEN
                    TMP = D( I )
                    D( I ) = D( J )
                    D( J ) = TMP
                    GO TO 90
                    END IF
                    IF( J-START > ENDD-J-1 ) THEN
                    STKPNT = STKPNT + 1
                    STACK( 1, STKPNT ) = START
                    STACK( 2, STKPNT ) = J
                    STKPNT = STKPNT + 1
                    STACK( 1, STKPNT ) = J + 1
                    STACK( 2, STKPNT ) = ENDD
                    ELSE
                    STKPNT = STKPNT + 1
                    STACK( 1, STKPNT ) = J + 1
                    STACK( 2, STKPNT ) = ENDD
                    STKPNT = STKPNT + 1
                    STACK( 1, STKPNT ) = START
                    STACK( 2, STKPNT ) = J
                    END IF
                END IF
            END IF
            IF( STKPNT > 0 ) &
                GO TO 10
            RETURN
        !
        !     End of DLASRT
        !
        END SUBROUTINE
    
    !==================================================================================================
        DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          NORM
            INTEGER            N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   D( * ), E( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE, ZERO
            PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            I
            DOUBLE PRECISION   ANORM, SCALE, SUM
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME, DISNAN
            EXTERNAL           LSAME, DISNAN
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DLASSQ
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          ABS, SQRT
        !     ..
        !     .. Executable Statements ..
        !
            IF( N <= 0 ) THEN
                ANORM = ZERO
            ELSE IF( LSAME( NORM, 'M' ) ) THEN
        !
        !        Find max(abs(A(i,j))).
        !
                ANORM = ABS( D( N ) )
                DO 10 I = 1, N - 1
                    SUM = ABS( D( I ) )
                    IF( ANORM  <  SUM .OR. DISNAN( SUM ) ) ANORM = SUM
                    SUM = ABS( E( I ) )
                    IF( ANORM  <  SUM .OR. DISNAN( SUM ) ) ANORM = SUM
        10    CONTINUE
            ELSE IF( LSAME( NORM, 'O' ) .OR. NORM == '1' .OR. &
                    LSAME( NORM, 'I' ) ) THEN
        !
        !        Find norm1(A).
        !
                IF( N == 1 ) THEN
                    ANORM = ABS( D( 1 ) )
                ELSE
                    ANORM = ABS( D( 1 ) )+ABS( E( 1 ) )
                    SUM = ABS( E( N-1 ) )+ABS( D( N ) )
                    IF( ANORM  <  SUM .OR. DISNAN( SUM ) ) ANORM = SUM
                    DO 20 I = 2, N - 1
                    SUM = ABS( D( I ) )+ABS( E( I ) )+ABS( E( I-1 ) )
                    IF( ANORM  <  SUM .OR. DISNAN( SUM ) ) ANORM = SUM
        20       CONTINUE
                END IF
            ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
        !
        !        Find normF(A).
        !
                SCALE = ZERO
                SUM = ONE
                IF( N > 1 ) THEN
                    CALL DLASSQ( N-1, E, 1, SCALE, SUM )
                    SUM = 2*SUM
                END IF
                CALL DLASSQ( N, D, 1, SCALE, SUM )
                ANORM = SCALE*SQRT( SUM )
            END IF
        !
            DLANST = ANORM
            RETURN
        !
        !     End of DLANST
        !
        END FUNCTION
    
    !==================================================================================================
        SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION   A, B, C, RT1, RT2
        !     ..
        !
        ! =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE
            PARAMETER          ( ONE = 1.0D0 )
            DOUBLE PRECISION   TWO
            PARAMETER          ( TWO = 2.0D0 )
            DOUBLE PRECISION   ZERO
            PARAMETER          ( ZERO = 0.0D0 )
            DOUBLE PRECISION   HALF
            PARAMETER          ( HALF = 0.5D0 )
        !     ..
        !     .. Local Scalars ..
            DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          ABS, SQRT
        !     ..
        !     .. Executable Statements ..
        !
        !     Compute the eigenvalues
        !
            SM = A + C
            DF = A - C
            ADF = ABS( DF )
            TB = B + B
            AB = ABS( TB )
            IF( ABS( A ) > ABS( C ) ) THEN
                ACMX = A
                ACMN = C
            ELSE
                ACMX = C
                ACMN = A
            END IF
            IF( ADF > AB ) THEN
                RT = ADF*SQRT( ONE+( AB / ADF )**2 )
            ELSE IF( ADF < AB ) THEN
                RT = AB*SQRT( ONE+( ADF / AB )**2 )
            ELSE
        !
        !        Includes case AB=ADF=0
        !
                RT = AB*SQRT( TWO )
            END IF
            IF( SM < ZERO ) THEN
                RT1 = HALF*( SM-RT )
        !
        !        Order of execution important.
        !        To get fully accurate smaller eigenvalue,
        !        next line needs to be executed in higher precision.
        !
                RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
            ELSE IF( SM > ZERO ) THEN
                RT1 = HALF*( SM+RT )
        !
        !        Order of execution important.
        !        To get fully accurate smaller eigenvalue,
        !        next line needs to be executed in higher precision.
        !
                RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
            ELSE
        !
        !        Includes case RT1 = RT2 = 0
        !
                RT1 = HALF*RT
                RT2 = -HALF*RT
            END IF
            RETURN
        !
        !     End of DLAE2
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
        !     ..
        !
        ! =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE
            PARAMETER          ( ONE = 1.0D0 )
            DOUBLE PRECISION   TWO
            PARAMETER          ( TWO = 2.0D0 )
            DOUBLE PRECISION   ZERO
            PARAMETER          ( ZERO = 0.0D0 )
            DOUBLE PRECISION   HALF
            PARAMETER          ( HALF = 0.5D0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            SGN1, SGN2
            DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM, &
                               TB, TN
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          ABS, SQRT
        !     ..
        !     .. Executable Statements ..
        !
        !     Compute the eigenvalues
        !
            SM = A + C
            DF = A - C
            ADF = ABS( DF )
            TB = B + B
            AB = ABS( TB )
            IF( ABS( A ) > ABS( C ) ) THEN
                ACMX = A
                ACMN = C
            ELSE
                ACMX = C
                ACMN = A
            END IF
            IF( ADF > AB ) THEN
                RT = ADF*SQRT( ONE+( AB / ADF )**2 )
            ELSE IF( ADF < AB ) THEN
                RT = AB*SQRT( ONE+( ADF / AB )**2 )
            ELSE
        !
        !        Includes case AB=ADF=0
        !
                RT = AB*SQRT( TWO )
            END IF
            IF( SM < ZERO ) THEN
                RT1 = HALF*( SM-RT )
                SGN1 = -1
        !
        !        Order of execution important.
        !        To get fully accurate smaller eigenvalue,
        !        next line needs to be executed in higher precision.
        !
                RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
            ELSE IF( SM > ZERO ) THEN
                RT1 = HALF*( SM+RT )
                SGN1 = 1
        !
        !        Order of execution important.
        !        To get fully accurate smaller eigenvalue,
        !        next line needs to be executed in higher precision.
        !
                RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
            ELSE
        !
        !        Includes case RT1 = RT2 = 0
        !
                RT1 = HALF*RT
                RT2 = -HALF*RT
                SGN1 = 1
            END IF
        !
        !     Compute the eigenvector
        !
            IF( DF >= ZERO ) THEN
                CS = DF + RT
                SGN2 = 1
            ELSE
                CS = DF - RT
                SGN2 = -1
            END IF
            ACS = ABS( CS )
            IF( ACS > AB ) THEN
                CT = -TB / CS
                SN1 = ONE / SQRT( ONE+CT*CT )
                CS1 = CT*SN1
            ELSE
                IF( AB == ZERO ) THEN
                    CS1 = ONE
                    SN1 = ZERO
                ELSE
                    TN = -CS / TB
                    CS1 = ONE / SQRT( ONE+TN*TN )
                    SN1 = TN*CS1
                END IF
            END IF
            IF( SGN1 == SGN2 ) THEN
                TN = CS1
                CS1 = -SN1
                SN1 = TN
            END IF
            RETURN
        !
        !     End of DLAEV2
        !
        END SUBROUTINE
    
    !==================================================================================================
        DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION   X, Y
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ZERO
            PARAMETER          ( ZERO = 0.0D0 )
            DOUBLE PRECISION   ONE
            PARAMETER          ( ONE = 1.0D0 )
        !     ..
        !     .. Local Scalars ..
            DOUBLE PRECISION   W, XABS, YABS, Z, HUGEVAL
            LOGICAL            X_IS_NAN, Y_IS_NAN
        !     ..
        !     .. External Functions ..
            LOGICAL            DISNAN
            EXTERNAL           DISNAN
        !     ..
        !     .. External Subroutines ..
            DOUBLE PRECISION   DLAMCH
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          ABS, MAX, MIN, SQRT
        !     ..
        !     .. Executable Statements ..
        !
            X_IS_NAN = DISNAN( X )
            Y_IS_NAN = DISNAN( Y )
            IF ( X_IS_NAN ) DLAPY2 = X
            IF ( Y_IS_NAN ) DLAPY2 = Y
            HUGEVAL = DLAMCH( 'Overflow' )
        !
            IF ( .NOT.( X_IS_NAN.OR.Y_IS_NAN ) ) THEN
                XABS = ABS( X )
                YABS = ABS( Y )
                W = MAX( XABS, YABS )
                Z = MIN( XABS, YABS )
                IF( Z == ZERO .OR. W > HUGEVAL ) THEN
                    DLAPY2 = W
                ELSE
                    DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
                END IF
            END IF
            RETURN
        !
        !     End of DLAPY2
        !
        END FUNCTION
    
    !==================================================================================================
        LOGICAL FUNCTION DISNAN( DIN )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION, INTENT(IN) :: DIN
        !     ..
        !
        !  =====================================================================
        !
        !  .. External Functions ..
            LOGICAL DLAISNAN
            EXTERNAL DLAISNAN
        !  ..
        !  .. Executable Statements ..
            DISNAN = DLAISNAN(DIN,DIN)
            RETURN
        END FUNCTION
    
    !==================================================================================================
        LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
        !     ..
        !
        !  =====================================================================
        !
        !  .. Executable Statements ..
            DLAISNAN = (DIN1 /= DIN2)
            RETURN
        END FUNCTION
    
    !==================================================================================================
        SUBROUTINE DLARTG( F, G, CS, SN, R )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine (version 3.2) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2006
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION   CS, F, G, R, SN
        !     ..
        !
        !  Purpose
        !  =======
        !
        !  DLARTG generate a plane rotation so that
        !
        !     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
        !     [ -SN  CS  ]     [ G ]     [ 0 ]
        !
        !  This is a slower, more accurate version of the BLAS1 routine DROTG,
        !  with the following other differences:
        !     F and G are unchanged on return.
        !     If G=0, then CS=1 and SN=0.
        !     If F=0 and (G  /=  0), then CS=0 and SN=1 without doing any
        !        floating point operations (saves work in DBDSQR when
        !        there are zeros on the diagonal).
        !
        !  If F exceeds G in magnitude, CS will be positive.
        !
        !  Arguments
        !  =========
        !
        !  F       (input) DOUBLE PRECISION
        !          The first component of vector to be rotated.
        !
        !  G       (input) DOUBLE PRECISION
        !          The second component of vector to be rotated.
        !
        !  CS      (output) DOUBLE PRECISION
        !          The cosine of the rotation.
        !
        !  SN      (output) DOUBLE PRECISION
        !          The sine of the rotation.
        !
        !  R       (output) DOUBLE PRECISION
        !          The nonzero component of the rotated vector.
        !
        !  This version has a few statements commented out for thread safety
        !  (machine parameters are computed on each entry). 10 feb 03, SJH.
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ZERO
            PARAMETER          ( ZERO = 0.0D0 )
            DOUBLE PRECISION   ONE
            PARAMETER          ( ONE = 1.0D0 )
            DOUBLE PRECISION   TWO
            PARAMETER          ( TWO = 2.0D0 )
        !     ..
        !     .. Local Scalars ..
        !     LOGICAL            FIRST
            INTEGER            COUNT, I
            DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
        !     ..
        !     .. External Functions ..
            DOUBLE PRECISION   DLAMCH
            EXTERNAL           DLAMCH
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          ABS, INT, LOG, MAX, SQRT
        !     ..
        !     .. Save statement ..
        !     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
        !     ..
        !     .. Data statements ..
        !     DATA               FIRST / .TRUE. /
        !     ..
        !     .. Executable Statements ..
        !
        !     IF( FIRST ) THEN
                SAFMIN = DLAMCH( 'S' )
                EPS = DLAMCH( 'E' )
                SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / &
                        LOG( DLAMCH( 'B' ) ) / TWO )
                SAFMX2 = ONE / SAFMN2
        !        FIRST = .FALSE.
        !     END IF
            IF( G == ZERO ) THEN
                CS = ONE
                SN = ZERO
                R = F
            ELSE IF( F == ZERO ) THEN
                CS = ZERO
                SN = ONE
                R = G
            ELSE
                F1 = F
                G1 = G
                SCALE = MAX( ABS( F1 ), ABS( G1 ) )
                IF( SCALE >= SAFMX2 ) THEN
                    COUNT = 0
        10       CONTINUE
                    COUNT = COUNT + 1
                    F1 = F1*SAFMN2
                    G1 = G1*SAFMN2
                    SCALE = MAX( ABS( F1 ), ABS( G1 ) )
                    IF( SCALE >= SAFMX2 ) &
                        GO TO 10
                    R = SQRT( F1**2+G1**2 )
                    CS = F1 / R
                    SN = G1 / R
                    DO 20 I = 1, COUNT
                    R = R*SAFMX2
        20       CONTINUE
                ELSE IF( SCALE <= SAFMN2 ) THEN
                    COUNT = 0
        30       CONTINUE
                    COUNT = COUNT + 1
                    F1 = F1*SAFMX2
                    G1 = G1*SAFMX2
                    SCALE = MAX( ABS( F1 ), ABS( G1 ) )
                    IF( SCALE <= SAFMN2 ) &
                        GO TO 30
                    R = SQRT( F1**2+G1**2 )
                    CS = F1 / R
                    SN = G1 / R
                    DO 40 I = 1, COUNT
                    R = R*SAFMN2
        40       CONTINUE
                ELSE
                    R = SQRT( F1**2+G1**2 )
                    CS = F1 / R
                    SN = G1 / R
                END IF
                IF( ABS( F ) > ABS( G ) .AND. CS < ZERO ) THEN
                    CS = -CS
                    SN = -SN
                    R = -R
                END IF
            END IF
            RETURN
        !
        !     End of DLARTG
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine (version 3.2) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2006
        !
        !     .. Scalar Arguments ..
            INTEGER            INCX, N
            DOUBLE PRECISION   SCALE, SUMSQ
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   X( * )
        !     ..
        !
        !  Purpose
        !  =======
        !
        !  DLASSQ  returns the values  scl  and  smsq  such that
        !
        !     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
        !
        !  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
        !  assumed to be non-negative and  scl  returns the value
        !
        !     scl = max( scale, abs( x( i ) ) ).
        !
        !  scale and sumsq must be supplied in SCALE and SUMSQ and
        !  scl and smsq are overwritten on SCALE and SUMSQ respectively.
        !
        !  The routine makes only one pass through the vector x.
        !
        !  Arguments
        !  =========
        !
        !  N       (input) INTEGER
        !          The number of elements to be used from the vector X.
        !
        !  X       (input) DOUBLE PRECISION array, dimension (N)
        !          The vector for which a scaled sum of squares is computed.
        !             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
        !
        !  INCX    (input) INTEGER
        !          The increment between successive values of the vector X.
        !          INCX > 0.
        !
        !  SCALE   (input/output) DOUBLE PRECISION
        !          On entry, the value  scale  in the equation above.
        !          On exit, SCALE is overwritten with  scl , the scaling factor
        !          for the sum of squares.
        !
        !  SUMSQ   (input/output) DOUBLE PRECISION
        !          On entry, the value  sumsq  in the equation above.
        !          On exit, SUMSQ is overwritten with  smsq , the basic sum of
        !          squares from which  scl  has been factored out.
        !
        ! =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ZERO
            PARAMETER          ( ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            IX
            DOUBLE PRECISION   ABSXI
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          ABS
        !     ..
        !     .. Executable Statements ..
        !
            IF( N > 0 ) THEN
                DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
                    IF( X( IX ) /= ZERO ) THEN
                    ABSXI = ABS( X( IX ) )
                    IF( SCALE < ABSXI ) THEN
                        SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                        SCALE = ABSXI
                    ELSE
                        SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
                    END IF
                    END IF
        10    CONTINUE
            END IF
            RETURN
        !
        !     End of DLASSQ
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    !==================================================================================================
        !
        !  -- Reference BLAS level3 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION ALPHA
            INTEGER LDA,LDB,M,N
            CHARACTER DIAG,SIDE,TRANSA,UPLO
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION A(LDA,*),B(LDB,*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. External Functions ..
            LOGICAL LSAME
            EXTERNAL LSAME
        !     ..
        !     .. External Subroutines ..
            EXTERNAL XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MAX
        !     ..
        !     .. Local Scalars ..
            DOUBLE PRECISION TEMP
            INTEGER I,INFO,J,K,NROWA
            LOGICAL LSIDE,NOUNIT,UPPER
        !     ..
        !     .. Parameters ..
            DOUBLE PRECISION ONE,ZERO
            PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
        !     ..
        !
        !     Test the input parameters.
        !
            LSIDE = LSAME(SIDE,'L')
            IF (LSIDE) THEN
                NROWA = M
            ELSE
                NROWA = N
            END IF
            NOUNIT = LSAME(DIAG,'N')
            UPPER = LSAME(UPLO,'U')
        !
            INFO = 0
            IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
                INFO = 1
            ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
                INFO = 2
            ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
                    (.NOT.LSAME(TRANSA,'T')) .AND. &
                    (.NOT.LSAME(TRANSA,'C'))) THEN
                INFO = 3
            ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
                INFO = 4
            ELSE IF (M < 0) THEN
                INFO = 5
            ELSE IF (N < 0) THEN
                INFO = 6
            ELSE IF (LDA < MAX(1,NROWA)) THEN
                INFO = 9
            ELSE IF (LDB < MAX(1,M)) THEN
                INFO = 11
            END IF
            IF (INFO /= 0) THEN
                CALL XERBLA('DTRMM ',INFO)
                RETURN
            END IF
        !
        !     Quick return if possible.
        !
            IF (M == 0 .OR. N == 0) RETURN
        !
        !     And when  alpha == zero.
        !
            IF (ALPHA == ZERO) THEN
                DO 20 J = 1,N
                    DO 10 I = 1,M
                        B(I,J) = ZERO
        10         CONTINUE
        20     CONTINUE
                RETURN
            END IF
        !
        !     Start the operations.
        !
            IF (LSIDE) THEN
                IF (LSAME(TRANSA,'N')) THEN
        !
        !           Form  B := alpha*A*B.
        !
                    IF (UPPER) THEN
                        DO 50 J = 1,N
                            DO 40 K = 1,M
                                IF (B(K,J) /= ZERO) THEN
                                    TEMP = ALPHA*B(K,J)
                                    DO 30 I = 1,K - 1
                                        B(I,J) = B(I,J) + TEMP*A(I,K)
        30                         CONTINUE
                                    IF (NOUNIT) TEMP = TEMP*A(K,K)
                                    B(K,J) = TEMP
                                END IF
        40                 CONTINUE
        50             CONTINUE
                    ELSE
                        DO 80 J = 1,N
                            DO 70 K = M,1,-1
                                IF (B(K,J) /= ZERO) THEN
                                    TEMP = ALPHA*B(K,J)
                                    B(K,J) = TEMP
                                    IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                                    DO 60 I = K + 1,M
                                        B(I,J) = B(I,J) + TEMP*A(I,K)
        60                         CONTINUE
                                END IF
        70                 CONTINUE
        80             CONTINUE
                    END IF
                ELSE
        !
        !           Form  B := alpha*A**T*B.
        !
                    IF (UPPER) THEN
                        DO 110 J = 1,N
                            DO 100 I = M,1,-1
                                TEMP = B(I,J)
                                IF (NOUNIT) TEMP = TEMP*A(I,I)
                                DO 90 K = 1,I - 1
                                    TEMP = TEMP + A(K,I)*B(K,J)
        90                     CONTINUE
                                B(I,J) = ALPHA*TEMP
        100                 CONTINUE
        110             CONTINUE
                    ELSE
                        DO 140 J = 1,N
                            DO 130 I = 1,M
                                TEMP = B(I,J)
                                IF (NOUNIT) TEMP = TEMP*A(I,I)
                                DO 120 K = I + 1,M
                                    TEMP = TEMP + A(K,I)*B(K,J)
        120                     CONTINUE
                                B(I,J) = ALPHA*TEMP
        130                 CONTINUE
        140             CONTINUE
                    END IF
                END IF
            ELSE
                IF (LSAME(TRANSA,'N')) THEN
        !
        !           Form  B := alpha*B*A.
        !
                    IF (UPPER) THEN
                        DO 180 J = N,1,-1
                            TEMP = ALPHA
                            IF (NOUNIT) TEMP = TEMP*A(J,J)
                            DO 150 I = 1,M
                                B(I,J) = TEMP*B(I,J)
        150                 CONTINUE
                            DO 170 K = 1,J - 1
                                IF (A(K,J) /= ZERO) THEN
                                    TEMP = ALPHA*A(K,J)
                                    DO 160 I = 1,M
                                        B(I,J) = B(I,J) + TEMP*B(I,K)
        160                         CONTINUE
                                END IF
        170                 CONTINUE
        180             CONTINUE
                    ELSE
                        DO 220 J = 1,N
                            TEMP = ALPHA
                            IF (NOUNIT) TEMP = TEMP*A(J,J)
                            DO 190 I = 1,M
                                B(I,J) = TEMP*B(I,J)
        190                 CONTINUE
                            DO 210 K = J + 1,N
                                IF (A(K,J) /= ZERO) THEN
                                    TEMP = ALPHA*A(K,J)
                                    DO 200 I = 1,M
                                        B(I,J) = B(I,J) + TEMP*B(I,K)
        200                         CONTINUE
                                END IF
        210                 CONTINUE
        220             CONTINUE
                    END IF
                ELSE
        !
        !           Form  B := alpha*B*A**T.
        !
                    IF (UPPER) THEN
                        DO 260 K = 1,N
                            DO 240 J = 1,K - 1
                                IF (A(J,K) /= ZERO) THEN
                                    TEMP = ALPHA*A(J,K)
                                    DO 230 I = 1,M
                                        B(I,J) = B(I,J) + TEMP*B(I,K)
        230                         CONTINUE
                                END IF
        240                 CONTINUE
                            TEMP = ALPHA
                            IF (NOUNIT) TEMP = TEMP*A(K,K)
                            IF (TEMP /= ONE) THEN
                                DO 250 I = 1,M
                                    B(I,K) = TEMP*B(I,K)
        250                     CONTINUE
                            END IF
        260             CONTINUE
                    ELSE
                        DO 300 K = N,1,-1
                            DO 280 J = K + 1,N
                                IF (A(J,K) /= ZERO) THEN
                                    TEMP = ALPHA*A(J,K)
                                    DO 270 I = 1,M
                                        B(I,J) = B(I,J) + TEMP*B(I,K)
        270                         CONTINUE
                                END IF
        280                 CONTINUE
                            TEMP = ALPHA
                            IF (NOUNIT) TEMP = TEMP*A(K,K)
                            IF (TEMP /= ONE) THEN
                                DO 290 I = 1,M
                                    B(I,K) = TEMP*B(I,K)
        290                     CONTINUE
                            END IF
        300             CONTINUE
                    END IF
                END IF
            END IF
        !
            RETURN
        !
        !     End of DTRMM
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          UPLO
            INTEGER            LDA, LDW, N, NB
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ZERO, ONE, HALF
            PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            I, IW
            DOUBLE PRECISION   ALPHA
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DAXPY, DGEMV, DLARFG, DSCAL, DSYMV
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            DOUBLE PRECISION   DDOT
            EXTERNAL           LSAME, DDOT
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Quick return if possible
        !
            IF( N <= 0 ) &
                RETURN
        !
            IF( LSAME( UPLO, 'U' ) ) THEN
        !
        !        Reduce last NB columns of upper triangle
        !
                DO 10 I = N, N - NB + 1, -1
                    IW = I - N + NB
                    IF( I < N ) THEN
        !
        !              Update A(1:i,i)
        !
                    CALL DGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ), &
                                LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
                    CALL DGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ), &
                                LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
                    END IF
                    IF( I > 1 ) THEN
        !
        !              Generate elementary reflector H(i) to annihilate
        !              A(1:i-2,i)
        !
                    CALL DLARFG( I-1, A( I-1, I ), A( 1, I ), 1, TAU( I-1 ) )
                    E( I-1 ) = A( I-1, I )
                    A( I-1, I ) = ONE
        !
        !              Compute W(1:i-1,i)
        !
                    CALL DSYMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1, &
                                ZERO, W( 1, IW ), 1 )
                    IF( I < N ) THEN
                        CALL DGEMV( 'Transpose', I-1, N-I, ONE, W( 1, IW+1 ), &
                                    LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                        CALL DGEMV( 'No transpose', I-1, N-I, -ONE, &
                                    A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE, &
                                    W( 1, IW ), 1 )
                        CALL DGEMV( 'Transpose', I-1, N-I, ONE, A( 1, I+1 ), &
                                    LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                        CALL DGEMV( 'No transpose', I-1, N-I, -ONE, &
                                    W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE, &
                                    W( 1, IW ), 1 )
                    END IF
                    CALL DSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
                    ALPHA = -HALF*TAU( I-1 )*DDOT( I-1, W( 1, IW ), 1, &
                                A( 1, I ), 1 )
                    CALL DAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
                    END IF
        !
        10    CONTINUE
            ELSE
        !
        !        Reduce first NB columns of lower triangle
        !
                DO 20 I = 1, NB
        !
        !           Update A(i:n,i)
        !
                    CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ), &
                                LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
                    CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ), &
                                LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
                    IF( I < N ) THEN
        !
        !              Generate elementary reflector H(i) to annihilate
        !              A(i+2:n,i)
        !
                    CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1, &
                                TAU( I ) )
                    E( I ) = A( I+1, I )
                    A( I+1, I ) = ONE
        !
        !              Compute W(i+1:n,i)
        !
                    CALL DSYMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA, &
                                A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
                    CALL DGEMV( 'Transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW, &
                                A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
                    CALL DGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ), &
                                LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
                    CALL DGEMV( 'Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, &
                                A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
                    CALL DGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ), &
                                LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
                    CALL DSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
                    ALPHA = -HALF*TAU( I )*DDOT( N-I, W( I+1, I ), 1, &
                                A( I+1, I ), 1 )
                    CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
                    END IF
        !
        20    CONTINUE
            END IF
        !
            RETURN
        !
        !     End of DLATRD
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
    !==================================================================================================
        !
        !  -- LAPACK computational routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          DIAG, UPLO
            INTEGER            INFO, LDA, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE
            PARAMETER          ( ONE = 1.0D+0 )
        !     ..
        !     .. Local Scalars ..
            LOGICAL            NOUNIT, UPPER
            INTEGER            J
            DOUBLE PRECISION   AJJ
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            EXTERNAL           LSAME
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DSCAL, DTRMV, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MAX
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters.
        !
            INFO = 0
            UPPER = LSAME( UPLO, 'U' )
            NOUNIT = LSAME( DIAG, 'N' )
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
                INFO = -1
            ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
                INFO = -2
            ELSE IF( N < 0 ) THEN
                INFO = -3
            ELSE IF( LDA < MAX( 1, N ) ) THEN
                INFO = -5
            END IF
            IF( INFO /= 0 ) THEN
                CALL XERBLA( 'DTRTI2', -INFO )
                RETURN
            END IF
        !
            IF( UPPER ) THEN
        !
        !        Compute inverse of upper triangular matrix.
        !
                DO 10 J = 1, N
                    IF( NOUNIT ) THEN
                    A( J, J ) = ONE / A( J, J )
                    AJJ = -A( J, J )
                    ELSE
                    AJJ = -ONE
                    END IF
        !
        !           Compute elements 1:j-1 of j-th column.
        !
                    CALL DTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA, &
                                A( 1, J ), 1 )
                    CALL DSCAL( J-1, AJJ, A( 1, J ), 1 )
        10    CONTINUE
            ELSE
        !
        !        Compute inverse of lower triangular matrix.
        !
                DO 20 J = N, 1, -1
                    IF( NOUNIT ) THEN
                    A( J, J ) = ONE / A( J, J )
                    AJJ = -A( J, J )
                    ELSE
                    AJJ = -ONE
                    END IF
                    IF( J < N ) THEN
        !
        !              Compute elements j+1:n of j-th column.
        !
                    CALL DTRMV( 'Lower', 'No transpose', DIAG, N-J, &
                                A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
                    CALL DSCAL( N-J, AJJ, A( J+1, J ), 1 )
                    END IF
        20    CONTINUE
            END IF
        !
            RETURN
        !
        !     End of DTRTI2
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    !==================================================================================================
        !
        !  -- Reference BLAS level3 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION ALPHA,BETA
            INTEGER K,LDA,LDB,LDC,N
            CHARACTER TRANS,UPLO
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. External Functions ..
            LOGICAL LSAME
            EXTERNAL LSAME
        !     ..
        !     .. External Subroutines ..
            EXTERNAL XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MAX
        !     ..
        !     .. Local Scalars ..
            DOUBLE PRECISION TEMP1,TEMP2
            INTEGER I,INFO,J,L,NROWA
            LOGICAL UPPER
        !     ..
        !     .. Parameters ..
            DOUBLE PRECISION ONE,ZERO
            PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
        !     ..
        !
        !     Test the input parameters.
        !
            IF (LSAME(TRANS,'N')) THEN
                NROWA = N
            ELSE
                NROWA = K
            END IF
            UPPER = LSAME(UPLO,'U')
        !
            INFO = 0
            IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
                INFO = 1
            ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND. &
                    (.NOT.LSAME(TRANS,'T')) .AND. &
                    (.NOT.LSAME(TRANS,'C'))) THEN
                INFO = 2
            ELSE IF (N < 0) THEN
                INFO = 3
            ELSE IF (K < 0) THEN
                INFO = 4
            ELSE IF (LDA < MAX(1,NROWA)) THEN
                INFO = 7
            ELSE IF (LDB < MAX(1,NROWA)) THEN
                INFO = 9
            ELSE IF (LDC < MAX(1,N)) THEN
                INFO = 12
            END IF
            IF (INFO /= 0) THEN
                CALL XERBLA('DSYR2K',INFO)
                RETURN
            END IF
        !
        !     Quick return if possible.
        !
            IF ((N == 0) .OR. (((ALPHA == ZERO).OR. &
                (K == 0)).AND. (BETA == ONE))) RETURN
        !
        !     And when  alpha == zero.
        !
            IF (ALPHA == ZERO) THEN
                IF (UPPER) THEN
                    IF (BETA == ZERO) THEN
                        DO 20 J = 1,N
                            DO 10 I = 1,J
                                C(I,J) = ZERO
        10                 CONTINUE
        20             CONTINUE
                    ELSE
                        DO 40 J = 1,N
                            DO 30 I = 1,J
                                C(I,J) = BETA*C(I,J)
        30                 CONTINUE
        40             CONTINUE
                    END IF
                ELSE
                    IF (BETA == ZERO) THEN
                        DO 60 J = 1,N
                            DO 50 I = J,N
                                C(I,J) = ZERO
        50                 CONTINUE
        60             CONTINUE
                    ELSE
                        DO 80 J = 1,N
                            DO 70 I = J,N
                                C(I,J) = BETA*C(I,J)
        70                 CONTINUE
        80             CONTINUE
                    END IF
                END IF
                RETURN
            END IF
        !
        !     Start the operations.
        !
            IF (LSAME(TRANS,'N')) THEN
        !
        !        Form  C := alpha*A*B**T + alpha*B*A**T + C.
        !
                IF (UPPER) THEN
                    DO 130 J = 1,N
                        IF (BETA == ZERO) THEN
                            DO 90 I = 1,J
                                C(I,J) = ZERO
        90                 CONTINUE
                        ELSE IF (BETA /= ONE) THEN
                            DO 100 I = 1,J
                                C(I,J) = BETA*C(I,J)
        100                 CONTINUE
                        END IF
                        DO 120 L = 1,K
                            IF ((A(J,L) /= ZERO) .OR. (B(J,L) /= ZERO)) THEN
                                TEMP1 = ALPHA*B(J,L)
                                TEMP2 = ALPHA*A(J,L)
                                DO 110 I = 1,J
                                    C(I,J) = C(I,J) + A(I,L)*TEMP1 + &
                                            B(I,L)*TEMP2
        110                     CONTINUE
                            END IF
        120             CONTINUE
        130         CONTINUE
                ELSE
                    DO 180 J = 1,N
                        IF (BETA == ZERO) THEN
                            DO 140 I = J,N
                                C(I,J) = ZERO
        140                 CONTINUE
                        ELSE IF (BETA /= ONE) THEN
                            DO 150 I = J,N
                                C(I,J) = BETA*C(I,J)
        150                 CONTINUE
                        END IF
                        DO 170 L = 1,K
                            IF ((A(J,L) /= ZERO) .OR. (B(J,L) /= ZERO)) THEN
                                TEMP1 = ALPHA*B(J,L)
                                TEMP2 = ALPHA*A(J,L)
                                DO 160 I = J,N
                                    C(I,J) = C(I,J) + A(I,L)*TEMP1 + &
                                            B(I,L)*TEMP2
        160                     CONTINUE
                            END IF
        170             CONTINUE
        180         CONTINUE
                END IF
            ELSE
        !
        !        Form  C := alpha*A**T*B + alpha*B**T*A + C.
        !
                IF (UPPER) THEN
                    DO 210 J = 1,N
                        DO 200 I = 1,J
                            TEMP1 = ZERO
                            TEMP2 = ZERO
                            DO 190 L = 1,K
                                TEMP1 = TEMP1 + A(L,I)*B(L,J)
                                TEMP2 = TEMP2 + B(L,I)*A(L,J)
        190                 CONTINUE
                            IF (BETA == ZERO) THEN
                                C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                            ELSE
                                C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + &
                                        ALPHA*TEMP2
                            END IF
        200             CONTINUE
        210         CONTINUE
                ELSE
                    DO 240 J = 1,N
                        DO 230 I = J,N
                            TEMP1 = ZERO
                            TEMP2 = ZERO
                            DO 220 L = 1,K
                                TEMP1 = TEMP1 + A(L,I)*B(L,J)
                                TEMP2 = TEMP2 + B(L,I)*A(L,J)
        220                 CONTINUE
                            IF (BETA == ZERO) THEN
                                C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                            ELSE
                                C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + &
                                        ALPHA*TEMP2
                            END IF
        230             CONTINUE
        240         CONTINUE
                END IF
            END IF
        !
            RETURN
        !
        !     End of DSYR2K
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
    !==================================================================================================
        !
        !  -- LAPACK computational routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          UPLO
            INTEGER            INFO, LDA, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE, ZERO, HALF
            PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, &
                               HALF = 1.0D0 / 2.0D0 )
        !     ..
        !     .. Local Scalars ..
            LOGICAL            UPPER
            INTEGER            I
            DOUBLE PRECISION   ALPHA, TAUI
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DAXPY, DLARFG, DSYMV, DSYR2, XERBLA
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            DOUBLE PRECISION   DDOT
            EXTERNAL           LSAME, DDOT
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MAX, MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters
        !
            INFO = 0
            UPPER = LSAME( UPLO, 'U' )
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
                INFO = -1
            ELSE IF( N < 0 ) THEN
                INFO = -2
            ELSE IF( LDA < MAX( 1, N ) ) THEN
                INFO = -4
            END IF
            IF( INFO /= 0 ) THEN
                CALL XERBLA( 'DSYTD2', -INFO )
                RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( N <= 0 ) &
                RETURN
        !
            IF( UPPER ) THEN
        !
        !        Reduce the upper triangle of A
        !
                DO 10 I = N - 1, 1, -1
        !
        !           Generate elementary reflector H(i) = I - tau * v * v**T
        !           to annihilate A(1:i-1,i+1)
        !
                    CALL DLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI )
                    E( I ) = A( I, I+1 )
        !
                    IF( TAUI /= ZERO ) THEN
        !
        !              Apply H(i) from both sides to A(1:i,1:i)
        !
                    A( I, I+1 ) = ONE
        !
        !              Compute  x := tau * A * v  storing x in TAU(1:i)
        !
                    CALL DSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO, &
                                TAU, 1 )
        !
        !              Compute  w := x - 1/2 * tau * (x**T * v) * v
        !
                    ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, A( 1, I+1 ), 1 )
                    CALL DAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
        !
        !              Apply the transformation as a rank-2 update:
        !                 A := A - v * w**T - w * v**T
        !
                    CALL DSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A, &
                                LDA )
        !
                    A( I, I+1 ) = E( I )
                    END IF
                    D( I+1 ) = A( I+1, I+1 )
                    TAU( I ) = TAUI
        10    CONTINUE
                D( 1 ) = A( 1, 1 )
            ELSE
        !
        !        Reduce the lower triangle of A
        !
                DO 20 I = 1, N - 1
        !
        !           Generate elementary reflector H(i) = I - tau * v * v**T
        !           to annihilate A(i+2:n,i)
        !
                    CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1, &
                                TAUI )
                    E( I ) = A( I+1, I )
        !
                    IF( TAUI /= ZERO ) THEN
        !
        !              Apply H(i) from both sides to A(i+1:n,i+1:n)
        !
                    A( I+1, I ) = ONE
        !
        !              Compute  x := tau * A * v  storing y in TAU(i:n-1)
        !
                    CALL DSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, &
                                A( I+1, I ), 1, ZERO, TAU( I ), 1 )
        !
        !              Compute  w := x - 1/2 * tau * (x**T * v) * v
        !
                    ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, A( I+1, I ), &
                            1 )
                    CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
        !
        !              Apply the transformation as a rank-2 update:
        !                 A := A - v * w**T - w * v**T
        !
                    CALL DSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1, &
                                A( I+1, I+1 ), LDA )
        !
                    A( I+1, I ) = E( I )
                    END IF
                    D( I ) = A( I, I )
                    TAU( I ) = TAUI
        20    CONTINUE
                D( N ) = A( N, N )
            END IF
        !
            RETURN
        !
        !     End of DSYTD2
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
    !==================================================================================================
        !
        !  -- Reference BLAS level1 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION DA
            INTEGER INCX,INCY,N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION DX(*),DY(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
            INTEGER I,IX,IY,M,MP1
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MOD
        !     ..
            IF (N <= 0) RETURN
            IF (DA == 0.0d0) RETURN
            IF (INCX == 1 .AND. INCY == 1) THEN
        !
        !        code for both increments equal to 1
        !
        !
        !        clean-up loop
        !
            M = MOD(N,4)
            IF (M /= 0) THEN
                DO I = 1,M
                    DY(I) = DY(I) + DA*DX(I)
                END DO
            END IF
            IF (N < 4) RETURN
            MP1 = M + 1
            DO I = MP1,N,4
                DY(I) = DY(I) + DA*DX(I)
                DY(I+1) = DY(I+1) + DA*DX(I+1)
                DY(I+2) = DY(I+2) + DA*DX(I+2)
                DY(I+3) = DY(I+3) + DA*DX(I+3)
            END DO
            ELSE
        !
        !        code for unequal increments or equal increments
        !          not equal to 1
        !
            IX = 1
            IY = 1
            IF (INCX < 0) IX = (-N+1)*INCX + 1
            IF (INCY < 0) IY = (-N+1)*INCY + 1
            DO I = 1,N
                DY(IY) = DY(IY) + DA*DX(IX)
                IX = IX + INCX
                IY = IY + INCY
            END DO
            END IF
            RETURN
        !
        !     End of DAXPY
        !
        END SUBROUTINE
    
    !==================================================================================================
        DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
    !==================================================================================================
        !
        !  -- Reference BLAS level1 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER INCX,INCY,N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION DX(*),DY(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
            DOUBLE PRECISION DTEMP
            INTEGER I,IX,IY,M,MP1
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MOD
        !     ..
            DDOT = 0.0d0
            DTEMP = 0.0d0
            IF (N <= 0) RETURN
            IF (INCX == 1 .AND. INCY == 1) THEN
        !
        !        code for both increments equal to 1
        !
        !
        !        clean-up loop
        !
                M = MOD(N,5)
                IF (M /= 0) THEN
                    DO I = 1,M
                    DTEMP = DTEMP + DX(I)*DY(I)
                    END DO
                    IF (N < 5) THEN
                    DDOT=DTEMP
                    RETURN
                    END IF
                END IF
                MP1 = M + 1
                DO I = MP1,N,5
                DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) + &
                        DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
                END DO
            ELSE
        !
        !        code for unequal increments or equal increments
        !          not equal to 1
        !
                IX = 1
                IY = 1
                IF (INCX < 0) IX = (-N+1)*INCX + 1
                IF (INCY < 0) IY = (-N+1)*INCY + 1
                DO I = 1,N
                    DTEMP = DTEMP + DX(IX)*DY(IY)
                    IX = IX + INCX
                    IY = IY + INCY
                END DO
            END IF
            DDOT = DTEMP
            RETURN
        !
        !     End of DDOT
        !
        END FUNCTION
    
    !==================================================================================================
        SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    !==================================================================================================
        !
        !  -- Reference BLAS level2 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION ALPHA,BETA
            INTEGER INCX,INCY,LDA,N
            CHARACTER UPLO
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION A(LDA,*),X(*),Y(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION ONE,ZERO
            PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
        !     ..
        !     .. Local Scalars ..
            DOUBLE PRECISION TEMP1,TEMP2
            INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
        !     ..
        !     .. External Functions ..
            LOGICAL LSAME
            EXTERNAL LSAME
        !     ..
        !     .. External Subroutines ..
            EXTERNAL XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MAX
        !     ..
        !
        !     Test the input parameters.
        !
            INFO = 0
            IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
                INFO = 1
            ELSE IF (N < 0) THEN
                INFO = 2
            ELSE IF (LDA < MAX(1,N)) THEN
                INFO = 5
            ELSE IF (INCX == 0) THEN
                INFO = 7
            ELSE IF (INCY == 0) THEN
                INFO = 10
            END IF
            IF (INFO /= 0) THEN
                CALL XERBLA('DSYMV ',INFO)
                RETURN
            END IF
        !
        !     Quick return if possible.
        !
            IF ((N == 0) .OR. ((ALPHA == ZERO).AND. (BETA == ONE))) RETURN
        !
        !     Set up the start points in  X  and  Y.
        !
            IF (INCX > 0) THEN
                KX = 1
            ELSE
                KX = 1 - (N-1)*INCX
            END IF
            IF (INCY > 0) THEN
                KY = 1
            ELSE
                KY = 1 - (N-1)*INCY
            END IF
        !
        !     Start the operations. In this version the elements of A are
        !     accessed sequentially with one pass through the triangular part
        !     of A.
        !
        !     First form  y := beta*y.
        !
            IF (BETA /= ONE) THEN
                IF (INCY == 1) THEN
                    IF (BETA == ZERO) THEN
                        DO 10 I = 1,N
                            Y(I) = ZERO
        10             CONTINUE
                    ELSE
                        DO 20 I = 1,N
                            Y(I) = BETA*Y(I)
        20             CONTINUE
                    END IF
                ELSE
                    IY = KY
                    IF (BETA == ZERO) THEN
                        DO 30 I = 1,N
                            Y(IY) = ZERO
                            IY = IY + INCY
        30             CONTINUE
                    ELSE
                        DO 40 I = 1,N
                            Y(IY) = BETA*Y(IY)
                            IY = IY + INCY
        40             CONTINUE
                    END IF
                END IF
            END IF
            IF (ALPHA == ZERO) RETURN
            IF (LSAME(UPLO,'U')) THEN
        !
        !        Form  y  when A is stored in upper triangle.
        !
                IF ((INCX == 1) .AND. (INCY == 1)) THEN
                    DO 60 J = 1,N
                        TEMP1 = ALPHA*X(J)
                        TEMP2 = ZERO
                        DO 50 I = 1,J - 1
                            Y(I) = Y(I) + TEMP1*A(I,J)
                            TEMP2 = TEMP2 + A(I,J)*X(I)
        50             CONTINUE
                        Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2
        60         CONTINUE
                ELSE
                    JX = KX
                    JY = KY
                    DO 80 J = 1,N
                        TEMP1 = ALPHA*X(JX)
                        TEMP2 = ZERO
                        IX = KX
                        IY = KY
                        DO 70 I = 1,J - 1
                            Y(IY) = Y(IY) + TEMP1*A(I,J)
                            TEMP2 = TEMP2 + A(I,J)*X(IX)
                            IX = IX + INCX
                            IY = IY + INCY
        70             CONTINUE
                        Y(JY) = Y(JY) + TEMP1*A(J,J) + ALPHA*TEMP2
                        JX = JX + INCX
                        JY = JY + INCY
        80         CONTINUE
                END IF
            ELSE
        !
        !        Form  y  when A is stored in lower triangle.
        !
                IF ((INCX == 1) .AND. (INCY == 1)) THEN
                    DO 100 J = 1,N
                        TEMP1 = ALPHA*X(J)
                        TEMP2 = ZERO
                        Y(J) = Y(J) + TEMP1*A(J,J)
                        DO 90 I = J + 1,N
                            Y(I) = Y(I) + TEMP1*A(I,J)
                            TEMP2 = TEMP2 + A(I,J)*X(I)
        90             CONTINUE
                        Y(J) = Y(J) + ALPHA*TEMP2
        100         CONTINUE
                ELSE
                    JX = KX
                    JY = KY
                    DO 120 J = 1,N
                        TEMP1 = ALPHA*X(JX)
                        TEMP2 = ZERO
                        Y(JY) = Y(JY) + TEMP1*A(J,J)
                        IX = JX
                        IY = JY
                        DO 110 I = J + 1,N
                            IX = IX + INCX
                            IY = IY + INCY
                            Y(IY) = Y(IY) + TEMP1*A(I,J)
                            TEMP2 = TEMP2 + A(I,J)*X(IX)
        110             CONTINUE
                        Y(JY) = Y(JY) + ALPHA*TEMP2
                        JX = JX + INCX
                        JY = JY + INCY
        120         CONTINUE
                END IF
            END IF
        !
            RETURN
        !
        !     End of DSYMV
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER            INCX, N
            DOUBLE PRECISION   ALPHA, TAU
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   X( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE, ZERO
            PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            J, KNT
            DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
        !     ..
        !     .. External Functions ..
            DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
            EXTERNAL           DLAMCH, DLAPY2, DNRM2
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          ABS, SIGN
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DSCAL
        !     ..
        !     .. Executable Statements ..
        !
            IF( N <= 1 ) THEN
                TAU = ZERO
                RETURN
            END IF
        !
            XNORM = DNRM2( N-1, X, INCX )
        !
            IF( XNORM == ZERO ) THEN
        !
        !        H  =  I
        !
                TAU = ZERO
            ELSE
        !
        !        general case
        !
                BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
                SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
                KNT = 0
                IF( ABS( BETA ) < SAFMIN ) THEN
        !
        !           XNORM, BETA may be inaccurate; scale X and recompute them
        !
                    RSAFMN = ONE / SAFMIN
        10       CONTINUE
                    KNT = KNT + 1
                    CALL DSCAL( N-1, RSAFMN, X, INCX )
                    BETA = BETA*RSAFMN
                    ALPHA = ALPHA*RSAFMN
                    IF( (ABS( BETA ) < SAFMIN) .AND. (KNT  <  20) ) &
                        GO TO 10
        !
        !           New BETA is at most 1, at least SAFMIN
        !
                    XNORM = DNRM2( N-1, X, INCX )
                    BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
                END IF
                TAU = ( BETA-ALPHA ) / BETA
                CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
        !
        !        If ALPHA is subnormal, it may lose relative accuracy
        !
                DO 20 J = 1, KNT
                    BETA = BETA*SAFMIN
        20      CONTINUE
                ALPHA = BETA
            END IF
        !
            RETURN
        !
        !     End of DLARFG
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    !==================================================================================================
        !
        !  -- Reference BLAS level2 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION ALPHA
            INTEGER INCX,INCY,LDA,N
            CHARACTER UPLO
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION A(LDA,*),X(*),Y(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION ZERO
            PARAMETER (ZERO=0.0D+0)
        !     ..
        !     .. Local Scalars ..
            DOUBLE PRECISION TEMP1,TEMP2
            INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
        !     ..
        !     .. External Functions ..
            LOGICAL LSAME
            EXTERNAL LSAME
        !     ..
        !     .. External Subroutines ..
            EXTERNAL XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MAX
        !     ..
        !
        !     Test the input parameters.
        !
            INFO = 0
            IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
                INFO = 1
            ELSE IF (N < 0) THEN
                INFO = 2
            ELSE IF (INCX == 0) THEN
                INFO = 5
            ELSE IF (INCY == 0) THEN
                INFO = 7
            ELSE IF (LDA < MAX(1,N)) THEN
                INFO = 9
            END IF
            IF (INFO /= 0) THEN
                CALL XERBLA('DSYR2 ',INFO)
                RETURN
            END IF
        !
        !     Quick return if possible.
        !
            IF ((N == 0) .OR. (ALPHA == ZERO)) RETURN
        !
        !     Set up the start points in X and Y if the increments are not both
        !     unity.
        !
            IF ((INCX /= 1) .OR. (INCY /= 1)) THEN
                IF (INCX > 0) THEN
                    KX = 1
                ELSE
                    KX = 1 - (N-1)*INCX
                END IF
                IF (INCY > 0) THEN
                    KY = 1
                ELSE
                    KY = 1 - (N-1)*INCY
                END IF
                JX = KX
                JY = KY
            END IF
        !
        !     Start the operations. In this version the elements of A are
        !     accessed sequentially with one pass through the triangular part
        !     of A.
        !
            IF (LSAME(UPLO,'U')) THEN
        !
        !        Form  A  when A is stored in the upper triangle.
        !
                IF ((INCX == 1) .AND. (INCY == 1)) THEN
                    DO 20 J = 1,N
                        IF ((X(J) /= ZERO) .OR. (Y(J) /= ZERO)) THEN
                            TEMP1 = ALPHA*Y(J)
                            TEMP2 = ALPHA*X(J)
                            DO 10 I = 1,J
                                A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
        10                 CONTINUE
                        END IF
        20         CONTINUE
                ELSE
                    DO 40 J = 1,N
                        IF ((X(JX) /= ZERO) .OR. (Y(JY) /= ZERO)) THEN
                            TEMP1 = ALPHA*Y(JY)
                            TEMP2 = ALPHA*X(JX)
                            IX = KX
                            IY = KY
                            DO 30 I = 1,J
                                A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                                IX = IX + INCX
                                IY = IY + INCY
        30                 CONTINUE
                        END IF
                        JX = JX + INCX
                        JY = JY + INCY
        40         CONTINUE
                END IF
            ELSE
        !
        !        Form  A  when A is stored in the lower triangle.
        !
                IF ((INCX == 1) .AND. (INCY == 1)) THEN
                    DO 60 J = 1,N
                        IF ((X(J) /= ZERO) .OR. (Y(J) /= ZERO)) THEN
                            TEMP1 = ALPHA*Y(J)
                            TEMP2 = ALPHA*X(J)
                            DO 50 I = J,N
                                A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
        50                 CONTINUE
                        END IF
        60         CONTINUE
                ELSE
                    DO 80 J = 1,N
                        IF ((X(JX) /= ZERO) .OR. (Y(JY) /= ZERO)) THEN
                            TEMP1 = ALPHA*Y(JY)
                            TEMP2 = ALPHA*X(JX)
                            IX = JX
                            IY = JY
                            DO 70 I = J,N
                                A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                                IX = IX + INCX
                                IY = IY + INCY
        70                 CONTINUE
                        END IF
                        JX = JX + INCX
                        JY = JY + INCY
        80         CONTINUE
                END IF
            END IF
        !
            RETURN
        !
        !     End of DSYR2
        !
        END SUBROUTINE
    
    !  =====================================================================
        function DNRM2( n, x, incx ) 
            integer, parameter :: wp = kind(1.d0)
            real(wp) :: DNRM2
        !
        !  -- Reference BLAS level1 routine (version 3.9.1) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     March 2021
        !
        !  .. Constants ..
            real(wp), parameter :: zero = 0.0_wp
            real(wp), parameter :: one  = 1.0_wp
            real(wp), parameter :: maxN = huge(0.0_wp)
        !  ..
        !  .. Blue's scaling constants ..
            real(wp), parameter :: tsml = real(radix(0._wp), wp)**ceiling( &
                (minexponent(0._wp) - 1) * 0.5_wp)
            real(wp), parameter :: tbig = real(radix(0._wp), wp)**floor( &
                (maxexponent(0._wp) - digits(0._wp) + 1) * 0.5_wp)
            real(wp), parameter :: ssml = real(radix(0._wp), wp)**( - floor( &
                (minexponent(0._wp) - digits(0._wp)) * 0.5_wp))
            real(wp), parameter :: sbig = real(radix(0._wp), wp)**( - ceiling( &
                (maxexponent(0._wp) + digits(0._wp) - 1) * 0.5_wp))
        !  ..
        !  .. Scalar Arguments ..
            integer :: incx, n
        !  ..
        !  .. Array Arguments ..
            real(wp) :: x(*)
        !  ..
        !  .. Local Scalars ..
            integer :: i, ix
            logical :: notbig
            real(wp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
        !
        !  Quick return if possible
        !
            DNRM2 = zero
            if( n <= 0 ) return
        !
            scl = one
            sumsq = zero
        !
        !  Compute the sum of squares in 3 accumulators:
        !     abig -- sums of squares scaled down to avoid overflow
        !     asml -- sums of squares scaled up to avoid underflow
        !     amed -- sums of squares that do not require scaling
        !  The thresholds and multipliers are
        !     tbig -- values bigger than this are scaled down by sbig
        !     tsml -- values smaller than this are scaled up by ssml
        !
            notbig = .true.
            asml = zero
            amed = zero
            abig = zero
            ix = 1
            if( incx < 0 ) ix = 1 - (n-1)*incx
            do i = 1, n
                ax = abs(x(ix))
                if (ax > tbig) then
                    abig = abig + (ax*sbig)**2
                    notbig = .false.
                else if (ax < tsml) then
                    if (notbig) asml = asml + (ax*ssml)**2
                else
                    amed = amed + ax**2
                end if
                ix = ix + incx
            end do
        !
        !  Combine abig and amed or amed and asml if more than one
        !  accumulator was used.
        !
            if (abig > zero) then
        !
        !     Combine abig and amed if abig > 0.
        !
                if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
                    abig = abig + (amed*sbig)*sbig
                end if
                scl = one / sbig
                sumsq = abig
            else if (asml > zero) then
        !
        !     Combine amed and asml if asml > 0.
        !
                if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
                    amed = sqrt(amed)
                    asml = sqrt(asml) / ssml
                    if (asml > amed) then
                        ymin = amed
                        ymax = asml
                    else
                        ymin = asml
                        ymax = amed
                    end if
                    scl = one
                    sumsq = ymax**2*( one + (ymin/ymax)**2 )
                else
                    scl = one / ssml
                    sumsq = asml
                end if
            else
        !
        !     Otherwise all values are mid-range
        !
                scl = one
                sumsq = amed
            end if
            DNRM2 = scl*sqrt( sumsq )
            return
        end function
    
    !==================================================================================================
        SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    !==================================================================================================
        !
        !  -- Reference BLAS level2 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER INCX,LDA,N
            CHARACTER DIAG,TRANS,UPLO
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION A(LDA,*),X(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION ZERO
            PARAMETER (ZERO=0.0D+0)
        !     ..
        !     .. Local Scalars ..
            DOUBLE PRECISION TEMP
            INTEGER I,INFO,IX,J,JX,KX
            LOGICAL NOUNIT
        !     ..
        !     .. External Functions ..
            LOGICAL LSAME
            EXTERNAL LSAME
        !     ..
        !     .. External Subroutines ..
            EXTERNAL XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MAX
        !     ..
        !
        !     Test the input parameters.
        !
            INFO = 0
            IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
                INFO = 1
            ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
                    .NOT.LSAME(TRANS,'C')) THEN
                INFO = 2
            ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
                INFO = 3
            ELSE IF (N < 0) THEN
                INFO = 4
            ELSE IF (LDA < MAX(1,N)) THEN
                INFO = 6
            ELSE IF (INCX == 0) THEN
                INFO = 8
            END IF
            IF (INFO /= 0) THEN
                CALL XERBLA('DTRMV ',INFO)
                RETURN
            END IF
        !
        !     Quick return if possible.
        !
            IF (N == 0) RETURN
        !
            NOUNIT = LSAME(DIAG,'N')
        !
        !     Set up the start point in X if the increment is not unity. This
        !     will be  ( N - 1 )*INCX  too small for descending loops.
        !
            IF (INCX <= 0) THEN
                KX = 1 - (N-1)*INCX
            ELSE IF (INCX /= 1) THEN
                KX = 1
            END IF
        !
        !     Start the operations. In this version the elements of A are
        !     accessed sequentially with one pass through A.
        !
            IF (LSAME(TRANS,'N')) THEN
        !
        !        Form  x := A*x.
        !
                IF (LSAME(UPLO,'U')) THEN
                    IF (INCX == 1) THEN
                        DO 20 J = 1,N
                            IF (X(J) /= ZERO) THEN
                                TEMP = X(J)
                                DO 10 I = 1,J - 1
                                    X(I) = X(I) + TEMP*A(I,J)
        10                     CONTINUE
                                IF (NOUNIT) X(J) = X(J)*A(J,J)
                            END IF
        20             CONTINUE
                    ELSE
                        JX = KX
                        DO 40 J = 1,N
                            IF (X(JX) /= ZERO) THEN
                                TEMP = X(JX)
                                IX = KX
                                DO 30 I = 1,J - 1
                                    X(IX) = X(IX) + TEMP*A(I,J)
                                    IX = IX + INCX
        30                     CONTINUE
                                IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                            END IF
                            JX = JX + INCX
        40             CONTINUE
                    END IF
                ELSE
                    IF (INCX == 1) THEN
                        DO 60 J = N,1,-1
                            IF (X(J) /= ZERO) THEN
                                TEMP = X(J)
                                DO 50 I = N,J + 1,-1
                                    X(I) = X(I) + TEMP*A(I,J)
        50                     CONTINUE
                                IF (NOUNIT) X(J) = X(J)*A(J,J)
                            END IF
        60             CONTINUE
                    ELSE
                        KX = KX + (N-1)*INCX
                        JX = KX
                        DO 80 J = N,1,-1
                            IF (X(JX) /= ZERO) THEN
                                TEMP = X(JX)
                                IX = KX
                                DO 70 I = N,J + 1,-1
                                    X(IX) = X(IX) + TEMP*A(I,J)
                                    IX = IX - INCX
        70                     CONTINUE
                                IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                            END IF
                            JX = JX - INCX
        80             CONTINUE
                    END IF
                END IF
            ELSE
        !
        !        Form  x := A**T*x.
        !
                IF (LSAME(UPLO,'U')) THEN
                    IF (INCX == 1) THEN
                        DO 100 J = N,1,-1
                            TEMP = X(J)
                            IF (NOUNIT) TEMP = TEMP*A(J,J)
                            DO 90 I = J - 1,1,-1
                                TEMP = TEMP + A(I,J)*X(I)
        90                 CONTINUE
                            X(J) = TEMP
        100             CONTINUE
                    ELSE
                        JX = KX + (N-1)*INCX
                        DO 120 J = N,1,-1
                            TEMP = X(JX)
                            IX = JX
                            IF (NOUNIT) TEMP = TEMP*A(J,J)
                            DO 110 I = J - 1,1,-1
                                IX = IX - INCX
                                TEMP = TEMP + A(I,J)*X(IX)
        110                 CONTINUE
                            X(JX) = TEMP
                            JX = JX - INCX
        120             CONTINUE
                    END IF
                ELSE
                    IF (INCX == 1) THEN
                        DO 140 J = 1,N
                            TEMP = X(J)
                            IF (NOUNIT) TEMP = TEMP*A(J,J)
                            DO 130 I = J + 1,N
                                TEMP = TEMP + A(I,J)*X(I)
        130                 CONTINUE
                            X(J) = TEMP
        140             CONTINUE
                    ELSE
                        JX = KX
                        DO 160 J = 1,N
                            TEMP = X(JX)
                            IX = JX
                            IF (NOUNIT) TEMP = TEMP*A(J,J)
                            DO 150 I = J + 1,N
                                IX = IX + INCX
                                TEMP = TEMP + A(I,J)*X(IX)
        150                 CONTINUE
                            X(JX) = TEMP
                            JX = JX + INCX
        160             CONTINUE
                    END IF
                END IF
            END IF
        !
            RETURN
        !
        !     End of DTRMV
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
    !==================================================================================================
        !
        !  -- LAPACK computational routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER            INFO, K, LDA, LWORK, M, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ZERO
            PARAMETER          ( ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            LOGICAL            LQUERY
            INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT, &
                               NB, NBMIN, NX
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DLARFB, DLARFT, DORG2L, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MAX, MIN
        !     ..
        !     .. External Functions ..
            INTEGER            ILAENV
            EXTERNAL           ILAENV
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
            INFO = 0
            LQUERY = ( LWORK == -1 )
            IF( M < 0 ) THEN
                INFO = -1
            ELSE IF( N < 0 .OR. N > M ) THEN
                INFO = -2
            ELSE IF( K < 0 .OR. K > N ) THEN
                INFO = -3
            ELSE IF( LDA < MAX( 1, M ) ) THEN
                INFO = -5
            END IF
        !
            IF( INFO == 0 ) THEN
                IF( N == 0 ) THEN
                    LWKOPT = 1
                ELSE
                    NB = ILAENV( 1, 'DORGQL', ' ', M, N, K, -1 )
                    LWKOPT = N*NB
                END IF
                WORK( 1 ) = LWKOPT
        !
                IF( LWORK < MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
                    INFO = -8
                END IF
            END IF
        !
            IF( INFO /= 0 ) THEN
                CALL XERBLA( 'DORGQL', -INFO )
                RETURN
            ELSE IF( LQUERY ) THEN
                RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( N <= 0 ) THEN
                RETURN
            END IF
        !
            NBMIN = 2
            NX = 0
            IWS = N
            IF( NB > 1 .AND. NB < K ) THEN
        !
        !        Determine when to cross over from blocked to unblocked code.
        !
                NX = MAX( 0, ILAENV( 3, 'DORGQL', ' ', M, N, K, -1 ) )
                IF( NX < K ) THEN
        !
        !           Determine if workspace is large enough for blocked code.
        !
                    LDWORK = N
                    IWS = LDWORK*NB
                    IF( LWORK < IWS ) THEN
        !
        !              Not enough workspace to use optimal NB:  reduce NB and
        !              determine the minimum value of NB.
        !
                    NB = LWORK / LDWORK
                    NBMIN = MAX( 2, ILAENV( 2, 'DORGQL', ' ', M, N, K, -1 ) )
                    END IF
                END IF
            END IF
        !
            IF( NB >= NBMIN .AND. NB < K .AND. NX < K ) THEN
        !
        !        Use blocked code after the first block.
        !        The last kk columns are handled by the block method.
        !
                KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
        !
        !        Set A(m-kk+1:m,1:n-kk) to zero.
        !
                DO 20 J = 1, N - KK
                    DO 10 I = M - KK + 1, M
                    A( I, J ) = ZERO
        10       CONTINUE
        20    CONTINUE
            ELSE
                KK = 0
            END IF
        !
        !     Use unblocked code for the first or only block.
        !
            CALL DORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
        !
            IF( KK > 0 ) THEN
        !
        !        Use blocked code
        !
                DO 50 I = K - KK + 1, K, NB
                    IB = MIN( NB, K-I+1 )
                    IF( N-K+I > 1 ) THEN
        !
        !              Form the triangular factor of the block reflector
        !              H = H(i+ib-1) . . . H(i+1) H(i)
        !
                    CALL DLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB, &
                                A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
        !
        !              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
        !
                    CALL DLARFB( 'Left', 'No transpose', 'Backward', &
                                'Columnwise', M-K+I+IB-1, N-K+I-1, IB, &
                                A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA, &
                                WORK( IB+1 ), LDWORK )
                    END IF
        !
        !           Apply H to rows 1:m-k+i+ib-1 of current block
        !
                    CALL DORG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA, &
                                TAU( I ), WORK, IINFO )
        !
        !           Set rows m-k+i+ib:m of current block to zero
        !
                    DO 40 J = N - K + I, N - K + I + IB - 1
                    DO 30 L = M - K + I + IB, M
                        A( L, J ) = ZERO
        30          CONTINUE
        40       CONTINUE
        50    CONTINUE
            END IF
        !
            WORK( 1 ) = IWS
            RETURN
        !
        !     End of DORGQL
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
    !==================================================================================================
        !
        !  -- LAPACK computational routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER            INFO, K, LDA, LWORK, M, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ZERO
            PARAMETER          ( ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            LOGICAL            LQUERY
            INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, &
                               LWKOPT, NB, NBMIN, NX
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MAX, MIN
        !     ..
        !     .. External Functions ..
            INTEGER            ILAENV
            EXTERNAL           ILAENV
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
            INFO = 0
            NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
            LWKOPT = MAX( 1, N )*NB
            WORK( 1 ) = LWKOPT
            LQUERY = ( LWORK == -1 )
            IF( M < 0 ) THEN
            INFO = -1
            ELSE IF( N < 0 .OR. N > M ) THEN
            INFO = -2
            ELSE IF( K < 0 .OR. K > N ) THEN
            INFO = -3
            ELSE IF( LDA < MAX( 1, M ) ) THEN
            INFO = -5
            ELSE IF( LWORK < MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
            END IF
            IF( INFO /= 0 ) THEN
            CALL XERBLA( 'DORGQR', -INFO )
            RETURN
            ELSE IF( LQUERY ) THEN
            RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( N <= 0 ) THEN
            WORK( 1 ) = 1
            RETURN
            END IF
        !
            NBMIN = 2
            NX = 0
            IWS = N
            IF( NB > 1 .AND. NB < K ) THEN
        !
        !        Determine when to cross over from blocked to unblocked code.
        !
            NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
            IF( NX < K ) THEN
        !
        !           Determine if workspace is large enough for blocked code.
        !
                LDWORK = N
                IWS = LDWORK*NB
                IF( LWORK < IWS ) THEN
        !
        !              Not enough workspace to use optimal NB:  reduce NB and
        !              determine the minimum value of NB.
        !
                    NB = LWORK / LDWORK
                    NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
                END IF
            END IF
            END IF
        !
            IF( NB >= NBMIN .AND. NB < K .AND. NX < K ) THEN
        !
        !        Use blocked code after the last block.
        !        The first kk columns are handled by the block method.
        !
            KI = ( ( K-NX-1 ) / NB )*NB
            KK = MIN( K, KI+NB )
        !
        !        Set A(1:kk,kk+1:n) to zero.
        !
            DO 20 J = KK + 1, N
                DO 10 I = 1, KK
                    A( I, J ) = ZERO
        10       CONTINUE
        20    CONTINUE
            ELSE
            KK = 0
            END IF
        !
        !     Use unblocked code for the last or only block.
        !
            IF( KK < N ) &
                CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, &
                            TAU( KK+1 ), WORK, IINFO )
        !
            IF( KK > 0 ) THEN
        !
        !        Use blocked code
        !
            DO 50 I = KI + 1, 1, -NB
                IB = MIN( NB, K-I+1 )
                IF( I+IB <= N ) THEN
        !
        !              Form the triangular factor of the block reflector
        !              H = H(i) H(i+1) . . . H(i+ib-1)
        !
                    CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                                A( I, I ), LDA, TAU( I ), WORK, LDWORK )
        !
        !              Apply H to A(i:m,i+ib:n) from the left
        !
                    CALL DLARFB( 'Left', 'No transpose', 'Forward', &
                               'Columnwise', M-I+1, N-I-IB+1, IB, &
                                A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                                LDA, WORK( IB+1 ), LDWORK )
                END IF
        !
        !           Apply H to rows i:m of current block
        !
                CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, &
                            IINFO )
        !
        !           Set rows 1:i-1 of current block to zero
        !
                DO 40 J = I, I + IB - 1
                    DO 30 L = 1, I - 1
                        A( L, J ) = ZERO
        30          CONTINUE
        40       CONTINUE
        50    CONTINUE
            END IF
        !
            WORK( 1 ) = IWS
            RETURN
        !
        !     End of DORGQR
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
                                T, LDT, C, LDC, WORK, LDWORK )
    !==================================================================================================
    
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          DIRECT, SIDE, STOREV, TRANS
            INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ), &
                                WORK( LDWORK, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE
            PARAMETER          ( ONE = 1.0D+0 )
        !     ..
        !     .. Local Scalars ..
            CHARACTER          TRANST
            INTEGER            I, J
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            EXTERNAL           LSAME
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DCOPY, DGEMM, DTRMM
        !     ..
        !     .. Executable Statements ..
        !
        !     Quick return if possible
        !
            IF( M <= 0 .OR. N <= 0 ) &
                RETURN
        !
            IF( LSAME( TRANS, 'N' ) ) THEN
                TRANST = 'T'
            ELSE
                TRANST = 'N'
            END IF
        !
            IF( LSAME( STOREV, 'C' ) ) THEN
        !
                IF( LSAME( DIRECT, 'F' ) ) THEN
        !
        !           Let  V =  ( V1 )    (first K rows)
        !                     ( V2 )
        !           where  V1  is unit lower triangular.
        !
                    IF( LSAME( SIDE, 'L' ) ) THEN
        !
        !              Form  H * C  or  H**T * C  where  C = ( C1 )
        !                                                    ( C2 )
        !
        !              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
        !
        !              W := C1**T
        !
                    DO 10 J = 1, K
                        CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
        10          CONTINUE
        !
        !              W := W * V1
        !
                    CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                                K, ONE, V, LDV, WORK, LDWORK )
                    IF( M > K ) THEN
        !
        !                 W := W + C2**T * V2
        !
                        CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                                    ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, &
                                    ONE, WORK, LDWORK )
                    END IF
        !
        !              W := W * T**T  or  W * T
        !
                    CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                                ONE, T, LDT, WORK, LDWORK )
        !
        !              C := C - V * W**T
        !
                    IF( M > K ) THEN
        !
        !                 C2 := C2 - V2 * W**T
        !
                        CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                                    -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, &
                                    C( K+1, 1 ), LDC )
                    END IF
        !
        !              W := W * V1**T
        !
                    CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                                ONE, V, LDV, WORK, LDWORK )
        !
        !              C1 := C1 - W**T
        !
                    DO 30 J = 1, K
                        DO 20 I = 1, N
                            C( J, I ) = C( J, I ) - WORK( I, J )
        20             CONTINUE
        30          CONTINUE
        !
                    ELSE IF( LSAME( SIDE, 'R' ) ) THEN
        !
        !              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
        !
        !              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
        !
        !              W := C1
        !
                    DO 40 J = 1, K
                        CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
        40          CONTINUE
        !
        !              W := W * V1
        !
                    CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                                K, ONE, V, LDV, WORK, LDWORK )
                    IF( N > K ) THEN
        !
        !                 W := W + C2 * V2
        !
                        CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                                    ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                                    ONE, WORK, LDWORK )
                    END IF
        !
        !              W := W * T  or  W * T**T
        !
                    CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                                ONE, T, LDT, WORK, LDWORK )
        !
        !              C := C - W * V**T
        !
                    IF( N > K ) THEN
        !
        !                 C2 := C2 - W * V2**T
        !
                        CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                                    -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                                    C( 1, K+1 ), LDC )
                    END IF
        !
        !              W := W * V1**T
        !
                    CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                                ONE, V, LDV, WORK, LDWORK )
        !
        !              C1 := C1 - W
        !
                    DO 60 J = 1, K
                        DO 50 I = 1, M
                            C( I, J ) = C( I, J ) - WORK( I, J )
        50             CONTINUE
        60          CONTINUE
                    END IF
        !
                ELSE
        !
        !           Let  V =  ( V1 )
        !                     ( V2 )    (last K rows)
        !           where  V2  is unit upper triangular.
        !
                    IF( LSAME( SIDE, 'L' ) ) THEN
        !
        !              Form  H * C  or  H**T * C  where  C = ( C1 )
        !                                                    ( C2 )
        !
        !              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
        !
        !              W := C2**T
        !
                    DO 70 J = 1, K
                        CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
        70          CONTINUE
        !
        !              W := W * V2
        !
                    CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                                K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
                    IF( M > K ) THEN
        !
        !                 W := W + C1**T * V1
        !
                        CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                                    ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
                    END IF
        !
        !              W := W * T**T  or  W * T
        !
                    CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                                ONE, T, LDT, WORK, LDWORK )
        !
        !              C := C - V * W**T
        !
                    IF( M > K ) THEN
        !
        !                 C1 := C1 - V1 * W**T
        !
                        CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                                    -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
                    END IF
        !
        !              W := W * V2**T
        !
                    CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                                ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
        !
        !              C2 := C2 - W**T
        !
                    DO 90 J = 1, K
                        DO 80 I = 1, N
                            C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
        80             CONTINUE
        90          CONTINUE
        !
                    ELSE IF( LSAME( SIDE, 'R' ) ) THEN
        !
        !              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
        !
        !              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
        !
        !              W := C2
        !
                    DO 100 J = 1, K
                        CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
        100          CONTINUE
        !
        !              W := W * V2
        !
                    CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                                K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
                    IF( N > K ) THEN
        !
        !                 W := W + C1 * V1
        !
                        CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                                    ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
                    END IF
        !
        !              W := W * T  or  W * T**T
        !
                    CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                                ONE, T, LDT, WORK, LDWORK )
        !
        !              C := C - W * V**T
        !
                    IF( N > K ) THEN
        !
        !                 C1 := C1 - W * V1**T
        !
                        CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                                    -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
                    END IF
        !
        !              W := W * V2**T
        !
                    CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                                ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
        !
        !              C2 := C2 - W
        !
                    DO 120 J = 1, K
                        DO 110 I = 1, M
                            C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
        110             CONTINUE
        120          CONTINUE
                    END IF
                END IF
        !
            ELSE IF( LSAME( STOREV, 'R' ) ) THEN
        !
                IF( LSAME( DIRECT, 'F' ) ) THEN
        !
        !           Let  V =  ( V1  V2 )    (V1: first K columns)
        !           where  V1  is unit upper triangular.
        !
                    IF( LSAME( SIDE, 'L' ) ) THEN
        !
        !              Form  H * C  or  H**T * C  where  C = ( C1 )
        !                                                    ( C2 )
        !
        !              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
        !
        !              W := C1**T
        !
                    DO 130 J = 1, K
                        CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
        130          CONTINUE
        !
        !              W := W * V1**T
        !
                    CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                                ONE, V, LDV, WORK, LDWORK )
                    IF( M > K ) THEN
        !
        !                 W := W + C2**T * V2**T
        !
                        CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                                    C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, &
                                    WORK, LDWORK )
                    END IF
        !
        !              W := W * T**T  or  W * T
        !
                    CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                                ONE, T, LDT, WORK, LDWORK )
        !
        !              C := C - V**T * W**T
        !
                    IF( M > K ) THEN
        !
        !                 C2 := C2 - V2**T * W**T
        !
                        CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                                    V( 1, K+1 ), LDV, WORK, LDWORK, ONE, &
                                    C( K+1, 1 ), LDC )
                    END IF
        !
        !              W := W * V1
        !
                    CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                                K, ONE, V, LDV, WORK, LDWORK )
        !
        !              C1 := C1 - W**T
        !
                    DO 150 J = 1, K
                        DO 140 I = 1, N
                            C( J, I ) = C( J, I ) - WORK( I, J )
        140             CONTINUE
        150          CONTINUE
        !
                    ELSE IF( LSAME( SIDE, 'R' ) ) THEN
        !
        !              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
        !
        !              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
        !
        !              W := C1
        !
                    DO 160 J = 1, K
                        CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
        160          CONTINUE
        !
        !              W := W * V1**T
        !
                    CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                                ONE, V, LDV, WORK, LDWORK )
                    IF( N > K ) THEN
        !
        !                 W := W + C2 * V2**T
        !
                        CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                                    ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, &
                                    ONE, WORK, LDWORK )
                    END IF
        !
        !              W := W * T  or  W * T**T
        !
                    CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                                ONE, T, LDT, WORK, LDWORK )
        !
        !              C := C - W * V
        !
                    IF( N > K ) THEN
        !
        !                 C2 := C2 - W * V2
        !
                        CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                                    -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
                                    C( 1, K+1 ), LDC )
                    END IF
        !
        !              W := W * V1
        !
                    CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                                K, ONE, V, LDV, WORK, LDWORK )
        !
        !              C1 := C1 - W
        !
                    DO 180 J = 1, K
                        DO 170 I = 1, M
                            C( I, J ) = C( I, J ) - WORK( I, J )
        170             CONTINUE
        180          CONTINUE
        !
                    END IF
        !
                ELSE
        !
        !           Let  V =  ( V1  V2 )    (V2: last K columns)
        !           where  V2  is unit lower triangular.
        !
                    IF( LSAME( SIDE, 'L' ) ) THEN
        !
        !              Form  H * C  or  H**T * C  where  C = ( C1 )
        !                                                    ( C2 )
        !
        !              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
        !
        !              W := C2**T
        !
                    DO 190 J = 1, K
                        CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
        190          CONTINUE
        !
        !              W := W * V2**T
        !
                    CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                                ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
                    IF( M > K ) THEN
        !
        !                 W := W + C1**T * V1**T
        !
                        CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                                    C, LDC, V, LDV, ONE, WORK, LDWORK )
                    END IF
        !
        !              W := W * T**T  or  W * T
        !
                    CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                                ONE, T, LDT, WORK, LDWORK )
        !
        !              C := C - V**T * W**T
        !
                    IF( M > K ) THEN
        !
        !                 C1 := C1 - V1**T * W**T
        !
                        CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                                    V, LDV, WORK, LDWORK, ONE, C, LDC )
                    END IF
        !
        !              W := W * V2
        !
                    CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                                K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
        !
        !              C2 := C2 - W**T
        !
                    DO 210 J = 1, K
                        DO 200 I = 1, N
                            C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
        200             CONTINUE
        210          CONTINUE
        !
                    ELSE IF( LSAME( SIDE, 'R' ) ) THEN
        !
        !              Form  C * H  or  C * H'  where  C = ( C1  C2 )
        !
        !              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
        !
        !              W := C2
        !
                    DO 220 J = 1, K
                        CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
        220          CONTINUE
        !
        !              W := W * V2**T
        !
                    CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                                ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
                    IF( N > K ) THEN
        !
        !                 W := W + C1 * V1**T
        !
                        CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                                    ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
                    END IF
        !
        !              W := W * T  or  W * T**T
        !
                    CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                                ONE, T, LDT, WORK, LDWORK )
        !
        !              C := C - W * V
        !
                    IF( N > K ) THEN
        !
        !                 C1 := C1 - W * V1
        !
                        CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                                    -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
                    END IF
        !
        !              W := W * V2
        !
                    CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                                K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
        !
        !              C1 := C1 - W
        !
                    DO 240 J = 1, K
                        DO 230 I = 1, M
                            C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
        230             CONTINUE
        240          CONTINUE
        !
                    END IF
        !
                END IF
            END IF
        !
            RETURN
        !
        !     End of DLARFB
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
    !==================================================================================================
        !
        !  -- Reference BLAS level1 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER INCX,INCY,N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION DX(*),DY(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
            INTEGER I,IX,IY,M,MP1
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MOD
        !     ..
            IF (N <= 0) RETURN
            IF (INCX == 1 .AND. INCY == 1) THEN
        !
        !        code for both increments equal to 1
        !
        !
        !        clean-up loop
        !
                M = MOD(N,7)
                IF (M /= 0) THEN
                    DO I = 1,M
                    DY(I) = DX(I)
                    END DO
                    IF (N < 7) RETURN
                END IF
                MP1 = M + 1
                DO I = MP1,N,7
                    DY(I) = DX(I)
                    DY(I+1) = DX(I+1)
                    DY(I+2) = DX(I+2)
                    DY(I+3) = DX(I+3)
                    DY(I+4) = DX(I+4)
                    DY(I+5) = DX(I+5)
                    DY(I+6) = DX(I+6)
                END DO
            ELSE
        !
        !        code for unequal increments or equal increments
        !          not equal to 1
        !
                IX = 1
                IY = 1
                IF (INCX < 0) IX = (-N+1)*INCX + 1
                IF (INCY < 0) IY = (-N+1)*INCY + 1
                DO I = 1,N
                    DY(IY) = DX(IX)
                    IX = IX + INCX
                    IY = IY + INCY
                END DO
            END IF
            RETURN
        !
        !     End of DCOPY
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
    !==================================================================================================
        !
        !  -- LAPACK computational routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER            INFO, K, LDA, M, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE, ZERO
            PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            I, II, J, L
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DLARF, DSCAL, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MAX
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
            INFO = 0
            IF( M < 0 ) THEN
                INFO = -1
            ELSE IF( N < 0 .OR. N > M ) THEN
                INFO = -2
            ELSE IF( K < 0 .OR. K > N ) THEN
                INFO = -3
            ELSE IF( LDA < MAX( 1, M ) ) THEN
                INFO = -5
            END IF
            IF( INFO /= 0 ) THEN
                CALL XERBLA( 'DORG2L', -INFO )
                RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( N <= 0 ) &
                RETURN
        !
        !     Initialise columns 1:n-k to columns of the unit matrix
        !
            DO 20 J = 1, N - K
                DO 10 L = 1, M
                    A( L, J ) = ZERO
        10    CONTINUE
                A( M-N+J, J ) = ONE
        20 CONTINUE
        !
            DO 40 I = 1, K
                II = N - K + I
        !
        !        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
        !
                A( M-N+II, II ) = ONE
                CALL DLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A, &
                            LDA, WORK )
                CALL DSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
                A( M-N+II, II ) = ONE - TAU( I )
        !
        !        Set A(m-k+i+1:m,n-k+i) to zero
        !
                DO 30 L = M - N + II + 1, M
                    A( L, II ) = ZERO
        30    CONTINUE
        40 CONTINUE
            RETURN
        !
        !     End of DORG2L
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
    !==================================================================================================
        !
        !  -- LAPACK computational routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER            INFO, K, LDA, M, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE, ZERO
            PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            I, J, L
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DLARF, DSCAL, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MAX
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
            INFO = 0
            IF( M < 0 ) THEN
                INFO = -1
            ELSE IF( N < 0 .OR. N > M ) THEN
                INFO = -2
            ELSE IF( K < 0 .OR. K > N ) THEN
                INFO = -3
            ELSE IF( LDA < MAX( 1, M ) ) THEN
                INFO = -5
            END IF
            IF( INFO /= 0 ) THEN
                CALL XERBLA( 'DORG2R', -INFO )
                RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( N <= 0 ) &
                RETURN
        !
        !     Initialise columns k+1:n to columns of the unit matrix
        !
            DO 20 J = K + 1, N
                DO 10 L = 1, M
                    A( L, J ) = ZERO
        10    CONTINUE
                A( J, J ) = ONE
        20 CONTINUE
        !
            DO 40 I = K, 1, -1
        !
        !        Apply H(i) to A(i:m,i:n) from the left
        !
                IF( I < N ) THEN
                    A( I, I ) = ONE
                    CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                                A( I, I+1 ), LDA, WORK )
                END IF
                IF( I < M ) &
                    CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
                A( I, I ) = ONE - TAU( I )
        !
        !        Set A(1:i-1,i) to zero
        !
                DO 30 L = 1, I - 1
                    A( L, I ) = ZERO
        30    CONTINUE
        40 CONTINUE
            RETURN
        !
        !     End of DORG2R
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          DIRECT, STOREV
            INTEGER            K, LDT, LDV, N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE, ZERO
            PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER            I, J, PREVLASTV, LASTV
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DGEMV, DTRMV
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            EXTERNAL           LSAME
        !     ..
        !     .. Executable Statements ..
        !
        !     Quick return if possible
        !
            IF( N == 0 ) &
                RETURN
        !
            IF( LSAME( DIRECT, 'F' ) ) THEN
                PREVLASTV = N
                DO I = 1, K
                    PREVLASTV = MAX( I, PREVLASTV )
                    IF( TAU( I ) == ZERO ) THEN
        !
        !              H(i)  =  I
        !
                    DO J = 1, I
                        T( J, I ) = ZERO
                    END DO
                    ELSE
        !
        !              general case
        !
                    IF( LSAME( STOREV, 'C' ) ) THEN
        !                 Skip any trailing zeros.
                        DO LASTV = N, I+1, -1
                            IF( V( LASTV, I ) /= ZERO ) EXIT
                        END DO
                        DO J = 1, I-1
                            T( J, I ) = -TAU( I ) * V( I , J )
                        END DO
                        J = MIN( LASTV, PREVLASTV )
        !
        !                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)
        !
                        CALL DGEMV( 'Transpose', J-I, I-1, -TAU( I ), &
                                    V( I+1, 1 ), LDV, V( I+1, I ), 1, ONE, &
                                    T( 1, I ), 1 )
                    ELSE
        !                 Skip any trailing zeros.
                        DO LASTV = N, I+1, -1
                            IF( V( I, LASTV ) /= ZERO ) EXIT
                        END DO
                        DO J = 1, I-1
                            T( J, I ) = -TAU( I ) * V( J , I )
                        END DO
                        J = MIN( LASTV, PREVLASTV )
        !
        !                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T
        !
                        CALL DGEMV( 'No transpose', I-1, J-I, -TAU( I ), &
                                    V( 1, I+1 ), LDV, V( I, I+1 ), LDV, ONE, &
                                    T( 1, I ), 1 )
                    END IF
        !
        !              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
        !
                    CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
                                LDT, T( 1, I ), 1 )
                    T( I, I ) = TAU( I )
                    IF( I > 1 ) THEN
                        PREVLASTV = MAX( PREVLASTV, LASTV )
                    ELSE
                        PREVLASTV = LASTV
                    END IF
                    END IF
                END DO
            ELSE
                PREVLASTV = 1
                DO I = K, 1, -1
                    IF( TAU( I ) == ZERO ) THEN
        !
        !              H(i)  =  I
        !
                    DO J = I, K
                        T( J, I ) = ZERO
                    END DO
                    ELSE
        !
        !              general case
        !
                    IF( I < K ) THEN
                        IF( LSAME( STOREV, 'C' ) ) THEN
        !                    Skip any leading zeros.
                            DO LASTV = 1, I-1
                                IF( V( LASTV, I ) /= ZERO ) EXIT
                            END DO
                            DO J = I+1, K
                                T( J, I ) = -TAU( I ) * V( N-K+I , J )
                            END DO
                            J = MAX( LASTV, PREVLASTV )
        !
        !                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)
        !
                            CALL DGEMV( 'Transpose', N-K+I-J, K-I, -TAU( I ), &
                                        V( J, I+1 ), LDV, V( J, I ), 1, ONE, &
                                        T( I+1, I ), 1 )
                        ELSE
        !                    Skip any leading zeros.
                            DO LASTV = 1, I-1
                                IF( V( I, LASTV ) /= ZERO ) EXIT
                            END DO
                            DO J = I+1, K
                                T( J, I ) = -TAU( I ) * V( J, N-K+I )
                            END DO
                            J = MAX( LASTV, PREVLASTV )
        !
        !                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T
        !
                            CALL DGEMV( 'No transpose', K-I, N-K+I-J, &
                                -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV, &
                                ONE, T( I+1, I ), 1 )
                        END IF
        !
        !                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
        !
                        CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
                                    T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
                        IF( I > 1 ) THEN
                            PREVLASTV = MIN( PREVLASTV, LASTV )
                        ELSE
                            PREVLASTV = LASTV
                        END IF
                    END IF
                    T( I, I ) = TAU( I )
                    END IF
                END DO
            END IF
            RETURN
        !
        !     End of DLARFT
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            CHARACTER          SIDE
            INTEGER            INCV, LDC, M, N
            DOUBLE PRECISION   TAU
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION   ONE, ZERO
            PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            LOGICAL            APPLYLEFT
            INTEGER            I, LASTV, LASTC
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           DGEMV, DGER
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            INTEGER            ILADLR, ILADLC
            EXTERNAL           LSAME, ILADLR, ILADLC
        !     ..
        !     .. Executable Statements ..
        !
            APPLYLEFT = LSAME( SIDE, 'L' )
            LASTV = 0
            LASTC = 0
            IF( TAU /= ZERO ) THEN
        !     Set up variables for scanning V.  LASTV begins pointing to the end
        !     of V.
                IF( APPLYLEFT ) THEN
                    LASTV = M
                ELSE
                    LASTV = N
                END IF
                IF( INCV > 0 ) THEN
                    I = 1 + (LASTV-1) * INCV
                ELSE
                    I = 1
                END IF
        !     Look for the last non-zero row in V.
                DO WHILE( LASTV > 0 .AND. V( I ) == ZERO )
                    LASTV = LASTV - 1
                    I = I - INCV
                END DO
                IF( APPLYLEFT ) THEN
        !     Scan for the last non-zero column in C(1:lastv,:).
                    LASTC = ILADLC(LASTV, N, C, LDC)
                ELSE
        !     Scan for the last non-zero row in C(:,1:lastv).
                    LASTC = ILADLR(M, LASTV, C, LDC)
                END IF
            END IF
        !     Note that lastc == 0 renders the BLAS operations null; no special
        !     case is needed at this level.
            IF( APPLYLEFT ) THEN
        !
        !        Form  H * C
        !
                IF( LASTV > 0 ) THEN
        !
        !           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
        !
                    CALL DGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV, &
                                ZERO, WORK, 1 )
        !
        !           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
        !
                    CALL DGER( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
                END IF
            ELSE
        !
        !        Form  C * H
        !
                IF( LASTV > 0 ) THEN
        !
        !           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
        !
                    CALL DGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC, &
                                V, INCV, ZERO, WORK, 1 )
        !
        !           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
        !
                    CALL DGER( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
                END IF
            END IF
            RETURN
        !
        !     End of DLARF
        !
        END SUBROUTINE
    
    !==================================================================================================
        SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
    !==================================================================================================
        !
        !  -- Reference BLAS level2 routine --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            DOUBLE PRECISION ALPHA
            INTEGER INCX,INCY,LDA,M,N
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION A(LDA,*),X(*),Y(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION ZERO
            PARAMETER (ZERO=0.0D+0)
        !     ..
        !     .. Local Scalars ..
            DOUBLE PRECISION TEMP
            INTEGER I,INFO,IX,J,JY,KX
        !     ..
        !     .. External Subroutines ..
            EXTERNAL XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC MAX
        !     ..
        !
        !     Test the input parameters.
        !
            INFO = 0
            IF (M < 0) THEN
                INFO = 1
            ELSE IF (N < 0) THEN
                INFO = 2
            ELSE IF (INCX == 0) THEN
                INFO = 5
            ELSE IF (INCY == 0) THEN
                INFO = 7
            ELSE IF (LDA < MAX(1,M)) THEN
                INFO = 9
            END IF
            IF (INFO /= 0) THEN
                CALL XERBLA('DGER  ',INFO)
                RETURN
            END IF
        !
        !     Quick return if possible.
        !
            IF ((M == 0) .OR. (N == 0) .OR. (ALPHA == ZERO)) RETURN
        !
        !     Start the operations. In this version the elements of A are
        !     accessed sequentially with one pass through A.
        !
            IF (INCY > 0) THEN
                JY = 1
            ELSE
                JY = 1 - (N-1)*INCY
            END IF
            IF (INCX == 1) THEN
                DO 20 J = 1,N
                    IF (Y(JY) /= ZERO) THEN
                        TEMP = ALPHA*Y(JY)
                        DO 10 I = 1,M
                            A(I,J) = A(I,J) + X(I)*TEMP
        10             CONTINUE
                    END IF
                    JY = JY + INCY
        20     CONTINUE
            ELSE
                IF (INCX > 0) THEN
                    KX = 1
                ELSE
                    KX = 1 - (M-1)*INCX
                END IF
                DO 40 J = 1,N
                    IF (Y(JY) /= ZERO) THEN
                        TEMP = ALPHA*Y(JY)
                        IX = KX
                        DO 30 I = 1,M
                            A(I,J) = A(I,J) + X(IX)*TEMP
                            IX = IX + INCX
        30             CONTINUE
                    END IF
                    JY = JY + INCY
        40     CONTINUE
            END IF
        !
            RETURN
        !
        !     End of DGER
        !
        END SUBROUTINE
    
    !==================================================================================================
        INTEGER FUNCTION ILADLC( M, N, A, LDA )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER            M, N, LDA
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION ZERO
            PARAMETER ( ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER I
        !     ..
        !     .. Executable Statements ..
        !
        !     Quick test for the common case where one corner is non-zero.
            IF( N == 0 ) THEN
                ILADLC = N
            ELSE IF( A(1, N) /= ZERO .OR. A(M, N) /= ZERO ) THEN
                ILADLC = N
            ELSE
        !     Now scan each column from the end, returning with the first non-zero.
                DO ILADLC = N, 1, -1
                    DO I = 1, M
                    IF( A(I, ILADLC) /= ZERO ) RETURN
                    END DO
                END DO
            END IF
            RETURN
        END FUNCTION
    
    !==================================================================================================
        INTEGER FUNCTION ILADLR( M, N, A, LDA )
    !==================================================================================================
        !
        !  -- LAPACK auxiliary routine --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !
        !     .. Scalar Arguments ..
            INTEGER            M, N, LDA
        !     ..
        !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            DOUBLE PRECISION ZERO
            PARAMETER ( ZERO = 0.0D+0 )
        !     ..
        !     .. Local Scalars ..
            INTEGER I, J
        !     ..
        !     .. Executable Statements ..
        !
        !     Quick test for the common case where one corner is non-zero.
            IF( M == 0 ) THEN
                ILADLR = M
            ELSE IF( A(M, 1) /= ZERO .OR. A(M, N) /= ZERO ) THEN
                ILADLR = M
            ELSE
        !     Scan up each column tracking the last zero row seen.
                ILADLR = 0
                DO J = 1, N
                    I=M
                    DO WHILE((A(MAX(I,1),J) == ZERO).AND.(I >= 1))
                    I=I-1
                    ENDDO
                    ILADLR = MAX( ILADLR, I )
                END DO
            END IF
            RETURN
        END FUNCTION