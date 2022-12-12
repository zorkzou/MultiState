      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JB, NB
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEMM, DPOTRF2, DSYRK, DTRSM, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRF', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 )
     $   RETURN
C
C     Determine the block size for this environment.
C
      NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
C
C        Use unblocked code.
C
         CALL DPOTRF2( UPLO, N, A, LDA, INFO )
      ELSE
C
C        Use blocked code.
C
         IF( UPPER ) THEN
C
C           Compute the Cholesky factorization A = U**T*U.
C
            DO 10 J = 1, N, NB
C
C              Update and factorize the current diagonal block and test
C              for non-positive-definiteness.
C
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Upper', 'Transpose', JB, J-1, -ONE,
     $                     A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTRF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
C
C                 Compute the current block row.
C
                  CALL DGEMM( 'Transpose', 'No transpose', JB, N-J-JB+1,
     $                        J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ),
     $                        LDA, ONE, A( J, J+JB ), LDA )
                  CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit',
     $                        JB, N-J-JB+1, ONE, A( J, J ), LDA,
     $                        A( J, J+JB ), LDA )
               END IF
   10       CONTINUE
C
         ELSE
C
C           Compute the Cholesky factorization A = L*L**T.
C
            DO 20 J = 1, N, NB
C
C              Update and factorize the current diagonal block and test
C              for non-positive-definiteness.
C
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Lower', 'No transpose', JB, J-1, -ONE,
     $                     A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTRF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
C
C                 Compute the current block column.
C
                  CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,
     $                        J-1, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ),
     $                        LDA, ONE, A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'Transpose', 'Non-unit',
     $                        N-J-JB+1, JB, ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
C
   30 CONTINUE
      INFO = INFO + J - 1
C
   40 CONTINUE
      RETURN
C
C     End of DPOTRF
C
      END

      SUBROUTINE DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D0, HALF = 0.5D0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K
      DOUBLE PRECISION   AKK, BKK, CT
C     ..
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, DSYR2, DTRMV, DTRSV, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYGS2', -INFO )
         RETURN
      END IF
C
      IF( ITYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
C
C           Compute inv(U**T)*A*inv(U)
C
            DO 10 K = 1, N
C
C              Update the upper triangle of A(k:n,k:n)
C
               AKK = A( K, K )
               BKK = B( K, K )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL DSCAL( N-K, ONE / BKK, A( K, K+1 ), LDA )
                  CT = -HALF*AKK
                  CALL DAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ),
     $                        LDA )
                  CALL DSYR2( UPLO, N-K, -ONE, A( K, K+1 ), LDA,
     $                        B( K, K+1 ), LDB, A( K+1, K+1 ), LDA )
                  CALL DAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ),
     $                        LDA )
                  CALL DTRSV( UPLO, 'Transpose', 'Non-unit', N-K,
     $                        B( K+1, K+1 ), LDB, A( K, K+1 ), LDA )
               END IF
   10       CONTINUE
         ELSE
C
C           Compute inv(L)*A*inv(L**T)
C
            DO 20 K = 1, N
C
C              Update the lower triangle of A(k:n,k:n)
C
               AKK = A( K, K )
               BKK = B( K, K )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL DSCAL( N-K, ONE / BKK, A( K+1, K ), 1 )
                  CT = -HALF*AKK
                  CALL DAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL DSYR2( UPLO, N-K, -ONE, A( K+1, K ), 1,
     $                        B( K+1, K ), 1, A( K+1, K+1 ), LDA )
                  CALL DAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL DTRSV( UPLO, 'No transpose', 'Non-unit', N-K,
     $                        B( K+1, K+1 ), LDB, A( K+1, K ), 1 )
               END IF
   20       CONTINUE
         END IF
      ELSE
         IF( UPPER ) THEN
C
C           Compute U*A*U**T
C
            DO 30 K = 1, N
C
C              Update the upper triangle of A(1:k,1:k)
C
               AKK = A( K, K )
               BKK = B( K, K )
               CALL DTRMV( UPLO, 'No transpose', 'Non-unit', K-1, B,
     $                     LDB, A( 1, K ), 1 )
               CT = HALF*AKK
               CALL DAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL DSYR2( UPLO, K-1, ONE, A( 1, K ), 1, B( 1, K ), 1,
     $                     A, LDA )
               CALL DAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL DSCAL( K-1, BKK, A( 1, K ), 1 )
               A( K, K ) = AKK*BKK**2
   30       CONTINUE
         ELSE
C
C           Compute L**T *A*L
C
            DO 40 K = 1, N
C
C              Update the lower triangle of A(1:k,1:k)
C
               AKK = A( K, K )
               BKK = B( K, K )
               CALL DTRMV( UPLO, 'Transpose', 'Non-unit', K-1, B, LDB,
     $                     A( K, 1 ), LDA )
               CT = HALF*AKK
               CALL DAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL DSYR2( UPLO, K-1, ONE, A( K, 1 ), LDA, B( K, 1 ),
     $                     LDB, A, LDA )
               CALL DAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL DSCAL( K-1, BKK, A( K, 1 ), LDA )
               A( K, K ) = AKK*BKK**2
   40       CONTINUE
         END IF
      END IF
      RETURN
C
C     End of DSYGS2
C
      END

      SUBROUTINE DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D0, HALF = 0.5D0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, KB, NB
C     ..
C     .. External Subroutines ..
      EXTERNAL           DSYGS2, DSYMM, DSYR2K, DTRMM, DTRSM, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYGST', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 )
     $   RETURN
C
C     Determine the block size for this environment.
C
      NB = ILAENV( 1, 'DSYGST', UPLO, N, -1, -1, -1 )
C
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
C
C        Use unblocked code
C
         CALL DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      ELSE
C
C        Use blocked code
C
         IF( ITYPE.EQ.1 ) THEN
            IF( UPPER ) THEN
C
C              Compute inv(U**T)*A*inv(U)
C
               DO 10 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
C
C                 Update the upper triangle of A(k:n,k:n)
C
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL DTRSM( 'Left', UPLO, 'Transpose', 'Non-unit',
     $                           KB, N-K-KB+1, ONE, B( K, K ), LDB,
     $                           A( K, K+KB ), LDA )
                     CALL DSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF,
     $                           A( K, K ), LDA, B( K, K+KB ), LDB, ONE,
     $                           A( K, K+KB ), LDA )
                     CALL DSYR2K( UPLO, 'Transpose', N-K-KB+1, KB, -ONE,
     $                            A( K, K+KB ), LDA, B( K, K+KB ), LDB,
     $                            ONE, A( K+KB, K+KB ), LDA )
                     CALL DSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF,
     $                           A( K, K ), LDA, B( K, K+KB ), LDB, ONE,
     $                           A( K, K+KB ), LDA )
                     CALL DTRSM( 'Right', UPLO, 'No transpose',
     $                           'Non-unit', KB, N-K-KB+1, ONE,
     $                           B( K+KB, K+KB ), LDB, A( K, K+KB ),
     $                           LDA )
                  END IF
   10          CONTINUE
            ELSE
C
C              Compute inv(L)*A*inv(L**T)
C
               DO 20 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
C
C                 Update the lower triangle of A(k:n,k:n)
C
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL DTRSM( 'Right', UPLO, 'Transpose', 'Non-unit',
     $                           N-K-KB+1, KB, ONE, B( K, K ), LDB,
     $                           A( K+KB, K ), LDA )
                     CALL DSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF,
     $                           A( K, K ), LDA, B( K+KB, K ), LDB, ONE,
     $                           A( K+KB, K ), LDA )
                     CALL DSYR2K( UPLO, 'No transpose', N-K-KB+1, KB,
     $                            -ONE, A( K+KB, K ), LDA, B( K+KB, K ),
     $                            LDB, ONE, A( K+KB, K+KB ), LDA )
                     CALL DSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF,
     $                           A( K, K ), LDA, B( K+KB, K ), LDB, ONE,
     $                           A( K+KB, K ), LDA )
                     CALL DTRSM( 'Left', UPLO, 'No transpose',
     $                           'Non-unit', N-K-KB+1, KB, ONE,
     $                           B( K+KB, K+KB ), LDB, A( K+KB, K ),
     $                           LDA )
                  END IF
   20          CONTINUE
            END IF
         ELSE
            IF( UPPER ) THEN
C
C              Compute U*A*U**T
C
               DO 30 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
C
C                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
C
                  CALL DTRMM( 'Left', UPLO, 'No transpose', 'Non-unit',
     $                        K-1, KB, ONE, B, LDB, A( 1, K ), LDA )
                  CALL DSYMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ),
     $                        LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA )
                  CALL DSYR2K( UPLO, 'No transpose', K-1, KB, ONE,
     $                         A( 1, K ), LDA, B( 1, K ), LDB, ONE, A,
     $                         LDA )
                  CALL DSYMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ),
     $                        LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA )
                  CALL DTRMM( 'Right', UPLO, 'Transpose', 'Non-unit',
     $                        K-1, KB, ONE, B( K, K ), LDB, A( 1, K ),
     $                        LDA )
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
   30          CONTINUE
            ELSE
C
C              Compute L**T*A*L
C
               DO 40 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
C
C                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
C
                  CALL DTRMM( 'Right', UPLO, 'No transpose', 'Non-unit',
     $                        KB, K-1, ONE, B, LDB, A( K, 1 ), LDA )
                  CALL DSYMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ),
     $                        LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA )
                  CALL DSYR2K( UPLO, 'Transpose', K-1, KB, ONE,
     $                         A( K, 1 ), LDA, B( K, 1 ), LDB, ONE, A,
     $                         LDA )
                  CALL DSYMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ),
     $                        LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA )
                  CALL DTRMM( 'Left', UPLO, 'Transpose', 'Non-unit', KB,
     $                        K-1, ONE, B( K, K ), LDB, A( K, 1 ), LDA )
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
C
C     End of DSYGST
C
      END

      SUBROUTINE DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
     $                  LWORK, INFO )
C
C  -- LAPACK driver routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER, WANTZ
      CHARACTER          TRANS
      INTEGER            LWKMIN, LWKOPT, NB, NEIG
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
C     ..
C     .. External Subroutines ..
      EXTERNAL           DPOTRF, DSYEV, DSYGST, DTRMM, DTRSM, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
C
      INFO = 0
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
C
      IF( INFO.EQ.0 ) THEN
         LWKMIN = MAX( 1, 3*N - 1 )
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( LWKMIN, ( NB + 2 )*N )
         WORK( 1 ) = LWKOPT
C
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -11
         END IF
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYGV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 )
     $   RETURN
C
C     Form a Cholesky factorization of B.
C
      CALL DPOTRF( UPLO, N, B, LDB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
C
C     Transform problem to standard eigenvalue problem and solve.
C
      CALL DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      CALL DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
C
      IF( WANTZ ) THEN
C
C        Backtransform eigenvectors to the original problem.
C
         NEIG = N
         IF( INFO.GT.0 )
     $      NEIG = INFO - 1
         IF( ITYPE.EQ.1 .OR. ITYPE.EQ.2 ) THEN
C
C           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
C           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
C
            IF( UPPER ) THEN
               TRANS = 'N'
            ELSE
               TRANS = 'T'
            END IF
C
            CALL DTRSM( 'Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE,
     $                  B, LDB, A, LDA )
C
         ELSE IF( ITYPE.EQ.3 ) THEN
C
C           For B*A*x=(lambda)*x;
C           backtransform eigenvectors: x = L*y or U**T*y
C
            IF( UPPER ) THEN
               TRANS = 'T'
            ELSE
               TRANS = 'N'
            END IF
C
            CALL DTRMM( 'Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE,
     $                  B, LDB, A, LDA )
         END IF
      END IF
C
      WORK( 1 ) = LWKOPT
      RETURN
C
C     End of DSYGV
C
      END

      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
C
C  -- LAPACK driver routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE,
     $                   LLWORK, LWKOPT, NB
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANSY
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLASCL, DORGTR, DSCAL, DSTEQR, DSTERF, DSYTRD,
     $                   XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
C
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
C
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, ( NB+2 )*N )
         WORK( 1 ) = LWKOPT
C
         IF( LWORK.LT.MAX( 1, 3*N-1 ) .AND. .NOT.LQUERY )
     $      INFO = -8
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 ) THEN
         RETURN
      END IF
C
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         WORK( 1 ) = 2
         IF( WANTZ )
     $      A( 1, 1 ) = ONE
         RETURN
      END IF
C
C     Get machine constants.
C
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
C
C     Scale matrix to allowable range, if necessary.
C
      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 )
     $   CALL DLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
C
C     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
C
      INDE = 1
      INDTAU = INDE + N
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      CALL DSYTRD( UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ),
     $             WORK( INDWRK ), LLWORK, IINFO )
C
C     For eigenvalues only, call DSTERF.  For eigenvectors, first call
C     DORGTR to generate the orthogonal matrix, then call DSTEQR.
C
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      ELSE
         CALL DORGTR( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ),
     $                LLWORK, IINFO )
         CALL DSTEQR( JOBZ, N, W, WORK( INDE ), A, LDA, WORK( INDTAU ),
     $                INFO )
      END IF
C
C     If matrix was scaled, then rescale eigenvalues appropriately.
C
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
C
C     Set WORK(1) to optimal workspace size.
C
      WORK( 1 ) = LWKOPT
C
      RETURN
C
C     End of DSYEV
C
      END

      SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ),
     $                   WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLATRD, DSYR2K, DSYTD2, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
C
      IF( INFO.EQ.0 ) THEN
C
C        Determine the block size.
C
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
C
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
C
C        Determine when to cross over from blocked to unblocked code
C        (last block is always handled by unblocked code).
C
         NX = MAX( NB, ILAENV( 3, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
C
C              Not enough workspace to use optimal NB:  determine the
C              minimum value of NB, and reduce NB or force use of
C              unblocked code by setting NX = N.
C
               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = ILAENV( 2, 'DSYTRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN )
     $            NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
C
      IF( UPPER ) THEN
C
C        Reduce the upper triangle of A.
C        Columns 1:kk are handled by the unblocked method.
C
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
C
C           Reduce columns i:i+nb-1 to tridiagonal form and form the
C           matrix W which is needed to update the unreduced part of
C           the matrix
C
            CALL DLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK,
     $                   LDWORK )
C
C           Update the unreduced submatrix A(1:i-1,1:i-1), using an
C           update of the form:  A := A - V*W**T - W*V**T
C
            CALL DSYR2K( UPLO, 'No transpose', I-1, NB, -ONE, A( 1, I ),
     $                   LDA, WORK, LDWORK, ONE, A, LDA )
C
C           Copy superdiagonal elements back into A, and diagonal
C           elements into D
C
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
C
C        Use unblocked code to reduce the last or only block
C
         CALL DSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
C
C        Reduce the lower triangle of A
C
         DO 40 I = 1, N - NX, NB
C
C           Reduce columns i:i+nb-1 to tridiagonal form and form the
C           matrix W which is needed to update the unreduced part of
C           the matrix
C
            CALL DLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ),
     $                   TAU( I ), WORK, LDWORK )
C
C           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
C           an update of the form:  A := A - V*W**T - W*V**T
C
            CALL DSYR2K( UPLO, 'No transpose', N-I-NB+1, NB, -ONE,
     $                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE,
     $                   A( I+NB, I+NB ), LDA )
C
C           Copy subdiagonal elements back into A, and diagonal
C           elements into D
C
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
C
C        Use unblocked code to reduce the last or only block
C
         CALL DSYTD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ),
     $                TAU( I ), IINFO )
      END IF
C
      WORK( 1 ) = LWKOPT
      RETURN
C
C     End of DSYTRD
C
      END

      SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0,
     $                   HALF = 1.0D0 / 2.0D0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      DOUBLE PRECISION   ALPHA, TAUI
C     ..
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DLARFG, DSYMV, DSYR2, XERBLA
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTD2', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.LE.0 )
     $   RETURN
C
      IF( UPPER ) THEN
C
C        Reduce the upper triangle of A
C
         DO 10 I = N - 1, 1, -1
C
C           Generate elementary reflector H(i) = I - tau * v * v**T
C           to annihilate A(1:i-1,i+1)
C
            CALL DLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI )
            E( I ) = A( I, I+1 )
C
            IF( TAUI.NE.ZERO ) THEN
C
C              Apply H(i) from both sides to A(1:i,1:i)
C
               A( I, I+1 ) = ONE
C
C              Compute  x := tau * A * v  storing x in TAU(1:i)
C
               CALL DSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO,
     $                     TAU, 1 )
C
C              Compute  w := x - 1/2 * tau * (x**T * v) * v
C
               ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL DAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
C
C              Apply the transformation as a rank-2 update:
C                 A := A - v * w**T - w * v**T
C
               CALL DSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A,
     $                     LDA )
C
               A( I, I+1 ) = E( I )
            END IF
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
C
C        Reduce the lower triangle of A
C
         DO 20 I = 1, N - 1
C
C           Generate elementary reflector H(i) = I - tau * v * v**T
C           to annihilate A(i+2:n,i)
C
            CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                   TAUI )
            E( I ) = A( I+1, I )
C
            IF( TAUI.NE.ZERO ) THEN
C
C              Apply H(i) from both sides to A(i+1:n,i+1:n)
C
               A( I+1, I ) = ONE
C
C              Compute  x := tau * A * v  storing y in TAU(i:n-1)
C
               CALL DSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, TAU( I ), 1 )
C
C              Compute  w := x - 1/2 * tau * (x**T * v) * v
C
               ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, A( I+1, I ),
     $                 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
C
C              Apply the transformation as a rank-2 update:
C                 A := A - v * w**T - w * v**T
C
               CALL DSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1,
     $                     A( I+1, I+1 ), LDA )
C
               A( I+1, I ) = E( I )
            END IF
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
C
      RETURN
C
C     End of DSYTD2
C
      END

      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
C     ..
C     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH, DISNAN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
C
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
C
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO .OR. DISNAN(CFROM) ) THEN
         INFO = -4
      ELSE IF( DISNAN(CTO) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
C
C     Get machine parameters
C
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
C
      CFROMC = CFROM
      CTOC = CTO
C
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      IF( CFROM1.EQ.CFROMC ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = .TRUE.
         CTO1 = CTOC
      ELSE
         CTO1 = CTOC / BIGNUM
         IF( CTO1.EQ.CTOC ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            MUL = CTOC
            DONE = .TRUE.
            CFROMC = ONE
         ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
            MUL = SMLNUM
            DONE = .FALSE.
            CFROMC = CFROM1
         ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
            MUL = BIGNUM
            DONE = .FALSE.
            CTOC = CTO1
         ELSE
            MUL = CTOC / CFROMC
            DONE = .TRUE.
         END IF
      END IF
C
      IF( ITYPE.EQ.0 ) THEN
C
C        Full matrix
C
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
C
      ELSE IF( ITYPE.EQ.1 ) THEN
C
C        Lower triangular matrix
C
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
C
      ELSE IF( ITYPE.EQ.2 ) THEN
C
C        Upper triangular matrix
C
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
C
      ELSE IF( ITYPE.EQ.3 ) THEN
C
C        Upper Hessenberg matrix
C
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
C
      ELSE IF( ITYPE.EQ.4 ) THEN
C
C        Lower half of a symmetric band matrix
C
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
C
      ELSE IF( ITYPE.EQ.5 ) THEN
C
C        Upper half of a symmetric band matrix
C
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
C
      ELSE IF( ITYPE.EQ.6 ) THEN
C
C        Band matrix
C
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
C
      END IF
C
      IF( .NOT.DONE )
     $   GO TO 10
C
      RETURN
C
C     End of DLASCL
C
      END

      DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
C     ..
C
C =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLASSQ
C     ..
C     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
C     ..
C     .. Executable Statements ..
C
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
C
C        Find max(abs(A(i,j))).
C
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J
                  SUM = ABS( A( I, J ) )
                  IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               DO 30 I = J, N
                  SUM = ABS( A( I, J ) )
                  IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
C
C        Find normI(A) ( = norm1(A), since A is symmetric).
C
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
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
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
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
C
C        Find normF(A).
C
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
C
      DLANSY = VALUE
      RETURN
C
C     End of DLANSY
C
      END

      DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, SCALE, SUM
C     ..
C     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLASSQ
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
C     ..
C     .. Executable Statements ..
C
      IF( N.LE.0 ) THEN
         ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
C
C        Find max(abs(A(i,j))).
C
         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            SUM = ABS( D( I ) )
            IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
            SUM = ABS( E( I ) )
            IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
   10    CONTINUE
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.
     $         LSAME( NORM, 'I' ) ) THEN
C
C        Find norm1(A).
C
         IF( N.EQ.1 ) THEN
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = ABS( D( 1 ) )+ABS( E( 1 ) )
            SUM = ABS( E( N-1 ) )+ABS( D( N ) )
            IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
            DO 20 I = 2, N - 1
               SUM = ABS( D( I ) )+ABS( E( I ) )+ABS( E( I-1 ) )
               IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
   20       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
C
C        Find normF(A).
C
         SCALE = ZERO
         SUM = ONE
         IF( N.GT.1 ) THEN
            CALL DLASSQ( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         CALL DLASSQ( N, D, 1, SCALE, SUM )
         ANORM = SCALE*SQRT( SUM )
      END IF
C
      DLANST = ANORM
      RETURN
C
C     End of DLANST
C
      END

      SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, IW
      DOUBLE PRECISION   ALPHA
C     ..
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DGEMV, DLARFG, DSCAL, DSYMV
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF( N.LE.0 )
     $   RETURN
C
      IF( LSAME( UPLO, 'U' ) ) THEN
C
C        Reduce last NB columns of upper triangle
C
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
C
C              Update A(1:i,i)
C
               CALL DGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ),
     $                     LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL DGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ),
     $                     LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
            END IF
            IF( I.GT.1 ) THEN
C
C              Generate elementary reflector H(i) to annihilate
C              A(1:i-2,i)
C
               CALL DLARFG( I-1, A( I-1, I ), A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = A( I-1, I )
               A( I-1, I ) = ONE
C
C              Compute W(1:i-1,i)
C
               CALL DSYMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1,
     $                     ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL DGEMV( 'Transpose', I-1, N-I, ONE, W( 1, IW+1 ),
     $                        LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
                  CALL DGEMV( 'Transpose', I-1, N-I, ONE, A( 1, I+1 ),
     $                        LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
               END IF
               CALL DSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*DDOT( I-1, W( 1, IW ), 1,
     $                 A( 1, I ), 1 )
               CALL DAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
C
   10    CONTINUE
      ELSE
C
C        Reduce first NB columns of lower triangle
C
         DO 20 I = 1, NB
C
C           Update A(i:n,i)
C
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ),
     $                  LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ),
     $                  LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            IF( I.LT.N ) THEN
C
C              Generate elementary reflector H(i) to annihilate
C              A(i+2:n,i)
C
               CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                      TAU( I ) )
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE
C
C              Compute W(i+1:n,i)
C
               CALL DSYMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ),
     $                     LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*DDOT( N-I, W( I+1, I ), 1,
     $                 A( I+1, I ), 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
C
   20    CONTINUE
      END IF
C
      RETURN
C
C     End of DLATRD
C
      END

      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEMV, DTRMV
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF( N.EQ.0 )
     $   RETURN
C
      IF( LSAME( DIRECT, 'F' ) ) THEN
         PREVLASTV = N
         DO I = 1, K
            PREVLASTV = MAX( I, PREVLASTV )
            IF( TAU( I ).EQ.ZERO ) THEN
C
C              H(i)  =  I
C
               DO J = 1, I
                  T( J, I ) = ZERO
               END DO
            ELSE
C
C              general case
C
               IF( LSAME( STOREV, 'C' ) ) THEN
C                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( LASTV, I ).NE.ZERO ) EXIT
                  END DO
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * V( I , J )
                  END DO
                  J = MIN( LASTV, PREVLASTV )
C
C                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)
C
                  CALL DGEMV( 'Transpose', J-I, I-1, -TAU( I ),
     $                        V( I+1, 1 ), LDV, V( I+1, I ), 1, ONE,
     $                        T( 1, I ), 1 )
               ELSE
C                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( I, LASTV ).NE.ZERO ) EXIT
                  END DO
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * V( J , I )
                  END DO
                  J = MIN( LASTV, PREVLASTV )
C
C                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T
C
                  CALL DGEMV( 'No transpose', I-1, J-I, -TAU( I ),
     $                        V( 1, I+1 ), LDV, V( I, I+1 ), LDV, ONE,
     $                        T( 1, I ), 1 )
               END IF
C
C              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
C
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
               IF( I.GT.1 ) THEN
                  PREVLASTV = MAX( PREVLASTV, LASTV )
               ELSE
                  PREVLASTV = LASTV
               END IF
            END IF
         END DO
      ELSE
         PREVLASTV = 1
         DO I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
C
C              H(i)  =  I
C
               DO J = I, K
                  T( J, I ) = ZERO
               END DO
            ELSE
C
C              general case
C
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
C                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( LASTV, I ).NE.ZERO ) EXIT
                     END DO
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * V( N-K+I , J )
                     END DO
                     J = MAX( LASTV, PREVLASTV )
C
C                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)
C
                     CALL DGEMV( 'Transpose', N-K+I-J, K-I, -TAU( I ),
     $                           V( J, I+1 ), LDV, V( J, I ), 1, ONE,
     $                           T( I+1, I ), 1 )
                  ELSE
C                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( I, LASTV ).NE.ZERO ) EXIT
                     END DO
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * V( J, N-K+I )
                     END DO
                     J = MAX( LASTV, PREVLASTV )
C
C                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T
C
                     CALL DGEMV( 'No transpose', K-I, N-K+I-J,
     $                    -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV,
     $                    ONE, T( I+1, I ), 1 )
                  END IF
C
C                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
C
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
                  IF( I.GT.1 ) THEN
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
C
C     End of DLARFT
C
      END

      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
C
C  -- LAPACK auxiliary routine (version 3.8.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     November 2017
C
C     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
C     ..
C     .. External Subroutines ..
      EXTERNAL           DSCAL
C     ..
C     .. Executable Statements ..
C
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
C
      XNORM = DNRM2( N-1, X, INCX )
C
      IF( XNORM.EQ.ZERO ) THEN
C
C        H  =  I
C
         TAU = ZERO
      ELSE
C
C        general case
C
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
C
C           XNORM, BETA may be inaccurate; scale X and recompute them
C
            RSAFMN = ONE / SAFMIN
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) )
     $         GO TO 10
C
C           New BETA is at most 1, at least SAFMIN
C
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         END IF
         TAU = ( BETA-ALPHA ) / BETA
         CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
C
C        If ALPHA is subnormal, it may lose relative accuracy
C
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
C
      RETURN
C
C     End of DLARFG
C
      END

      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
C     ..
C     .. External Functions ..
      LOGICAL            DISNAN
      EXTERNAL           DISNAN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
C
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            ABSXI = ABS( X( IX ) )
            IF( ABSXI.GT.ZERO.OR.DISNAN( ABSXI ) ) THEN
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
C
C     End of DLASSQ
C
      END

      SUBROUTINE DLASRT( ID, N, D, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
C     ..
C     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
C     ..
C     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASRT', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.LE.1 )
     $   RETURN
C
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
C
C        Do Insertion sort on D( START:ENDD )
C
         IF( DIR.EQ.0 ) THEN
C
C           Sort into decreasing order
C
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
C
         ELSE
C
C           Sort into increasing order
C
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
C
         END IF
C
      ELSE IF( ENDD-START.GT.SELECT ) THEN
C
C        Partition D( START:ENDD ) and stack parts, largest one first
C
C        Choose partition entry as median of 3
C
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
C
         IF( DIR.EQ.0 ) THEN
C
C           Sort into decreasing order
C
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX )
     $         GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX )
     $         GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
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
C
C           Sort into increasing order
C
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX )
     $         GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX )
     $         GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
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
      IF( STKPNT.GT.0 )
     $   GO TO 10
      RETURN
C
C     End of DLASRT
C
      END

      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
C     ..
C
C =====================================================================
C
C     .. Local Scalars ..
      INTEGER            I, J
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
C
      IF( LSAME( UPLO, 'U' ) ) THEN
C
C        Set the strictly upper triangular or trapezoidal part of the
C        array to ALPHA.
C
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
C
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
C
C        Set the strictly lower triangular or trapezoidal part of the
C        array to ALPHA.
C
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
C
      ELSE
C
C        Set the leading m-by-n submatrix to ALPHA.
C
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
C
C     Set the first min(M,N) diagonal elements to BETA.
C
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
C
      RETURN
C
C     End of DLASET
C
      END

      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2013
C
C     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
C     ..
C     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DTRMM
C     ..
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
C
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
C
      IF( LSAME( STOREV, 'C' ) ) THEN
C
         IF( LSAME( DIRECT, 'F' ) ) THEN
C
C           Let  V =  ( V1 )    (first K rows)
C                     ( V2 )
C           where  V1  is unit lower triangular.
C
            IF( LSAME( SIDE, 'L' ) ) THEN
C
C              Form  H * C  or  H**T * C  where  C = ( C1 )
C                                                    ( C2 )
C
C              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
C
C              W := C1**T
C
               DO 10 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
C
C              W := W * V1
C
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
C
C                 W := W + C2**T * V2
C
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
C
C              W := W * T**T  or  W * T
C
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - V * W**T
C
               IF( M.GT.K ) THEN
C
C                 C2 := C2 - V2 * W**T
C
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
C
C              W := W * V1**T
C
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
C
C              C1 := C1 - W**T
C
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
C
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
C
C              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
C
C              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
C
C              W := C1
C
               DO 40 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
C
C              W := W * V1
C
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
C
C                 W := W + C2 * V2
C
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
C
C              W := W * T  or  W * T**T
C
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - W * V**T
C
               IF( N.GT.K ) THEN
C
C                 C2 := C2 - W * V2**T
C
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
C
C              W := W * V1**T
C
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
C
C              C1 := C1 - W
C
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
C
         ELSE
C
C           Let  V =  ( V1 )
C                     ( V2 )    (last K rows)
C           where  V2  is unit upper triangular.
C
            IF( LSAME( SIDE, 'L' ) ) THEN
C
C              Form  H * C  or  H**T * C  where  C = ( C1 )
C                                                    ( C2 )
C
C              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
C
C              W := C2**T
C
               DO 70 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
C
C              W := W * V2
C
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
C
C                 W := W + C1**T * V1
C
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
C
C              W := W * T**T  or  W * T
C
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - V * W**T
C
               IF( M.GT.K ) THEN
C
C                 C1 := C1 - V1 * W**T
C
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
C
C              W := W * V2**T
C
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
C
C              C2 := C2 - W**T
C
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
C
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
C
C              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
C
C              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
C
C              W := C2
C
               DO 100 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
C
C              W := W * V2
C
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
C
C                 W := W + C1 * V1
C
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
C
C              W := W * T  or  W * T**T
C
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - W * V**T
C
               IF( N.GT.K ) THEN
C
C                 C1 := C1 - W * V1**T
C
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
C
C              W := W * V2**T
C
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
C
C              C2 := C2 - W
C
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
C
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
C
         IF( LSAME( DIRECT, 'F' ) ) THEN
C
C           Let  V =  ( V1  V2 )    (V1: first K columns)
C           where  V1  is unit upper triangular.
C
            IF( LSAME( SIDE, 'L' ) ) THEN
C
C              Form  H * C  or  H**T * C  where  C = ( C1 )
C                                                    ( C2 )
C
C              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
C
C              W := C1**T
C
               DO 130 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
C
C              W := W * V1**T
C
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
C
C                 W := W + C2**T * V2**T
C
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
     $                        WORK, LDWORK )
               END IF
C
C              W := W * T**T  or  W * T
C
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - V**T * W**T
C
               IF( M.GT.K ) THEN
C
C                 C2 := C2 - V2**T * W**T
C
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
C
C              W := W * V1
C
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
C
C              C1 := C1 - W**T
C
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
C
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
C
C              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
C
C              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
C
C              W := C1
C
               DO 160 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
C
C              W := W * V1**T
C
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
C
C                 W := W + C2 * V2**T
C
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
C
C              W := W * T  or  W * T**T
C
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - W * V
C
               IF( N.GT.K ) THEN
C
C                 C2 := C2 - W * V2
C
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
C
C              W := W * V1
C
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
C
C              C1 := C1 - W
C
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
C
            END IF
C
         ELSE
C
C           Let  V =  ( V1  V2 )    (V2: last K columns)
C           where  V2  is unit lower triangular.
C
            IF( LSAME( SIDE, 'L' ) ) THEN
C
C              Form  H * C  or  H**T * C  where  C = ( C1 )
C                                                    ( C2 )
C
C              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
C
C              W := C2**T
C
               DO 190 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
C
C              W := W * V2**T
C
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
C
C                 W := W + C1**T * V1**T
C
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
C
C              W := W * T**T  or  W * T
C
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - V**T * W**T
C
               IF( M.GT.K ) THEN
C
C                 C1 := C1 - V1**T * W**T
C
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
C
C              W := W * V2
C
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
C
C              C2 := C2 - W**T
C
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
C
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
C
C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
C
C              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
C
C              W := C2
C
               DO 220 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
C
C              W := W * V2**T
C
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
C
C                 W := W + C1 * V1**T
C
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
C
C              W := W * T  or  W * T**T
C
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - W * V
C
               IF( N.GT.K ) THEN
C
C                 C1 := C1 - W * V1
C
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
C
C              W := W * V2
C
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
C
C              C1 := C1 - W
C
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
C
            END IF
C
         END IF
      END IF
C
      RETURN
C
C     End of DLARFB
C
      END

      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILADLR, ILADLC
      EXTERNAL           LSAME, ILADLR, ILADLC
C     ..
C     .. Executable Statements ..
C
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
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
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
      IF( APPLYLEFT ) THEN
C
C        Form  H * C
C
         IF( LASTV.GT.0 ) THEN
C
C           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
C
            CALL DGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV,
     $           ZERO, WORK, 1 )
C
C           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
C
            CALL DGER( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
C
C        Form  C * H
C
         IF( LASTV.GT.0 ) THEN
C
C           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
C
            CALL DGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC,
     $           V, INCV, ZERO, WORK, 1 )
C
C           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
C
            CALL DGER( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
C
C     End of DLARF
C
      END

      SUBROUTINE DLARTG( F, G, CS, SN, R )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
C     ..
C     .. Local Scalars ..
C     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
C     ..
C     .. Save statement ..
C     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
C     ..
C     .. Data statements ..
C     DATA               FIRST / .TRUE. /
C     ..
C     .. Executable Statements ..
C
C     IF( FIRST ) THEN
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
C        FIRST = .FALSE.
C     END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 30
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
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
C
C     End of DLARTG
C
      END

      LOGICAL FUNCTION DISNAN( DIN )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2017
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN
C     ..
C
C  =====================================================================
C
C  .. External Functions ..
      LOGICAL DLAISNAN
      EXTERNAL DLAISNAN
C  ..
C  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END

      INTEGER FUNCTION ILADLC( M, N, A, LDA )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            M, N, LDA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Executable Statements ..
C
C     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 ) THEN
         ILADLC = N
      ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILADLC = N
      ELSE
C     Now scan each column from the end, returning with the first non-zero.
         DO ILADLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILADLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END

      INTEGER FUNCTION ILADLR( M, N, A, LDA )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            M, N, LDA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER I, J
C     ..
C     .. Executable Statements ..
C
C     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILADLR = M
      ELSE IF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILADLR = M
      ELSE
C     Scan up each column tracking the last zero row seen.
         ILADLR = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ILADLR = MAX( ILADLR, I )
         END DO
      END IF
      RETURN
      END

      LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2017
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
C     ..
C
C  =====================================================================
C
C  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END

      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
C
C  -- LAPACK auxiliary routine (version 3.8.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     November 2017
C
C     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
C     ..
C
C  =====================================================================
C
C     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME, TWOSTAGE
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*16
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
C     ..
C     .. External Functions ..
      INTEGER            IEEECK, IPARMQ, IPARAM2STAGE
      EXTERNAL           IEEECK, IPARMQ, IPARAM2STAGE
C     ..
C     .. Executable Statements ..
C
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160)ISPEC
C
C     Invalid value for ISPEC
C
      ILAENV = -1
      RETURN
C
   10 CONTINUE
C
C     Convert NAME to upper case if the first character is lower case.
C
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
C
C        ASCII character set
C
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
C
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
C
C        EBCDIC character set
C
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $             I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
C
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
C
C        Prime machines:  ASCII+128
C
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
C
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
      TWOSTAGE = LEN( SUBNAM ).GE.11
     $           .AND. SUBNAM( 11: 11 ).EQ.'2'
C
      GO TO ( 50, 60, 70 )ISPEC
C
   50 CONTINUE
C
C     ISPEC = 1:  block size
C
C     In these examples, separate code is provided for setting NB for
C     real and complex.  We assume that NB will take the same value in
C     single or double precision.
C
      NB = 1
C
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'QR ') THEN
            IF( N3 .EQ. 1) THEN
               IF( SNAME ) THEN
C     M*N
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               ELSE
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
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
         ELSE IF( C3.EQ.'LQ ') THEN
            IF( N3 .EQ. 2) THEN
               IF( SNAME ) THEN
C     M*N
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               ELSE
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
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
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
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
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( TWOSTAGE ) THEN
               NB = 192
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF ( C3.EQ.'EVC' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      ELSE IF( C2.EQ.'GG' ) THEN
         NB = 32
         IF( C3.EQ.'HD3' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         END IF
      END IF
      ILAENV = NB
      RETURN
C
   60 CONTINUE
C
C     ISPEC = 2:  minimum block size
C
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'GG' ) THEN
         NBMIN = 2
         IF( C3.EQ.'HD3' ) THEN
            NBMIN = 2
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
C
   70 CONTINUE
C
C     ISPEC = 3:  crossover point
C
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'GG' ) THEN
         NX = 128
         IF( C3.EQ.'HD3' ) THEN
            NX = 128
         END IF
      END IF
      ILAENV = NX
      RETURN
C
   80 CONTINUE
C
C     ISPEC = 4:  number of shifts (used by xHSEQR)
C
      ILAENV = 6
      RETURN
C
   90 CONTINUE
C
C     ISPEC = 5:  minimum column dimension (not used)
C
      ILAENV = 2
      RETURN
C
  100 CONTINUE
C
C     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
C
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
C
  110 CONTINUE
C
C     ISPEC = 7:  number of processors (not used)
C
      ILAENV = 1
      RETURN
C
  120 CONTINUE
C
C     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
C
      ILAENV = 50
      RETURN
C
  130 CONTINUE
C
C     ISPEC = 9:  maximum size of the subproblems at the bottom of the
C                 computation tree in the divide-and-conquer algorithm
C                 (used by xGELSD and xGESDD)
C
      ILAENV = 25
      RETURN
C
  140 CONTINUE
C
C     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
C
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
C
  150 CONTINUE
C
C     ISPEC = 11: infinity arithmetic can be trusted not to trap
C
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
C
  160 CONTINUE
C
C     12 <= ISPEC <= 16: xHSEQR or related subroutines.
C
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
C
C     End of ILAENV
C
      END

      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2017
C
C     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
C
C  ================================================================
C     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,
     $                   ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14,
     $                   NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
C     ..
C     .. Local Scalars ..
      INTEGER            NH, NS
      INTEGER            I, IC, IZ
      CHARACTER          SUBNAM*6
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
C     ..
C     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR.
     $    ( ISPEC.EQ.IACC22 ) ) THEN
C
C        ==== Set the number simultaneous shifts ====
C
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 )
     $      NS = 4
         IF( NH.GE.60 )
     $      NS = 10
         IF( NH.GE.150 )
     $      NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 )
     $      NS = 64
         IF( NH.GE.3000 )
     $      NS = 128
         IF( NH.GE.6000 )
     $      NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
C
      IF( ISPEC.EQ.INMIN ) THEN
C
C
C        ===== Matrices of order smaller than NMIN get sent
C        .     to xLAHQR, the classic double shift algorithm.
C        .     This must be at least 11. ====
C
         IPARMQ = NMIN
C
      ELSE IF( ISPEC.EQ.INIBL ) THEN
C
C        ==== INIBL: skip a multi-shift qr iteration and
C        .    whenever aggressive early deflation finds
C        .    at least (NIBBLE*(window size)/100) deflations. ====
C
         IPARMQ = NIBBLE
C
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
C
C        ==== NSHFTS: The number of simultaneous shifts =====
C
         IPARMQ = NS
C
      ELSE IF( ISPEC.EQ.INWIN ) THEN
C
C        ==== NW: deflation window size.  ====
C
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
C
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
C
C        ==== IACC22: Whether to accumulate reflections
C        .     before updating the far-from-diagonal elements
C        .     and whether to use 2-by-2 block structure while
C        .     doing it.  A small amount of work could be saved
C        .     by making this choice dependent also upon the
C        .     NH=IHI-ILO+1.
C
C
C        Convert NAME to upper case if the first character is lower case.
C
         IPARMQ = 0
         SUBNAM = NAME
         IC = ICHAR( SUBNAM( 1: 1 ) )
         IZ = ICHAR( 'Z' )
         IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
C
C           ASCII character set
C
            IF( IC.GE.97 .AND. IC.LE.122 ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( IC.GE.97 .AND. IC.LE.122 )
     $               SUBNAM( I: I ) = CHAR( IC-32 )
               END DO
            END IF
C
         ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
C
C           EBCDIC character set
C
            IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $          ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $          ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC+64 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $                ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $                ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $                I ) = CHAR( IC+64 )
               END DO
            END IF
C
         ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
C
C           Prime machines:  ASCII+128
C
            IF( IC.GE.225 .AND. IC.LE.250 ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( IC.GE.225 .AND. IC.LE.250 )
     $               SUBNAM( I: I ) = CHAR( IC-32 )
               END DO
            END IF
         END IF
C
         IF( SUBNAM( 2:6 ).EQ.'GGHRD' .OR.
     $       SUBNAM( 2:6 ).EQ.'GGHD3' ) THEN
            IPARMQ = 1
            IF( NH.GE.K22MIN )
     $         IPARMQ = 2
         ELSE IF ( SUBNAM( 4:6 ).EQ.'EXC' ) THEN
            IF( NH.GE.KACMIN )
     $         IPARMQ = 1
            IF( NH.GE.K22MIN )
     $         IPARMQ = 2
         ELSE IF ( SUBNAM( 2:6 ).EQ.'HSEQR' .OR.
     $             SUBNAM( 2:5 ).EQ.'LAQR' ) THEN
            IF( NS.GE.KACMIN )
     $         IPARMQ = 1
            IF( NS.GE.K22MIN )
     $         IPARMQ = 2
         END IF
C
      ELSE
C        ===== invalid value of ispec =====
         IPARMQ = -1
C
      END IF
C
C     ==== End of IPARMQ ====
C
      END

      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
C     ..
C
C  =====================================================================
C
C     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
C     ..
C     .. Executable Statements ..
      IEEECK = 1
C
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
C
C
C
C
C     Return if we were only asked to check infinity arithmetic
C
      IF( ISPEC.EQ.0 )
     $   RETURN
C
      NAN1 = POSINF + NEGINF
C
      NAN2 = POSINF / NEGINF
C
      NAN3 = POSINF / POSINF
C
      NAN4 = POSINF*ZERO
C
      NAN5 = NEGINF*NEGZRO
C
      NAN6 = NAN5*ZERO
C
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
C
      RETURN
      END

      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
C
C        Form  P * A
C
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
C
C        Form A * P**T
C
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
C
      RETURN
C
C     End of DLASR
C
      END

      SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, J, LWKOPT, NB
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
C     ..
C     .. External Subroutines ..
      EXTERNAL           DORGQL, DORGQR, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N-1 ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            NB = ILAENV( 1, 'DORGQL', ' ', N-1, N-1, N-1, -1 )
         ELSE
            NB = ILAENV( 1, 'DORGQR', ' ', N-1, N-1, N-1, -1 )
         END IF
         LWKOPT = MAX( 1, N-1 )*NB
         WORK( 1 ) = LWKOPT
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
C
      IF( UPPER ) THEN
C
C        Q was determined by a call to DSYTRD with UPLO = 'U'
C
C        Shift the vectors which define the elementary reflectors one
C        column to the left, and set the last row and column of Q to
C        those of the unit matrix
C
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
C
C        Generate Q(1:n-1,1:n-1)
C
         CALL DORGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
C
      ELSE
C
C        Q was determined by a call to DSYTRD with UPLO = 'L'.
C
C        Shift the vectors which define the elementary reflectors one
C        column to the right, and set the first row and column of Q to
C        those of the unit matrix
C
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
         IF( N.GT.1 ) THEN
C
C           Generate Q(2:n,2:n)
C
            CALL DORGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,
     $                   LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
C
C     End of DORGTR
C
      END

      SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,
     $                   NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASET, DLASR,
     $                   DLASRT, DSWAP, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
C
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEQR', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 )
     $   RETURN
C
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
C
C     Determine the unit roundoff and over/underflow thresholds.
C
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
C
C     Compute the eigenvalues and eigenvectors of the tridiagonal
C     matrix.
C
      IF( ICOMPZ.EQ.2 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
C
      NMAXIT = N*MAXIT
      JTOT = 0
C
C     Determine where the matrix splits and choose QL or QR iteration
C     for each block, according to whether top or bottom diagonal
C     element is smaller.
C
      L1 = 1
      NM1 = N - 1
C
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
C
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
C
C     Scale submatrix in rows and columns L to LEND
C
      ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
C
C     Choose between QL and QR iteration
C
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
C
      IF( LEND.GT.L ) THEN
C
C        QL Iteration
C
C        Look for small subdiagonal element.
C
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GO TO 60
   50       CONTINUE
         END IF
C
         M = LEND
C
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
C
C        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
C        to compute its eigensystem.
C
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL DLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF
C
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
C
C        Form shift.
C
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
C
         S = ONE
         C = ONE
         P = ZERO
C
C        Inner loop
C
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
C
C           If eigenvectors are desired, then save rotations.
C
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
C
   70    CONTINUE
C
C        If eigenvectors are desired, then apply saved rotations.
C
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL DLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
C
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
C
C        Eigenvalue found.
C
   80    CONTINUE
         D( L ) = P
C
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140
C
      ELSE
C
C        QR Iteration
C
C        Look for small superdiagonal element.
C
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GO TO 110
  100       CONTINUE
         END IF
C
         M = LEND
C
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
C
C        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
C        to compute its eigensystem.
C
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL DLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF
C
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
C
C        Form shift.
C
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
C
         S = ONE
         C = ONE
         P = ZERO
C
C        Inner loop
C
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
C
C           If eigenvectors are desired, then save rotations.
C
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
C
  120    CONTINUE
C
C        If eigenvectors are desired, then apply saved rotations.
C
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL DLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
C
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
C
C        Eigenvalue found.
C
  130    CONTINUE
         D( L ) = P
C
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140
C
      END IF
C
C     Undo scaling if necessary
C
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      END IF
C
C     Check for no convergence to an eigenvalue after a total
C     of N*MAXIT iterations.
C
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 150 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  150 CONTINUE
      GO TO 190
C
C     Order eigenvalues and eigenvectors.
C
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
C
C        Use Quick Sort
C
         CALL DLASRT( 'I', N, D, INFO )
C
      ELSE
C
C        Use Selection Sort to minimize swaps of eigenvectors
C
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
C
  190 CONTINUE
      RETURN
C
C     End of DSTEQR
C
      END

      SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT,
     $                   NB, NBMIN, NX
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2L, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            NB = ILAENV( 1, 'DORGQL', ' ', M, N, K, -1 )
            LWKOPT = N*NB
         END IF
         WORK( 1 ) = LWKOPT
C
         IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
         END IF
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         RETURN
      END IF
C
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         NX = MAX( 0, ILAENV( 3, 'DORGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
C
C              Not enough workspace to use optimal NB:  reduce NB and
C              determine the minimum value of NB.
C
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
C
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
C
C        Use blocked code after the first block.
C        The last kk columns are handled by the block method.
C
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
C
C        Set A(m-kk+1:m,1:n-kk) to zero.
C
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
C
C     Use unblocked code for the first or only block.
C
      CALL DORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
C
      IF( KK.GT.0 ) THEN
C
C        Use blocked code
C
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
C
C              Form the triangular factor of the block reflector
C              H = H(i+ib-1) . . . H(i+1) H(i)
C
               CALL DLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB,
     $                      A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
C
C              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
C
               CALL DLARFB( 'Left', 'No transpose', 'Backward',
     $                      'Columnwise', M-K+I+IB-1, N-K+I-1, IB,
     $                      A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA,
     $                      WORK( IB+1 ), LDWORK )
            END IF
C
C           Apply H to rows 1:m-k+i+ib-1 of current block
C
            CALL DORG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA,
     $                   TAU( I ), WORK, IINFO )
C
C           Set rows m-k+i+ib:m of current block to zero
C
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
C
      WORK( 1 ) = IWS
      RETURN
C
C     End of DORGQL
C
      END

      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
C
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
C
C              Not enough workspace to use optimal NB:  reduce NB and
C              determine the minimum value of NB.
C
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
C
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
C
C        Use blocked code after the last block.
C        The first kk columns are handled by the block method.
C
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
C
C        Set A(1:kk,kk+1:n) to zero.
C
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
C
C     Use unblocked code for the last or only block.
C
      IF( KK.LT.N )
     $   CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
C
      IF( KK.GT.0 ) THEN
C
C        Use blocked code
C
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
C
C              Form the triangular factor of the block reflector
C              H = H(i) H(i+1) . . . H(i+ib-1)
C
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
C
C              Apply H to A(i:m,i+ib:n) from the left
C
               CALL DLARFB( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
C
C           Apply H to rows i:m of current block
C
            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
C
C           Set rows 1:i-1 of current block to zero
C
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
C
      WORK( 1 ) = IWS
      RETURN
C
C     End of DORGQR
C
      END

      SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, II, J, L
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2L', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.LE.0 )
     $   RETURN
C
C     Initialise columns 1:n-k to columns of the unit matrix
C
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
C
      DO 40 I = 1, K
         II = N - K + I
C
C        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
C
         A( M-N+II, II ) = ONE
         CALL DLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A,
     $               LDA, WORK )
         CALL DSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
C
C        Set A(m-k+i+1:m,n-k+i) to zero
C
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
C
C     End of DORG2L
C
      END

      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, J, L
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2R', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.LE.0 )
     $   RETURN
C
C     Initialise columns k+1:n to columns of the unit matrix
C
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
C
      DO 40 I = K, 1, -1
C
C        Apply H(i) to A(i:m,i:n) from the left
C
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
C
C        Set A(1:i-1,i) to zero
C
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
C
C     End of DORG2R
C
      END

      SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, RT1, RT2
C     ..
C
C =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
C     ..
C     .. Executable Statements ..
C
C     Compute the eigenvalues
C
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
C
C        Includes case AB=ADF=0
C
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
C
C        Order of execution important.
C        To get fully accurate smaller eigenvalue,
C        next line needs to be executed in higher precision.
C
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
C
C        Order of execution important.
C        To get fully accurate smaller eigenvalue,
C        next line needs to be executed in higher precision.
C
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
C
C        Includes case RT1 = RT2 = 0
C
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
C
C     End of DLAE2
C
      END

      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
C     ..
C
C =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
C     ..
C     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
C     ..
C     .. Executable Statements ..
C
C     Compute the eigenvalues
C
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
C
C        Includes case AB=ADF=0
C
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
C
C        Order of execution important.
C        To get fully accurate smaller eigenvalue,
C        next line needs to be executed in higher precision.
C
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
C
C        Order of execution important.
C        To get fully accurate smaller eigenvalue,
C        next line needs to be executed in higher precision.
C
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
C
C        Includes case RT1 = RT2 = 0
C
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
C
C     Compute the eigenvector
C
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
C
C     End of DLAEV2
C
      END

      SUBROUTINE DSTERF( N, D, E, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INFO, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M,
     $                   NMAXIT
      DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC,
     $                   OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN,
     $                   SIGMA, SSFMAX, SSFMIN, RMAX
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           DLAMCH, DLANST, DLAPY2
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
C
C     Quick return if possible
C
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DSTERF', -INFO )
         RETURN
      END IF
      IF( N.LE.1 )
     $   RETURN
C
C     Determine the unit roundoff for this environment.
C
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
      RMAX = DLAMCH( 'O' )
C
C     Compute the eigenvalues of the tridiagonal matrix.
C
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
C
C     Determine where the matrix splits and choose QL or QR iteration
C     for each block, according to whether top or bottom diagonal
C     element is smaller.
C
      L1 = 1
C
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 170
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      DO 20 M = L1, N - 1
         IF( ABS( E( M ) ).LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $       1 ) ) ) )*EPS ) THEN
            E( M ) = ZERO
            GO TO 30
         END IF
   20 CONTINUE
      M = N
C
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
C
C     Scale submatrix in rows and columns L to LEND
C
      ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( (ANORM.GT.SSFMAX) ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
C
      DO 40 I = L, LEND - 1
         E( I ) = E( I )**2
   40 CONTINUE
C
C     Choose between QL and QR iteration
C
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
C
      IF( LEND.GE.L ) THEN
C
C        QL Iteration
C
C        Look for small subdiagonal element.
C
   50    CONTINUE
         IF( L.NE.LEND ) THEN
            DO 60 M = L, LEND - 1
               IF( ABS( E( M ) ).LE.EPS2*ABS( D( M )*D( M+1 ) ) )
     $            GO TO 70
   60       CONTINUE
         END IF
         M = LEND
C
   70    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 90
C
C        If remaining matrix is 2 by 2, use DLAE2 to compute its
C        eigenvalues.
C
         IF( M.EQ.L+1 ) THEN
            RTE = SQRT( E( L ) )
            CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 50
            GO TO 150
         END IF
C
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
C
C        Form shift.
C
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
C
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
C
C        Inner loop
C
         DO 80 I = M - 1, L, -1
            BB = E( I )
            R = P + BB
            IF( I.NE.M-1 )
     $         E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
C
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
C
C        Eigenvalue found.
C
   90    CONTINUE
         D( L ) = P
C
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 50
         GO TO 150
C
      ELSE
C
C        QR Iteration
C
C        Look for small superdiagonal element.
C
  100    CONTINUE
         DO 110 M = L, LEND + 1, -1
            IF( ABS( E( M-1 ) ).LE.EPS2*ABS( D( M )*D( M-1 ) ) )
     $         GO TO 120
  110    CONTINUE
         M = LEND
C
  120    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 140
C
C        If remaining matrix is 2 by 2, use DLAE2 to compute its
C        eigenvalues.
C
         IF( M.EQ.L-1 ) THEN
            RTE = SQRT( E( L-1 ) )
            CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 100
            GO TO 150
         END IF
C
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
C
C        Form shift.
C
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
C
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
C
C        Inner loop
C
         DO 130 I = M, L - 1
            BB = E( I )
            R = P + BB
            IF( I.NE.M )
     $         E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
C
         E( L-1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
C
C        Eigenvalue found.
C
  140    CONTINUE
         D( L ) = P
C
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 100
         GO TO 150
C
      END IF
C
C     Undo scaling if necessary
C
  150 CONTINUE
      IF( ISCALE.EQ.1 )
     $   CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
      IF( ISCALE.EQ.2 )
     $   CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
C
C     Check for no convergence to an eigenvalue after a total
C     of N*MAXIT iterations.
C
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 160 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  160 CONTINUE
      GO TO 180
C
C     Sort eigenvalues in increasing order.
C
  170 CONTINUE
      CALL DLASRT( 'I', N, D, INFO )
C
  180 CONTINUE
      RETURN
C
C     End of DSTERF
C
      END

      SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
     $                  INFO )
C
C  -- LAPACK driver routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )
C     ..
C     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE,
     $                   LLWORK, LWKOPT, NB
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANHE
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANHE
C     ..
C     .. External Subroutines ..
      EXTERNAL           DSCAL, DSTERF, XERBLA, ZHETRD, ZLASCL, ZSTEQR,
     $                   ZUNGTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
C
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
C
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, ( NB+1 )*N )
         WORK( 1 ) = LWKOPT
C
         IF( LWORK.LT.MAX( 1, 2*N-1 ) .AND. .NOT.LQUERY )
     $      INFO = -8
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHEEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 ) THEN
         RETURN
      END IF
C
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         WORK( 1 ) = 1
         IF( WANTZ )
     $      A( 1, 1 ) = CONE
         RETURN
      END IF
C
C     Get machine constants.
C
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
C
C     Scale matrix to allowable range, if necessary.
C
      ANRM = ZLANHE( 'M', UPLO, N, A, LDA, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 )
     $   CALL ZLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
C
C     Call ZHETRD to reduce Hermitian matrix to tridiagonal form.
C
      INDE = 1
      INDTAU = 1
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      CALL ZHETRD( UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ),
     $             WORK( INDWRK ), LLWORK, IINFO )
C
C     For eigenvalues only, call DSTERF.  For eigenvectors, first call
C     ZUNGTR to generate the unitary matrix, then call ZSTEQR.
C
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL ZUNGTR( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ),
     $                LLWORK, IINFO )
         INDWRK = INDE + N
         CALL ZSTEQR( JOBZ, N, W, RWORK( INDE ), A, LDA,
     $                RWORK( INDWRK ), INFO )
      END IF
C
C     If matrix was scaled, then rescale eigenvalues appropriately.
C
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
C
C     Set WORK(1) to optimal complex workspace size.
C
      WORK( 1 ) = LWKOPT
C
      RETURN
C
C     End of ZHEEV
C
      END

      SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZHER2K, ZHETD2, ZLATRD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
C
      IF( INFO.EQ.0 ) THEN
C
C        Determine the block size.
C
         NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHETRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
C
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
C
C        Determine when to cross over from blocked to unblocked code
C        (last block is always handled by unblocked code).
C
         NX = MAX( NB, ILAENV( 3, 'ZHETRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
C
C              Not enough workspace to use optimal NB:  determine the
C              minimum value of NB, and reduce NB or force use of
C              unblocked code by setting NX = N.
C
               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = ILAENV( 2, 'ZHETRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN )
     $            NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
C
      IF( UPPER ) THEN
C
C        Reduce the upper triangle of A.
C        Columns 1:kk are handled by the unblocked method.
C
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
C
C           Reduce columns i:i+nb-1 to tridiagonal form and form the
C           matrix W which is needed to update the unreduced part of
C           the matrix
C
            CALL ZLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK,
     $                   LDWORK )
C
C           Update the unreduced submatrix A(1:i-1,1:i-1), using an
C           update of the form:  A := A - V*W**H - W*V**H
C
            CALL ZHER2K( UPLO, 'No transpose', I-1, NB, -CONE,
     $                   A( 1, I ), LDA, WORK, LDWORK, ONE, A, LDA )
C
C           Copy superdiagonal elements back into A, and diagonal
C           elements into D
C
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
C
C        Use unblocked code to reduce the last or only block
C
         CALL ZHETD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
C
C        Reduce the lower triangle of A
C
         DO 40 I = 1, N - NX, NB
C
C           Reduce columns i:i+nb-1 to tridiagonal form and form the
C           matrix W which is needed to update the unreduced part of
C           the matrix
C
            CALL ZLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ),
     $                   TAU( I ), WORK, LDWORK )
C
C           Update the unreduced submatrix A(i+nb:n,i+nb:n), using
C           an update of the form:  A := A - V*W**H - W*V**H
C
            CALL ZHER2K( UPLO, 'No transpose', N-I-NB+1, NB, -CONE,
     $                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE,
     $                   A( I+NB, I+NB ), LDA )
C
C           Copy subdiagonal elements back into A, and diagonal
C           elements into D
C
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
C
C        Use unblocked code to reduce the last or only block
C
         CALL ZHETD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ),
     $                TAU( I ), IINFO )
      END IF
C
      WORK( 1 ) = LWKOPT
      RETURN
C
C     End of ZHETRD
C
      END

      SUBROUTINE ZUNGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
C     ..
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, J, LWKOPT, NB
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNGQL, ZUNGQR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N-1 ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            NB = ILAENV( 1, 'ZUNGQL', ' ', N-1, N-1, N-1, -1 )
         ELSE
            NB = ILAENV( 1, 'ZUNGQR', ' ', N-1, N-1, N-1, -1 )
         END IF
         LWKOPT = MAX( 1, N-1 )*NB
         WORK( 1 ) = LWKOPT
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
C
      IF( UPPER ) THEN
C
C        Q was determined by a call to ZHETRD with UPLO = 'U'
C
C        Shift the vectors which define the elementary reflectors one
C        column to the left, and set the last row and column of Q to
C        those of the unit matrix
C
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
C
C        Generate Q(1:n-1,1:n-1)
C
         CALL ZUNGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
C
      ELSE
C
C        Q was determined by a call to ZHETRD with UPLO = 'L'.
C
C        Shift the vectors which define the elementary reflectors one
C        column to the right, and set the first row and column of Q to
C        those of the unit matrix
C
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
         IF( N.GT.1 ) THEN
C
C           Generate Q(2:n,2:n)
C
            CALL ZUNGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,
     $                   LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
C
C     End of ZUNGTR
C
      END

      SUBROUTINE ZHETD2( UPLO, N, A, LDA, D, E, TAU, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( LDA, * ), TAU( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ONE, ZERO, HALF
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   HALF = ( 0.5D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      COMPLEX*16         ALPHA, TAUI
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZHEMV, ZHER2, ZLARFG
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      UPPER = LSAME( UPLO, 'U')
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHETD2', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.LE.0 )
     $   RETURN
C
      IF( UPPER ) THEN
C
C        Reduce the upper triangle of A
C
         A( N, N ) = DBLE( A( N, N ) )
         DO 10 I = N - 1, 1, -1
C
C           Generate elementary reflector H(i) = I - tau * v * v**H
C           to annihilate A(1:i-1,i+1)
C
            ALPHA = A( I, I+1 )
            CALL ZLARFG( I, ALPHA, A( 1, I+1 ), 1, TAUI )
            E( I ) = ALPHA
C
            IF( TAUI.NE.ZERO ) THEN
C
C              Apply H(i) from both sides to A(1:i,1:i)
C
               A( I, I+1 ) = ONE
C
C              Compute  x := tau * A * v  storing x in TAU(1:i)
C
               CALL ZHEMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO,
     $                     TAU, 1 )
C
C              Compute  w := x - 1/2 * tau * (x**H * v) * v
C
               ALPHA = -HALF*TAUI*ZDOTC( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL ZAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
C
C              Apply the transformation as a rank-2 update:
C                 A := A - v * w**H - w * v**H
C
               CALL ZHER2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A,
     $                     LDA )
C
            ELSE
               A( I, I ) = DBLE( A( I, I ) )
            END IF
            A( I, I+1 ) = E( I )
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
C
C        Reduce the lower triangle of A
C
         A( 1, 1 ) = DBLE( A( 1, 1 ) )
         DO 20 I = 1, N - 1
C
C           Generate elementary reflector H(i) = I - tau * v * v**H
C           to annihilate A(i+2:n,i)
C
            ALPHA = A( I+1, I )
            CALL ZLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAUI )
            E( I ) = ALPHA
C
            IF( TAUI.NE.ZERO ) THEN
C
C              Apply H(i) from both sides to A(i+1:n,i+1:n)
C
               A( I+1, I ) = ONE
C
C              Compute  x := tau * A * v  storing y in TAU(i:n-1)
C
               CALL ZHEMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, TAU( I ), 1 )
C
C              Compute  w := x - 1/2 * tau * (x**H * v) * v
C
               ALPHA = -HALF*TAUI*ZDOTC( N-I, TAU( I ), 1, A( I+1, I ),
     $                 1 )
               CALL ZAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
C
C              Apply the transformation as a rank-2 update:
C                 A := A - v * w**H - w * v**H
C
               CALL ZHER2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1,
     $                     A( I+1, I+1 ), LDA )
C
            ELSE
               A( I+1, I+1 ) = DBLE( A( I+1, I+1 ) )
            END IF
            A( I+1, I ) = E( I )
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
C
      RETURN
C
C     End of ZHETD2
C
      END

      SUBROUTINE ZLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   E( * )
      COMPLEX*16         A( LDA, * ), TAU( * ), W( LDW, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ZERO, ONE, HALF
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ),
     $                   HALF = ( 0.5D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      INTEGER            I, IW
      COMPLEX*16         ALPHA
C     ..
C     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZGEMV, ZHEMV, ZLACGV, ZLARFG, ZSCAL
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN
C     ..
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF( N.LE.0 )
     $   RETURN
C
      IF( LSAME( UPLO, 'U' ) ) THEN
C
C        Reduce last NB columns of upper triangle
C
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
C
C              Update A(1:i,i)
C
               A( I, I ) = DBLE( A( I, I ) )
               CALL ZLACGV( N-I, W( I, IW+1 ), LDW )
               CALL ZGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ),
     $                     LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL ZLACGV( N-I, W( I, IW+1 ), LDW )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               CALL ZGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ),
     $                     LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               A( I, I ) = DBLE( A( I, I ) )
            END IF
            IF( I.GT.1 ) THEN
C
C              Generate elementary reflector H(i) to annihilate
C              A(1:i-2,i)
C
               ALPHA = A( I-1, I )
               CALL ZLARFG( I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = ALPHA
               A( I-1, I ) = ONE
C
C              Compute W(1:i-1,i)
C
               CALL ZHEMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1,
     $                     ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL ZGEMV( 'Conjugate transpose', I-1, N-I, ONE,
     $                        W( 1, IW+1 ), LDW, A( 1, I ), 1, ZERO,
     $                        W( I+1, IW ), 1 )
                  CALL ZGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
                  CALL ZGEMV( 'Conjugate transpose', I-1, N-I, ONE,
     $                        A( 1, I+1 ), LDA, A( 1, I ), 1, ZERO,
     $                        W( I+1, IW ), 1 )
                  CALL ZGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
               END IF
               CALL ZSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*ZDOTC( I-1, W( 1, IW ), 1,
     $                 A( 1, I ), 1 )
               CALL ZAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
C
   10    CONTINUE
      ELSE
C
C        Reduce first NB columns of lower triangle
C
         DO 20 I = 1, NB
C
C           Update A(i:n,i)
C
            A( I, I ) = DBLE( A( I, I ) )
            CALL ZLACGV( I-1, W( I, 1 ), LDW )
            CALL ZGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ),
     $                  LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL ZLACGV( I-1, W( I, 1 ), LDW )
            CALL ZLACGV( I-1, A( I, 1 ), LDA )
            CALL ZGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ),
     $                  LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            CALL ZLACGV( I-1, A( I, 1 ), LDA )
            A( I, I ) = DBLE( A( I, I ) )
            IF( I.LT.N ) THEN
C
C              Generate elementary reflector H(i) to annihilate
C              A(i+2:n,i)
C
               ALPHA = A( I+1, I )
               CALL ZLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1,
     $                      TAU( I ) )
               E( I ) = ALPHA
               A( I+1, I ) = ONE
C
C              Compute W(i+1:n,i)
C
               CALL ZHEMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', N-I, I-1, ONE,
     $                     W( I+1, 1 ), LDW, A( I+1, I ), 1, ZERO,
     $                     W( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', N-I, I-1, ONE,
     $                     A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO,
     $                     W( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ),
     $                     LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL ZSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*ZDOTC( N-I, W( I+1, I ), 1,
     $                 A( I+1, I ), 1 )
               CALL ZAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
C
   20    CONTINUE
      END IF
C
      RETURN
C
C     End of ZLATRD
C
      END

      SUBROUTINE ZUNGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
C     ..
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT,
     $                   NB, NBMIN, NX
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNG2L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
C
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            NB = ILAENV( 1, 'ZUNGQL', ' ', M, N, K, -1 )
            LWKOPT = N*NB
         END IF
         WORK( 1 ) = LWKOPT
C
         IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
         END IF
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         RETURN
      END IF
C
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         NX = MAX( 0, ILAENV( 3, 'ZUNGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
C
C              Not enough workspace to use optimal NB:  reduce NB and
C              determine the minimum value of NB.
C
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZUNGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
C
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
C
C        Use blocked code after the first block.
C        The last kk columns are handled by the block method.
C
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
C
C        Set A(m-kk+1:m,1:n-kk) to zero.
C
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
C
C     Use unblocked code for the first or only block.
C
      CALL ZUNG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
C
      IF( KK.GT.0 ) THEN
C
C        Use blocked code
C
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
C
C              Form the triangular factor of the block reflector
C              H = H(i+ib-1) . . . H(i+1) H(i)
C
               CALL ZLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB,
     $                      A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
C
C              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
C
               CALL ZLARFB( 'Left', 'No transpose', 'Backward',
     $                      'Columnwise', M-K+I+IB-1, N-K+I-1, IB,
     $                      A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA,
     $                      WORK( IB+1 ), LDWORK )
            END IF
C
C           Apply H to rows 1:m-k+i+ib-1 of current block
C
            CALL ZUNG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA,
     $                   TAU( I ), WORK, IINFO )
C
C           Set rows m-k+i+ib:m of current block to zero
C
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
C
      WORK( 1 ) = IWS
      RETURN
C
C     End of ZUNGQL
C
      END

      DOUBLE PRECISION FUNCTION ZLANHE( NORM, UPLO, N, A, LDA, WORK )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         A( LDA, * )
C     ..
C
C =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
C     ..
C     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
C     ..
C     .. External Subroutines ..
      EXTERNAL           ZLASSQ
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, SQRT
C     ..
C     .. Executable Statements ..
C
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
C
C        Find max(abs(A(i,j))).
C
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J - 1
                  SUM = ABS( A( I, J ) )
                  IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   10          CONTINUE
               SUM = ABS( DBLE( A( J, J ) ) )
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               SUM = ABS( DBLE( A( J, J ) ) )
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
               DO 30 I = J + 1, N
                  SUM = ABS( A( I, J ) )
                  IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
C
C        Find normI(A) ( = norm1(A), since A is hermitian).
C
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( DBLE( A( J, J ) ) )
   60       CONTINUE
            DO 70 I = 1, N
               SUM = WORK( I )
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( DBLE( A( J, J ) ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
C
C        Find normF(A).
C
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL ZLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL ZLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         DO 130 I = 1, N
            IF( DBLE( A( I, I ) ).NE.ZERO ) THEN
               ABSA = ABS( DBLE( A( I, I ) ) )
               IF( SCALE.LT.ABSA ) THEN
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + ( ABSA / SCALE )**2
               END IF
            END IF
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
C
      ZLANHE = VALUE
      RETURN
C
C     End of ZLANHE
C
      END

      INTEGER FUNCTION IPARAM2STAGE( ISPEC, NAME, OPTS, 
     $                              NI, NBI, IBI, NXI )
      IMPLICIT NONE
C
C  -- LAPACK auxiliary routine (version 3.8.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, NI, NBI, IBI, NXI
C
C  ================================================================
C     ..
C     .. Local Scalars ..
      INTEGER            I, IC, IZ, KD, IB, LHOUS, LWORK, NTHREADS,
     $                   FACTOPTNB, QROPTNB, LQOPTNB
      LOGICAL            RPREC, CPREC
      CHARACTER          PREC*1, ALGO*3, STAG*5, SUBNAM*12, VECT*1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, MAX
C     ..
C     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
C     ..
C     .. Executable Statements ..
C
C     Invalid value for ISPEC
C
      IF( (ISPEC.LT.17).OR.(ISPEC.GT.21) ) THEN
          IPARAM2STAGE = -1
          RETURN
      ENDIF
C
C     Get the number of threads
C      
      NTHREADS = 1
C      WRITE(*,*) 'IPARAM VOICI NTHREADS ISPEC ',NTHREADS, ISPEC
C
      IF( ISPEC .NE. 19 ) THEN
C
C        Convert NAME to upper case if the first character is lower case.
C
         IPARAM2STAGE = -1
         SUBNAM = NAME
         IC = ICHAR( SUBNAM( 1: 1 ) )
         IZ = ICHAR( 'Z' )
         IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
C
C           ASCII character set
C
            IF( IC.GE.97 .AND. IC.LE.122 ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               DO 100 I = 2, 12
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( IC.GE.97 .AND. IC.LE.122 )
     $               SUBNAM( I: I ) = CHAR( IC-32 )
  100          CONTINUE
            END IF
C
         ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
C
C           EBCDIC character set
C
            IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $          ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $          ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC+64 )
               DO 110 I = 2, 12
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $                ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $                ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $                I ) = CHAR( IC+64 )
  110          CONTINUE
            END IF
C
         ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
C
C           Prime machines:  ASCII+128
C
            IF( IC.GE.225 .AND. IC.LE.250 ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               DO 120 I = 2, 12
                 IC = ICHAR( SUBNAM( I: I ) )
                 IF( IC.GE.225 .AND. IC.LE.250 )
     $             SUBNAM( I: I ) = CHAR( IC-32 )
  120          CONTINUE
            END IF
         END IF
C
         PREC  = SUBNAM( 1: 1 )
         ALGO  = SUBNAM( 4: 6 )
         STAG  = SUBNAM( 8:12 )
         RPREC = PREC.EQ.'S' .OR. PREC.EQ.'D'
         CPREC = PREC.EQ.'C' .OR. PREC.EQ.'Z'
C
C        Invalid value for PRECISION
C      
         IF( .NOT.( RPREC .OR. CPREC ) ) THEN
             IPARAM2STAGE = -1
             RETURN
         ENDIF
      ENDIF
C      WRITE(*,*),'RPREC,CPREC ',RPREC,CPREC,
C     $           '   ALGO ',ALGO,'    STAGE ',STAG
C      
C
      IF (( ISPEC .EQ. 17 ) .OR. ( ISPEC .EQ. 18 )) THEN 
C
C     ISPEC = 17, 18:  block size KD, IB
C     Could be also dependent from N but for now it
C     depend only on sequential or parallel
C
         IF( NTHREADS.GT.4 ) THEN
            IF( CPREC ) THEN
               KD = 128
               IB = 32
            ELSE
               KD = 160
               IB = 40
            ENDIF
         ELSE IF( NTHREADS.GT.1 ) THEN
            IF( CPREC ) THEN
               KD = 64
               IB = 32
            ELSE
               KD = 64
               IB = 32
            ENDIF
         ELSE
            IF( CPREC ) THEN
               KD = 16
               IB = 16
            ELSE
               KD = 32
               IB = 16
            ENDIF
         ENDIF
         IF( ISPEC.EQ.17 ) IPARAM2STAGE = KD
         IF( ISPEC.EQ.18 ) IPARAM2STAGE = IB
C
      ELSE IF ( ISPEC .EQ. 19 ) THEN
C
C     ISPEC = 19:  
C     LHOUS length of the Houselholder representation
C     matrix (V,T) of the second stage. should be >= 1.
C
C     Will add the VECT OPTION HERE next release
         VECT  = OPTS(1:1)
         IF( VECT.EQ.'N' ) THEN
            LHOUS = MAX( 1, 4*NI )
         ELSE
C           This is not correct, it need to call the ALGO and the stage2
            LHOUS = MAX( 1, 4*NI ) + IBI
         ENDIF
         IF( LHOUS.GE.0 ) THEN
            IPARAM2STAGE = LHOUS
         ELSE
            IPARAM2STAGE = -1
         ENDIF
C
      ELSE IF ( ISPEC .EQ. 20 ) THEN
C
C     ISPEC = 20: (21 for future use)  
C     LWORK length of the workspace for 
C     either or both stages for TRD and BRD. should be >= 1.
C     TRD:
C     TRD_stage 1: = LT + LW + LS1 + LS2
C                  = LDT*KD + N*KD + N*MAX(KD,FACTOPTNB) + LDS2*KD 
C                    where LDT=LDS2=KD
C                  = N*KD + N*max(KD,FACTOPTNB) + 2*KD*KD
C     TRD_stage 2: = (2NB+1)*N + KD*NTHREADS
C     TRD_both   : = max(stage1,stage2) + AB ( AB=(KD+1)*N )
C                  = N*KD + N*max(KD+1,FACTOPTNB) 
C                    + max(2*KD*KD, KD*NTHREADS) 
C                    + (KD+1)*N
         LWORK        = -1
         SUBNAM(1:1)  = PREC
         SUBNAM(2:6)  = 'GEQRF'
         QROPTNB      = ILAENV( 1, SUBNAM, ' ', NI, NBI, -1, -1 )
         SUBNAM(2:6)  = 'GELQF'
         LQOPTNB      = ILAENV( 1, SUBNAM, ' ', NBI, NI, -1, -1 )
C        Could be QR or LQ for TRD and the max for BRD
         FACTOPTNB    = MAX(QROPTNB, LQOPTNB)
         IF( ALGO.EQ.'TRD' ) THEN
            IF( STAG.EQ.'2STAG' ) THEN
               LWORK = NI*NBI + NI*MAX(NBI+1,FACTOPTNB) 
     $              + MAX(2*NBI*NBI, NBI*NTHREADS) 
     $              + (NBI+1)*NI
            ELSE IF( (STAG.EQ.'HE2HB').OR.(STAG.EQ.'SY2SB') ) THEN
               LWORK = NI*NBI + NI*MAX(NBI,FACTOPTNB) + 2*NBI*NBI
            ELSE IF( (STAG.EQ.'HB2ST').OR.(STAG.EQ.'SB2ST') ) THEN
               LWORK = (2*NBI+1)*NI + NBI*NTHREADS
            ENDIF
         ELSE IF( ALGO.EQ.'BRD' ) THEN
            IF( STAG.EQ.'2STAG' ) THEN
               LWORK = 2*NI*NBI + NI*MAX(NBI+1,FACTOPTNB) 
     $              + MAX(2*NBI*NBI, NBI*NTHREADS) 
     $              + (NBI+1)*NI
            ELSE IF( STAG.EQ.'GE2GB' ) THEN
               LWORK = NI*NBI + NI*MAX(NBI,FACTOPTNB) + 2*NBI*NBI
            ELSE IF( STAG.EQ.'GB2BD' ) THEN
               LWORK = (3*NBI+1)*NI + NBI*NTHREADS
            ENDIF
         ENDIF
         LWORK = MAX ( 1, LWORK )

         IF( LWORK.GT.0 ) THEN
            IPARAM2STAGE = LWORK
         ELSE
            IPARAM2STAGE = -1
         ENDIF
C
      ELSE IF ( ISPEC .EQ. 21 ) THEN
C
C     ISPEC = 21 for future use 
         IPARAM2STAGE = NXI
      ENDIF
C
C     ==== End of IPARAM2STAGE ====
C
      END

      RECURSIVE SUBROUTINE DPOTRF2( UPLO, N, A, LDA, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            N1, N2, IINFO
C     ..
C     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
C     ..
C     .. External Subroutines ..
      EXTERNAL           DSYRK, DTRSM, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRF2', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 )
     $   RETURN
C
C     N=1 case
C
      IF( N.EQ.1 ) THEN
C
C        Test for non-positive-definiteness
C
         IF( A( 1, 1 ).LE.ZERO.OR.DISNAN( A( 1, 1 ) ) ) THEN
            INFO = 1
            RETURN
         END IF
C
C        Factor
C
         A( 1, 1 ) = SQRT( A( 1, 1 ) )
C
C     Use recursive code
C
      ELSE
         N1 = N/2
         N2 = N-N1
C
C        Factor A11
C
         CALL DPOTRF2( UPLO, N1, A( 1, 1 ), LDA, IINFO )
         IF ( IINFO.NE.0 ) THEN
            INFO = IINFO
            RETURN
         END IF
C
C        Compute the Cholesky factorization A = U**T*U
C
         IF( UPPER ) THEN
C
C           Update and scale A12
C
            CALL DTRSM( 'L', 'U', 'T', 'N', N1, N2, ONE,
     $                  A( 1, 1 ), LDA, A( 1, N1+1 ), LDA )
C
C           Update and factor A22
C
            CALL DSYRK( UPLO, 'T', N2, N1, -ONE, A( 1, N1+1 ), LDA,
     $                  ONE, A( N1+1, N1+1 ), LDA )
            CALL DPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
            IF ( IINFO.NE.0 ) THEN
               INFO = IINFO + N1
               RETURN
            END IF
C
C        Compute the Cholesky factorization A = L*L**T
C
         ELSE
C
C           Update and scale A21
C
            CALL DTRSM( 'R', 'L', 'T', 'N', N2, N1, ONE,
     $                  A( 1, 1 ), LDA, A( N1+1, 1 ), LDA )
C
C           Update and factor A22
C
            CALL DSYRK( UPLO, 'N', N2, N1, -ONE, A( N1+1, 1 ), LDA,
     $                  ONE, A( N1+1, N1+1 ), LDA )
            CALL DPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
            IF ( IINFO.NE.0 ) THEN
               INFO = IINFO + N1
               RETURN
            END IF
         END IF
      END IF
      RETURN
C
C     End of DPOTRF2
C
      END

      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          CMACH
C     ..
C
C =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, RMACH
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT,
     $                   MINEXPONENT, RADIX, TINY
C     ..
C     .. Executable Statements ..
C
C
C     Assume rounding, not chopping. Always.
C
      RND = ONE
C
      IF( ONE.EQ.RND ) THEN
         EPS = EPSILON(ZERO) * 0.5
      ELSE
         EPS = EPSILON(ZERO)
      END IF
C
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         SFMIN = TINY(ZERO)
         SMALL = ONE / HUGE(ZERO)
         IF( SMALL.GE.SFMIN ) THEN
C
C           Use SMALL plus a bit, to avoid the possibility of rounding
C           causing overflow when computing  1/sfmin.
C
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
C
      DLAMCH = RMACH
      RETURN
C
C     End of DLAMCH
C
      END

      SUBROUTINE ZLACGV( N, X, INCX )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INCX, N
C     ..
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C =====================================================================
C
C     .. Local Scalars ..
      INTEGER            I, IOFF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
C     ..
C     .. Executable Statements ..
C
      IF( INCX.EQ.1 ) THEN
         DO 10 I = 1, N
            X( I ) = DCONJG( X( I ) )
   10    CONTINUE
      ELSE
         IOFF = 1
         IF( INCX.LT.0 )
     $      IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = DCONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      END IF
      RETURN
C
C     End of ZLACGV
C
      END

      SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
C     ..
C     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
C     ..
C     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZGERC
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAZLR, ILAZLC
      EXTERNAL           LSAME, ILAZLR, ILAZLC
C     ..
C     .. Executable Statements ..
C
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
C     Set up variables for scanning V.  LASTV begins pointing to the end
C     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
C     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
C     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILAZLC(LASTV, N, C, LDC)
         ELSE
C     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILAZLR(M, LASTV, C, LDC)
         END IF
      END IF
C     Note that lastc.eq.0 renders the BLAS operations null; no special
C     case is needed at this level.
      IF( APPLYLEFT ) THEN
C
C        Form  H * C
C
         IF( LASTV.GT.0 ) THEN
C
C           w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1)
C
            CALL ZGEMV( 'Conjugate transpose', LASTV, LASTC, ONE,
     $           C, LDC, V, INCV, ZERO, WORK, 1 )
C
C           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H
C
            CALL ZGERC( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
C
C        Form  C * H
C
         IF( LASTV.GT.0 ) THEN
C
C           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
C
            CALL ZGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC,
     $           V, INCV, ZERO, WORK, 1 )
C
C           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H
C
            CALL ZGERC( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
C
C     End of ZLARF
C
      END

      SUBROUTINE ZLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2013
C
C     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
C     ..
C     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZGEMM, ZLACGV, ZTRMM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
C     ..
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
C
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'C'
      ELSE
         TRANST = 'N'
      END IF
C
      IF( LSAME( STOREV, 'C' ) ) THEN
C
         IF( LSAME( DIRECT, 'F' ) ) THEN
C
C           Let  V =  ( V1 )    (first K rows)
C                     ( V2 )
C           where  V1  is unit lower triangular.
C
            IF( LSAME( SIDE, 'L' ) ) THEN
C
C              Form  H * C  or  H**H * C  where  C = ( C1 )
C                                                    ( C2 )
C
C              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
C
C              W := C1**H
C
               DO 10 J = 1, K
                  CALL ZCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
   10          CONTINUE
C
C              W := W * V1
C
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
C
C                 W := W + C2**H * V2
C
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', N,
     $                        K, M-K, ONE, C( K+1, 1 ), LDC,
     $                        V( K+1, 1 ), LDV, ONE, WORK, LDWORK )
               END IF
C
C              W := W * T**H  or  W * T
C
               CALL ZTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - V * W**H
C
               IF( M.GT.K ) THEN
C
C                 C2 := C2 - V2 * W**H
C
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose',
     $                        M-K, N, K, -ONE, V( K+1, 1 ), LDV, WORK,
     $                        LDWORK, ONE, C( K+1, 1 ), LDC )
               END IF
C
C              W := W * V1**H
C
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
C
C              C1 := C1 - W**H
C
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - DCONJG( WORK( I, J ) )
   20             CONTINUE
   30          CONTINUE
C
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
C
C              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
C
C              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
C
C              W := C1
C
               DO 40 J = 1, K
                  CALL ZCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
C
C              W := W * V1
C
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
C
C                 W := W + C2 * V2
C
                  CALL ZGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
C
C              W := W * T  or  W * T**H
C
               CALL ZTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - W * V**H
C
               IF( N.GT.K ) THEN
C
C                 C2 := C2 - W * V2**H
C
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M,
     $                        N-K, K, -ONE, WORK, LDWORK, V( K+1, 1 ),
     $                        LDV, ONE, C( 1, K+1 ), LDC )
               END IF
C
C              W := W * V1**H
C
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
C
C              C1 := C1 - W
C
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
C
         ELSE
C
C           Let  V =  ( V1 )
C                     ( V2 )    (last K rows)
C           where  V2  is unit upper triangular.
C
            IF( LSAME( SIDE, 'L' ) ) THEN
C
C              Form  H * C  or  H**H * C  where  C = ( C1 )
C                                                    ( C2 )
C
C              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
C
C              W := C2**H
C
               DO 70 J = 1, K
                  CALL ZCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
   70          CONTINUE
C
C              W := W * V2
C
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
C
C                 W := W + C1**H * V1
C
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', N,
     $                        K, M-K, ONE, C, LDC, V, LDV, ONE, WORK,
     $                        LDWORK )
               END IF
C
C              W := W * T**H  or  W * T
C
               CALL ZTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - V * W**H
C
               IF( M.GT.K ) THEN
C
C                 C1 := C1 - V1 * W**H
C
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose',
     $                        M-K, N, K, -ONE, V, LDV, WORK, LDWORK,
     $                        ONE, C, LDC )
               END IF
C
C              W := W * V2**H
C
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V( M-K+1, 1 ), LDV, WORK,
     $                     LDWORK )
C
C              C2 := C2 - W**H
C
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) -
     $                               DCONJG( WORK( I, J ) )
   80             CONTINUE
   90          CONTINUE
C
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
C
C              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
C
C              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
C
C              W := C2
C
               DO 100 J = 1, K
                  CALL ZCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
C
C              W := W * V2
C
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
C
C                 W := W + C1 * V1
C
                  CALL ZGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
C
C              W := W * T  or  W * T**H
C
               CALL ZTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - W * V**H
C
               IF( N.GT.K ) THEN
C
C                 C1 := C1 - W * V1**H
C
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M,
     $                        N-K, K, -ONE, WORK, LDWORK, V, LDV, ONE,
     $                        C, LDC )
               END IF
C
C              W := W * V2**H
C
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V( N-K+1, 1 ), LDV, WORK,
     $                     LDWORK )
C
C              C2 := C2 - W
C
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
C
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
C
         IF( LSAME( DIRECT, 'F' ) ) THEN
C
C           Let  V =  ( V1  V2 )    (V1: first K columns)
C           where  V1  is unit upper triangular.
C
            IF( LSAME( SIDE, 'L' ) ) THEN
C
C              Form  H * C  or  H**H * C  where  C = ( C1 )
C                                                    ( C2 )
C
C              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
C
C              W := C1**H
C
               DO 130 J = 1, K
                  CALL ZCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
  130          CONTINUE
C
C              W := W * V1**H
C
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
C
C                 W := W + C2**H * V2**H
C
                  CALL ZGEMM( 'Conjugate transpose',
     $                        'Conjugate transpose', N, K, M-K, ONE,
     $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
     $                        WORK, LDWORK )
               END IF
C
C              W := W * T**H  or  W * T
C
               CALL ZTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - V**H * W**H
C
               IF( M.GT.K ) THEN
C
C                 C2 := C2 - V2**H * W**H
C
                  CALL ZGEMM( 'Conjugate transpose',
     $                        'Conjugate transpose', M-K, N, K, -ONE,
     $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
C
C              W := W * V1
C
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
C
C              C1 := C1 - W**H
C
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - DCONJG( WORK( I, J ) )
  140             CONTINUE
  150          CONTINUE
C
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
C
C              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
C
C              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
C
C              W := C1
C
               DO 160 J = 1, K
                  CALL ZCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
C
C              W := W * V1**H
C
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
C
C                 W := W + C2 * V2**H
C
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M,
     $                        K, N-K, ONE, C( 1, K+1 ), LDC,
     $                        V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
               END IF
C
C              W := W * T  or  W * T**H
C
               CALL ZTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - W * V
C
               IF( N.GT.K ) THEN
C
C                 C2 := C2 - W * V2
C
                  CALL ZGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
C
C              W := W * V1
C
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
C
C              C1 := C1 - W
C
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
C
            END IF
C
         ELSE
C
C           Let  V =  ( V1  V2 )    (V2: last K columns)
C           where  V2  is unit lower triangular.
C
            IF( LSAME( SIDE, 'L' ) ) THEN
C
C              Form  H * C  or  H**H * C  where  C = ( C1 )
C                                                    ( C2 )
C
C              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
C
C              W := C2**H
C
               DO 190 J = 1, K
                  CALL ZCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
  190          CONTINUE
C
C              W := W * V2**H
C
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V( 1, M-K+1 ), LDV, WORK,
     $                     LDWORK )
               IF( M.GT.K ) THEN
C
C                 W := W + C1**H * V1**H
C
                  CALL ZGEMM( 'Conjugate transpose',
     $                        'Conjugate transpose', N, K, M-K, ONE, C,
     $                        LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
C
C              W := W * T**H  or  W * T
C
               CALL ZTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - V**H * W**H
C
               IF( M.GT.K ) THEN
C
C                 C1 := C1 - V1**H * W**H
C
                  CALL ZGEMM( 'Conjugate transpose',
     $                        'Conjugate transpose', M-K, N, K, -ONE, V,
     $                        LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
C
C              W := W * V2
C
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
C
C              C2 := C2 - W**H
C
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) -
     $                               DCONJG( WORK( I, J ) )
  200             CONTINUE
  210          CONTINUE
C
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
C
C              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
C
C              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
C
C              W := C2
C
               DO 220 J = 1, K
                  CALL ZCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
C
C              W := W * V2**H
C
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V( 1, N-K+1 ), LDV, WORK,
     $                     LDWORK )
               IF( N.GT.K ) THEN
C
C                 W := W + C1 * V1**H
C
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M,
     $                        K, N-K, ONE, C, LDC, V, LDV, ONE, WORK,
     $                        LDWORK )
               END IF
C
C              W := W * T  or  W * T**H
C
               CALL ZTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
C
C              C := C - W * V
C
               IF( N.GT.K ) THEN
C
C                 C1 := C1 - W * V1
C
                  CALL ZGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
C
C              W := W * V2
C
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
C
C              C1 := C1 - W
C
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
C
            END IF
C
         END IF
      END IF
C
      RETURN
C
C     End of ZLARFB
C
      END

      SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
C
C  -- LAPACK auxiliary routine (version 3.8.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     November 2017
C
C     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*16         ALPHA, TAU
C     ..
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY3, DZNRM2
      COMPLEX*16         ZLADIV
      EXTERNAL           DLAMCH, DLAPY3, DZNRM2, ZLADIV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, SIGN
C     ..
C     .. External Subroutines ..
      EXTERNAL           ZDSCAL, ZSCAL
C     ..
C     .. Executable Statements ..
C
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
C
      XNORM = DZNRM2( N-1, X, INCX )
      ALPHR = DBLE( ALPHA )
      ALPHI = DIMAG( ALPHA )
C
      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
C
C        H  =  I
C
         TAU = ZERO
      ELSE
C
C        general case
C
         BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         RSAFMN = ONE / SAFMIN
C
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
C
C           XNORM, BETA may be inaccurate; scale X and recompute them
C
   10       CONTINUE
            KNT = KNT + 1
            CALL ZDSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) )
     $         GO TO 10
C
C           New BETA is at most 1, at least SAFMIN
C
            XNORM = DZNRM2( N-1, X, INCX )
            ALPHA = DCMPLX( ALPHR, ALPHI )
            BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         END IF
         TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
         ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
         CALL ZSCAL( N-1, ALPHA, X, INCX )
C
C        If ALPHA is subnormal, it may lose relative accuracy
C
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
C
      RETURN
C
C     End of ZLARFG
C
      END

      SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
C     ..
C     .. Array Arguments ..
      COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
C     ..
C     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZTRMV, ZGEMM
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF( N.EQ.0 )
     $   RETURN
C
      IF( LSAME( DIRECT, 'F' ) ) THEN
         PREVLASTV = N
         DO I = 1, K
            PREVLASTV = MAX( PREVLASTV, I )
            IF( TAU( I ).EQ.ZERO ) THEN
C
C              H(i)  =  I
C
               DO J = 1, I
                  T( J, I ) = ZERO
               END DO
            ELSE
C
C              general case
C
               IF( LSAME( STOREV, 'C' ) ) THEN
C                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( LASTV, I ).NE.ZERO ) EXIT
                  END DO
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * CONJG( V( I , J ) )
                  END DO
                  J = MIN( LASTV, PREVLASTV )
C
C                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**H * V(i:j,i)
C
                  CALL ZGEMV( 'Conjugate transpose', J-I, I-1,
     $                        -TAU( I ), V( I+1, 1 ), LDV,
     $                        V( I+1, I ), 1, ONE, T( 1, I ), 1 )
               ELSE
C                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( I, LASTV ).NE.ZERO ) EXIT
                  END DO
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * V( J , I )
                  END DO
                  J = MIN( LASTV, PREVLASTV )
C
C                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**H
C
                  CALL ZGEMM( 'N', 'C', I-1, 1, J-I, -TAU( I ),
     $                        V( 1, I+1 ), LDV, V( I, I+1 ), LDV,
     $                        ONE, T( 1, I ), LDT )
               END IF
C
C              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
C
               CALL ZTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
               IF( I.GT.1 ) THEN
                  PREVLASTV = MAX( PREVLASTV, LASTV )
               ELSE
                  PREVLASTV = LASTV
               END IF
             END IF
         END DO
      ELSE
         PREVLASTV = 1
         DO I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
C
C              H(i)  =  I
C
               DO J = I, K
                  T( J, I ) = ZERO
               END DO
            ELSE
C
C              general case
C
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
C                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( LASTV, I ).NE.ZERO ) EXIT
                     END DO
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * CONJG( V( N-K+I , J ) )
                     END DO
                     J = MAX( LASTV, PREVLASTV )
C
C                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**H * V(j:n-k+i,i)
C
                     CALL ZGEMV( 'Conjugate transpose', N-K+I-J, K-I,
     $                           -TAU( I ), V( J, I+1 ), LDV, V( J, I ),
     $                           1, ONE, T( I+1, I ), 1 )
                  ELSE
C                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( I, LASTV ).NE.ZERO ) EXIT
                     END DO
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * V( J, N-K+I )
                     END DO
                     J = MAX( LASTV, PREVLASTV )
C
C                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**H
C
                     CALL ZGEMM( 'N', 'C', K-I, 1, N-K+I-J, -TAU( I ),
     $                           V( I+1, J ), LDV, V( I, J ), LDV,
     $                           ONE, T( I+1, I ), LDT )
                  END IF
C
C                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
C
                  CALL ZTRMV( 'Lower', 'No transpose', 'Non-unit', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
                  IF( I.GT.1 ) THEN
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
C
C     End of ZLARFT
C
      END

      COMPLEX*16     FUNCTION ZLADIV( X, Y )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      COMPLEX*16         X, Y
C     ..
C
C  =====================================================================
C
C     .. Local Scalars ..
      DOUBLE PRECISION   ZI, ZR
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLADIV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DIMAG
C     ..
C     .. Executable Statements ..
C
      CALL DLADIV( DBLE( X ), DIMAG( X ), DBLE( Y ), DIMAG( Y ), ZR,
     $             ZI )
      ZLADIV = DCMPLX( ZR, ZI )
C
      RETURN
C
C     End of ZLADIV
C
      END

      SUBROUTINE ZLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
C     ..
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
C     ..
C     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH, DISNAN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
C
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
C
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO .OR. DISNAN(CFROM) ) THEN
         INFO = -4
      ELSE IF( DISNAN(CTO) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLASCL', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
C
C     Get machine parameters
C
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
C
      CFROMC = CFROM
      CTOC = CTO
C
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      IF( CFROM1.EQ.CFROMC ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = .TRUE.
         CTO1 = CTOC
      ELSE
         CTO1 = CTOC / BIGNUM
         IF( CTO1.EQ.CTOC ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            MUL = CTOC
            DONE = .TRUE.
            CFROMC = ONE
         ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
            MUL = SMLNUM
            DONE = .FALSE.
            CFROMC = CFROM1
         ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
            MUL = BIGNUM
            DONE = .FALSE.
            CTOC = CTO1
         ELSE
            MUL = CTOC / CFROMC
            DONE = .TRUE.
         END IF
      END IF
C
      IF( ITYPE.EQ.0 ) THEN
C
C        Full matrix
C
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
C
      ELSE IF( ITYPE.EQ.1 ) THEN
C
C        Lower triangular matrix
C
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
C
      ELSE IF( ITYPE.EQ.2 ) THEN
C
C        Upper triangular matrix
C
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
C
      ELSE IF( ITYPE.EQ.3 ) THEN
C
C        Upper Hessenberg matrix
C
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
C
      ELSE IF( ITYPE.EQ.4 ) THEN
C
C        Lower half of a symmetric band matrix
C
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
C
      ELSE IF( ITYPE.EQ.5 ) THEN
C
C        Upper half of a symmetric band matrix
C
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
C
      ELSE IF( ITYPE.EQ.6 ) THEN
C
C        Band matrix
C
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
C
      END IF
C
      IF( .NOT.DONE )
     $   GO TO 10
C
      RETURN
C
C     End of ZLASCL
C
      END

      SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
C     ..
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   TEMP1
C     ..
C     .. External Functions ..
      LOGICAL            DISNAN
      EXTERNAL           DISNAN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG
C     ..
C     .. Executable Statements ..
C
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            TEMP1 = ABS( DBLE( X( IX ) ) )
            IF( TEMP1.GT.ZERO.OR.DISNAN( TEMP1 ) ) THEN
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
            TEMP1 = ABS( DIMAG( X( IX ) ) )
            IF( TEMP1.GT.ZERO.OR.DISNAN( TEMP1 ) ) THEN
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
C
      RETURN
C
C     End of ZLASSQ
C
      END

      SUBROUTINE ZUNG2L( M, N, K, A, LDA, TAU, WORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
C     ..
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      INTEGER            I, II, J, L
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2L', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.LE.0 )
     $   RETURN
C
C     Initialise columns 1:n-k to columns of the unit matrix
C
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
C
      DO 40 I = 1, K
         II = N - K + I
C
C        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
C
         A( M-N+II, II ) = ONE
         CALL ZLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A,
     $               LDA, WORK )
         CALL ZSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
C
C        Set A(m-k+i+1:m,n-k+i) to zero
C
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
C
C     End of ZUNG2L
C
      END

      SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
C     ..
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNG2R
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      NB = ILAENV( 1, 'ZUNGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
C
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         NX = MAX( 0, ILAENV( 3, 'ZUNGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
C
C           Determine if workspace is large enough for blocked code.
C
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
C
C              Not enough workspace to use optimal NB:  reduce NB and
C              determine the minimum value of NB.
C
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZUNGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
C
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
C
C        Use blocked code after the last block.
C        The first kk columns are handled by the block method.
C
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
C
C        Set A(1:kk,kk+1:n) to zero.
C
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
C
C     Use unblocked code for the last or only block.
C
      IF( KK.LT.N )
     $   CALL ZUNG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
C
      IF( KK.GT.0 ) THEN
C
C        Use blocked code
C
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
C
C              Form the triangular factor of the block reflector
C              H = H(i) H(i+1) . . . H(i+ib-1)
C
               CALL ZLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
C
C              Apply H to A(i:m,i+ib:n) from the left
C
               CALL ZLARFB( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
C
C           Apply H to rows i:m of current block
C
            CALL ZUNG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
C
C           Set rows 1:i-1 of current block to zero
C
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
C
      WORK( 1 ) = IWS
      RETURN
C
C     End of ZUNGQR
C
      END

      SUBROUTINE ZUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
C     ..
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      INTEGER            I, J, L
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments
C
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2R', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.LE.0 )
     $   RETURN
C
C     Initialise columns k+1:n to columns of the unit matrix
C
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
C
      DO 40 I = K, 1, -1
C
C        Apply H(i) to A(i:m,i:n) from the left
C
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL ZLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL ZSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
C
C        Set A(1:i-1,i) to zero
C
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
C
C     End of ZUNG2R
C
      END

      INTEGER FUNCTION ILAZLR( M, N, A, LDA )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            M, N, LDA
C     ..
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16       ZERO
      PARAMETER ( ZERO = (0.0D+0, 0.0D+0) )
C     ..
C     .. Local Scalars ..
      INTEGER I, J
C     ..
C     .. Executable Statements ..
C
C     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILAZLR = M
      ELSE IF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILAZLR = M
      ELSE
C     Scan up each column tracking the last zero row seen.
         ILAZLR = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ILAZLR = MAX( ILAZLR, I )
         END DO
      END IF
      RETURN
      END

      INTEGER FUNCTION ILAZLC( M, N, A, LDA )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            M, N, LDA
C     ..
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16       ZERO
      PARAMETER ( ZERO = (0.0D+0, 0.0D+0) )
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Executable Statements ..
C
C     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 ) THEN
         ILAZLC = N
      ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILAZLC = N
      ELSE
C     Now scan each column from the end, returning with the first non-zero.
         DO ILAZLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILAZLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END

      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2017
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
      LOGICAL            X_IS_NAN, Y_IS_NAN
C     ..
C     .. External Functions ..
      LOGICAL            DISNAN
      EXTERNAL           DISNAN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
C     ..
C     .. Executable Statements ..
C
      X_IS_NAN = DISNAN( X )
      Y_IS_NAN = DISNAN( Y )
      IF ( X_IS_NAN ) DLAPY2 = X
      IF ( Y_IS_NAN ) DLAPY2 = Y
C
      IF ( .NOT.( X_IS_NAN.OR.Y_IS_NAN ) ) THEN
         XABS = ABS( X )
         YABS = ABS( Y )
         W = MAX( XABS, YABS )
         Z = MIN( XABS, YABS )
         IF( Z.EQ.ZERO ) THEN
            DLAPY2 = W
         ELSE
            DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
         END IF
      END IF
      RETURN
C
C     End of DLAPY2
C
      END

      DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y, Z
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, ZABS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
C     ..
C     .. Executable Statements ..
C
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
C     W can be zero for max(0,nan,0)
C     adding all three entries together will make sure
C     NaN will not disappear.
         DLAPY3 =  XABS + YABS + ZABS
      ELSE
         DLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+
     $            ( ZABS / W )**2 )
      END IF
      RETURN
C
C     End of DLAPY3
C
      END

      SUBROUTINE DLADIV( A, B, C, D, P, Q )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     January 2013
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   BS
      PARAMETER          ( BS = 2.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
C
C     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLADIV1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
C     ..
C     .. Executable Statements ..
C
      AA = A
      BB = B
      CC = C
      DD = D
      AB = MAX( ABS(A), ABS(B) )
      CD = MAX( ABS(C), ABS(D) )
      S = 1.0D0

      OV = DLAMCH( 'Overflow threshold' )
      UN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Epsilon' )
      BE = BS / (EPS*EPS)

      IF( AB >= HALF*OV ) THEN
         AA = HALF * AA
         BB = HALF * BB
         S  = TWO * S
      END IF
      IF( CD >= HALF*OV ) THEN
         CC = HALF * CC
         DD = HALF * DD
         S  = HALF * S
      END IF
      IF( AB <= UN*BS/EPS ) THEN
         AA = AA * BE
         BB = BB * BE
         S  = S / BE
      END IF
      IF( CD <= UN*BS/EPS ) THEN
         CC = CC * BE
         DD = DD * BE
         S  = S * BE
      END IF
      IF( ABS( D ).LE.ABS( C ) ) THEN
         CALL DLADIV1(AA, BB, CC, DD, P, Q)
      ELSE
         CALL DLADIV1(BB, AA, DD, CC, P, Q)
         Q = -Q
      END IF
      P = P * S
      Q = Q * S
C
      RETURN
C
C     End of DLADIV
C
      END

      SUBROUTINE DLADIV1( A, B, C, D, P, Q )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     January 2013
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
C
C     .. Local Scalars ..
      DOUBLE PRECISION   R, T
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLADIV2
      EXTERNAL           DLADIV2
C     ..
C     .. Executable Statements ..
C
      R = D / C
      T = ONE / (C + D * R)
      P = DLADIV2(A, B, C, D, R, T)
      A = -A
      Q = DLADIV2(B, A, C, D, R, T)
C
      RETURN
C
C     End of DLADIV1
C
      END

      DOUBLE PRECISION FUNCTION DLADIV2( A, B, C, D, R, T )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     January 2013
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, R, T
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
C
C     .. Local Scalars ..
      DOUBLE PRECISION   BR
C     ..
C     .. Executable Statements ..
C
      IF( R.NE.ZERO ) THEN
         BR = B * R
         IF( BR.NE.ZERO ) THEN
            DLADIV2 = (A + BR) * T
         ELSE
            DLADIV2 = A * T + (B * T) * R
         END IF
      ELSE
         DLADIV2 = (A + D * (B / C)) * T
      END IF
C
      RETURN
C
C     End of DLADIV12
C
      END

      SUBROUTINE ZSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
     $           ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
     $           LIWORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE
      INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, M, N
      DOUBLE PRECISION ABSTOL, VL, VU
C     ..
C     .. Array Arguments ..
      INTEGER            ISUPPZ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
      COMPLEX*16         Z( LDZ, * )
C     ..
C
C  =====================================================================
C
C     .. Local Scalars ..
      LOGICAL TRYRAC
C     ..
C     .. External Subroutines ..
      EXTERNAL ZSTEMR
C     ..
C     .. Executable Statements ..
      INFO = 0
      TRYRAC = .FALSE.

      CALL ZSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
     $                   M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK, LWORK,
     $                   IWORK, LIWORK, INFO )
C
C     End of ZSTEGR
C
      END

      SUBROUTINE ZSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
     $                   M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC, WORK, LWORK,
     $                   IWORK, LIWORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE
      LOGICAL            TRYRAC
      INTEGER            IL, INFO, IU, LDZ, NZC, LIWORK, LWORK, M, N
      DOUBLE PRECISION VL, VU
C     ..
C     .. Array Arguments ..
      INTEGER            ISUPPZ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
      COMPLEX*16         Z( LDZ, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, FOUR, MINRGP
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0,
     $                     FOUR = 4.0D0,
     $                     MINRGP = 1.0D-3 )
C     ..
C     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LQUERY, VALEIG, WANTZ, ZQUERY
      INTEGER            I, IBEGIN, IEND, IFIRST, IIL, IINDBL, IINDW,
     $                   IINDWK, IINFO, IINSPL, IIU, ILAST, IN, INDD,
     $                   INDE2, INDERR, INDGP, INDGRS, INDWRK, ITMP,
     $                   ITMP2, J, JBLK, JJ, LIWMIN, LWMIN, NSPLIT,
     $                   NZCMIN, OFFSET, WBEGIN, WEND
      DOUBLE PRECISION   BIGNUM, CS, EPS, PIVMIN, R1, R2, RMAX, RMIN,
     $                   RTOL1, RTOL2, SAFMIN, SCALE, SMLNUM, SN,
     $                   THRESH, TMP, TNRM, WL, WU
C     ..
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, DLAMCH, DLANST
C     ..
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAE2, DLAEV2, DLARRC, DLARRE, DLARRJ,
     $                   DLARRR, DLASRT, DSCAL, XERBLA, ZLARRV, ZSWAP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT


C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
C
      LQUERY = ( ( LWORK.EQ.-1 ).OR.( LIWORK.EQ.-1 ) )
      ZQUERY = ( NZC.EQ.-1 )

C     DSTEMR needs WORK of size 6*N, IWORK of size 3*N.
C     In addition, DLARRE needs WORK of size 6*N, IWORK of size 5*N.
C     Furthermore, ZLARRV needs WORK of size 12*N, IWORK of size 7*N.
      IF( WANTZ ) THEN
         LWMIN = 18*N
         LIWMIN = 10*N
      ELSE
C        need less workspace if only the eigenvalues are wanted
         LWMIN = 12*N
         LIWMIN = 8*N
      ENDIF

      WL = ZERO
      WU = ZERO
      IIL = 0
      IIU = 0
      NSPLIT = 0

      IF( VALEIG ) THEN
C        We do not reference VL, VU in the cases RANGE = 'I','A'
C        The interval (WL, WU] contains all the wanted eigenvalues.
C        It is either given by the user or computed in DLARRE.
         WL = VL
         WU = VU
      ELSEIF( INDEIG ) THEN
C        We do not reference IL, IU in the cases RANGE = 'V','A'
         IIL = IL
         IIU = IU
      ENDIF
C
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( VALEIG .AND. N.GT.0 .AND. WU.LE.WL ) THEN
         INFO = -7
      ELSE IF( INDEIG .AND. ( IIL.LT.1 .OR. IIL.GT.N ) ) THEN
         INFO = -8
      ELSE IF( INDEIG .AND. ( IIU.LT.IIL .OR. IIU.GT.N ) ) THEN
         INFO = -9
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -13
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -17
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -19
      END IF
C
C     Get machine constants.
C
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
C
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
C
         IF( WANTZ .AND. ALLEIG ) THEN
            NZCMIN = N
         ELSE IF( WANTZ .AND. VALEIG ) THEN
            CALL DLARRC( 'T', N, VL, VU, D, E, SAFMIN,
     $                            NZCMIN, ITMP, ITMP2, INFO )
         ELSE IF( WANTZ .AND. INDEIG ) THEN
            NZCMIN = IIU-IIL+1
         ELSE
C           WANTZ .EQ. FALSE.
            NZCMIN = 0
         ENDIF
         IF( ZQUERY .AND. INFO.EQ.0 ) THEN
            Z( 1,1 ) = NZCMIN
         ELSE IF( NZC.LT.NZCMIN .AND. .NOT.ZQUERY ) THEN
            INFO = -14
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
C
         CALL XERBLA( 'ZSTEMR', -INFO )
C
         RETURN
      ELSE IF( LQUERY .OR. ZQUERY ) THEN
         RETURN
      END IF
C
C     Handle N = 0, 1, and 2 cases immediately
C
      M = 0
      IF( N.EQ.0 )
     $   RETURN
C
      IF( N.EQ.1 ) THEN
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = D( 1 )
         ELSE
            IF( WL.LT.D( 1 ) .AND. WU.GE.D( 1 ) ) THEN
               M = 1
               W( 1 ) = D( 1 )
            END IF
         END IF
         IF( WANTZ.AND.(.NOT.ZQUERY) ) THEN
            Z( 1, 1 ) = ONE
            ISUPPZ(1) = 1
            ISUPPZ(2) = 1
         END IF
         RETURN
      END IF
C
      IF( N.EQ.2 ) THEN
         IF( .NOT.WANTZ ) THEN
            CALL DLAE2( D(1), E(1), D(2), R1, R2 )
         ELSE IF( WANTZ.AND.(.NOT.ZQUERY) ) THEN
            CALL DLAEV2( D(1), E(1), D(2), R1, R2, CS, SN )
         END IF
         IF( ALLEIG.OR.
     $      (VALEIG.AND.(R2.GT.WL).AND.
     $                  (R2.LE.WU)).OR.
     $      (INDEIG.AND.(IIL.EQ.1)) ) THEN
            M = M+1
            W( M ) = R2
            IF( WANTZ.AND.(.NOT.ZQUERY) ) THEN
               Z( 1, M ) = -SN
               Z( 2, M ) = CS
C              Note: At most one of SN and CS can be zero.
               IF (SN.NE.ZERO) THEN
                  IF (CS.NE.ZERO) THEN
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 2
                  ELSE
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 1
                  END IF
               ELSE
                  ISUPPZ(2*M-1) = 2
                  ISUPPZ(2*M) = 2
               END IF
            ENDIF
         ENDIF
         IF( ALLEIG.OR.
     $      (VALEIG.AND.(R1.GT.WL).AND.
     $                  (R1.LE.WU)).OR.
     $      (INDEIG.AND.(IIU.EQ.2)) ) THEN
            M = M+1
            W( M ) = R1
            IF( WANTZ.AND.(.NOT.ZQUERY) ) THEN
               Z( 1, M ) = CS
               Z( 2, M ) = SN
C              Note: At most one of SN and CS can be zero.
               IF (SN.NE.ZERO) THEN
                  IF (CS.NE.ZERO) THEN
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 2
                  ELSE
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 1
                  END IF
               ELSE
                  ISUPPZ(2*M-1) = 2
                  ISUPPZ(2*M) = 2
               END IF
            ENDIF
         ENDIF
      ELSE

C        Continue with general N

         INDGRS = 1
         INDERR = 2*N + 1
         INDGP = 3*N + 1
         INDD = 4*N + 1
         INDE2 = 5*N + 1
         INDWRK = 6*N + 1
C
         IINSPL = 1
         IINDBL = N + 1
         IINDW = 2*N + 1
         IINDWK = 3*N + 1
C
C        Scale matrix to allowable range, if necessary.
C        The allowable range is related to the PIVMIN parameter; see the
C        comments in DLARRD.  The preference for scaling small values
C        up is heuristic; we expect users' matrices not to be close to the
C        RMAX threshold.
C
         SCALE = ONE
         TNRM = DLANST( 'M', N, D, E )
         IF( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) THEN
            SCALE = RMIN / TNRM
         ELSE IF( TNRM.GT.RMAX ) THEN
            SCALE = RMAX / TNRM
         END IF
         IF( SCALE.NE.ONE ) THEN
            CALL DSCAL( N, SCALE, D, 1 )
            CALL DSCAL( N-1, SCALE, E, 1 )
            TNRM = TNRM*SCALE
            IF( VALEIG ) THEN
C              If eigenvalues in interval have to be found,
C              scale (WL, WU] accordingly
               WL = WL*SCALE
               WU = WU*SCALE
            ENDIF
         END IF
C
C        Compute the desired eigenvalues of the tridiagonal after splitting
C        into smaller subblocks if the corresponding off-diagonal elements
C        are small
C        THRESH is the splitting parameter for DLARRE
C        A negative THRESH forces the old splitting criterion based on the
C        size of the off-diagonal. A positive THRESH switches to splitting
C        which preserves relative accuracy.
C
         IF( TRYRAC ) THEN
C           Test whether the matrix warrants the more expensive relative approach.
            CALL DLARRR( N, D, E, IINFO )
         ELSE
C           The user does not care about relative accurately eigenvalues
            IINFO = -1
         ENDIF
C        Set the splitting criterion
         IF (IINFO.EQ.0) THEN
            THRESH = EPS
         ELSE
            THRESH = -EPS
C           relative accuracy is desired but T does not guarantee it
            TRYRAC = .FALSE.
         ENDIF
C
         IF( TRYRAC ) THEN
C           Copy original diagonal, needed to guarantee relative accuracy
            CALL DCOPY(N,D,1,WORK(INDD),1)
         ENDIF
C        Store the squares of the offdiagonal values of T
         DO 5 J = 1, N-1
            WORK( INDE2+J-1 ) = E(J)**2
 5    CONTINUE

C        Set the tolerance parameters for bisection
         IF( .NOT.WANTZ ) THEN
C           DLARRE computes the eigenvalues to full precision.
            RTOL1 = FOUR * EPS
            RTOL2 = FOUR * EPS
         ELSE
C           DLARRE computes the eigenvalues to less than full precision.
C           ZLARRV will refine the eigenvalue approximations, and we only
C           need less accurate initial bisection in DLARRE.
C           Note: these settings do only affect the subset case and DLARRE
            RTOL1 = SQRT(EPS)
            RTOL2 = MAX( SQRT(EPS)*5.0D-3, FOUR * EPS )
         ENDIF
         CALL DLARRE( RANGE, N, WL, WU, IIL, IIU, D, E,
     $             WORK(INDE2), RTOL1, RTOL2, THRESH, NSPLIT,
     $             IWORK( IINSPL ), M, W, WORK( INDERR ),
     $             WORK( INDGP ), IWORK( IINDBL ),
     $             IWORK( IINDW ), WORK( INDGRS ), PIVMIN,
     $             WORK( INDWRK ), IWORK( IINDWK ), IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 10 + ABS( IINFO )
            RETURN
         END IF
C        Note that if RANGE .NE. 'V', DLARRE computes bounds on the desired
C        part of the spectrum. All desired eigenvalues are contained in
C        (WL,WU]


         IF( WANTZ ) THEN
C
C           Compute the desired eigenvectors corresponding to the computed
C           eigenvalues
C
            CALL ZLARRV( N, WL, WU, D, E,
     $                PIVMIN, IWORK( IINSPL ), M,
     $                1, M, MINRGP, RTOL1, RTOL2,
     $                W, WORK( INDERR ), WORK( INDGP ), IWORK( IINDBL ),
     $                IWORK( IINDW ), WORK( INDGRS ), Z, LDZ,
     $                ISUPPZ, WORK( INDWRK ), IWORK( IINDWK ), IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = 20 + ABS( IINFO )
               RETURN
            END IF
         ELSE
C           DLARRE computes eigenvalues of the (shifted) root representation
C           ZLARRV returns the eigenvalues of the unshifted matrix.
C           However, if the eigenvectors are not desired by the user, we need
C           to apply the corresponding shifts from DLARRE to obtain the
C           eigenvalues of the original matrix.
            DO 20 J = 1, M
               ITMP = IWORK( IINDBL+J-1 )
               W( J ) = W( J ) + E( IWORK( IINSPL+ITMP-1 ) )
 20      CONTINUE
         END IF
C

         IF ( TRYRAC ) THEN
C           Refine computed eigenvalues so that they are relatively accurate
C           with respect to the original matrix T.
            IBEGIN = 1
            WBEGIN = 1
            DO 39  JBLK = 1, IWORK( IINDBL+M-1 )
               IEND = IWORK( IINSPL+JBLK-1 )
               IN = IEND - IBEGIN + 1
               WEND = WBEGIN - 1
C              check if any eigenvalues have to be refined in this block
 36         CONTINUE
               IF( WEND.LT.M ) THEN
                  IF( IWORK( IINDBL+WEND ).EQ.JBLK ) THEN
                     WEND = WEND + 1
                     GO TO 36
                  END IF
               END IF
               IF( WEND.LT.WBEGIN ) THEN
                  IBEGIN = IEND + 1
                  GO TO 39
               END IF

               OFFSET = IWORK(IINDW+WBEGIN-1)-1
               IFIRST = IWORK(IINDW+WBEGIN-1)
               ILAST = IWORK(IINDW+WEND-1)
               RTOL2 = FOUR * EPS
               CALL DLARRJ( IN,
     $                   WORK(INDD+IBEGIN-1), WORK(INDE2+IBEGIN-1),
     $                   IFIRST, ILAST, RTOL2, OFFSET, W(WBEGIN),
     $                   WORK( INDERR+WBEGIN-1 ),
     $                   WORK( INDWRK ), IWORK( IINDWK ), PIVMIN,
     $                   TNRM, IINFO )
               IBEGIN = IEND + 1
               WBEGIN = WEND + 1
 39      CONTINUE
         ENDIF
C
C        If matrix was scaled, then rescale eigenvalues appropriately.
C
         IF( SCALE.NE.ONE ) THEN
            CALL DSCAL( M, ONE / SCALE, W, 1 )
         END IF
      END IF
C
C     If eigenvalues are not in increasing order, then sort them,
C     possibly along with eigenvectors.
C
      IF( NSPLIT.GT.1 .OR. N.EQ.2 ) THEN
         IF( .NOT. WANTZ ) THEN
            CALL DLASRT( 'I', M, W, IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = 3
               RETURN
            END IF
         ELSE
            DO 60 J = 1, M - 1
               I = 0
               TMP = W( J )
               DO 50 JJ = J + 1, M
                  IF( W( JJ ).LT.TMP ) THEN
                     I = JJ
                     TMP = W( JJ )
                  END IF
 50            CONTINUE
               IF( I.NE.0 ) THEN
                  W( I ) = W( J )
                  W( J ) = TMP
                  IF( WANTZ ) THEN
                     CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
                     ITMP = ISUPPZ( 2*I-1 )
                     ISUPPZ( 2*I-1 ) = ISUPPZ( 2*J-1 )
                     ISUPPZ( 2*J-1 ) = ITMP
                     ITMP = ISUPPZ( 2*I )
                     ISUPPZ( 2*I ) = ISUPPZ( 2*J )
                     ISUPPZ( 2*J ) = ITMP
                  END IF
               END IF
 60         CONTINUE
         END IF
      ENDIF
C
C
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
C
C     End of ZSTEMR
C
      END

      SUBROUTINE ZLARRV( N, VL, VU, D, L, PIVMIN,
     $                   ISPLIT, M, DOL, DOU, MINRGP,
     $                   RTOL1, RTOL2, W, WERR, WGAP,
     $                   IBLOCK, INDEXW, GERS, Z, LDZ, ISUPPZ,
     $                   WORK, IWORK, INFO )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      INTEGER            DOL, DOU, INFO, LDZ, M, N
      DOUBLE PRECISION   MINRGP, PIVMIN, RTOL1, RTOL2, VL, VU
C     ..
C     .. Array Arguments ..
      INTEGER            IBLOCK( * ), INDEXW( * ), ISPLIT( * ),
     $                   ISUPPZ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), GERS( * ), L( * ), W( * ), WERR( * ),
     $                   WGAP( * ), WORK( * )
      COMPLEX*16        Z( LDZ, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      INTEGER            MAXITR
      PARAMETER          ( MAXITR = 10 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ) )
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, HALF
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0,
     $                     TWO = 2.0D0, THREE = 3.0D0,
     $                     FOUR = 4.0D0, HALF = 0.5D0)
C     ..
C     .. Local Scalars ..
      LOGICAL            ESKIP, NEEDBS, STP2II, TRYRQC, USEDBS, USEDRQ
      INTEGER            DONE, I, IBEGIN, IDONE, IEND, II, IINDC1,
     $                   IINDC2, IINDR, IINDWK, IINFO, IM, IN, INDEIG,
     $                   INDLD, INDLLD, INDWRK, ISUPMN, ISUPMX, ITER,
     $                   ITMP1, J, JBLK, K, MINIWSIZE, MINWSIZE, NCLUS,
     $                   NDEPTH, NEGCNT, NEWCLS, NEWFST, NEWFTT, NEWLST,
     $                   NEWSIZ, OFFSET, OLDCLS, OLDFST, OLDIEN, OLDLST,
     $                   OLDNCL, P, PARITY, Q, WBEGIN, WEND, WINDEX,
     $                   WINDMN, WINDPL, ZFROM, ZTO, ZUSEDL, ZUSEDU,
     $                   ZUSEDW
      INTEGER            INDIN1, INDIN2
      DOUBLE PRECISION   BSTRES, BSTW, EPS, FUDGE, GAP, GAPTOL, GL, GU,
     $                   LAMBDA, LEFT, LGAP, MINGMA, NRMINV, RESID,
     $                   RGAP, RIGHT, RQCORR, RQTOL, SAVGAP, SGNDEF,
     $                   SIGMA, SPDIAM, SSIGMA, TAU, TMP, TOL, ZTZ
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DLARRB, DLARRF, ZDSCAL, ZLAR1V,
     $                   ZLASET
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS, DBLE, MAX, MIN
      INTRINSIC DCMPLX
C     ..
C     .. Executable Statements ..
C     ..

      INFO = 0
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         RETURN
      END IF
C
C     The first N entries of WORK are reserved for the eigenvalues
      INDLD = N+1
      INDLLD= 2*N+1
      INDIN1 = 3*N + 1
      INDIN2 = 4*N + 1
      INDWRK = 5*N + 1
      MINWSIZE = 12 * N

      DO 5 I= 1,MINWSIZE
         WORK( I ) = ZERO
 5    CONTINUE

C     IWORK(IINDR+1:IINDR+N) hold the twist indices R for the
C     factorization used to compute the FP vector
      IINDR = 0
C     IWORK(IINDC1+1:IINC2+N) are used to store the clusters of the current
C     layer and the one above.
      IINDC1 = N
      IINDC2 = 2*N
      IINDWK = 3*N + 1

      MINIWSIZE = 7 * N
      DO 10 I= 1,MINIWSIZE
         IWORK( I ) = 0
 10   CONTINUE

      ZUSEDL = 1
      IF(DOL.GT.1) THEN
C        Set lower bound for use of Z
         ZUSEDL = DOL-1
      ENDIF
      ZUSEDU = M
      IF(DOU.LT.M) THEN
C        Set lower bound for use of Z
         ZUSEDU = DOU+1
      ENDIF
C     The width of the part of Z that is used
      ZUSEDW = ZUSEDU - ZUSEDL + 1


      CALL ZLASET( 'Full', N, ZUSEDW, CZERO, CZERO,
     $                    Z(1,ZUSEDL), LDZ )

      EPS = DLAMCH( 'Precision' )
      RQTOL = TWO * EPS
C
C     Set expert flags for standard code.
      TRYRQC = .TRUE.

      IF((DOL.EQ.1).AND.(DOU.EQ.M)) THEN
      ELSE
C        Only selected eigenpairs are computed. Since the other evalues
C        are not refined by RQ iteration, bisection has to compute to full
C        accuracy.
         RTOL1 = FOUR * EPS
         RTOL2 = FOUR * EPS
      ENDIF

C     The entries WBEGIN:WEND in W, WERR, WGAP correspond to the
C     desired eigenvalues. The support of the nonzero eigenvector
C     entries is contained in the interval IBEGIN:IEND.
C     Remark that if k eigenpairs are desired, then the eigenvectors
C     are stored in k contiguous columns of Z.

C     DONE is the number of eigenvectors already computed
      DONE = 0
      IBEGIN = 1
      WBEGIN = 1
      DO 170 JBLK = 1, IBLOCK( M )
         IEND = ISPLIT( JBLK )
         SIGMA = L( IEND )
C        Find the eigenvectors of the submatrix indexed IBEGIN
C        through IEND.
         WEND = WBEGIN - 1
 15      CONTINUE
         IF( WEND.LT.M ) THEN
            IF( IBLOCK( WEND+1 ).EQ.JBLK ) THEN
               WEND = WEND + 1
               GO TO 15
            END IF
         END IF
         IF( WEND.LT.WBEGIN ) THEN
            IBEGIN = IEND + 1
            GO TO 170
         ELSEIF( (WEND.LT.DOL).OR.(WBEGIN.GT.DOU) ) THEN
            IBEGIN = IEND + 1
            WBEGIN = WEND + 1
            GO TO 170
         END IF

C        Find local spectral diameter of the block
         GL = GERS( 2*IBEGIN-1 )
         GU = GERS( 2*IBEGIN )
         DO 20 I = IBEGIN+1 , IEND
            GL = MIN( GERS( 2*I-1 ), GL )
            GU = MAX( GERS( 2*I ), GU )
 20      CONTINUE
         SPDIAM = GU - GL

C        OLDIEN is the last index of the previous block
         OLDIEN = IBEGIN - 1
C        Calculate the size of the current block
         IN = IEND - IBEGIN + 1
C        The number of eigenvalues in the current block
         IM = WEND - WBEGIN + 1

C        This is for a 1x1 block
         IF( IBEGIN.EQ.IEND ) THEN
            DONE = DONE+1
            Z( IBEGIN, WBEGIN ) = DCMPLX( ONE, ZERO )
            ISUPPZ( 2*WBEGIN-1 ) = IBEGIN
            ISUPPZ( 2*WBEGIN ) = IBEGIN
            W( WBEGIN ) = W( WBEGIN ) + SIGMA
            WORK( WBEGIN ) = W( WBEGIN )
            IBEGIN = IEND + 1
            WBEGIN = WBEGIN + 1
            GO TO 170
         END IF

C        The desired (shifted) eigenvalues are stored in W(WBEGIN:WEND)
C        Note that these can be approximations, in this case, the corresp.
C        entries of WERR give the size of the uncertainty interval.
C        The eigenvalue approximations will be refined when necessary as
C        high relative accuracy is required for the computation of the
C        corresponding eigenvectors.
         CALL DCOPY( IM, W( WBEGIN ), 1,
     $                   WORK( WBEGIN ), 1 )

C        We store in W the eigenvalue approximations w.r.t. the original
C        matrix T.
         DO 30 I=1,IM
            W(WBEGIN+I-1) = W(WBEGIN+I-1)+SIGMA
 30      CONTINUE


C        NDEPTH is the current depth of the representation tree
         NDEPTH = 0
C        PARITY is either 1 or 0
         PARITY = 1
C        NCLUS is the number of clusters for the next level of the
C        representation tree, we start with NCLUS = 1 for the root
         NCLUS = 1
         IWORK( IINDC1+1 ) = 1
         IWORK( IINDC1+2 ) = IM

C        IDONE is the number of eigenvectors already computed in the current
C        block
         IDONE = 0
C        loop while( IDONE.LT.IM )
C        generate the representation tree for the current block and
C        compute the eigenvectors
   40    CONTINUE
         IF( IDONE.LT.IM ) THEN
C           This is a crude protection against infinitely deep trees
            IF( NDEPTH.GT.M ) THEN
               INFO = -2
               RETURN
            ENDIF
C           breadth first processing of the current level of the representation
C           tree: OLDNCL = number of clusters on current level
            OLDNCL = NCLUS
C           reset NCLUS to count the number of child clusters
            NCLUS = 0
C
            PARITY = 1 - PARITY
            IF( PARITY.EQ.0 ) THEN
               OLDCLS = IINDC1
               NEWCLS = IINDC2
            ELSE
               OLDCLS = IINDC2
               NEWCLS = IINDC1
            END IF
C           Process the clusters on the current level
            DO 150 I = 1, OLDNCL
               J = OLDCLS + 2*I
C              OLDFST, OLDLST = first, last index of current cluster.
C                               cluster indices start with 1 and are relative
C                               to WBEGIN when accessing W, WGAP, WERR, Z
               OLDFST = IWORK( J-1 )
               OLDLST = IWORK( J )
               IF( NDEPTH.GT.0 ) THEN
C                 Retrieve relatively robust representation (RRR) of cluster
C                 that has been computed at the previous level
C                 The RRR is stored in Z and overwritten once the eigenvectors
C                 have been computed or when the cluster is refined

                  IF((DOL.EQ.1).AND.(DOU.EQ.M)) THEN
C                    Get representation from location of the leftmost evalue
C                    of the cluster
                     J = WBEGIN + OLDFST - 1
                  ELSE
                     IF(WBEGIN+OLDFST-1.LT.DOL) THEN
C                       Get representation from the left end of Z array
                        J = DOL - 1
                     ELSEIF(WBEGIN+OLDFST-1.GT.DOU) THEN
C                       Get representation from the right end of Z array
                        J = DOU
                     ELSE
                        J = WBEGIN + OLDFST - 1
                     ENDIF
                  ENDIF
                  DO 45 K = 1, IN - 1
                     D( IBEGIN+K-1 ) = DBLE( Z( IBEGIN+K-1,
     $                                 J ) )
                     L( IBEGIN+K-1 ) = DBLE( Z( IBEGIN+K-1,
     $                                 J+1 ) )
   45             CONTINUE
                  D( IEND ) = DBLE( Z( IEND, J ) )
                  SIGMA = DBLE( Z( IEND, J+1 ) )

C                 Set the corresponding entries in Z to zero
                  CALL ZLASET( 'Full', IN, 2, CZERO, CZERO,
     $                         Z( IBEGIN, J), LDZ )
               END IF

C              Compute DL and DLL of current RRR
               DO 50 J = IBEGIN, IEND-1
                  TMP = D( J )*L( J )
                  WORK( INDLD-1+J ) = TMP
                  WORK( INDLLD-1+J ) = TMP*L( J )
   50          CONTINUE

               IF( NDEPTH.GT.0 ) THEN
C                 P and Q are index of the first and last eigenvalue to compute
C                 within the current block
                  P = INDEXW( WBEGIN-1+OLDFST )
                  Q = INDEXW( WBEGIN-1+OLDLST )
C                 Offset for the arrays WORK, WGAP and WERR, i.e., the P-OFFSET
C                 through the Q-OFFSET elements of these arrays are to be used.
C                  OFFSET = P-OLDFST
                  OFFSET = INDEXW( WBEGIN ) - 1
C                 perform limited bisection (if necessary) to get approximate
C                 eigenvalues to the precision needed.
                  CALL DLARRB( IN, D( IBEGIN ),
     $                         WORK(INDLLD+IBEGIN-1),
     $                         P, Q, RTOL1, RTOL2, OFFSET,
     $                         WORK(WBEGIN),WGAP(WBEGIN),WERR(WBEGIN),
     $                         WORK( INDWRK ), IWORK( IINDWK ),
     $                         PIVMIN, SPDIAM, IN, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     INFO = -1
                     RETURN
                  ENDIF
C                 We also recompute the extremal gaps. W holds all eigenvalues
C                 of the unshifted matrix and must be used for computation
C                 of WGAP, the entries of WORK might stem from RRRs with
C                 different shifts. The gaps from WBEGIN-1+OLDFST to
C                 WBEGIN-1+OLDLST are correctly computed in DLARRB.
C                 However, we only allow the gaps to become greater since
C                 this is what should happen when we decrease WERR
                  IF( OLDFST.GT.1) THEN
                     WGAP( WBEGIN+OLDFST-2 ) =
     $             MAX(WGAP(WBEGIN+OLDFST-2),
     $                 W(WBEGIN+OLDFST-1)-WERR(WBEGIN+OLDFST-1)
     $                 - W(WBEGIN+OLDFST-2)-WERR(WBEGIN+OLDFST-2) )
                  ENDIF
                  IF( WBEGIN + OLDLST -1 .LT. WEND ) THEN
                     WGAP( WBEGIN+OLDLST-1 ) =
     $               MAX(WGAP(WBEGIN+OLDLST-1),
     $                   W(WBEGIN+OLDLST)-WERR(WBEGIN+OLDLST)
     $                   - W(WBEGIN+OLDLST-1)-WERR(WBEGIN+OLDLST-1) )
                  ENDIF
C                 Each time the eigenvalues in WORK get refined, we store
C                 the newly found approximation with all shifts applied in W
                  DO 53 J=OLDFST,OLDLST
                     W(WBEGIN+J-1) = WORK(WBEGIN+J-1)+SIGMA
 53               CONTINUE
               END IF

C              Process the current node.
               NEWFST = OLDFST
               DO 140 J = OLDFST, OLDLST
                  IF( J.EQ.OLDLST ) THEN
C                    we are at the right end of the cluster, this is also the
C                    boundary of the child cluster
                     NEWLST = J
                  ELSE IF ( WGAP( WBEGIN + J -1).GE.
     $                    MINRGP* ABS( WORK(WBEGIN + J -1) ) ) THEN
C                    the right relative gap is big enough, the child cluster
C                    (NEWFST,..,NEWLST) is well separated from the following
                     NEWLST = J
                   ELSE
C                    inside a child cluster, the relative gap is not
C                    big enough.
                     GOTO 140
                  END IF

C                 Compute size of child cluster found
                  NEWSIZ = NEWLST - NEWFST + 1

C                 NEWFTT is the place in Z where the new RRR or the computed
C                 eigenvector is to be stored
                  IF((DOL.EQ.1).AND.(DOU.EQ.M)) THEN
C                    Store representation at location of the leftmost evalue
C                    of the cluster
                     NEWFTT = WBEGIN + NEWFST - 1
                  ELSE
                     IF(WBEGIN+NEWFST-1.LT.DOL) THEN
C                       Store representation at the left end of Z array
                        NEWFTT = DOL - 1
                     ELSEIF(WBEGIN+NEWFST-1.GT.DOU) THEN
C                       Store representation at the right end of Z array
                        NEWFTT = DOU
                     ELSE
                        NEWFTT = WBEGIN + NEWFST - 1
                     ENDIF
                  ENDIF

                  IF( NEWSIZ.GT.1) THEN
C
C                    Current child is not a singleton but a cluster.
C                    Compute and store new representation of child.
C
C
C                    Compute left and right cluster gap.
C
C                    LGAP and RGAP are not computed from WORK because
C                    the eigenvalue approximations may stem from RRRs
C                    different shifts. However, W hold all eigenvalues
C                    of the unshifted matrix. Still, the entries in WGAP
C                    have to be computed from WORK since the entries
C                    in W might be of the same order so that gaps are not
C                    exhibited correctly for very close eigenvalues.
                     IF( NEWFST.EQ.1 ) THEN
                        LGAP = MAX( ZERO,
     $                       W(WBEGIN)-WERR(WBEGIN) - VL )
                    ELSE
                        LGAP = WGAP( WBEGIN+NEWFST-2 )
                     ENDIF
                     RGAP = WGAP( WBEGIN+NEWLST-1 )
C
C                    Compute left- and rightmost eigenvalue of child
C                    to high precision in order to shift as close
C                    as possible and obtain as large relative gaps
C                    as possible
C
                     DO 55 K =1,2
                        IF(K.EQ.1) THEN
                           P = INDEXW( WBEGIN-1+NEWFST )
                        ELSE
                           P = INDEXW( WBEGIN-1+NEWLST )
                        ENDIF
                        OFFSET = INDEXW( WBEGIN ) - 1
                        CALL DLARRB( IN, D(IBEGIN),
     $                       WORK( INDLLD+IBEGIN-1 ),P,P,
     $                       RQTOL, RQTOL, OFFSET,
     $                       WORK(WBEGIN),WGAP(WBEGIN),
     $                       WERR(WBEGIN),WORK( INDWRK ),
     $                       IWORK( IINDWK ), PIVMIN, SPDIAM,
     $                       IN, IINFO )
 55                  CONTINUE
C
                     IF((WBEGIN+NEWLST-1.LT.DOL).OR.
     $                  (WBEGIN+NEWFST-1.GT.DOU)) THEN
C                       if the cluster contains no desired eigenvalues
C                       skip the computation of that branch of the rep. tree
C
C                       We could skip before the refinement of the extremal
C                       eigenvalues of the child, but then the representation
C                       tree could be different from the one when nothing is
C                       skipped. For this reason we skip at this place.
                        IDONE = IDONE + NEWLST - NEWFST + 1
                        GOTO 139
                     ENDIF
C
C                    Compute RRR of child cluster.
C                    Note that the new RRR is stored in Z
C
C                    DLARRF needs LWORK = 2*N
                     CALL DLARRF( IN, D( IBEGIN ), L( IBEGIN ),
     $                         WORK(INDLD+IBEGIN-1),
     $                         NEWFST, NEWLST, WORK(WBEGIN),
     $                         WGAP(WBEGIN), WERR(WBEGIN),
     $                         SPDIAM, LGAP, RGAP, PIVMIN, TAU,
     $                         WORK( INDIN1 ), WORK( INDIN2 ),
     $                         WORK( INDWRK ), IINFO )
C                    In the complex case, DLARRF cannot write
C                    the new RRR directly into Z and needs an intermediate
C                    workspace
                     DO 56 K = 1, IN-1
                        Z( IBEGIN+K-1, NEWFTT ) =
     $                     DCMPLX( WORK( INDIN1+K-1 ), ZERO )
                        Z( IBEGIN+K-1, NEWFTT+1 ) =
     $                     DCMPLX( WORK( INDIN2+K-1 ), ZERO )
   56                CONTINUE
                     Z( IEND, NEWFTT ) =
     $                  DCMPLX( WORK( INDIN1+IN-1 ), ZERO )
                     IF( IINFO.EQ.0 ) THEN
C                       a new RRR for the cluster was found by DLARRF
C                       update shift and store it
                        SSIGMA = SIGMA + TAU
                        Z( IEND, NEWFTT+1 ) = DCMPLX( SSIGMA, ZERO )
C                       WORK() are the midpoints and WERR() the semi-width
C                       Note that the entries in W are unchanged.
                        DO 116 K = NEWFST, NEWLST
                           FUDGE =
     $                          THREE*EPS*ABS(WORK(WBEGIN+K-1))
                           WORK( WBEGIN + K - 1 ) =
     $                          WORK( WBEGIN + K - 1) - TAU
                           FUDGE = FUDGE +
     $                          FOUR*EPS*ABS(WORK(WBEGIN+K-1))
C                          Fudge errors
                           WERR( WBEGIN + K - 1 ) =
     $                          WERR( WBEGIN + K - 1 ) + FUDGE
C                          Gaps are not fudged. Provided that WERR is small
C                          when eigenvalues are close, a zero gap indicates
C                          that a new representation is needed for resolving
C                          the cluster. A fudge could lead to a wrong decision
C                          of judging eigenvalues 'separated' which in
C                          reality are not. This could have a negative impact
C                          on the orthogonality of the computed eigenvectors.
 116                    CONTINUE

                        NCLUS = NCLUS + 1
                        K = NEWCLS + 2*NCLUS
                        IWORK( K-1 ) = NEWFST
                        IWORK( K ) = NEWLST
                     ELSE
                        INFO = -2
                        RETURN
                     ENDIF
                  ELSE
C
C                    Compute eigenvector of singleton
C
                     ITER = 0
C
                     TOL = FOUR * LOG(DBLE(IN)) * EPS
C
                     K = NEWFST
                     WINDEX = WBEGIN + K - 1
                     WINDMN = MAX(WINDEX - 1,1)
                     WINDPL = MIN(WINDEX + 1,M)
                     LAMBDA = WORK( WINDEX )
                     DONE = DONE + 1
C                    Check if eigenvector computation is to be skipped
                     IF((WINDEX.LT.DOL).OR.
     $                  (WINDEX.GT.DOU)) THEN
                        ESKIP = .TRUE.
                        GOTO 125
                     ELSE
                        ESKIP = .FALSE.
                     ENDIF
                     LEFT = WORK( WINDEX ) - WERR( WINDEX )
                     RIGHT = WORK( WINDEX ) + WERR( WINDEX )
                     INDEIG = INDEXW( WINDEX )
C                    Note that since we compute the eigenpairs for a child,
C                    all eigenvalue approximations are w.r.t the same shift.
C                    In this case, the entries in WORK should be used for
C                    computing the gaps since they exhibit even very small
C                    differences in the eigenvalues, as opposed to the
C                    entries in W which might "look" the same.

                     IF( K .EQ. 1) THEN
C                       In the case RANGE='I' and with not much initial
C                       accuracy in LAMBDA and VL, the formula
C                       LGAP = MAX( ZERO, (SIGMA - VL) + LAMBDA )
C                       can lead to an overestimation of the left gap and
C                       thus to inadequately early RQI 'convergence'.
C                       Prevent this by forcing a small left gap.
                        LGAP = EPS*MAX(ABS(LEFT),ABS(RIGHT))
                     ELSE
                        LGAP = WGAP(WINDMN)
                     ENDIF
                     IF( K .EQ. IM) THEN
C                       In the case RANGE='I' and with not much initial
C                       accuracy in LAMBDA and VU, the formula
C                       can lead to an overestimation of the right gap and
C                       thus to inadequately early RQI 'convergence'.
C                       Prevent this by forcing a small right gap.
                        RGAP = EPS*MAX(ABS(LEFT),ABS(RIGHT))
                     ELSE
                        RGAP = WGAP(WINDEX)
                     ENDIF
                     GAP = MIN( LGAP, RGAP )
                     IF(( K .EQ. 1).OR.(K .EQ. IM)) THEN
C                       The eigenvector support can become wrong
C                       because significant entries could be cut off due to a
C                       large GAPTOL parameter in LAR1V. Prevent this.
                        GAPTOL = ZERO
                     ELSE
                        GAPTOL = GAP * EPS
                     ENDIF
                     ISUPMN = IN
                     ISUPMX = 1
C                    Update WGAP so that it holds the minimum gap
C                    to the left or the right. This is crucial in the
C                    case where bisection is used to ensure that the
C                    eigenvalue is refined up to the required precision.
C                    The correct value is restored afterwards.
                     SAVGAP = WGAP(WINDEX)
                     WGAP(WINDEX) = GAP
C                    We want to use the Rayleigh Quotient Correction
C                    as often as possible since it converges quadratically
C                    when we are close enough to the desired eigenvalue.
C                    However, the Rayleigh Quotient can have the wrong sign
C                    and lead us away from the desired eigenvalue. In this
C                    case, the best we can do is to use bisection.
                     USEDBS = .FALSE.
                     USEDRQ = .FALSE.
C                    Bisection is initially turned off unless it is forced
                     NEEDBS =  .NOT.TRYRQC
 120                 CONTINUE
C                    Check if bisection should be used to refine eigenvalue
                     IF(NEEDBS) THEN
C                       Take the bisection as new iterate
                        USEDBS = .TRUE.
                        ITMP1 = IWORK( IINDR+WINDEX )
                        OFFSET = INDEXW( WBEGIN ) - 1
                        CALL DLARRB( IN, D(IBEGIN),
     $                       WORK(INDLLD+IBEGIN-1),INDEIG,INDEIG,
     $                       ZERO, TWO*EPS, OFFSET,
     $                       WORK(WBEGIN),WGAP(WBEGIN),
     $                       WERR(WBEGIN),WORK( INDWRK ),
     $                       IWORK( IINDWK ), PIVMIN, SPDIAM,
     $                       ITMP1, IINFO )
                        IF( IINFO.NE.0 ) THEN
                           INFO = -3
                           RETURN
                        ENDIF
                        LAMBDA = WORK( WINDEX )
C                       Reset twist index from inaccurate LAMBDA to
C                       force computation of true MINGMA
                        IWORK( IINDR+WINDEX ) = 0
                     ENDIF
C                    Given LAMBDA, compute the eigenvector.
                     CALL ZLAR1V( IN, 1, IN, LAMBDA, D( IBEGIN ),
     $                    L( IBEGIN ), WORK(INDLD+IBEGIN-1),
     $                    WORK(INDLLD+IBEGIN-1),
     $                    PIVMIN, GAPTOL, Z( IBEGIN, WINDEX ),
     $                    .NOT.USEDBS, NEGCNT, ZTZ, MINGMA,
     $                    IWORK( IINDR+WINDEX ), ISUPPZ( 2*WINDEX-1 ),
     $                    NRMINV, RESID, RQCORR, WORK( INDWRK ) )
                     IF(ITER .EQ. 0) THEN
                        BSTRES = RESID
                        BSTW = LAMBDA
                     ELSEIF(RESID.LT.BSTRES) THEN
                        BSTRES = RESID
                        BSTW = LAMBDA
                     ENDIF
                     ISUPMN = MIN(ISUPMN,ISUPPZ( 2*WINDEX-1 ))
                     ISUPMX = MAX(ISUPMX,ISUPPZ( 2*WINDEX ))
                     ITER = ITER + 1

C                    sin alpha <= |resid|/gap
C                    Note that both the residual and the gap are
C                    proportional to the matrix, so ||T|| doesn't play
C                    a role in the quotient

C
C                    Convergence test for Rayleigh-Quotient iteration
C                    (omitted when Bisection has been used)
C
                     IF( RESID.GT.TOL*GAP .AND. ABS( RQCORR ).GT.
     $                    RQTOL*ABS( LAMBDA ) .AND. .NOT. USEDBS)
     $                    THEN
C                       We need to check that the RQCORR update doesn't
C                       move the eigenvalue away from the desired one and
C                       towards a neighbor. -> protection with bisection
                        IF(INDEIG.LE.NEGCNT) THEN
C                          The wanted eigenvalue lies to the left
                           SGNDEF = -ONE
                        ELSE
C                          The wanted eigenvalue lies to the right
                           SGNDEF = ONE
                        ENDIF
C                       We only use the RQCORR if it improves the
C                       the iterate reasonably.
                        IF( ( RQCORR*SGNDEF.GE.ZERO )
     $                       .AND.( LAMBDA + RQCORR.LE. RIGHT)
     $                       .AND.( LAMBDA + RQCORR.GE. LEFT)
     $                       ) THEN
                           USEDRQ = .TRUE.
C                          Store new midpoint of bisection interval in WORK
                           IF(SGNDEF.EQ.ONE) THEN
C                             The current LAMBDA is on the left of the true
C                             eigenvalue
                              LEFT = LAMBDA
C                             We prefer to assume that the error estimate
C                             is correct. We could make the interval not
C                             as a bracket but to be modified if the RQCORR
C                             chooses to. In this case, the RIGHT side should
C                             be modified as follows:
C                              RIGHT = MAX(RIGHT, LAMBDA + RQCORR)
                           ELSE
C                             The current LAMBDA is on the right of the true
C                             eigenvalue
                              RIGHT = LAMBDA
C                             See comment about assuming the error estimate is
C                             correct above.
C                              LEFT = MIN(LEFT, LAMBDA + RQCORR)
                           ENDIF
                           WORK( WINDEX ) =
     $                       HALF * (RIGHT + LEFT)
C                          Take RQCORR since it has the correct sign and
C                          improves the iterate reasonably
                           LAMBDA = LAMBDA + RQCORR
C                          Update width of error interval
                           WERR( WINDEX ) =
     $                             HALF * (RIGHT-LEFT)
                        ELSE
                           NEEDBS = .TRUE.
                        ENDIF
                        IF(RIGHT-LEFT.LT.RQTOL*ABS(LAMBDA)) THEN
C                             The eigenvalue is computed to bisection accuracy
C                             compute eigenvector and stop
                           USEDBS = .TRUE.
                           GOTO 120
                        ELSEIF( ITER.LT.MAXITR ) THEN
                           GOTO 120
                        ELSEIF( ITER.EQ.MAXITR ) THEN
                           NEEDBS = .TRUE.
                           GOTO 120
                        ELSE
                           INFO = 5
                           RETURN
                        END IF
                     ELSE
                        STP2II = .FALSE.
        IF(USEDRQ .AND. USEDBS .AND.
     $                     BSTRES.LE.RESID) THEN
                           LAMBDA = BSTW
                           STP2II = .TRUE.
                        ENDIF
                        IF (STP2II) THEN
C                          improve error angle by second step
                           CALL ZLAR1V( IN, 1, IN, LAMBDA,
     $                          D( IBEGIN ), L( IBEGIN ),
     $                          WORK(INDLD+IBEGIN-1),
     $                          WORK(INDLLD+IBEGIN-1),
     $                          PIVMIN, GAPTOL, Z( IBEGIN, WINDEX ),
     $                          .NOT.USEDBS, NEGCNT, ZTZ, MINGMA,
     $                          IWORK( IINDR+WINDEX ),
     $                          ISUPPZ( 2*WINDEX-1 ),
     $                          NRMINV, RESID, RQCORR, WORK( INDWRK ) )
                        ENDIF
                        WORK( WINDEX ) = LAMBDA
                     END IF
C
C                    Compute FP-vector support w.r.t. whole matrix
C
                     ISUPPZ( 2*WINDEX-1 ) = ISUPPZ( 2*WINDEX-1 )+OLDIEN
                     ISUPPZ( 2*WINDEX ) = ISUPPZ( 2*WINDEX )+OLDIEN
                     ZFROM = ISUPPZ( 2*WINDEX-1 )
                     ZTO = ISUPPZ( 2*WINDEX )
                     ISUPMN = ISUPMN + OLDIEN
                     ISUPMX = ISUPMX + OLDIEN
C                    Ensure vector is ok if support in the RQI has changed
                     IF(ISUPMN.LT.ZFROM) THEN
                        DO 122 II = ISUPMN,ZFROM-1
                           Z( II, WINDEX ) = ZERO
 122                    CONTINUE
                     ENDIF
                     IF(ISUPMX.GT.ZTO) THEN
                        DO 123 II = ZTO+1,ISUPMX
                           Z( II, WINDEX ) = ZERO
 123                    CONTINUE
                     ENDIF
                     CALL ZDSCAL( ZTO-ZFROM+1, NRMINV,
     $                       Z( ZFROM, WINDEX ), 1 )
 125                 CONTINUE
C                    Update W
                     W( WINDEX ) = LAMBDA+SIGMA
C                    Recompute the gaps on the left and right
C                    But only allow them to become larger and not
C                    smaller (which can only happen through "bad"
C                    cancellation and doesn't reflect the theory
C                    where the initial gaps are underestimated due
C                    to WERR being too crude.)
                     IF(.NOT.ESKIP) THEN
                        IF( K.GT.1) THEN
                           WGAP( WINDMN ) = MAX( WGAP(WINDMN),
     $                          W(WINDEX)-WERR(WINDEX)
     $                          - W(WINDMN)-WERR(WINDMN) )
                        ENDIF
                        IF( WINDEX.LT.WEND ) THEN
                           WGAP( WINDEX ) = MAX( SAVGAP,
     $                          W( WINDPL )-WERR( WINDPL )
     $                          - W( WINDEX )-WERR( WINDEX) )
                        ENDIF
                     ENDIF
                     IDONE = IDONE + 1
                  ENDIF
C                 here ends the code for the current child
C
 139              CONTINUE
C                 Proceed to any remaining child nodes
                  NEWFST = J + 1
 140           CONTINUE
 150        CONTINUE
            NDEPTH = NDEPTH + 1
            GO TO 40
         END IF
         IBEGIN = IEND + 1
         WBEGIN = WEND + 1
 170  CONTINUE
C

      RETURN
C
C     End of ZLARRV
C
      END

      SUBROUTINE ZLAR1V( N, B1, BN, LAMBDA, D, L, LD, LLD,
     $           PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ, MINGMA,
     $           R, ISUPPZ, NRMINV, RESID, RQCORR, WORK )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      LOGICAL            WANTNC
      INTEGER   B1, BN, N, NEGCNT, R
      DOUBLE PRECISION   GAPTOL, LAMBDA, MINGMA, NRMINV, PIVMIN, RESID,
     $                   RQCORR, ZTZ
C     ..
C     .. Array Arguments ..
      INTEGER            ISUPPZ( * )
      DOUBLE PRECISION   D( * ), L( * ), LD( * ), LLD( * ),
     $                  WORK( * )
      COMPLEX*16       Z( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )

C     ..
C     .. Local Scalars ..
      LOGICAL            SAWNAN1, SAWNAN2
      INTEGER            I, INDLPL, INDP, INDS, INDUMN, NEG1, NEG2, R1,
     $                   R2
      DOUBLE PRECISION   DMINUS, DPLUS, EPS, S, TMP
C     ..
C     .. External Functions ..
      LOGICAL DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DISNAN, DLAMCH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE
C     ..
C     .. Executable Statements ..
C
      EPS = DLAMCH( 'Precision' )


      IF( R.EQ.0 ) THEN
         R1 = B1
         R2 = BN
      ELSE
         R1 = R
         R2 = R
      END IF

C     Storage for LPLUS
      INDLPL = 0
C     Storage for UMINUS
      INDUMN = N
      INDS = 2*N + 1
      INDP = 3*N + 1

      IF( B1.EQ.1 ) THEN
         WORK( INDS ) = ZERO
      ELSE
         WORK( INDS+B1-1 ) = LLD( B1-1 )
      END IF

C
C     Compute the stationary transform (using the differential form)
C     until the index R2.
C
      SAWNAN1 = .FALSE.
      NEG1 = 0
      S = WORK( INDS+B1-1 ) - LAMBDA
      DO 50 I = B1, R1 - 1
         DPLUS = D( I ) + S
         WORK( INDLPL+I ) = LD( I ) / DPLUS
         IF(DPLUS.LT.ZERO) NEG1 = NEG1 + 1
         WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
         S = WORK( INDS+I ) - LAMBDA
 50   CONTINUE
      SAWNAN1 = DISNAN( S )
      IF( SAWNAN1 ) GOTO 60
      DO 51 I = R1, R2 - 1
         DPLUS = D( I ) + S
         WORK( INDLPL+I ) = LD( I ) / DPLUS
         WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
         S = WORK( INDS+I ) - LAMBDA
 51   CONTINUE
      SAWNAN1 = DISNAN( S )
C
 60   CONTINUE
      IF( SAWNAN1 ) THEN
C        Runs a slower version of the above loop if a NaN is detected
         NEG1 = 0
         S = WORK( INDS+B1-1 ) - LAMBDA
         DO 70 I = B1, R1 - 1
            DPLUS = D( I ) + S
            IF(ABS(DPLUS).LT.PIVMIN) DPLUS = -PIVMIN
            WORK( INDLPL+I ) = LD( I ) / DPLUS
            IF(DPLUS.LT.ZERO) NEG1 = NEG1 + 1
            WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
            IF( WORK( INDLPL+I ).EQ.ZERO )
     $                      WORK( INDS+I ) = LLD( I )
            S = WORK( INDS+I ) - LAMBDA
 70      CONTINUE
         DO 71 I = R1, R2 - 1
            DPLUS = D( I ) + S
            IF(ABS(DPLUS).LT.PIVMIN) DPLUS = -PIVMIN
            WORK( INDLPL+I ) = LD( I ) / DPLUS
            WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
            IF( WORK( INDLPL+I ).EQ.ZERO )
     $                      WORK( INDS+I ) = LLD( I )
            S = WORK( INDS+I ) - LAMBDA
 71      CONTINUE
      END IF
C
C     Compute the progressive transform (using the differential form)
C     until the index R1
C
      SAWNAN2 = .FALSE.
      NEG2 = 0
      WORK( INDP+BN-1 ) = D( BN ) - LAMBDA
      DO 80 I = BN - 1, R1, -1
         DMINUS = LLD( I ) + WORK( INDP+I )
         TMP = D( I ) / DMINUS
         IF(DMINUS.LT.ZERO) NEG2 = NEG2 + 1
         WORK( INDUMN+I ) = L( I )*TMP
         WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - LAMBDA
 80   CONTINUE
      TMP = WORK( INDP+R1-1 )
      SAWNAN2 = DISNAN( TMP )

      IF( SAWNAN2 ) THEN
C        Runs a slower version of the above loop if a NaN is detected
         NEG2 = 0
         DO 100 I = BN-1, R1, -1
            DMINUS = LLD( I ) + WORK( INDP+I )
            IF(ABS(DMINUS).LT.PIVMIN) DMINUS = -PIVMIN
            TMP = D( I ) / DMINUS
            IF(DMINUS.LT.ZERO) NEG2 = NEG2 + 1
            WORK( INDUMN+I ) = L( I )*TMP
            WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - LAMBDA
            IF( TMP.EQ.ZERO )
     $          WORK( INDP+I-1 ) = D( I ) - LAMBDA
 100     CONTINUE
      END IF
C
C     Find the index (from R1 to R2) of the largest (in magnitude)
C     diagonal element of the inverse
C
      MINGMA = WORK( INDS+R1-1 ) + WORK( INDP+R1-1 )
      IF( MINGMA.LT.ZERO ) NEG1 = NEG1 + 1
      IF( WANTNC ) THEN
         NEGCNT = NEG1 + NEG2
      ELSE
         NEGCNT = -1
      ENDIF
      IF( ABS(MINGMA).EQ.ZERO )
     $   MINGMA = EPS*WORK( INDS+R1-1 )
      R = R1
      DO 110 I = R1, R2 - 1
         TMP = WORK( INDS+I ) + WORK( INDP+I )
         IF( TMP.EQ.ZERO )
     $      TMP = EPS*WORK( INDS+I )
         IF( ABS( TMP ).LE.ABS( MINGMA ) ) THEN
            MINGMA = TMP
            R = I + 1
         END IF
 110  CONTINUE
C
C     Compute the FP vector: solve N^T v = e_r
C
      ISUPPZ( 1 ) = B1
      ISUPPZ( 2 ) = BN
      Z( R ) = CONE
      ZTZ = ONE
C
C     Compute the FP vector upwards from R
C
      IF( .NOT.SAWNAN1 .AND. .NOT.SAWNAN2 ) THEN
         DO 210 I = R-1, B1, -1
            Z( I ) = -( WORK( INDLPL+I )*Z( I+1 ) )
            IF( (ABS(Z(I))+ABS(Z(I+1)))* ABS(LD(I)).LT.GAPTOL )
     $           THEN
               Z( I ) = ZERO
               ISUPPZ( 1 ) = I + 1
               GOTO 220
            ENDIF
            ZTZ = ZTZ + DBLE( Z( I )*Z( I ) )
 210     CONTINUE
 220     CONTINUE
      ELSE
C        Run slower loop if NaN occurred.
         DO 230 I = R - 1, B1, -1
            IF( Z( I+1 ).EQ.ZERO ) THEN
               Z( I ) = -( LD( I+1 ) / LD( I ) )*Z( I+2 )
            ELSE
               Z( I ) = -( WORK( INDLPL+I )*Z( I+1 ) )
            END IF
            IF( (ABS(Z(I))+ABS(Z(I+1)))* ABS(LD(I)).LT.GAPTOL )
     $           THEN
               Z( I ) = ZERO
               ISUPPZ( 1 ) = I + 1
               GO TO 240
            END IF
            ZTZ = ZTZ + DBLE( Z( I )*Z( I ) )
 230     CONTINUE
 240     CONTINUE
      ENDIF

C     Compute the FP vector downwards from R in blocks of size BLKSIZ
      IF( .NOT.SAWNAN1 .AND. .NOT.SAWNAN2 ) THEN
         DO 250 I = R, BN-1
            Z( I+1 ) = -( WORK( INDUMN+I )*Z( I ) )
            IF( (ABS(Z(I))+ABS(Z(I+1)))* ABS(LD(I)).LT.GAPTOL )
     $         THEN
               Z( I+1 ) = ZERO
               ISUPPZ( 2 ) = I
               GO TO 260
            END IF
            ZTZ = ZTZ + DBLE( Z( I+1 )*Z( I+1 ) )
 250     CONTINUE
 260     CONTINUE
      ELSE
C        Run slower loop if NaN occurred.
         DO 270 I = R, BN - 1
            IF( Z( I ).EQ.ZERO ) THEN
               Z( I+1 ) = -( LD( I-1 ) / LD( I ) )*Z( I-1 )
            ELSE
               Z( I+1 ) = -( WORK( INDUMN+I )*Z( I ) )
            END IF
            IF( (ABS(Z(I))+ABS(Z(I+1)))* ABS(LD(I)).LT.GAPTOL )
     $           THEN
               Z( I+1 ) = ZERO
               ISUPPZ( 2 ) = I
               GO TO 280
            END IF
            ZTZ = ZTZ + DBLE( Z( I+1 )*Z( I+1 ) )
 270     CONTINUE
 280     CONTINUE
      END IF
C
C     Compute quantities for convergence test
C
      TMP = ONE / ZTZ
      NRMINV = SQRT( TMP )
      RESID = ABS( MINGMA )*NRMINV
      RQCORR = MINGMA*TMP
C
C
      RETURN
C
C     End of ZLAR1V
C
      END

      SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      COMPLEX*16         ALPHA, BETA
C     ..
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
C     ..
C
C  =====================================================================
C
C     .. Local Scalars ..
      INTEGER            I, J
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
C
      IF( LSAME( UPLO, 'U' ) ) THEN
C
C        Set the diagonal to BETA and the strictly upper triangular
C        part of the array to ALPHA.
C
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
         DO 30 I = 1, MIN( N, M )
            A( I, I ) = BETA
   30    CONTINUE
C
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
C
C        Set the diagonal to BETA and the strictly lower triangular
C        part of the array to ALPHA.
C
         DO 50 J = 1, MIN( M, N )
            DO 40 I = J + 1, M
               A( I, J ) = ALPHA
   40       CONTINUE
   50    CONTINUE
         DO 60 I = 1, MIN( N, M )
            A( I, I ) = BETA
   60    CONTINUE
C
      ELSE
C
C        Set the array to BETA on the diagonal and ALPHA on the
C        offdiagonal.
C
         DO 80 J = 1, N
            DO 70 I = 1, M
               A( I, J ) = ALPHA
   70       CONTINUE
   80    CONTINUE
         DO 90 I = 1, MIN( M, N )
            A( I, I ) = BETA
   90    CONTINUE
      END IF
C
      RETURN
C
C     End of ZLASET
C
      END

      SUBROUTINE DLARRB( N, D, LLD, IFIRST, ILAST, RTOL1,
     $                   RTOL2, OFFSET, W, WGAP, WERR, WORK, IWORK,
     $                   PIVMIN, SPDIAM, TWIST, INFO )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2017
C
C     .. Scalar Arguments ..
      INTEGER            IFIRST, ILAST, INFO, N, OFFSET, TWIST
      DOUBLE PRECISION   PIVMIN, RTOL1, RTOL2, SPDIAM
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), LLD( * ), W( * ),
     $                   WERR( * ), WGAP( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, TWO, HALF
      PARAMETER        ( ZERO = 0.0D0, TWO = 2.0D0,
     $                   HALF = 0.5D0 )
      INTEGER   MAXITR
C     ..
C     .. Local Scalars ..
      INTEGER            I, I1, II, IP, ITER, K, NEGCNT, NEXT, NINT,
     $                   OLNINT, PREV, R
      DOUBLE PRECISION   BACK, CVRGD, GAP, LEFT, LGAP, MID, MNWDTH,
     $                   RGAP, RIGHT, TMP, WIDTH
C     ..
C     .. External Functions ..
      INTEGER            DLANEG
      EXTERNAL           DLANEG
C
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     ..
C     .. Executable Statements ..
C
      INFO = 0
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         RETURN
      END IF
C
      MAXITR = INT( ( LOG( SPDIAM+PIVMIN )-LOG( PIVMIN ) ) /
     $           LOG( TWO ) ) + 2
      MNWDTH = TWO * PIVMIN
C
      R = TWIST
      IF((R.LT.1).OR.(R.GT.N)) R = N
C
C     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
C     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while
C     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 )
C     for an unconverged interval is set to the index of the next unconverged
C     interval, and is -1 or 0 for a converged interval. Thus a linked
C     list of unconverged intervals is set up.
C
      I1 = IFIRST
C     The number of unconverged intervals
      NINT = 0
C     The last unconverged interval found
      PREV = 0

      RGAP = WGAP( I1-OFFSET )
      DO 75 I = I1, ILAST
         K = 2*I
         II = I - OFFSET
         LEFT = W( II ) - WERR( II )
         RIGHT = W( II ) + WERR( II )
         LGAP = RGAP
         RGAP = WGAP( II )
         GAP = MIN( LGAP, RGAP )

C        Make sure that [LEFT,RIGHT] contains the desired eigenvalue
C        Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT
C
C        Do while( NEGCNT(LEFT).GT.I-1 )
C
         BACK = WERR( II )
 20      CONTINUE
         NEGCNT = DLANEG( N, D, LLD, LEFT, PIVMIN, R )
         IF( NEGCNT.GT.I-1 ) THEN
            LEFT = LEFT - BACK
            BACK = TWO*BACK
            GO TO 20
         END IF
C
C        Do while( NEGCNT(RIGHT).LT.I )
C        Compute negcount from dstqds facto L+D+L+^T = L D L^T - RIGHT
C
         BACK = WERR( II )
 50      CONTINUE

         NEGCNT = DLANEG( N, D, LLD, RIGHT, PIVMIN, R )
          IF( NEGCNT.LT.I ) THEN
             RIGHT = RIGHT + BACK
             BACK = TWO*BACK
             GO TO 50
          END IF
         WIDTH = HALF*ABS( LEFT - RIGHT )
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )
         CVRGD = MAX(RTOL1*GAP,RTOL2*TMP)
         IF( WIDTH.LE.CVRGD .OR. WIDTH.LE.MNWDTH ) THEN
C           This interval has already converged and does not need refinement.
C           (Note that the gaps might change through refining the
C            eigenvalues, however, they can only get bigger.)
C           Remove it from the list.
            IWORK( K-1 ) = -1
C           Make sure that I1 always points to the first unconverged interval
            IF((I.EQ.I1).AND.(I.LT.ILAST)) I1 = I + 1
            IF((PREV.GE.I1).AND.(I.LE.ILAST)) IWORK( 2*PREV-1 ) = I + 1
         ELSE
C           unconverged interval found
            PREV = I
            NINT = NINT + 1
            IWORK( K-1 ) = I + 1
            IWORK( K ) = NEGCNT
         END IF
         WORK( K-1 ) = LEFT
         WORK( K ) = RIGHT
 75   CONTINUE

C
C     Do while( NINT.GT.0 ), i.e. there are still unconverged intervals
C     and while (ITER.LT.MAXITR)
C
      ITER = 0
 80   CONTINUE
      PREV = I1 - 1
      I = I1
      OLNINT = NINT

      DO 100 IP = 1, OLNINT
         K = 2*I
         II = I - OFFSET
         RGAP = WGAP( II )
         LGAP = RGAP
         IF(II.GT.1) LGAP = WGAP( II-1 )
         GAP = MIN( LGAP, RGAP )
         NEXT = IWORK( K-1 )
         LEFT = WORK( K-1 )
         RIGHT = WORK( K )
         MID = HALF*( LEFT + RIGHT )

C        semiwidth of interval
         WIDTH = RIGHT - MID
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )
         CVRGD = MAX(RTOL1*GAP,RTOL2*TMP)
         IF( ( WIDTH.LE.CVRGD ) .OR. ( WIDTH.LE.MNWDTH ).OR.
     $       ( ITER.EQ.MAXITR ) )THEN
C           reduce number of unconverged intervals
            NINT = NINT - 1
C           Mark interval as converged.
            IWORK( K-1 ) = 0
            IF( I1.EQ.I ) THEN
               I1 = NEXT
            ELSE
C              Prev holds the last unconverged interval previously examined
               IF(PREV.GE.I1) IWORK( 2*PREV-1 ) = NEXT
            END IF
            I = NEXT
            GO TO 100
         END IF
         PREV = I
C
C        Perform one bisection step
C
         NEGCNT = DLANEG( N, D, LLD, MID, PIVMIN, R )
         IF( NEGCNT.LE.I-1 ) THEN
            WORK( K-1 ) = MID
         ELSE
            WORK( K ) = MID
         END IF
         I = NEXT
 100  CONTINUE
      ITER = ITER + 1
C     do another loop if there are still unconverged intervals
C     However, in the last iteration, all intervals are accepted
C     since this is the best we can do.
      IF( ( NINT.GT.0 ).AND.(ITER.LE.MAXITR) ) GO TO 80
C
C
C     At this point, all the intervals have converged
      DO 110 I = IFIRST, ILAST
         K = 2*I
         II = I - OFFSET
C        All intervals marked by '0' have been refined.
         IF( IWORK( K-1 ).EQ.0 ) THEN
            W( II ) = HALF*( WORK( K-1 )+WORK( K ) )
            WERR( II ) = WORK( K ) - W( II )
         END IF
 110  CONTINUE
C
      DO 111 I = IFIRST+1, ILAST
         K = 2*I
         II = I - OFFSET
         WGAP( II-1 ) = MAX( ZERO,
     $                     W(II) - WERR (II) - W( II-1 ) - WERR( II-1 ))
 111  CONTINUE

      RETURN
C
C     End of DLARRB
C
      END

      SUBROUTINE DLARRC( JOBT, N, VL, VU, D, E, PIVMIN,
     $                            EIGCNT, LCNT, RCNT, INFO )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      CHARACTER          JOBT
      INTEGER            EIGCNT, INFO, LCNT, N, RCNT
      DOUBLE PRECISION   PIVMIN, VL, VU
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I
      LOGICAL            MATT
      DOUBLE PRECISION   LPIVOT, RPIVOT, SL, SU, TMP, TMP2

C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. Executable Statements ..
C
      INFO = 0
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         RETURN
      END IF
C
      LCNT = 0
      RCNT = 0
      EIGCNT = 0
      MATT = LSAME( JOBT, 'T' )


      IF (MATT) THEN
C        Sturm sequence count on T
         LPIVOT = D( 1 ) - VL
         RPIVOT = D( 1 ) - VU
         IF( LPIVOT.LE.ZERO ) THEN
            LCNT = LCNT + 1
         ENDIF
         IF( RPIVOT.LE.ZERO ) THEN
            RCNT = RCNT + 1
         ENDIF
         DO 10 I = 1, N-1
            TMP = E(I)**2
            LPIVOT = ( D( I+1 )-VL ) - TMP/LPIVOT
            RPIVOT = ( D( I+1 )-VU ) - TMP/RPIVOT
            IF( LPIVOT.LE.ZERO ) THEN
               LCNT = LCNT + 1
            ENDIF
            IF( RPIVOT.LE.ZERO ) THEN
               RCNT = RCNT + 1
            ENDIF
 10      CONTINUE
      ELSE
C        Sturm sequence count on L D L^T
         SL = -VL
         SU = -VU
         DO 20 I = 1, N - 1
            LPIVOT = D( I ) + SL
            RPIVOT = D( I ) + SU
            IF( LPIVOT.LE.ZERO ) THEN
               LCNT = LCNT + 1
            ENDIF
            IF( RPIVOT.LE.ZERO ) THEN
               RCNT = RCNT + 1
            ENDIF
            TMP = E(I) * D(I) * E(I)
C
            TMP2 = TMP / LPIVOT
            IF( TMP2.EQ.ZERO ) THEN
               SL =  TMP - VL
            ELSE
               SL = SL*TMP2 - VL
            END IF
C
            TMP2 = TMP / RPIVOT
            IF( TMP2.EQ.ZERO ) THEN
               SU =  TMP - VU
            ELSE
               SU = SU*TMP2 - VU
            END IF
 20      CONTINUE
         LPIVOT = D( N ) + SL
         RPIVOT = D( N ) + SU
         IF( LPIVOT.LE.ZERO ) THEN
            LCNT = LCNT + 1
         ENDIF
         IF( RPIVOT.LE.ZERO ) THEN
            RCNT = RCNT + 1
         ENDIF
      ENDIF
      EIGCNT = RCNT - LCNT

      RETURN
C
C     end of DLARRC
C
      END

      SUBROUTINE DLARRE( RANGE, N, VL, VU, IL, IU, D, E, E2,
     $                    RTOL1, RTOL2, SPLTOL, NSPLIT, ISPLIT, M,
     $                    W, WERR, WGAP, IBLOCK, INDEXW, GERS, PIVMIN,
     $                    WORK, IWORK, INFO )
C
C  -- LAPACK auxiliary routine (version 3.8.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      CHARACTER          RANGE
      INTEGER            IL, INFO, IU, M, N, NSPLIT
      DOUBLE PRECISION  PIVMIN, RTOL1, RTOL2, SPLTOL, VL, VU
C     ..
C     .. Array Arguments ..
      INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * ),
     $                   INDEXW( * )
      DOUBLE PRECISION   D( * ), E( * ), E2( * ), GERS( * ),
     $                   W( * ),WERR( * ), WGAP( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   FAC, FOUR, FOURTH, FUDGE, HALF, HNDRD,
     $                   MAXGROWTH, ONE, PERT, TWO, ZERO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0,
     $                     TWO = 2.0D0, FOUR=4.0D0,
     $                     HNDRD = 100.0D0,
     $                     PERT = 8.0D0,
     $                     HALF = ONE/TWO, FOURTH = ONE/FOUR, FAC= HALF,
     $                     MAXGROWTH = 64.0D0, FUDGE = 2.0D0 )
      INTEGER            MAXTRY, ALLRNG, INDRNG, VALRNG
      PARAMETER          ( MAXTRY = 6, ALLRNG = 1, INDRNG = 2,
     $                     VALRNG = 3 )
C     ..
C     .. Local Scalars ..
      LOGICAL            FORCEB, NOREP, USEDQD
      INTEGER            CNT, CNT1, CNT2, I, IBEGIN, IDUM, IEND, IINFO,
     $                   IN, INDL, INDU, IRANGE, J, JBLK, MB, MM,
     $                   WBEGIN, WEND
      DOUBLE PRECISION   AVGAP, BSRTOL, CLWDTH, DMAX, DPIVOT, EABS,
     $                   EMAX, EOLD, EPS, GL, GU, ISLEFT, ISRGHT, RTL,
     $                   RTOL, S1, S2, SAFMIN, SGNDEF, SIGMA, SPDIAM,
     $                   TAU, TMP, TMP1


C     ..
C     .. Local Arrays ..
      INTEGER            ISEED( 4 )
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION            DLAMCH
      EXTERNAL           DLAMCH, LSAME

C     ..
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DLARNV, DLARRA, DLARRB, DLARRC, DLARRD,
     $                   DLASQ2, DLARRK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN

C     ..
C     .. Executable Statements ..
C

      INFO = 0
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         RETURN
      END IF
C
C     Decode RANGE
C
      IF( LSAME( RANGE, 'A' ) ) THEN
         IRANGE = ALLRNG
      ELSE IF( LSAME( RANGE, 'V' ) ) THEN
         IRANGE = VALRNG
      ELSE IF( LSAME( RANGE, 'I' ) ) THEN
         IRANGE = INDRNG
      END IF

      M = 0

C     Get machine constants
      SAFMIN = DLAMCH( 'S' )
      EPS = DLAMCH( 'P' )

C     Set parameters
      RTL = SQRT(EPS)
      BSRTOL = SQRT(EPS)

C     Treat case of 1x1 matrix for quick return
      IF( N.EQ.1 ) THEN
         IF( (IRANGE.EQ.ALLRNG).OR.
     $       ((IRANGE.EQ.VALRNG).AND.(D(1).GT.VL).AND.(D(1).LE.VU)).OR.
     $       ((IRANGE.EQ.INDRNG).AND.(IL.EQ.1).AND.(IU.EQ.1)) ) THEN
            M = 1
            W(1) = D(1)
C           The computation error of the eigenvalue is zero
            WERR(1) = ZERO
            WGAP(1) = ZERO
            IBLOCK( 1 ) = 1
            INDEXW( 1 ) = 1
            GERS(1) = D( 1 )
            GERS(2) = D( 1 )
         ENDIF
C        store the shift for the initial RRR, which is zero in this case
         E(1) = ZERO
         RETURN
      END IF

C     General case: tridiagonal matrix of order > 1
C
C     Init WERR, WGAP. Compute Gerschgorin intervals and spectral diameter.
C     Compute maximum off-diagonal entry and pivmin.
      GL = D(1)
      GU = D(1)
      EOLD = ZERO
      EMAX = ZERO
      E(N) = ZERO
      DO 5 I = 1,N
         WERR(I) = ZERO
         WGAP(I) = ZERO
         EABS = ABS( E(I) )
         IF( EABS .GE. EMAX ) THEN
            EMAX = EABS
         END IF
         TMP1 = EABS + EOLD
         GERS( 2*I-1) = D(I) - TMP1
         GL =  MIN( GL, GERS( 2*I - 1))
         GERS( 2*I ) = D(I) + TMP1
         GU = MAX( GU, GERS(2*I) )
         EOLD  = EABS
 5    CONTINUE
C     The minimum pivot allowed in the Sturm sequence for T
      PIVMIN = SAFMIN * MAX( ONE, EMAX**2 )
C     Compute spectral diameter. The Gerschgorin bounds give an
C     estimate that is wrong by at most a factor of SQRT(2)
      SPDIAM = GU - GL

C     Compute splitting points
      CALL DLARRA( N, D, E, E2, SPLTOL, SPDIAM,
     $                    NSPLIT, ISPLIT, IINFO )

C     Can force use of bisection instead of faster DQDS.
C     Option left in the code for future multisection work.
      FORCEB = .FALSE.

C     Initialize USEDQD, DQDS should be used for ALLRNG unless someone
C     explicitly wants bisection.
      USEDQD = (( IRANGE.EQ.ALLRNG ) .AND. (.NOT.FORCEB))

      IF( (IRANGE.EQ.ALLRNG) .AND. (.NOT. FORCEB) ) THEN
C        Set interval [VL,VU] that contains all eigenvalues
         VL = GL
         VU = GU
      ELSE
C        We call DLARRD to find crude approximations to the eigenvalues
C        in the desired range. In case IRANGE = INDRNG, we also obtain the
C        interval (VL,VU] that contains all the wanted eigenvalues.
C        An interval [LEFT,RIGHT] has converged if
C        RIGHT-LEFT.LT.RTOL*MAX(ABS(LEFT),ABS(RIGHT))
C        DLARRD needs a WORK of size 4*N, IWORK of size 3*N
         CALL DLARRD( RANGE, 'B', N, VL, VU, IL, IU, GERS,
     $                    BSRTOL, D, E, E2, PIVMIN, NSPLIT, ISPLIT,
     $                    MM, W, WERR, VL, VU, IBLOCK, INDEXW,
     $                    WORK, IWORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = -1
            RETURN
         ENDIF
C        Make sure that the entries M+1 to N in W, WERR, IBLOCK, INDEXW are 0
         DO 14 I = MM+1,N
            W( I ) = ZERO
            WERR( I ) = ZERO
            IBLOCK( I ) = 0
            INDEXW( I ) = 0
 14      CONTINUE
      END IF


C**
C     Loop over unreduced blocks
      IBEGIN = 1
      WBEGIN = 1
      DO 170 JBLK = 1, NSPLIT
         IEND = ISPLIT( JBLK )
         IN = IEND - IBEGIN + 1

C        1 X 1 block
         IF( IN.EQ.1 ) THEN
            IF( (IRANGE.EQ.ALLRNG).OR.( (IRANGE.EQ.VALRNG).AND.
     $         ( D( IBEGIN ).GT.VL ).AND.( D( IBEGIN ).LE.VU ) )
     $        .OR. ( (IRANGE.EQ.INDRNG).AND.(IBLOCK(WBEGIN).EQ.JBLK))
     $        ) THEN
               M = M + 1
               W( M ) = D( IBEGIN )
               WERR(M) = ZERO
C              The gap for a single block doesn't matter for the later
C              algorithm and is assigned an arbitrary large value
               WGAP(M) = ZERO
               IBLOCK( M ) = JBLK
               INDEXW( M ) = 1
               WBEGIN = WBEGIN + 1
            ENDIF
C           E( IEND ) holds the shift for the initial RRR
            E( IEND ) = ZERO
            IBEGIN = IEND + 1
            GO TO 170
         END IF
C
C        Blocks of size larger than 1x1
C
C        E( IEND ) will hold the shift for the initial RRR, for now set it =0
         E( IEND ) = ZERO
C
C        Find local outer bounds GL,GU for the block
         GL = D(IBEGIN)
         GU = D(IBEGIN)
         DO 15 I = IBEGIN , IEND
            GL = MIN( GERS( 2*I-1 ), GL )
            GU = MAX( GERS( 2*I ), GU )
 15      CONTINUE
         SPDIAM = GU - GL

         IF(.NOT. ((IRANGE.EQ.ALLRNG).AND.(.NOT.FORCEB)) ) THEN
C           Count the number of eigenvalues in the current block.
            MB = 0
            DO 20 I = WBEGIN,MM
               IF( IBLOCK(I).EQ.JBLK ) THEN
                  MB = MB+1
               ELSE
                  GOTO 21
               ENDIF
 20         CONTINUE
 21         CONTINUE

            IF( MB.EQ.0) THEN
C              No eigenvalue in the current block lies in the desired range
C              E( IEND ) holds the shift for the initial RRR
               E( IEND ) = ZERO
               IBEGIN = IEND + 1
               GO TO 170
            ELSE

C              Decide whether dqds or bisection is more efficient
               USEDQD = ( (MB .GT. FAC*IN) .AND. (.NOT.FORCEB) )
               WEND = WBEGIN + MB - 1
C              Calculate gaps for the current block
C              In later stages, when representations for individual
C              eigenvalues are different, we use SIGMA = E( IEND ).
               SIGMA = ZERO
               DO 30 I = WBEGIN, WEND - 1
                  WGAP( I ) = MAX( ZERO,
     $                        W(I+1)-WERR(I+1) - (W(I)+WERR(I)) )
 30            CONTINUE
               WGAP( WEND ) = MAX( ZERO,
     $                     VU - SIGMA - (W( WEND )+WERR( WEND )))
C              Find local index of the first and last desired evalue.
               INDL = INDEXW(WBEGIN)
               INDU = INDEXW( WEND )
            ENDIF
         ENDIF
         IF(( (IRANGE.EQ.ALLRNG) .AND. (.NOT. FORCEB) ).OR.USEDQD) THEN
C           Case of DQDS
C           Find approximations to the extremal eigenvalues of the block
            CALL DLARRK( IN, 1, GL, GU, D(IBEGIN),
     $               E2(IBEGIN), PIVMIN, RTL, TMP, TMP1, IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = -1
               RETURN
            ENDIF
            ISLEFT = MAX(GL, TMP - TMP1
     $               - HNDRD * EPS* ABS(TMP - TMP1))

            CALL DLARRK( IN, IN, GL, GU, D(IBEGIN),
     $               E2(IBEGIN), PIVMIN, RTL, TMP, TMP1, IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = -1
               RETURN
            ENDIF
            ISRGHT = MIN(GU, TMP + TMP1
     $                 + HNDRD * EPS * ABS(TMP + TMP1))
C           Improve the estimate of the spectral diameter
            SPDIAM = ISRGHT - ISLEFT
         ELSE
C           Case of bisection
C           Find approximations to the wanted extremal eigenvalues
            ISLEFT = MAX(GL, W(WBEGIN) - WERR(WBEGIN)
     $                  - HNDRD * EPS*ABS(W(WBEGIN)- WERR(WBEGIN) ))
            ISRGHT = MIN(GU,W(WEND) + WERR(WEND)
     $                  + HNDRD * EPS * ABS(W(WEND)+ WERR(WEND)))
         ENDIF


C        Decide whether the base representation for the current block
C        L_JBLK D_JBLK L_JBLK^T = T_JBLK - sigma_JBLK I
C        should be on the left or the right end of the current block.
C        The strategy is to shift to the end which is "more populated"
C        Furthermore, decide whether to use DQDS for the computation of
C        the eigenvalue approximations at the end of DLARRE or bisection.
C        dqds is chosen if all eigenvalues are desired or the number of
C        eigenvalues to be computed is large compared to the blocksize.
         IF( ( IRANGE.EQ.ALLRNG ) .AND. (.NOT.FORCEB) ) THEN
C           If all the eigenvalues have to be computed, we use dqd
            USEDQD = .TRUE.
C           INDL is the local index of the first eigenvalue to compute
            INDL = 1
            INDU = IN
C           MB =  number of eigenvalues to compute
            MB = IN
            WEND = WBEGIN + MB - 1
C           Define 1/4 and 3/4 points of the spectrum
            S1 = ISLEFT + FOURTH * SPDIAM
            S2 = ISRGHT - FOURTH * SPDIAM
         ELSE
C           DLARRD has computed IBLOCK and INDEXW for each eigenvalue
C           approximation.
C           choose sigma
            IF( USEDQD ) THEN
               S1 = ISLEFT + FOURTH * SPDIAM
               S2 = ISRGHT - FOURTH * SPDIAM
            ELSE
               TMP = MIN(ISRGHT,VU) -  MAX(ISLEFT,VL)
               S1 =  MAX(ISLEFT,VL) + FOURTH * TMP
               S2 =  MIN(ISRGHT,VU) - FOURTH * TMP
            ENDIF
         ENDIF

C        Compute the negcount at the 1/4 and 3/4 points
         IF(MB.GT.1) THEN
            CALL DLARRC( 'T', IN, S1, S2, D(IBEGIN),
     $                    E(IBEGIN), PIVMIN, CNT, CNT1, CNT2, IINFO)
         ENDIF

         IF(MB.EQ.1) THEN
            SIGMA = GL
            SGNDEF = ONE
         ELSEIF( CNT1 - INDL .GE. INDU - CNT2 ) THEN
            IF( ( IRANGE.EQ.ALLRNG ) .AND. (.NOT.FORCEB) ) THEN
               SIGMA = MAX(ISLEFT,GL)
            ELSEIF( USEDQD ) THEN
C              use Gerschgorin bound as shift to get pos def matrix
C              for dqds
               SIGMA = ISLEFT
            ELSE
C              use approximation of the first desired eigenvalue of the
C              block as shift
               SIGMA = MAX(ISLEFT,VL)
            ENDIF
            SGNDEF = ONE
         ELSE
            IF( ( IRANGE.EQ.ALLRNG ) .AND. (.NOT.FORCEB) ) THEN
               SIGMA = MIN(ISRGHT,GU)
            ELSEIF( USEDQD ) THEN
C              use Gerschgorin bound as shift to get neg def matrix
C              for dqds
               SIGMA = ISRGHT
            ELSE
C              use approximation of the first desired eigenvalue of the
C              block as shift
               SIGMA = MIN(ISRGHT,VU)
            ENDIF
            SGNDEF = -ONE
         ENDIF


C        An initial SIGMA has been chosen that will be used for computing
C        T - SIGMA I = L D L^T
C        Define the increment TAU of the shift in case the initial shift
C        needs to be refined to obtain a factorization with not too much
C        element growth.
         IF( USEDQD ) THEN
C           The initial SIGMA was to the outer end of the spectrum
C           the matrix is definite and we need not retreat.
            TAU = SPDIAM*EPS*N + TWO*PIVMIN
            TAU = MAX( TAU,TWO*EPS*ABS(SIGMA) )
         ELSE
            IF(MB.GT.1) THEN
               CLWDTH = W(WEND) + WERR(WEND) - W(WBEGIN) - WERR(WBEGIN)
               AVGAP = ABS(CLWDTH / DBLE(WEND-WBEGIN))
               IF( SGNDEF.EQ.ONE ) THEN
                  TAU = HALF*MAX(WGAP(WBEGIN),AVGAP)
                  TAU = MAX(TAU,WERR(WBEGIN))
               ELSE
                  TAU = HALF*MAX(WGAP(WEND-1),AVGAP)
                  TAU = MAX(TAU,WERR(WEND))
               ENDIF
            ELSE
               TAU = WERR(WBEGIN)
            ENDIF
         ENDIF
C
         DO 80 IDUM = 1, MAXTRY
C           Compute L D L^T factorization of tridiagonal matrix T - sigma I.
C           Store D in WORK(1:IN), L in WORK(IN+1:2*IN), and reciprocals of
C           pivots in WORK(2*IN+1:3*IN)
            DPIVOT = D( IBEGIN ) - SIGMA
            WORK( 1 ) = DPIVOT
            DMAX = ABS( WORK(1) )
            J = IBEGIN
            DO 70 I = 1, IN - 1
               WORK( 2*IN+I ) = ONE / WORK( I )
               TMP = E( J )*WORK( 2*IN+I )
               WORK( IN+I ) = TMP
               DPIVOT = ( D( J+1 )-SIGMA ) - TMP*E( J )
               WORK( I+1 ) = DPIVOT
               DMAX = MAX( DMAX, ABS(DPIVOT) )
               J = J + 1
 70         CONTINUE
C           check for element growth
            IF( DMAX .GT. MAXGROWTH*SPDIAM ) THEN
               NOREP = .TRUE.
            ELSE
               NOREP = .FALSE.
            ENDIF
            IF( USEDQD .AND. .NOT.NOREP ) THEN
C              Ensure the definiteness of the representation
C              All entries of D (of L D L^T) must have the same sign
               DO 71 I = 1, IN
                  TMP = SGNDEF*WORK( I )
                  IF( TMP.LT.ZERO ) NOREP = .TRUE.
 71            CONTINUE
            ENDIF
            IF(NOREP) THEN
C              Note that in the case of IRANGE=ALLRNG, we use the Gerschgorin
C              shift which makes the matrix definite. So we should end up
C              here really only in the case of IRANGE = VALRNG or INDRNG.
               IF( IDUM.EQ.MAXTRY-1 ) THEN
                  IF( SGNDEF.EQ.ONE ) THEN
C                    The fudged Gerschgorin shift should succeed
                     SIGMA =
     $                    GL - FUDGE*SPDIAM*EPS*N - FUDGE*TWO*PIVMIN
                  ELSE
                     SIGMA =
     $                    GU + FUDGE*SPDIAM*EPS*N + FUDGE*TWO*PIVMIN
                  END IF
               ELSE
                  SIGMA = SIGMA - SGNDEF * TAU
                  TAU = TWO * TAU
               END IF
            ELSE
C              an initial RRR is found
               GO TO 83
            END IF
 80      CONTINUE
C        if the program reaches this point, no base representation could be
C        found in MAXTRY iterations.
         INFO = 2
         RETURN

 83      CONTINUE
C        At this point, we have found an initial base representation
C        T - SIGMA I = L D L^T with not too much element growth.
C        Store the shift.
         E( IEND ) = SIGMA
C        Store D and L.
         CALL DCOPY( IN, WORK, 1, D( IBEGIN ), 1 )
         CALL DCOPY( IN-1, WORK( IN+1 ), 1, E( IBEGIN ), 1 )


         IF(MB.GT.1 ) THEN
C
C           Perturb each entry of the base representation by a small
C           (but random) relative amount to overcome difficulties with
C           glued matrices.
C
            DO 122 I = 1, 4
               ISEED( I ) = 1
 122        CONTINUE

            CALL DLARNV(2, ISEED, 2*IN-1, WORK(1))
            DO 125 I = 1,IN-1
               D(IBEGIN+I-1) = D(IBEGIN+I-1)*(ONE+EPS*PERT*WORK(I))
               E(IBEGIN+I-1) = E(IBEGIN+I-1)*(ONE+EPS*PERT*WORK(IN+I))
 125        CONTINUE
            D(IEND) = D(IEND)*(ONE+EPS*FOUR*WORK(IN))
C
         ENDIF
C
C        Don't update the Gerschgorin intervals because keeping track
C        of the updates would be too much work in DLARRV.
C        We update W instead and use it to locate the proper Gerschgorin
C        intervals.

C        Compute the required eigenvalues of L D L' by bisection or dqds
         IF ( .NOT.USEDQD ) THEN
C           If DLARRD has been used, shift the eigenvalue approximations
C           according to their representation. This is necessary for
C           a uniform DLARRV since dqds computes eigenvalues of the
C           shifted representation. In DLARRV, W will always hold the
C           UNshifted eigenvalue approximation.
            DO 134 J=WBEGIN,WEND
               W(J) = W(J) - SIGMA
               WERR(J) = WERR(J) + ABS(W(J)) * EPS
 134        CONTINUE
C           call DLARRB to reduce eigenvalue error of the approximations
C           from DLARRD
            DO 135 I = IBEGIN, IEND-1
               WORK( I ) = D( I ) * E( I )**2
 135        CONTINUE
C           use bisection to find EV from INDL to INDU
            CALL DLARRB(IN, D(IBEGIN), WORK(IBEGIN),
     $                  INDL, INDU, RTOL1, RTOL2, INDL-1,
     $                  W(WBEGIN), WGAP(WBEGIN), WERR(WBEGIN),
     $                  WORK( 2*N+1 ), IWORK, PIVMIN, SPDIAM,
     $                  IN, IINFO )
            IF( IINFO .NE. 0 ) THEN
               INFO = -4
               RETURN
            END IF
C           DLARRB computes all gaps correctly except for the last one
C           Record distance to VU/GU
            WGAP( WEND ) = MAX( ZERO,
     $           ( VU-SIGMA ) - ( W( WEND ) + WERR( WEND ) ) )
            DO 138 I = INDL, INDU
               M = M + 1
               IBLOCK(M) = JBLK
               INDEXW(M) = I
 138        CONTINUE
         ELSE
C           Call dqds to get all eigs (and then possibly delete unwanted
C           eigenvalues).
C           Note that dqds finds the eigenvalues of the L D L^T representation
C           of T to high relative accuracy. High relative accuracy
C           might be lost when the shift of the RRR is subtracted to obtain
C           the eigenvalues of T. However, T is not guaranteed to define its
C           eigenvalues to high relative accuracy anyway.
C           Set RTOL to the order of the tolerance used in DLASQ2
C           This is an ESTIMATED error, the worst case bound is 4*N*EPS
C           which is usually too large and requires unnecessary work to be
C           done by bisection when computing the eigenvectors
            RTOL = LOG(DBLE(IN)) * FOUR * EPS
            J = IBEGIN
            DO 140 I = 1, IN - 1
               WORK( 2*I-1 ) = ABS( D( J ) )
               WORK( 2*I ) = E( J )*E( J )*WORK( 2*I-1 )
               J = J + 1
  140       CONTINUE
            WORK( 2*IN-1 ) = ABS( D( IEND ) )
            WORK( 2*IN ) = ZERO
            CALL DLASQ2( IN, WORK, IINFO )
            IF( IINFO .NE. 0 ) THEN
C              If IINFO = -5 then an index is part of a tight cluster
C              and should be changed. The index is in IWORK(1) and the
C              gap is in WORK(N+1)
               INFO = -5
               RETURN
            ELSE
C              Test that all eigenvalues are positive as expected
               DO 149 I = 1, IN
                  IF( WORK( I ).LT.ZERO ) THEN
                     INFO = -6
                     RETURN
                  ENDIF
 149           CONTINUE
            END IF
            IF( SGNDEF.GT.ZERO ) THEN
               DO 150 I = INDL, INDU
                  M = M + 1
                  W( M ) = WORK( IN-I+1 )
                  IBLOCK( M ) = JBLK
                  INDEXW( M ) = I
 150           CONTINUE
            ELSE
               DO 160 I = INDL, INDU
                  M = M + 1
                  W( M ) = -WORK( I )
                  IBLOCK( M ) = JBLK
                  INDEXW( M ) = I
 160           CONTINUE
            END IF

            DO 165 I = M - MB + 1, M
C              the value of RTOL below should be the tolerance in DLASQ2
               WERR( I ) = RTOL * ABS( W(I) )
 165        CONTINUE
            DO 166 I = M - MB + 1, M - 1
C              compute the right gap between the intervals
               WGAP( I ) = MAX( ZERO,
     $                          W(I+1)-WERR(I+1) - (W(I)+WERR(I)) )
 166        CONTINUE
            WGAP( M ) = MAX( ZERO,
     $           ( VU-SIGMA ) - ( W( M ) + WERR( M ) ) )
         END IF
C        proceed with next block
         IBEGIN = IEND + 1
         WBEGIN = WEND + 1
 170  CONTINUE
C

      RETURN
C
C     end of DLARRE
C
      END

      SUBROUTINE DLARRF( N, D, L, LD, CLSTRT, CLEND,
     $                   W, WGAP, WERR,
     $                   SPDIAM, CLGAPL, CLGAPR, PIVMIN, SIGMA,
     $                   DPLUS, LPLUS, WORK, INFO )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      INTEGER            CLSTRT, CLEND, INFO, N
      DOUBLE PRECISION   CLGAPL, CLGAPR, PIVMIN, SIGMA, SPDIAM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DPLUS( * ), L( * ), LD( * ),
     $          LPLUS( * ), W( * ), WGAP( * ), WERR( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   FOUR, MAXGROWTH1, MAXGROWTH2, ONE, QUART, TWO
      PARAMETER          ( ONE = 1.0D0, TWO = 2.0D0, FOUR = 4.0D0,
     $                     QUART = 0.25D0,
     $                     MAXGROWTH1 = 8.D0,
     $                     MAXGROWTH2 = 8.D0 )
C     ..
C     .. Local Scalars ..
      LOGICAL   DORRR1, FORCER, NOFAIL, SAWNAN1, SAWNAN2, TRYRRR1
      INTEGER            I, INDX, KTRY, KTRYMAX, SLEFT, SRIGHT, SHIFT
      PARAMETER          ( KTRYMAX = 1, SLEFT = 1, SRIGHT = 2 )
      DOUBLE PRECISION   AVGAP, BESTSHIFT, CLWDTH, EPS, FACT, FAIL,
     $                   FAIL2, GROWTHBOUND, LDELTA, LDMAX, LSIGMA,
     $                   MAX1, MAX2, MINGAP, OLDP, PROD, RDELTA, RDMAX,
     $                   RRR1, RRR2, RSIGMA, S, SMLGROWTH, TMP, ZNM2
C     ..
C     .. External Functions ..
      LOGICAL DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DISNAN, DLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL           DCOPY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
C
      INFO = 0
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         RETURN
      END IF
C
      FACT = DBLE(2**KTRYMAX)
      EPS = DLAMCH( 'Precision' )
      SHIFT = 0
      FORCER = .FALSE.


C     Note that we cannot guarantee that for any of the shifts tried,
C     the factorization has a small or even moderate element growth.
C     There could be Ritz values at both ends of the cluster and despite
C     backing off, there are examples where all factorizations tried
C     (in IEEE mode, allowing zero pivots & infinities) have INFINITE
C     element growth.
C     For this reason, we should use PIVMIN in this subroutine so that at
C     least the L D L^T factorization exists. It can be checked afterwards
C     whether the element growth caused bad residuals/orthogonality.

C     Decide whether the code should accept the best among all
C     representations despite large element growth or signal INFO=1
C     Setting NOFAIL to .FALSE. for quick fix for bug 113
      NOFAIL = .FALSE.
C

C     Compute the average gap length of the cluster
      CLWDTH = ABS(W(CLEND)-W(CLSTRT)) + WERR(CLEND) + WERR(CLSTRT)
      AVGAP = CLWDTH / DBLE(CLEND-CLSTRT)
      MINGAP = MIN(CLGAPL, CLGAPR)
C     Initial values for shifts to both ends of cluster
      LSIGMA = MIN(W( CLSTRT ),W( CLEND )) - WERR( CLSTRT )
      RSIGMA = MAX(W( CLSTRT ),W( CLEND )) + WERR( CLEND )

C     Use a small fudge to make sure that we really shift to the outside
      LSIGMA = LSIGMA - ABS(LSIGMA)* FOUR * EPS
      RSIGMA = RSIGMA + ABS(RSIGMA)* FOUR * EPS

C     Compute upper bounds for how much to back off the initial shifts
      LDMAX = QUART * MINGAP + TWO * PIVMIN
      RDMAX = QUART * MINGAP + TWO * PIVMIN

      LDELTA = MAX(AVGAP,WGAP( CLSTRT ))/FACT
      RDELTA = MAX(AVGAP,WGAP( CLEND-1 ))/FACT
C
C     Initialize the record of the best representation found
C
      S = DLAMCH( 'S' )
      SMLGROWTH = ONE / S
      FAIL = DBLE(N-1)*MINGAP/(SPDIAM*EPS)
      FAIL2 = DBLE(N-1)*MINGAP/(SPDIAM*SQRT(EPS))
      BESTSHIFT = LSIGMA
C
C     while (KTRY <= KTRYMAX)
      KTRY = 0
      GROWTHBOUND = MAXGROWTH1*SPDIAM

 5    CONTINUE
      SAWNAN1 = .FALSE.
      SAWNAN2 = .FALSE.
C     Ensure that we do not back off too much of the initial shifts
      LDELTA = MIN(LDMAX,LDELTA)
      RDELTA = MIN(RDMAX,RDELTA)

C     Compute the element growth when shifting to both ends of the cluster
C     accept the shift if there is no element growth at one of the two ends

C     Left end
      S = -LSIGMA
      DPLUS( 1 ) = D( 1 ) + S
      IF(ABS(DPLUS(1)).LT.PIVMIN) THEN
         DPLUS(1) = -PIVMIN
C        Need to set SAWNAN1 because refined RRR test should not be used
C        in this case
         SAWNAN1 = .TRUE.
      ENDIF
      MAX1 = ABS( DPLUS( 1 ) )
      DO 6 I = 1, N - 1
         LPLUS( I ) = LD( I ) / DPLUS( I )
         S = S*LPLUS( I )*L( I ) - LSIGMA
         DPLUS( I+1 ) = D( I+1 ) + S
         IF(ABS(DPLUS(I+1)).LT.PIVMIN) THEN
            DPLUS(I+1) = -PIVMIN
C           Need to set SAWNAN1 because refined RRR test should not be used
C           in this case
            SAWNAN1 = .TRUE.
         ENDIF
         MAX1 = MAX( MAX1,ABS(DPLUS(I+1)) )
 6    CONTINUE
      SAWNAN1 = SAWNAN1 .OR.  DISNAN( MAX1 )

      IF( FORCER .OR.
     $   (MAX1.LE.GROWTHBOUND .AND. .NOT.SAWNAN1 ) ) THEN
         SIGMA = LSIGMA
         SHIFT = SLEFT
         GOTO 100
      ENDIF

C     Right end
      S = -RSIGMA
      WORK( 1 ) = D( 1 ) + S
      IF(ABS(WORK(1)).LT.PIVMIN) THEN
         WORK(1) = -PIVMIN
C        Need to set SAWNAN2 because refined RRR test should not be used
C        in this case
         SAWNAN2 = .TRUE.
      ENDIF
      MAX2 = ABS( WORK( 1 ) )
      DO 7 I = 1, N - 1
         WORK( N+I ) = LD( I ) / WORK( I )
         S = S*WORK( N+I )*L( I ) - RSIGMA
         WORK( I+1 ) = D( I+1 ) + S
         IF(ABS(WORK(I+1)).LT.PIVMIN) THEN
            WORK(I+1) = -PIVMIN
C           Need to set SAWNAN2 because refined RRR test should not be used
C           in this case
            SAWNAN2 = .TRUE.
         ENDIF
         MAX2 = MAX( MAX2,ABS(WORK(I+1)) )
 7    CONTINUE
      SAWNAN2 = SAWNAN2 .OR.  DISNAN( MAX2 )

      IF( FORCER .OR.
     $   (MAX2.LE.GROWTHBOUND .AND. .NOT.SAWNAN2 ) ) THEN
         SIGMA = RSIGMA
         SHIFT = SRIGHT
         GOTO 100
      ENDIF
C     If we are at this point, both shifts led to too much element growth

C     Record the better of the two shifts (provided it didn't lead to NaN)
      IF(SAWNAN1.AND.SAWNAN2) THEN
C        both MAX1 and MAX2 are NaN
         GOTO 50
      ELSE
         IF( .NOT.SAWNAN1 ) THEN
            INDX = 1
            IF(MAX1.LE.SMLGROWTH) THEN
               SMLGROWTH = MAX1
               BESTSHIFT = LSIGMA
            ENDIF
         ENDIF
         IF( .NOT.SAWNAN2 ) THEN
            IF(SAWNAN1 .OR. MAX2.LE.MAX1) INDX = 2
            IF(MAX2.LE.SMLGROWTH) THEN
               SMLGROWTH = MAX2
               BESTSHIFT = RSIGMA
            ENDIF
         ENDIF
      ENDIF

C     If we are here, both the left and the right shift led to
C     element growth. If the element growth is moderate, then
C     we may still accept the representation, if it passes a
C     refined test for RRR. This test supposes that no NaN occurred.
C     Moreover, we use the refined RRR test only for isolated clusters.
      IF((CLWDTH.LT.MINGAP/DBLE(128)) .AND.
     $   (MIN(MAX1,MAX2).LT.FAIL2)
     $  .AND.(.NOT.SAWNAN1).AND.(.NOT.SAWNAN2)) THEN
         DORRR1 = .TRUE.
      ELSE
         DORRR1 = .FALSE.
      ENDIF
      TRYRRR1 = .TRUE.
      IF( TRYRRR1 .AND. DORRR1 ) THEN
      IF(INDX.EQ.1) THEN
         TMP = ABS( DPLUS( N ) )
         ZNM2 = ONE
         PROD = ONE
         OLDP = ONE
         DO 15 I = N-1, 1, -1
            IF( PROD .LE. EPS ) THEN
               PROD =
     $         ((DPLUS(I+1)*WORK(N+I+1))/(DPLUS(I)*WORK(N+I)))*OLDP
            ELSE
               PROD = PROD*ABS(WORK(N+I))
            END IF
            OLDP = PROD
            ZNM2 = ZNM2 + PROD**2
            TMP = MAX( TMP, ABS( DPLUS( I ) * PROD ))
 15      CONTINUE
         RRR1 = TMP/( SPDIAM * SQRT( ZNM2 ) )
         IF (RRR1.LE.MAXGROWTH2) THEN
            SIGMA = LSIGMA
            SHIFT = SLEFT
            GOTO 100
         ENDIF
      ELSE IF(INDX.EQ.2) THEN
         TMP = ABS( WORK( N ) )
         ZNM2 = ONE
         PROD = ONE
         OLDP = ONE
         DO 16 I = N-1, 1, -1
            IF( PROD .LE. EPS ) THEN
               PROD = ((WORK(I+1)*LPLUS(I+1))/(WORK(I)*LPLUS(I)))*OLDP
            ELSE
               PROD = PROD*ABS(LPLUS(I))
            END IF
            OLDP = PROD
            ZNM2 = ZNM2 + PROD**2
            TMP = MAX( TMP, ABS( WORK( I ) * PROD ))
 16      CONTINUE
         RRR2 = TMP/( SPDIAM * SQRT( ZNM2 ) )
         IF (RRR2.LE.MAXGROWTH2) THEN
            SIGMA = RSIGMA
            SHIFT = SRIGHT
            GOTO 100
         ENDIF
      END IF
      ENDIF

 50   CONTINUE

      IF (KTRY.LT.KTRYMAX) THEN
C        If we are here, both shifts failed also the RRR test.
C        Back off to the outside
         LSIGMA = MAX( LSIGMA - LDELTA,
     $     LSIGMA - LDMAX)
         RSIGMA = MIN( RSIGMA + RDELTA,
     $     RSIGMA + RDMAX )
         LDELTA = TWO * LDELTA
         RDELTA = TWO * RDELTA
         KTRY = KTRY + 1
         GOTO 5
      ELSE
C        None of the representations investigated satisfied our
C        criteria. Take the best one we found.
         IF((SMLGROWTH.LT.FAIL).OR.NOFAIL) THEN
            LSIGMA = BESTSHIFT
            RSIGMA = BESTSHIFT
            FORCER = .TRUE.
            GOTO 5
         ELSE
            INFO = 1
            RETURN
         ENDIF
      END IF

 100  CONTINUE
      IF (SHIFT.EQ.SLEFT) THEN
      ELSEIF (SHIFT.EQ.SRIGHT) THEN
C        store new L and D back into DPLUS, LPLUS
         CALL DCOPY( N, WORK, 1, DPLUS, 1 )
         CALL DCOPY( N-1, WORK(N+1), 1, LPLUS, 1 )
      ENDIF

      RETURN
C
C     End of DLARRF
C
      END

      SUBROUTINE DLARRJ( N, D, E2, IFIRST, ILAST,
     $                   RTOL, OFFSET, W, WERR, WORK, IWORK,
     $                   PIVMIN, SPDIAM, INFO )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2017
C
C     .. Scalar Arguments ..
      INTEGER            IFIRST, ILAST, INFO, N, OFFSET
      DOUBLE PRECISION   PIVMIN, RTOL, SPDIAM
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E2( * ), W( * ),
     $                   WERR( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   HALF = 0.5D0 )
      INTEGER   MAXITR
C     ..
C     .. Local Scalars ..
      INTEGER            CNT, I, I1, I2, II, ITER, J, K, NEXT, NINT,
     $                   OLNINT, P, PREV, SAVI1
      DOUBLE PRECISION   DPLUS, FAC, LEFT, MID, RIGHT, S, TMP, WIDTH
C
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
C     ..
C     .. Executable Statements ..
C
      INFO = 0
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         RETURN
      END IF
C
      MAXITR = INT( ( LOG( SPDIAM+PIVMIN )-LOG( PIVMIN ) ) /
     $           LOG( TWO ) ) + 2
C
C     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
C     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while
C     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 )
C     for an unconverged interval is set to the index of the next unconverged
C     interval, and is -1 or 0 for a converged interval. Thus a linked
C     list of unconverged intervals is set up.
C

      I1 = IFIRST
      I2 = ILAST
C     The number of unconverged intervals
      NINT = 0
C     The last unconverged interval found
      PREV = 0
      DO 75 I = I1, I2
         K = 2*I
         II = I - OFFSET
         LEFT = W( II ) - WERR( II )
         MID = W(II)
         RIGHT = W( II ) + WERR( II )
         WIDTH = RIGHT - MID
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )

C        The following test prevents the test of converged intervals
         IF( WIDTH.LT.RTOL*TMP ) THEN
C           This interval has already converged and does not need refinement.
C           (Note that the gaps might change through refining the
C            eigenvalues, however, they can only get bigger.)
C           Remove it from the list.
            IWORK( K-1 ) = -1
C           Make sure that I1 always points to the first unconverged interval
            IF((I.EQ.I1).AND.(I.LT.I2)) I1 = I + 1
            IF((PREV.GE.I1).AND.(I.LE.I2)) IWORK( 2*PREV-1 ) = I + 1
         ELSE
C           unconverged interval found
            PREV = I
C           Make sure that [LEFT,RIGHT] contains the desired eigenvalue
C
C           Do while( CNT(LEFT).GT.I-1 )
C
            FAC = ONE
 20         CONTINUE
            CNT = 0
            S = LEFT
            DPLUS = D( 1 ) - S
            IF( DPLUS.LT.ZERO ) CNT = CNT + 1
            DO 30 J = 2, N
               DPLUS = D( J ) - S - E2( J-1 )/DPLUS
               IF( DPLUS.LT.ZERO ) CNT = CNT + 1
 30         CONTINUE
            IF( CNT.GT.I-1 ) THEN
               LEFT = LEFT - WERR( II )*FAC
               FAC = TWO*FAC
               GO TO 20
            END IF
C
C           Do while( CNT(RIGHT).LT.I )
C
            FAC = ONE
 50         CONTINUE
            CNT = 0
            S = RIGHT
            DPLUS = D( 1 ) - S
            IF( DPLUS.LT.ZERO ) CNT = CNT + 1
            DO 60 J = 2, N
               DPLUS = D( J ) - S - E2( J-1 )/DPLUS
               IF( DPLUS.LT.ZERO ) CNT = CNT + 1
 60         CONTINUE
            IF( CNT.LT.I ) THEN
               RIGHT = RIGHT + WERR( II )*FAC
               FAC = TWO*FAC
               GO TO 50
            END IF
            NINT = NINT + 1
            IWORK( K-1 ) = I + 1
            IWORK( K ) = CNT
         END IF
         WORK( K-1 ) = LEFT
         WORK( K ) = RIGHT
 75   CONTINUE


      SAVI1 = I1
C
C     Do while( NINT.GT.0 ), i.e. there are still unconverged intervals
C     and while (ITER.LT.MAXITR)
C
      ITER = 0
 80   CONTINUE
      PREV = I1 - 1
      I = I1
      OLNINT = NINT

      DO 100 P = 1, OLNINT
         K = 2*I
         II = I - OFFSET
         NEXT = IWORK( K-1 )
         LEFT = WORK( K-1 )
         RIGHT = WORK( K )
         MID = HALF*( LEFT + RIGHT )

C        semiwidth of interval
         WIDTH = RIGHT - MID
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )

         IF( ( WIDTH.LT.RTOL*TMP ) .OR.
     $      (ITER.EQ.MAXITR) )THEN
C           reduce number of unconverged intervals
            NINT = NINT - 1
C           Mark interval as converged.
            IWORK( K-1 ) = 0
            IF( I1.EQ.I ) THEN
               I1 = NEXT
            ELSE
C              Prev holds the last unconverged interval previously examined
               IF(PREV.GE.I1) IWORK( 2*PREV-1 ) = NEXT
            END IF
            I = NEXT
            GO TO 100
         END IF
         PREV = I
C
C        Perform one bisection step
C
         CNT = 0
         S = MID
         DPLUS = D( 1 ) - S
         IF( DPLUS.LT.ZERO ) CNT = CNT + 1
         DO 90 J = 2, N
            DPLUS = D( J ) - S - E2( J-1 )/DPLUS
            IF( DPLUS.LT.ZERO ) CNT = CNT + 1
 90      CONTINUE
         IF( CNT.LE.I-1 ) THEN
            WORK( K-1 ) = MID
         ELSE
            WORK( K ) = MID
         END IF
         I = NEXT

 100  CONTINUE
      ITER = ITER + 1
C     do another loop if there are still unconverged intervals
C     However, in the last iteration, all intervals are accepted
C     since this is the best we can do.
      IF( ( NINT.GT.0 ).AND.(ITER.LE.MAXITR) ) GO TO 80
C
C
C     At this point, all the intervals have converged
      DO 110 I = SAVI1, ILAST
         K = 2*I
         II = I - OFFSET
C        All intervals marked by '0' have been refined.
         IF( IWORK( K-1 ).EQ.0 ) THEN
            W( II ) = HALF*( WORK( K-1 )+WORK( K ) )
            WERR( II ) = WORK( K ) - W( II )
         END IF
 110  CONTINUE
C

      RETURN
C
C     End of DLARRJ
C
      END

      SUBROUTINE DLARRR( N, D, E, INFO )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2017
C
C     .. Scalar Arguments ..
      INTEGER            N, INFO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
C     ..
C
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, RELCOND
      PARAMETER          ( ZERO = 0.0D0,
     $                     RELCOND = 0.999D0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I
      LOGICAL            YESREL
      DOUBLE PRECISION   EPS, SAFMIN, SMLNUM, RMIN, TMP, TMP2,
     $          OFFDIG, OFFDIG2

C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         INFO = 0
         RETURN
      END IF
C
C     As a default, do NOT go for relative-accuracy preserving computations.
      INFO = 1

      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      RMIN = SQRT( SMLNUM )

C     Tests for relative accuracy
C
C     Test for scaled diagonal dominance
C     Scale the diagonal entries to one and check whether the sum of the
C     off-diagonals is less than one
C
C     The sdd relative error bounds have a 1/(1- 2*x) factor in them,
C     x = max(OFFDIG + OFFDIG2), so when x is close to 1/2, no relative
C     accuracy is promised.  In the notation of the code fragment below,
C     1/(1 - (OFFDIG + OFFDIG2)) is the condition number.
C     We don't think it is worth going into "sdd mode" unless the relative
C     condition number is reasonable, not 1/macheps.
C     The threshold should be compatible with other thresholds used in the
C     code. We set  OFFDIG + OFFDIG2 <= .999 =: RELCOND, it corresponds
C     to losing at most 3 decimal digits: 1 / (1 - (OFFDIG + OFFDIG2)) <= 1000
C     instead of the current OFFDIG + OFFDIG2 < 1
C
      YESREL = .TRUE.
      OFFDIG = ZERO
      TMP = SQRT(ABS(D(1)))
      IF (TMP.LT.RMIN) YESREL = .FALSE.
      IF(.NOT.YESREL) GOTO 11
      DO 10 I = 2, N
         TMP2 = SQRT(ABS(D(I)))
         IF (TMP2.LT.RMIN) YESREL = .FALSE.
         IF(.NOT.YESREL) GOTO 11
         OFFDIG2 = ABS(E(I-1))/(TMP*TMP2)
         IF(OFFDIG+OFFDIG2.GE.RELCOND) YESREL = .FALSE.
         IF(.NOT.YESREL) GOTO 11
         TMP = TMP2
         OFFDIG = OFFDIG2
 10   CONTINUE
 11   CONTINUE

      IF( YESREL ) THEN
         INFO = 0
         RETURN
      ELSE
      ENDIF
C

C
C     *** MORE TO BE IMPLEMENTED ***
C

C
C     Test if the lower bidiagonal matrix L from T = L D L^T
C     (zero shift facto) is well conditioned
C

C
C     Test if the upper bidiagonal matrix U from T = U D U^T
C     (zero shift facto) is well conditioned.
C     In this case, the matrix needs to be flipped and, at the end
C     of the eigenvector computation, the flip needs to be applied
C     to the computed eigenvectors (and the support)
C

C
      RETURN
C
C     END OF DLARRR
C
      END

      SUBROUTINE DLARRA( N, D, E, E2, SPLTOL, TNRM,
     $                    NSPLIT, ISPLIT, INFO )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2017
C
C     .. Scalar Arguments ..
      INTEGER            INFO, N, NSPLIT
      DOUBLE PRECISION    SPLTOL, TNRM
C     ..
C     .. Array Arguments ..
      INTEGER            ISPLIT( * )
      DOUBLE PRECISION   D( * ), E( * ), E2( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   EABS, TMP1

C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
C
      INFO = 0
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         RETURN
      END IF
C
C     Compute splitting points
      NSPLIT = 1
      IF(SPLTOL.LT.ZERO) THEN
C        Criterion based on absolute off-diagonal value
         TMP1 = ABS(SPLTOL)* TNRM
         DO 9 I = 1, N-1
            EABS = ABS( E(I) )
            IF( EABS .LE. TMP1) THEN
               E(I) = ZERO
               E2(I) = ZERO
               ISPLIT( NSPLIT ) = I
               NSPLIT = NSPLIT + 1
            END IF
 9       CONTINUE
      ELSE
C        Criterion that guarantees relative accuracy
         DO 10 I = 1, N-1
            EABS = ABS( E(I) )
            IF( EABS .LE. SPLTOL * SQRT(ABS(D(I)))*SQRT(ABS(D(I+1))) )
     $      THEN
               E(I) = ZERO
               E2(I) = ZERO
               ISPLIT( NSPLIT ) = I
               NSPLIT = NSPLIT + 1
            END IF
 10      CONTINUE
      ENDIF
      ISPLIT( NSPLIT ) = N

      RETURN
C
C     End of DLARRA
C
      END

      SUBROUTINE DLARRD( RANGE, ORDER, N, VL, VU, IL, IU, GERS,
     $                    RELTOL, D, E, E2, PIVMIN, NSPLIT, ISPLIT,
     $                    M, W, WERR, WL, WU, IBLOCK, INDEXW,
     $                    WORK, IWORK, INFO )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      CHARACTER          ORDER, RANGE
      INTEGER            IL, INFO, IU, M, N, NSPLIT
      DOUBLE PRECISION    PIVMIN, RELTOL, VL, VU, WL, WU
C     ..
C     .. Array Arguments ..
      INTEGER            IBLOCK( * ), INDEXW( * ),
     $                   ISPLIT( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), E2( * ),
     $                   GERS( * ), W( * ), WERR( * ), WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF, FUDGE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0,
     $                     TWO = 2.0D0, HALF = ONE/TWO,
     $                     FUDGE = TWO )
      INTEGER   ALLRNG, VALRNG, INDRNG
      PARAMETER ( ALLRNG = 1, VALRNG = 2, INDRNG = 3 )
C     ..
C     .. Local Scalars ..
      LOGICAL            NCNVRG, TOOFEW
      INTEGER            I, IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO,
     $                   IM, IN, IOFF, IOUT, IRANGE, ITMAX, ITMP1,
     $                   ITMP2, IW, IWOFF, J, JBLK, JDISC, JE, JEE, NB,
     $                   NWL, NWU
      DOUBLE PRECISION   ATOLI, EPS, GL, GU, RTOLI, TMP1, TMP2,
     $                   TNORM, UFLOW, WKILL, WLU, WUL

C     ..
C     .. Local Arrays ..
      INTEGER            IDUMMA( 1 )
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, ILAENV, DLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLAEBZ
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN
C     ..
C     .. Executable Statements ..
C
      INFO = 0
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         RETURN
      END IF
C
C     Decode RANGE
C
      IF( LSAME( RANGE, 'A' ) ) THEN
         IRANGE = ALLRNG
      ELSE IF( LSAME( RANGE, 'V' ) ) THEN
         IRANGE = VALRNG
      ELSE IF( LSAME( RANGE, 'I' ) ) THEN
         IRANGE = INDRNG
      ELSE
         IRANGE = 0
      END IF
C
C     Check for Errors
C
      IF( IRANGE.LE.0 ) THEN
         INFO = -1
      ELSE IF( .NOT.(LSAME(ORDER,'B').OR.LSAME(ORDER,'E')) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( IRANGE.EQ.VALRNG ) THEN
         IF( VL.GE.VU )
     $      INFO = -5
      ELSE IF( IRANGE.EQ.INDRNG .AND.
     $        ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) ) THEN
         INFO = -6
      ELSE IF( IRANGE.EQ.INDRNG .AND.
     $        ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) ) THEN
         INFO = -7
      END IF
C
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF

C     Initialize error flags
      INFO = 0
      NCNVRG = .FALSE.
      TOOFEW = .FALSE.

C     Quick return if possible
      M = 0
      IF( N.EQ.0 ) RETURN

C     Simplification:
      IF( IRANGE.EQ.INDRNG .AND. IL.EQ.1 .AND. IU.EQ.N ) IRANGE = 1

C     Get machine constants
      EPS = DLAMCH( 'P' )
      UFLOW = DLAMCH( 'U' )


C     Special Case when N=1
C     Treat case of 1x1 matrix for quick return
      IF( N.EQ.1 ) THEN
         IF( (IRANGE.EQ.ALLRNG).OR.
     $       ((IRANGE.EQ.VALRNG).AND.(D(1).GT.VL).AND.(D(1).LE.VU)).OR.
     $       ((IRANGE.EQ.INDRNG).AND.(IL.EQ.1).AND.(IU.EQ.1)) ) THEN
            M = 1
            W(1) = D(1)
C           The computation error of the eigenvalue is zero
            WERR(1) = ZERO
            IBLOCK( 1 ) = 1
            INDEXW( 1 ) = 1
         ENDIF
         RETURN
      END IF

C     NB is the minimum vector length for vector bisection, or 0
C     if only scalar is to be done.
      NB = ILAENV( 1, 'DSTEBZ', ' ', N, -1, -1, -1 )
      IF( NB.LE.1 ) NB = 0

C     Find global spectral radius
      GL = D(1)
      GU = D(1)
      DO 5 I = 1,N
         GL =  MIN( GL, GERS( 2*I - 1))
         GU = MAX( GU, GERS(2*I) )
 5    CONTINUE
C     Compute global Gerschgorin bounds and spectral diameter
      TNORM = MAX( ABS( GL ), ABS( GU ) )
      GL = GL - FUDGE*TNORM*EPS*N - FUDGE*TWO*PIVMIN
      GU = GU + FUDGE*TNORM*EPS*N + FUDGE*TWO*PIVMIN
C     [JAN/28/2009] remove the line below since SPDIAM variable not use
C     SPDIAM = GU - GL
C     Input arguments for DLAEBZ:
C     The relative tolerance.  An interval (a,b] lies within
C     "relative tolerance" if  b-a < RELTOL*max(|a|,|b|),
      RTOLI = RELTOL
C     Set the absolute tolerance for interval convergence to zero to force
C     interval convergence based on relative size of the interval.
C     This is dangerous because intervals might not converge when RELTOL is
C     small. But at least a very small number should be selected so that for
C     strongly graded matrices, the code can get relatively accurate
C     eigenvalues.
      ATOLI = FUDGE*TWO*UFLOW + FUDGE*TWO*PIVMIN

      IF( IRANGE.EQ.INDRNG ) THEN

C        RANGE='I': Compute an interval containing eigenvalues
C        IL through IU. The initial interval [GL,GU] from the global
C        Gerschgorin bounds GL and GU is refined by DLAEBZ.
         ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) /
     $           LOG( TWO ) ) + 2
         WORK( N+1 ) = GL
         WORK( N+2 ) = GL
         WORK( N+3 ) = GU
         WORK( N+4 ) = GU
         WORK( N+5 ) = GL
         WORK( N+6 ) = GU
         IWORK( 1 ) = -1
         IWORK( 2 ) = -1
         IWORK( 3 ) = N + 1
         IWORK( 4 ) = N + 1
         IWORK( 5 ) = IL - 1
         IWORK( 6 ) = IU
C
         CALL DLAEBZ( 3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN,
     $         D, E, E2, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT,
     $                IWORK, W, IBLOCK, IINFO )
         IF( IINFO .NE. 0 ) THEN
            INFO = IINFO
            RETURN
         END IF
C        On exit, output intervals may not be ordered by ascending negcount
         IF( IWORK( 6 ).EQ.IU ) THEN
            WL = WORK( N+1 )
            WLU = WORK( N+3 )
            NWL = IWORK( 1 )
            WU = WORK( N+4 )
            WUL = WORK( N+2 )
            NWU = IWORK( 4 )
         ELSE
            WL = WORK( N+2 )
            WLU = WORK( N+4 )
            NWL = IWORK( 2 )
            WU = WORK( N+3 )
            WUL = WORK( N+1 )
            NWU = IWORK( 3 )
         END IF
C        On exit, the interval [WL, WLU] contains a value with negcount NWL,
C        and [WUL, WU] contains a value with negcount NWU.
         IF( NWL.LT.0 .OR. NWL.GE.N .OR. NWU.LT.1 .OR. NWU.GT.N ) THEN
            INFO = 4
            RETURN
         END IF

      ELSEIF( IRANGE.EQ.VALRNG ) THEN
         WL = VL
         WU = VU

      ELSEIF( IRANGE.EQ.ALLRNG ) THEN
         WL = GL
         WU = GU
      ENDIF



C     Find Eigenvalues -- Loop Over blocks and recompute NWL and NWU.
C     NWL accumulates the number of eigenvalues .le. WL,
C     NWU accumulates the number of eigenvalues .le. WU
      M = 0
      IEND = 0
      INFO = 0
      NWL = 0
      NWU = 0
C
      DO 70 JBLK = 1, NSPLIT
         IOFF = IEND
         IBEGIN = IOFF + 1
         IEND = ISPLIT( JBLK )
         IN = IEND - IOFF
C
         IF( IN.EQ.1 ) THEN
C           1x1 block
            IF( WL.GE.D( IBEGIN )-PIVMIN )
     $         NWL = NWL + 1
            IF( WU.GE.D( IBEGIN )-PIVMIN )
     $         NWU = NWU + 1
            IF( IRANGE.EQ.ALLRNG .OR.
     $           ( WL.LT.D( IBEGIN )-PIVMIN
     $             .AND. WU.GE. D( IBEGIN )-PIVMIN ) ) THEN
               M = M + 1
               W( M ) = D( IBEGIN )
               WERR(M) = ZERO
C              The gap for a single block doesn't matter for the later
C              algorithm and is assigned an arbitrary large value
               IBLOCK( M ) = JBLK
               INDEXW( M ) = 1
            END IF

C        Disabled 2x2 case because of a failure on the following matrix
C        RANGE = 'I', IL = IU = 4
C          Original Tridiagonal, d = [
C           -0.150102010615740E+00
C           -0.849897989384260E+00
C           -0.128208148052635E-15
C            0.128257718286320E-15
C          ];
C          e = [
C           -0.357171383266986E+00
C           -0.180411241501588E-15
C           -0.175152352710251E-15
C          ];
C
C         ELSE IF( IN.EQ.2 ) THEN
C*           2x2 block
C            DISC = SQRT( (HALF*(D(IBEGIN)-D(IEND)))**2 + E(IBEGIN)**2 )
C            TMP1 = HALF*(D(IBEGIN)+D(IEND))
C            L1 = TMP1 - DISC
C            IF( WL.GE. L1-PIVMIN )
C     $         NWL = NWL + 1
C            IF( WU.GE. L1-PIVMIN )
C     $         NWU = NWU + 1
C            IF( IRANGE.EQ.ALLRNG .OR. ( WL.LT.L1-PIVMIN .AND. WU.GE.
C     $          L1-PIVMIN ) ) THEN
C               M = M + 1
C               W( M ) = L1
C*              The uncertainty of eigenvalues of a 2x2 matrix is very small
C               WERR( M ) = EPS * ABS( W( M ) ) * TWO
C               IBLOCK( M ) = JBLK
C               INDEXW( M ) = 1
C            ENDIF
C            L2 = TMP1 + DISC
C            IF( WL.GE. L2-PIVMIN )
C     $         NWL = NWL + 1
C            IF( WU.GE. L2-PIVMIN )
C     $         NWU = NWU + 1
C            IF( IRANGE.EQ.ALLRNG .OR. ( WL.LT.L2-PIVMIN .AND. WU.GE.
C     $          L2-PIVMIN ) ) THEN
C               M = M + 1
C               W( M ) = L2
C*              The uncertainty of eigenvalues of a 2x2 matrix is very small
C               WERR( M ) = EPS * ABS( W( M ) ) * TWO
C               IBLOCK( M ) = JBLK
C               INDEXW( M ) = 2
C            ENDIF
         ELSE
C           General Case - block of size IN >= 2
C           Compute local Gerschgorin interval and use it as the initial
C           interval for DLAEBZ
            GU = D( IBEGIN )
            GL = D( IBEGIN )
            TMP1 = ZERO

            DO 40 J = IBEGIN, IEND
               GL =  MIN( GL, GERS( 2*J - 1))
               GU = MAX( GU, GERS(2*J) )
   40       CONTINUE
C           [JAN/28/2009]
C           change SPDIAM by TNORM in lines 2 and 3 thereafter
C           line 1: remove computation of SPDIAM (not useful anymore)
C           SPDIAM = GU - GL
C           GL = GL - FUDGE*SPDIAM*EPS*IN - FUDGE*PIVMIN
C           GU = GU + FUDGE*SPDIAM*EPS*IN + FUDGE*PIVMIN
            GL = GL - FUDGE*TNORM*EPS*IN - FUDGE*PIVMIN
            GU = GU + FUDGE*TNORM*EPS*IN + FUDGE*PIVMIN
C
            IF( IRANGE.GT.1 ) THEN
               IF( GU.LT.WL ) THEN
C                 the local block contains none of the wanted eigenvalues
                  NWL = NWL + IN
                  NWU = NWU + IN
                  GO TO 70
               END IF
C              refine search interval if possible, only range (WL,WU] matters
               GL = MAX( GL, WL )
               GU = MIN( GU, WU )
               IF( GL.GE.GU )
     $            GO TO 70
            END IF

C           Find negcount of initial interval boundaries GL and GU
            WORK( N+1 ) = GL
            WORK( N+IN+1 ) = GU
            CALL DLAEBZ( 1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
     $                   D( IBEGIN ), E( IBEGIN ), E2( IBEGIN ),
     $                   IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM,
     $                   IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
            IF( IINFO .NE. 0 ) THEN
               INFO = IINFO
               RETURN
            END IF
C
            NWL = NWL + IWORK( 1 )
            NWU = NWU + IWORK( IN+1 )
            IWOFF = M - IWORK( 1 )

C           Compute Eigenvalues
            ITMAX = INT( ( LOG( GU-GL+PIVMIN )-LOG( PIVMIN ) ) /
     $              LOG( TWO ) ) + 2
            CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
     $                   D( IBEGIN ), E( IBEGIN ), E2( IBEGIN ),
     $                   IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT,
     $                   IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
            IF( IINFO .NE. 0 ) THEN
               INFO = IINFO
               RETURN
            END IF
C
C           Copy eigenvalues into W and IBLOCK
C           Use -JBLK for block number for unconverged eigenvalues.
C           Loop over the number of output intervals from DLAEBZ
            DO 60 J = 1, IOUT
C              eigenvalue approximation is middle point of interval
               TMP1 = HALF*( WORK( J+N )+WORK( J+IN+N ) )
C              semi length of error interval
               TMP2 = HALF*ABS( WORK( J+N )-WORK( J+IN+N ) )
               IF( J.GT.IOUT-IINFO ) THEN
C                 Flag non-convergence.
                  NCNVRG = .TRUE.
                  IB = -JBLK
               ELSE
                  IB = JBLK
               END IF
               DO 50 JE = IWORK( J ) + 1 + IWOFF,
     $                 IWORK( J+IN ) + IWOFF
                  W( JE ) = TMP1
                  WERR( JE ) = TMP2
                  INDEXW( JE ) = JE - IWOFF
                  IBLOCK( JE ) = IB
   50          CONTINUE
   60       CONTINUE
C
            M = M + IM
         END IF
   70 CONTINUE

C     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
C     If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
      IF( IRANGE.EQ.INDRNG ) THEN
         IDISCL = IL - 1 - NWL
         IDISCU = NWU - IU
C
         IF( IDISCL.GT.0 ) THEN
            IM = 0
            DO 80 JE = 1, M
C              Remove some of the smallest eigenvalues from the left so that
C              at the end IDISCL =0. Move all eigenvalues up to the left.
               IF( W( JE ).LE.WLU .AND. IDISCL.GT.0 ) THEN
                  IDISCL = IDISCL - 1
               ELSE
                  IM = IM + 1
                  W( IM ) = W( JE )
                  WERR( IM ) = WERR( JE )
                  INDEXW( IM ) = INDEXW( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               END IF
 80         CONTINUE
            M = IM
         END IF
         IF( IDISCU.GT.0 ) THEN
C           Remove some of the largest eigenvalues from the right so that
C           at the end IDISCU =0. Move all eigenvalues up to the left.
            IM=M+1
            DO 81 JE = M, 1, -1
               IF( W( JE ).GE.WUL .AND. IDISCU.GT.0 ) THEN
                  IDISCU = IDISCU - 1
               ELSE
                  IM = IM - 1
                  W( IM ) = W( JE )
                  WERR( IM ) = WERR( JE )
                  INDEXW( IM ) = INDEXW( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               END IF
 81         CONTINUE
            JEE = 0
            DO 82 JE = IM, M
               JEE = JEE + 1
               W( JEE ) = W( JE )
               WERR( JEE ) = WERR( JE )
               INDEXW( JEE ) = INDEXW( JE )
               IBLOCK( JEE ) = IBLOCK( JE )
 82         CONTINUE
            M = M-IM+1
         END IF

         IF( IDISCL.GT.0 .OR. IDISCU.GT.0 ) THEN
C           Code to deal with effects of bad arithmetic. (If N(w) is
C           monotone non-decreasing, this should never happen.)
C           Some low eigenvalues to be discarded are not in (WL,WLU],
C           or high eigenvalues to be discarded are not in (WUL,WU]
C           so just kill off the smallest IDISCL/largest IDISCU
C           eigenvalues, by marking the corresponding IBLOCK = 0
            IF( IDISCL.GT.0 ) THEN
               WKILL = WU
               DO 100 JDISC = 1, IDISCL
                  IW = 0
                  DO 90 JE = 1, M
                     IF( IBLOCK( JE ).NE.0 .AND.
     $                    ( W( JE ).LT.WKILL .OR. IW.EQ.0 ) ) THEN
                        IW = JE
                        WKILL = W( JE )
                     END IF
 90               CONTINUE
                  IBLOCK( IW ) = 0
 100           CONTINUE
            END IF
            IF( IDISCU.GT.0 ) THEN
               WKILL = WL
               DO 120 JDISC = 1, IDISCU
                  IW = 0
                  DO 110 JE = 1, M
                     IF( IBLOCK( JE ).NE.0 .AND.
     $                    ( W( JE ).GE.WKILL .OR. IW.EQ.0 ) ) THEN
                        IW = JE
                        WKILL = W( JE )
                     END IF
 110              CONTINUE
                  IBLOCK( IW ) = 0
 120           CONTINUE
            END IF
C           Now erase all eigenvalues with IBLOCK set to zero
            IM = 0
            DO 130 JE = 1, M
               IF( IBLOCK( JE ).NE.0 ) THEN
                  IM = IM + 1
                  W( IM ) = W( JE )
                  WERR( IM ) = WERR( JE )
                  INDEXW( IM ) = INDEXW( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               END IF
 130        CONTINUE
            M = IM
         END IF
         IF( IDISCL.LT.0 .OR. IDISCU.LT.0 ) THEN
            TOOFEW = .TRUE.
         END IF
      END IF
C
      IF(( IRANGE.EQ.ALLRNG .AND. M.NE.N ).OR.
     $   ( IRANGE.EQ.INDRNG .AND. M.NE.IU-IL+1 ) ) THEN
         TOOFEW = .TRUE.
      END IF

C     If ORDER='B', do nothing the eigenvalues are already sorted by
C        block.
C     If ORDER='E', sort the eigenvalues from smallest to largest

      IF( LSAME(ORDER,'E') .AND. NSPLIT.GT.1 ) THEN
         DO 150 JE = 1, M - 1
            IE = 0
            TMP1 = W( JE )
            DO 140 J = JE + 1, M
               IF( W( J ).LT.TMP1 ) THEN
                  IE = J
                  TMP1 = W( J )
               END IF
  140       CONTINUE
            IF( IE.NE.0 ) THEN
               TMP2 = WERR( IE )
               ITMP1 = IBLOCK( IE )
               ITMP2 = INDEXW( IE )
               W( IE ) = W( JE )
               WERR( IE ) = WERR( JE )
               IBLOCK( IE ) = IBLOCK( JE )
               INDEXW( IE ) = INDEXW( JE )
               W( JE ) = TMP1
               WERR( JE ) = TMP2
               IBLOCK( JE ) = ITMP1
               INDEXW( JE ) = ITMP2
            END IF
  150    CONTINUE
      END IF
C
      INFO = 0
      IF( NCNVRG )
     $   INFO = INFO + 1
      IF( TOOFEW )
     $   INFO = INFO + 2
      RETURN
C
C     End of DLARRD
C
      END

      SUBROUTINE DLARRK( N, IW, GL, GU,
     $                    D, E2, PIVMIN, RELTOL, W, WERR, INFO)
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2017
C
C     .. Scalar Arguments ..
      INTEGER   INFO, IW, N
      DOUBLE PRECISION    PIVMIN, RELTOL, GL, GU, W, WERR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E2( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   FUDGE, HALF, TWO, ZERO
      PARAMETER          ( HALF = 0.5D0, TWO = 2.0D0,
     $                     FUDGE = TWO, ZERO = 0.0D0 )
C     ..
C     .. Local Scalars ..
      INTEGER   I, IT, ITMAX, NEGCNT
      DOUBLE PRECISION   ATOLI, EPS, LEFT, MID, RIGHT, RTOLI, TMP1,
     $                   TMP2, TNORM
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL   DLAMCH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX
C     ..
C     .. Executable Statements ..
C
C     Quick return if possible
C
      IF( N.LE.0 ) THEN
         INFO = 0
         RETURN
      END IF
C
C     Get machine constants
      EPS = DLAMCH( 'P' )

      TNORM = MAX( ABS( GL ), ABS( GU ) )
      RTOLI = RELTOL
      ATOLI = FUDGE*TWO*PIVMIN

      ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) /
     $           LOG( TWO ) ) + 2

      INFO = -1

      LEFT = GL - FUDGE*TNORM*EPS*N - FUDGE*TWO*PIVMIN
      RIGHT = GU + FUDGE*TNORM*EPS*N + FUDGE*TWO*PIVMIN
      IT = 0

 10   CONTINUE
C
C     Check if interval converged or maximum number of iterations reached
C
      TMP1 = ABS( RIGHT - LEFT )
      TMP2 = MAX( ABS(RIGHT), ABS(LEFT) )
      IF( TMP1.LT.MAX( ATOLI, PIVMIN, RTOLI*TMP2 ) ) THEN
         INFO = 0
         GOTO 30
      ENDIF
      IF(IT.GT.ITMAX)
     $   GOTO 30

C
C     Count number of negative pivots for mid-point
C
      IT = IT + 1
      MID = HALF * (LEFT + RIGHT)
      NEGCNT = 0
      TMP1 = D( 1 ) - MID
      IF( ABS( TMP1 ).LT.PIVMIN )
     $   TMP1 = -PIVMIN
      IF( TMP1.LE.ZERO )
     $   NEGCNT = NEGCNT + 1
C
      DO 20 I = 2, N
         TMP1 = D( I ) - E2( I-1 ) / TMP1 - MID
         IF( ABS( TMP1 ).LT.PIVMIN )
     $      TMP1 = -PIVMIN
         IF( TMP1.LE.ZERO )
     $      NEGCNT = NEGCNT + 1
 20   CONTINUE

      IF(NEGCNT.GE.IW) THEN
         RIGHT = MID
      ELSE
         LEFT = MID
      ENDIF
      GOTO 10

 30   CONTINUE
C
C     Converged or maximum number of iterations reached
C
      W = HALF * (LEFT + RIGHT)
      WERR = HALF * ABS( RIGHT - LEFT )

      RETURN
C
C     End of DLARRK
C
      END

      SUBROUTINE DLARUV( ISEED, N, X )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            N
C     ..
C     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( N )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      INTEGER            LV, IPW2
      DOUBLE PRECISION   R
      PARAMETER          ( LV = 128, IPW2 = 4096, R = ONE / IPW2 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, I1, I2, I3, I4, IT1, IT2, IT3, IT4, J
C     ..
C     .. Local Arrays ..
      INTEGER            MM( LV, 4 )
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN, MOD
C     ..
C     .. Data statements ..
      DATA               ( MM( 1, J ), J = 1, 4 ) / 494, 322, 2508,
     $                   2549 /
      DATA               ( MM( 2, J ), J = 1, 4 ) / 2637, 789, 3754,
     $                   1145 /
      DATA               ( MM( 3, J ), J = 1, 4 ) / 255, 1440, 1766,
     $                   2253 /
      DATA               ( MM( 4, J ), J = 1, 4 ) / 2008, 752, 3572,
     $                   305 /
      DATA               ( MM( 5, J ), J = 1, 4 ) / 1253, 2859, 2893,
     $                   3301 /
      DATA               ( MM( 6, J ), J = 1, 4 ) / 3344, 123, 307,
     $                   1065 /
      DATA               ( MM( 7, J ), J = 1, 4 ) / 4084, 1848, 1297,
     $                   3133 /
      DATA               ( MM( 8, J ), J = 1, 4 ) / 1739, 643, 3966,
     $                   2913 /
      DATA               ( MM( 9, J ), J = 1, 4 ) / 3143, 2405, 758,
     $                   3285 /
      DATA               ( MM( 10, J ), J = 1, 4 ) / 3468, 2638, 2598,
     $                   1241 /
      DATA               ( MM( 11, J ), J = 1, 4 ) / 688, 2344, 3406,
     $                   1197 /
      DATA               ( MM( 12, J ), J = 1, 4 ) / 1657, 46, 2922,
     $                   3729 /
      DATA               ( MM( 13, J ), J = 1, 4 ) / 1238, 3814, 1038,
     $                   2501 /
      DATA               ( MM( 14, J ), J = 1, 4 ) / 3166, 913, 2934,
     $                   1673 /
      DATA               ( MM( 15, J ), J = 1, 4 ) / 1292, 3649, 2091,
     $                   541 /
      DATA               ( MM( 16, J ), J = 1, 4 ) / 3422, 339, 2451,
     $                   2753 /
      DATA               ( MM( 17, J ), J = 1, 4 ) / 1270, 3808, 1580,
     $                   949 /
      DATA               ( MM( 18, J ), J = 1, 4 ) / 2016, 822, 1958,
     $                   2361 /
      DATA               ( MM( 19, J ), J = 1, 4 ) / 154, 2832, 2055,
     $                   1165 /
      DATA               ( MM( 20, J ), J = 1, 4 ) / 2862, 3078, 1507,
     $                   4081 /
      DATA               ( MM( 21, J ), J = 1, 4 ) / 697, 3633, 1078,
     $                   2725 /
      DATA               ( MM( 22, J ), J = 1, 4 ) / 1706, 2970, 3273,
     $                   3305 /
      DATA               ( MM( 23, J ), J = 1, 4 ) / 491, 637, 17,
     $                   3069 /
      DATA               ( MM( 24, J ), J = 1, 4 ) / 931, 2249, 854,
     $                   3617 /
      DATA               ( MM( 25, J ), J = 1, 4 ) / 1444, 2081, 2916,
     $                   3733 /
      DATA               ( MM( 26, J ), J = 1, 4 ) / 444, 4019, 3971,
     $                   409 /
      DATA               ( MM( 27, J ), J = 1, 4 ) / 3577, 1478, 2889,
     $                   2157 /
      DATA               ( MM( 28, J ), J = 1, 4 ) / 3944, 242, 3831,
     $                   1361 /
      DATA               ( MM( 29, J ), J = 1, 4 ) / 2184, 481, 2621,
     $                   3973 /
      DATA               ( MM( 30, J ), J = 1, 4 ) / 1661, 2075, 1541,
     $                   1865 /
      DATA               ( MM( 31, J ), J = 1, 4 ) / 3482, 4058, 893,
     $                   2525 /
      DATA               ( MM( 32, J ), J = 1, 4 ) / 657, 622, 736,
     $                   1409 /
      DATA               ( MM( 33, J ), J = 1, 4 ) / 3023, 3376, 3992,
     $                   3445 /
      DATA               ( MM( 34, J ), J = 1, 4 ) / 3618, 812, 787,
     $                   3577 /
      DATA               ( MM( 35, J ), J = 1, 4 ) / 1267, 234, 2125,
     $                   77 /
      DATA               ( MM( 36, J ), J = 1, 4 ) / 1828, 641, 2364,
     $                   3761 /
      DATA               ( MM( 37, J ), J = 1, 4 ) / 164, 4005, 2460,
     $                   2149 /
      DATA               ( MM( 38, J ), J = 1, 4 ) / 3798, 1122, 257,
     $                   1449 /
      DATA               ( MM( 39, J ), J = 1, 4 ) / 3087, 3135, 1574,
     $                   3005 /
      DATA               ( MM( 40, J ), J = 1, 4 ) / 2400, 2640, 3912,
     $                   225 /
      DATA               ( MM( 41, J ), J = 1, 4 ) / 2870, 2302, 1216,
     $                   85 /
      DATA               ( MM( 42, J ), J = 1, 4 ) / 3876, 40, 3248,
     $                   3673 /
      DATA               ( MM( 43, J ), J = 1, 4 ) / 1905, 1832, 3401,
     $                   3117 /
      DATA               ( MM( 44, J ), J = 1, 4 ) / 1593, 2247, 2124,
     $                   3089 /
      DATA               ( MM( 45, J ), J = 1, 4 ) / 1797, 2034, 2762,
     $                   1349 /
      DATA               ( MM( 46, J ), J = 1, 4 ) / 1234, 2637, 149,
     $                   2057 /
      DATA               ( MM( 47, J ), J = 1, 4 ) / 3460, 1287, 2245,
     $                   413 /
      DATA               ( MM( 48, J ), J = 1, 4 ) / 328, 1691, 166,
     $                   65 /
      DATA               ( MM( 49, J ), J = 1, 4 ) / 2861, 496, 466,
     $                   1845 /
      DATA               ( MM( 50, J ), J = 1, 4 ) / 1950, 1597, 4018,
     $                   697 /
      DATA               ( MM( 51, J ), J = 1, 4 ) / 617, 2394, 1399,
     $                   3085 /
      DATA               ( MM( 52, J ), J = 1, 4 ) / 2070, 2584, 190,
     $                   3441 /
      DATA               ( MM( 53, J ), J = 1, 4 ) / 3331, 1843, 2879,
     $                   1573 /
      DATA               ( MM( 54, J ), J = 1, 4 ) / 769, 336, 153,
     $                   3689 /
      DATA               ( MM( 55, J ), J = 1, 4 ) / 1558, 1472, 2320,
     $                   2941 /
      DATA               ( MM( 56, J ), J = 1, 4 ) / 2412, 2407, 18,
     $                   929 /
      DATA               ( MM( 57, J ), J = 1, 4 ) / 2800, 433, 712,
     $                   533 /
      DATA               ( MM( 58, J ), J = 1, 4 ) / 189, 2096, 2159,
     $                   2841 /
      DATA               ( MM( 59, J ), J = 1, 4 ) / 287, 1761, 2318,
     $                   4077 /
      DATA               ( MM( 60, J ), J = 1, 4 ) / 2045, 2810, 2091,
     $                   721 /
      DATA               ( MM( 61, J ), J = 1, 4 ) / 1227, 566, 3443,
     $                   2821 /
      DATA               ( MM( 62, J ), J = 1, 4 ) / 2838, 442, 1510,
     $                   2249 /
      DATA               ( MM( 63, J ), J = 1, 4 ) / 209, 41, 449,
     $                   2397 /
      DATA               ( MM( 64, J ), J = 1, 4 ) / 2770, 1238, 1956,
     $                   2817 /
      DATA               ( MM( 65, J ), J = 1, 4 ) / 3654, 1086, 2201,
     $                   245 /
      DATA               ( MM( 66, J ), J = 1, 4 ) / 3993, 603, 3137,
     $                   1913 /
      DATA               ( MM( 67, J ), J = 1, 4 ) / 192, 840, 3399,
     $                   1997 /
      DATA               ( MM( 68, J ), J = 1, 4 ) / 2253, 3168, 1321,
     $                   3121 /
      DATA               ( MM( 69, J ), J = 1, 4 ) / 3491, 1499, 2271,
     $                   997 /
      DATA               ( MM( 70, J ), J = 1, 4 ) / 2889, 1084, 3667,
     $                   1833 /
      DATA               ( MM( 71, J ), J = 1, 4 ) / 2857, 3438, 2703,
     $                   2877 /
      DATA               ( MM( 72, J ), J = 1, 4 ) / 2094, 2408, 629,
     $                   1633 /
      DATA               ( MM( 73, J ), J = 1, 4 ) / 1818, 1589, 2365,
     $                   981 /
      DATA               ( MM( 74, J ), J = 1, 4 ) / 688, 2391, 2431,
     $                   2009 /
      DATA               ( MM( 75, J ), J = 1, 4 ) / 1407, 288, 1113,
     $                   941 /
      DATA               ( MM( 76, J ), J = 1, 4 ) / 634, 26, 3922,
     $                   2449 /
      DATA               ( MM( 77, J ), J = 1, 4 ) / 3231, 512, 2554,
     $                   197 /
      DATA               ( MM( 78, J ), J = 1, 4 ) / 815, 1456, 184,
     $                   2441 /
      DATA               ( MM( 79, J ), J = 1, 4 ) / 3524, 171, 2099,
     $                   285 /
      DATA               ( MM( 80, J ), J = 1, 4 ) / 1914, 1677, 3228,
     $                   1473 /
      DATA               ( MM( 81, J ), J = 1, 4 ) / 516, 2657, 4012,
     $                   2741 /
      DATA               ( MM( 82, J ), J = 1, 4 ) / 164, 2270, 1921,
     $                   3129 /
      DATA               ( MM( 83, J ), J = 1, 4 ) / 303, 2587, 3452,
     $                   909 /
      DATA               ( MM( 84, J ), J = 1, 4 ) / 2144, 2961, 3901,
     $                   2801 /
      DATA               ( MM( 85, J ), J = 1, 4 ) / 3480, 1970, 572,
     $                   421 /
      DATA               ( MM( 86, J ), J = 1, 4 ) / 119, 1817, 3309,
     $                   4073 /
      DATA               ( MM( 87, J ), J = 1, 4 ) / 3357, 676, 3171,
     $                   2813 /
      DATA               ( MM( 88, J ), J = 1, 4 ) / 837, 1410, 817,
     $                   2337 /
      DATA               ( MM( 89, J ), J = 1, 4 ) / 2826, 3723, 3039,
     $                   1429 /
      DATA               ( MM( 90, J ), J = 1, 4 ) / 2332, 2803, 1696,
     $                   1177 /
      DATA               ( MM( 91, J ), J = 1, 4 ) / 2089, 3185, 1256,
     $                   1901 /
      DATA               ( MM( 92, J ), J = 1, 4 ) / 3780, 184, 3715,
     $                   81 /
      DATA               ( MM( 93, J ), J = 1, 4 ) / 1700, 663, 2077,
     $                   1669 /
      DATA               ( MM( 94, J ), J = 1, 4 ) / 3712, 499, 3019,
     $                   2633 /
      DATA               ( MM( 95, J ), J = 1, 4 ) / 150, 3784, 1497,
     $                   2269 /
      DATA               ( MM( 96, J ), J = 1, 4 ) / 2000, 1631, 1101,
     $                   129 /
      DATA               ( MM( 97, J ), J = 1, 4 ) / 3375, 1925, 717,
     $                   1141 /
      DATA               ( MM( 98, J ), J = 1, 4 ) / 1621, 3912, 51,
     $                   249 /
      DATA               ( MM( 99, J ), J = 1, 4 ) / 3090, 1398, 981,
     $                   3917 /
      DATA               ( MM( 100, J ), J = 1, 4 ) / 3765, 1349, 1978,
     $                   2481 /
      DATA               ( MM( 101, J ), J = 1, 4 ) / 1149, 1441, 1813,
     $                   3941 /
      DATA               ( MM( 102, J ), J = 1, 4 ) / 3146, 2224, 3881,
     $                   2217 /
      DATA               ( MM( 103, J ), J = 1, 4 ) / 33, 2411, 76,
     $                   2749 /
      DATA               ( MM( 104, J ), J = 1, 4 ) / 3082, 1907, 3846,
     $                   3041 /
      DATA               ( MM( 105, J ), J = 1, 4 ) / 2741, 3192, 3694,
     $                   1877 /
      DATA               ( MM( 106, J ), J = 1, 4 ) / 359, 2786, 1682,
     $                   345 /
      DATA               ( MM( 107, J ), J = 1, 4 ) / 3316, 382, 124,
     $                   2861 /
      DATA               ( MM( 108, J ), J = 1, 4 ) / 1749, 37, 1660,
     $                   1809 /
      DATA               ( MM( 109, J ), J = 1, 4 ) / 185, 759, 3997,
     $                   3141 /
      DATA               ( MM( 110, J ), J = 1, 4 ) / 2784, 2948, 479,
     $                   2825 /
      DATA               ( MM( 111, J ), J = 1, 4 ) / 2202, 1862, 1141,
     $                   157 /
      DATA               ( MM( 112, J ), J = 1, 4 ) / 2199, 3802, 886,
     $                   2881 /
      DATA               ( MM( 113, J ), J = 1, 4 ) / 1364, 2423, 3514,
     $                   3637 /
      DATA               ( MM( 114, J ), J = 1, 4 ) / 1244, 2051, 1301,
     $                   1465 /
      DATA               ( MM( 115, J ), J = 1, 4 ) / 2020, 2295, 3604,
     $                   2829 /
      DATA               ( MM( 116, J ), J = 1, 4 ) / 3160, 1332, 1888,
     $                   2161 /
      DATA               ( MM( 117, J ), J = 1, 4 ) / 2785, 1832, 1836,
     $                   3365 /
      DATA               ( MM( 118, J ), J = 1, 4 ) / 2772, 2405, 1990,
     $                   361 /
      DATA               ( MM( 119, J ), J = 1, 4 ) / 1217, 3638, 2058,
     $                   2685 /
      DATA               ( MM( 120, J ), J = 1, 4 ) / 1822, 3661, 692,
     $                   3745 /
      DATA               ( MM( 121, J ), J = 1, 4 ) / 1245, 327, 1194,
     $                   2325 /
      DATA               ( MM( 122, J ), J = 1, 4 ) / 2252, 3660, 20,
     $                   3609 /
      DATA               ( MM( 123, J ), J = 1, 4 ) / 3904, 716, 3285,
     $                   3821 /
      DATA               ( MM( 124, J ), J = 1, 4 ) / 2774, 1842, 2046,
     $                   3537 /
      DATA               ( MM( 125, J ), J = 1, 4 ) / 997, 3987, 2107,
     $                   517 /
      DATA               ( MM( 126, J ), J = 1, 4 ) / 2573, 1368, 3508,
     $                   3017 /
      DATA               ( MM( 127, J ), J = 1, 4 ) / 1148, 1848, 3525,
     $                   2141 /
      DATA               ( MM( 128, J ), J = 1, 4 ) / 545, 2366, 3801,
     $                   1537 /
C     ..
C     .. Executable Statements ..
C
      I1 = ISEED( 1 )
      I2 = ISEED( 2 )
      I3 = ISEED( 3 )
      I4 = ISEED( 4 )
C
      DO 10 I = 1, MIN( N, LV )
C
  20     CONTINUE
C
C        Multiply the seed by i-th power of the multiplier modulo 2**48
C
         IT4 = I4*MM( I, 4 )
         IT3 = IT4 / IPW2
         IT4 = IT4 - IPW2*IT3
         IT3 = IT3 + I3*MM( I, 4 ) + I4*MM( I, 3 )
         IT2 = IT3 / IPW2
         IT3 = IT3 - IPW2*IT2
         IT2 = IT2 + I2*MM( I, 4 ) + I3*MM( I, 3 ) + I4*MM( I, 2 )
         IT1 = IT2 / IPW2
         IT2 = IT2 - IPW2*IT1
         IT1 = IT1 + I1*MM( I, 4 ) + I2*MM( I, 3 ) + I3*MM( I, 2 ) +
     $         I4*MM( I, 1 )
         IT1 = MOD( IT1, IPW2 )
C
C        Convert 48-bit integer to a real number in the interval (0,1)
C
         X( I ) = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )+R*
     $            DBLE( IT4 ) ) ) )
C
         IF (X( I ).EQ.1.0D0) THEN
C           If a real number has n bits of precision, and the first
C           n bits of the 48-bit integer above happen to be all 1 (which
C           will occur about once every 2**n calls), then X( I ) will
C           be rounded to exactly 1.0.
C           Since X( I ) is not supposed to return exactly 0.0 or 1.0,
C           the statistically correct thing to do in this situation is
C           simply to iterate again.
C           N.B. the case X( I ) = 0.0 should not be possible.
            I1 = I1 + 2
            I2 = I2 + 2
            I3 = I3 + 2
            I4 = I4 + 2
            GOTO 20
         END IF
C
   10 CONTINUE
C
C     Return final value of seed
C
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
      RETURN
C
C     End of DLARUV
C
      END

      SUBROUTINE DLAEBZ( IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL,
     $                   RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT,
     $                   NAB, WORK, IWORK, INFO )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            IJOB, INFO, MINP, MMAX, MOUT, N, NBMIN, NITMAX
      DOUBLE PRECISION   ABSTOL, PIVMIN, RELTOL
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * ), NAB( MMAX, * ), NVAL( * )
      DOUBLE PRECISION   AB( MMAX, * ), C( * ), D( * ), E( * ), E2( * ),
     $                   WORK( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, TWO, HALF
      PARAMETER          ( ZERO = 0.0D0, TWO = 2.0D0,
     $                   HALF = 1.0D0 / TWO )
C     ..
C     .. Local Scalars ..
      INTEGER            ITMP1, ITMP2, J, JI, JIT, JP, KF, KFNEW, KL,
     $                   KLNEW
      DOUBLE PRECISION   TMP1, TMP2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Check for Errors
C
      INFO = 0
      IF( IJOB.LT.1 .OR. IJOB.GT.3 ) THEN
         INFO = -1
         RETURN
      END IF
C
C     Initialize NAB
C
      IF( IJOB.EQ.1 ) THEN
C
C        Compute the number of eigenvalues in the initial intervals.
C
         MOUT = 0
         DO 30 JI = 1, MINP
            DO 20 JP = 1, 2
               TMP1 = D( 1 ) - AB( JI, JP )
               IF( ABS( TMP1 ).LT.PIVMIN )
     $            TMP1 = -PIVMIN
               NAB( JI, JP ) = 0
               IF( TMP1.LE.ZERO )
     $            NAB( JI, JP ) = 1
C
               DO 10 J = 2, N
                  TMP1 = D( J ) - E2( J-1 ) / TMP1 - AB( JI, JP )
                  IF( ABS( TMP1 ).LT.PIVMIN )
     $               TMP1 = -PIVMIN
                  IF( TMP1.LE.ZERO )
     $               NAB( JI, JP ) = NAB( JI, JP ) + 1
   10          CONTINUE
   20       CONTINUE
            MOUT = MOUT + NAB( JI, 2 ) - NAB( JI, 1 )
   30    CONTINUE
         RETURN
      END IF
C
C     Initialize for loop
C
C     KF and KL have the following meaning:
C        Intervals 1,...,KF-1 have converged.
C        Intervals KF,...,KL  still need to be refined.
C
      KF = 1
      KL = MINP
C
C     If IJOB=2, initialize C.
C     If IJOB=3, use the user-supplied starting point.
C
      IF( IJOB.EQ.2 ) THEN
         DO 40 JI = 1, MINP
            C( JI ) = HALF*( AB( JI, 1 )+AB( JI, 2 ) )
   40    CONTINUE
      END IF
C
C     Iteration loop
C
      DO 130 JIT = 1, NITMAX
C
C        Loop over intervals
C
         IF( KL-KF+1.GE.NBMIN .AND. NBMIN.GT.0 ) THEN
C
C           Begin of Parallel Version of the loop
C
            DO 60 JI = KF, KL
C
C              Compute N(c), the number of eigenvalues less than c
C
               WORK( JI ) = D( 1 ) - C( JI )
               IWORK( JI ) = 0
               IF( WORK( JI ).LE.PIVMIN ) THEN
                  IWORK( JI ) = 1
                  WORK( JI ) = MIN( WORK( JI ), -PIVMIN )
               END IF
C
               DO 50 J = 2, N
                  WORK( JI ) = D( J ) - E2( J-1 ) / WORK( JI ) - C( JI )
                  IF( WORK( JI ).LE.PIVMIN ) THEN
                     IWORK( JI ) = IWORK( JI ) + 1
                     WORK( JI ) = MIN( WORK( JI ), -PIVMIN )
                  END IF
   50          CONTINUE
   60       CONTINUE
C
            IF( IJOB.LE.2 ) THEN
C
C              IJOB=2: Choose all intervals containing eigenvalues.
C
               KLNEW = KL
               DO 70 JI = KF, KL
C
C                 Insure that N(w) is monotone
C
                  IWORK( JI ) = MIN( NAB( JI, 2 ),
     $                          MAX( NAB( JI, 1 ), IWORK( JI ) ) )
C
C                 Update the Queue -- add intervals if both halves
C                 contain eigenvalues.
C
                  IF( IWORK( JI ).EQ.NAB( JI, 2 ) ) THEN
C
C                    No eigenvalue in the upper interval:
C                    just use the lower interval.
C
                     AB( JI, 2 ) = C( JI )
C
                  ELSE IF( IWORK( JI ).EQ.NAB( JI, 1 ) ) THEN
C
C                    No eigenvalue in the lower interval:
C                    just use the upper interval.
C
                     AB( JI, 1 ) = C( JI )
                  ELSE
                     KLNEW = KLNEW + 1
                     IF( KLNEW.LE.MMAX ) THEN
C
C                       Eigenvalue in both intervals -- add upper to
C                       queue.
C
                        AB( KLNEW, 2 ) = AB( JI, 2 )
                        NAB( KLNEW, 2 ) = NAB( JI, 2 )
                        AB( KLNEW, 1 ) = C( JI )
                        NAB( KLNEW, 1 ) = IWORK( JI )
                        AB( JI, 2 ) = C( JI )
                        NAB( JI, 2 ) = IWORK( JI )
                     ELSE
                        INFO = MMAX + 1
                     END IF
                  END IF
   70          CONTINUE
               IF( INFO.NE.0 )
     $            RETURN
               KL = KLNEW
            ELSE
C
C              IJOB=3: Binary search.  Keep only the interval containing
C                      w   s.t. N(w) = NVAL
C
               DO 80 JI = KF, KL
                  IF( IWORK( JI ).LE.NVAL( JI ) ) THEN
                     AB( JI, 1 ) = C( JI )
                     NAB( JI, 1 ) = IWORK( JI )
                  END IF
                  IF( IWORK( JI ).GE.NVAL( JI ) ) THEN
                     AB( JI, 2 ) = C( JI )
                     NAB( JI, 2 ) = IWORK( JI )
                  END IF
   80          CONTINUE
            END IF
C
         ELSE
C
C           End of Parallel Version of the loop
C
C           Begin of Serial Version of the loop
C
            KLNEW = KL
            DO 100 JI = KF, KL
C
C              Compute N(w), the number of eigenvalues less than w
C
               TMP1 = C( JI )
               TMP2 = D( 1 ) - TMP1
               ITMP1 = 0
               IF( TMP2.LE.PIVMIN ) THEN
                  ITMP1 = 1
                  TMP2 = MIN( TMP2, -PIVMIN )
               END IF
C
               DO 90 J = 2, N
                  TMP2 = D( J ) - E2( J-1 ) / TMP2 - TMP1
                  IF( TMP2.LE.PIVMIN ) THEN
                     ITMP1 = ITMP1 + 1
                     TMP2 = MIN( TMP2, -PIVMIN )
                  END IF
   90          CONTINUE
C
               IF( IJOB.LE.2 ) THEN
C
C                 IJOB=2: Choose all intervals containing eigenvalues.
C
C                 Insure that N(w) is monotone
C
                  ITMP1 = MIN( NAB( JI, 2 ),
     $                    MAX( NAB( JI, 1 ), ITMP1 ) )
C
C                 Update the Queue -- add intervals if both halves
C                 contain eigenvalues.
C
                  IF( ITMP1.EQ.NAB( JI, 2 ) ) THEN
C
C                    No eigenvalue in the upper interval:
C                    just use the lower interval.
C
                     AB( JI, 2 ) = TMP1
C
                  ELSE IF( ITMP1.EQ.NAB( JI, 1 ) ) THEN
C
C                    No eigenvalue in the lower interval:
C                    just use the upper interval.
C
                     AB( JI, 1 ) = TMP1
                  ELSE IF( KLNEW.LT.MMAX ) THEN
C
C                    Eigenvalue in both intervals -- add upper to queue.
C
                     KLNEW = KLNEW + 1
                     AB( KLNEW, 2 ) = AB( JI, 2 )
                     NAB( KLNEW, 2 ) = NAB( JI, 2 )
                     AB( KLNEW, 1 ) = TMP1
                     NAB( KLNEW, 1 ) = ITMP1
                     AB( JI, 2 ) = TMP1
                     NAB( JI, 2 ) = ITMP1
                  ELSE
                     INFO = MMAX + 1
                     RETURN
                  END IF
               ELSE
C
C                 IJOB=3: Binary search.  Keep only the interval
C                         containing  w  s.t. N(w) = NVAL
C
                  IF( ITMP1.LE.NVAL( JI ) ) THEN
                     AB( JI, 1 ) = TMP1
                     NAB( JI, 1 ) = ITMP1
                  END IF
                  IF( ITMP1.GE.NVAL( JI ) ) THEN
                     AB( JI, 2 ) = TMP1
                     NAB( JI, 2 ) = ITMP1
                  END IF
               END IF
  100       CONTINUE
            KL = KLNEW
C
         END IF
C
C        Check for convergence
C
         KFNEW = KF
         DO 110 JI = KF, KL
            TMP1 = ABS( AB( JI, 2 )-AB( JI, 1 ) )
            TMP2 = MAX( ABS( AB( JI, 2 ) ), ABS( AB( JI, 1 ) ) )
            IF( TMP1.LT.MAX( ABSTOL, PIVMIN, RELTOL*TMP2 ) .OR.
     $          NAB( JI, 1 ).GE.NAB( JI, 2 ) ) THEN
C
C              Converged -- Swap with position KFNEW,
C                           then increment KFNEW
C
               IF( JI.GT.KFNEW ) THEN
                  TMP1 = AB( JI, 1 )
                  TMP2 = AB( JI, 2 )
                  ITMP1 = NAB( JI, 1 )
                  ITMP2 = NAB( JI, 2 )
                  AB( JI, 1 ) = AB( KFNEW, 1 )
                  AB( JI, 2 ) = AB( KFNEW, 2 )
                  NAB( JI, 1 ) = NAB( KFNEW, 1 )
                  NAB( JI, 2 ) = NAB( KFNEW, 2 )
                  AB( KFNEW, 1 ) = TMP1
                  AB( KFNEW, 2 ) = TMP2
                  NAB( KFNEW, 1 ) = ITMP1
                  NAB( KFNEW, 2 ) = ITMP2
                  IF( IJOB.EQ.3 ) THEN
                     ITMP1 = NVAL( JI )
                     NVAL( JI ) = NVAL( KFNEW )
                     NVAL( KFNEW ) = ITMP1
                  END IF
               END IF
               KFNEW = KFNEW + 1
            END IF
  110    CONTINUE
         KF = KFNEW
C
C        Choose Midpoints
C
         DO 120 JI = KF, KL
            C( JI ) = HALF*( AB( JI, 1 )+AB( JI, 2 ) )
  120    CONTINUE
C
C        If no more intervals to refine, quit.
C
         IF( KF.GT.KL )
     $      GO TO 140
  130 CONTINUE
C
C     Converged
C
  140 CONTINUE
      INFO = MAX( KL+1-KF, 0 )
      MOUT = KL
C
      RETURN
C
C     End of DLAEBZ
C
      END

      INTEGER FUNCTION DLANEG( N, D, LLD, SIGMA, PIVMIN, R )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            N, R
      DOUBLE PRECISION   PIVMIN, SIGMA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), LLD( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0 )
C     Some architectures propagate Infinities and NaNs very slowly, so
C     the code computes counts in BLKLEN chunks.  Then a NaN can
C     propagate at most BLKLEN columns before being detected.  This is
C     not a general tuning parameter; it needs only to be just large
C     enough that the overhead is tiny in common cases.
      INTEGER BLKLEN
      PARAMETER ( BLKLEN = 128 )
C     ..
C     .. Local Scalars ..
      INTEGER            BJ, J, NEG1, NEG2, NEGCNT
      DOUBLE PRECISION   BSAV, DMINUS, DPLUS, GAMMA, P, T, TMP
      LOGICAL SAWNAN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN, MAX
C     ..
C     .. External Functions ..
      LOGICAL DISNAN
      EXTERNAL DISNAN
C     ..
C     .. Executable Statements ..

      NEGCNT = 0

C     I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
      T = -SIGMA
      DO 210 BJ = 1, R-1, BLKLEN
         NEG1 = 0
         BSAV = T
         DO 21 J = BJ, MIN(BJ+BLKLEN-1, R-1)
            DPLUS = D( J ) + T
            IF( DPLUS.LT.ZERO ) NEG1 = NEG1 + 1
            TMP = T / DPLUS
            T = TMP * LLD( J ) - SIGMA
 21      CONTINUE
         SAWNAN = DISNAN( T )
C     Run a slower version of the above loop if a NaN is detected.
C     A NaN should occur only with a zero pivot after an infinite
C     pivot.  In that case, substituting 1 for T/DPLUS is the
C     correct limit.
         IF( SAWNAN ) THEN
            NEG1 = 0
            T = BSAV
            DO 22 J = BJ, MIN(BJ+BLKLEN-1, R-1)
               DPLUS = D( J ) + T
               IF( DPLUS.LT.ZERO ) NEG1 = NEG1 + 1
               TMP = T / DPLUS
               IF (DISNAN(TMP)) TMP = ONE
               T = TMP * LLD(J) - SIGMA
 22         CONTINUE
         END IF
         NEGCNT = NEGCNT + NEG1
 210  CONTINUE
C
C     II) lower part: L D L^T - SIGMA I = U- D- U-^T
      P = D( N ) - SIGMA
      DO 230 BJ = N-1, R, -BLKLEN
         NEG2 = 0
         BSAV = P
         DO 23 J = BJ, MAX(BJ-BLKLEN+1, R), -1
            DMINUS = LLD( J ) + P
            IF( DMINUS.LT.ZERO ) NEG2 = NEG2 + 1
            TMP = P / DMINUS
            P = TMP * D( J ) - SIGMA
 23      CONTINUE
         SAWNAN = DISNAN( P )
C     As above, run a slower version that substitutes 1 for Inf/Inf.
C
         IF( SAWNAN ) THEN
            NEG2 = 0
            P = BSAV
            DO 24 J = BJ, MAX(BJ-BLKLEN+1, R), -1
               DMINUS = LLD( J ) + P
               IF( DMINUS.LT.ZERO ) NEG2 = NEG2 + 1
               TMP = P / DMINUS
               IF (DISNAN(TMP)) TMP = ONE
               P = TMP * D(J) - SIGMA
 24         CONTINUE
         END IF
         NEGCNT = NEGCNT + NEG2
 230  CONTINUE
C
C     III) Twist index
C       T was shifted by SIGMA initially.
      GAMMA = (T + SIGMA) + P
      IF( GAMMA.LT.ZERO ) NEGCNT = NEGCNT+1

      DLANEG = NEGCNT
      END

      SUBROUTINE DLARNV( IDIST, ISEED, N, X )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            IDIST, N
C     ..
C     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
      INTEGER            LV
      PARAMETER          ( LV = 128 )
      DOUBLE PRECISION   TWOPI
      PARAMETER          ( TWOPI = 6.2831853071795864769252867663D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, IL, IL2, IV
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION   U( LV )
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          COS, LOG, MIN, SQRT
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLARUV
C     ..
C     .. Executable Statements ..
C
      DO 40 IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )
         IF( IDIST.EQ.3 ) THEN
            IL2 = 2*IL
         ELSE
            IL2 = IL
         END IF
C
C        Call DLARUV to generate IL2 numbers from a uniform (0,1)
C        distribution (IL2 <= LV)
C
         CALL DLARUV( ISEED, IL2, U )
C
         IF( IDIST.EQ.1 ) THEN
C
C           Copy generated numbers
C
            DO 10 I = 1, IL
               X( IV+I-1 ) = U( I )
   10       CONTINUE
         ELSE IF( IDIST.EQ.2 ) THEN
C
C           Convert generated numbers to uniform (-1,1) distribution
C
            DO 20 I = 1, IL
               X( IV+I-1 ) = TWO*U( I ) - ONE
   20       CONTINUE
         ELSE IF( IDIST.EQ.3 ) THEN
C
C           Convert generated numbers to normal (0,1) distribution
C
            DO 30 I = 1, IL
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )*
     $                       COS( TWOPI*U( 2*I ) )
   30       CONTINUE
         END IF
   40 CONTINUE
      RETURN
C
C     End of DLARNV
C
      END

      SUBROUTINE DLASQ2( N, Z, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INFO, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   CBIAS
      PARAMETER          ( CBIAS = 1.50D0 )
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO, FOUR, HUNDRD
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0,
     $                     TWO = 2.0D0, FOUR = 4.0D0, HUNDRD = 100.0D0 )
C     ..
C     .. Local Scalars ..
      LOGICAL            IEEE
      INTEGER            I0, I1, I4, IINFO, IPN4, ITER, IWHILA, IWHILB,
     $                   K, KMIN, N0, N1, NBIG, NDIV, NFAIL, PP, SPLT,
     $                   TTYPE
      DOUBLE PRECISION   D, DEE, DEEMIN, DESIG, DMIN, DMIN1, DMIN2, DN,
     $                   DN1, DN2, E, EMAX, EMIN, EPS, G, OLDEMN, QMAX,
     $                   QMIN, S, SAFMIN, SIGMA, T, TAU, TEMP, TOL,
     $                   TOL2, TRACE, ZMAX, TEMPE, TEMPQ
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLASQ3, DLASRT, XERBLA
C     ..
C     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, ILAENV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input arguments.
C     (in case DLASQ2 is not called by DLASQ1)
C
      INFO = 0
      EPS = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      TOL = EPS*HUNDRD
      TOL2 = TOL**2
C
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DLASQ2', 1 )
         RETURN
      ELSE IF( N.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
C
C        1-by-1 case.
C
         IF( Z( 1 ).LT.ZERO ) THEN
            INFO = -201
            CALL XERBLA( 'DLASQ2', 2 )
         END IF
         RETURN
      ELSE IF( N.EQ.2 ) THEN
C
C        2-by-2 case.
C
         IF( Z( 2 ).LT.ZERO .OR. Z( 3 ).LT.ZERO ) THEN
            INFO = -2
            CALL XERBLA( 'DLASQ2', 2 )
            RETURN
         ELSE IF( Z( 3 ).GT.Z( 1 ) ) THEN
            D = Z( 3 )
            Z( 3 ) = Z( 1 )
            Z( 1 ) = D
         END IF
         Z( 5 ) = Z( 1 ) + Z( 2 ) + Z( 3 )
         IF( Z( 2 ).GT.Z( 3 )*TOL2 ) THEN
            T = HALF*( ( Z( 1 )-Z( 3 ) )+Z( 2 ) )
            S = Z( 3 )*( Z( 2 ) / T )
            IF( S.LE.T ) THEN
               S = Z( 3 )*( Z( 2 ) / ( T*( ONE+SQRT( ONE+S / T ) ) ) )
            ELSE
               S = Z( 3 )*( Z( 2 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
            END IF
            T = Z( 1 ) + ( S+Z( 2 ) )
            Z( 3 ) = Z( 3 )*( Z( 1 ) / T )
            Z( 1 ) = T
         END IF
         Z( 2 ) = Z( 3 )
         Z( 6 ) = Z( 2 ) + Z( 1 )
         RETURN
      END IF
C
C     Check for negative data and compute sums of q's and e's.
C
      Z( 2*N ) = ZERO
      EMIN = Z( 2 )
      QMAX = ZERO
      ZMAX = ZERO
      D = ZERO
      E = ZERO
C
      DO 10 K = 1, 2*( N-1 ), 2
         IF( Z( K ).LT.ZERO ) THEN
            INFO = -( 200+K )
            CALL XERBLA( 'DLASQ2', 2 )
            RETURN
         ELSE IF( Z( K+1 ).LT.ZERO ) THEN
            INFO = -( 200+K+1 )
            CALL XERBLA( 'DLASQ2', 2 )
            RETURN
         END IF
         D = D + Z( K )
         E = E + Z( K+1 )
         QMAX = MAX( QMAX, Z( K ) )
         EMIN = MIN( EMIN, Z( K+1 ) )
         ZMAX = MAX( QMAX, ZMAX, Z( K+1 ) )
   10 CONTINUE
      IF( Z( 2*N-1 ).LT.ZERO ) THEN
         INFO = -( 200+2*N-1 )
         CALL XERBLA( 'DLASQ2', 2 )
         RETURN
      END IF
      D = D + Z( 2*N-1 )
      QMAX = MAX( QMAX, Z( 2*N-1 ) )
      ZMAX = MAX( QMAX, ZMAX )
C
C     Check for diagonality.
C
      IF( E.EQ.ZERO ) THEN
         DO 20 K = 2, N
            Z( K ) = Z( 2*K-1 )
   20    CONTINUE
         CALL DLASRT( 'D', N, Z, IINFO )
         Z( 2*N-1 ) = D
         RETURN
      END IF
C
      TRACE = D + E
C
C     Check for zero data.
C
      IF( TRACE.EQ.ZERO ) THEN
         Z( 2*N-1 ) = ZERO
         RETURN
      END IF
C
C     Check whether the machine is IEEE conformable.
C
      IEEE = ILAENV( 10, 'DLASQ2', 'N', 1, 2, 3, 4 ).EQ.1 .AND.
     $       ILAENV( 11, 'DLASQ2', 'N', 1, 2, 3, 4 ).EQ.1
C
C     Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
C
      DO 30 K = 2*N, 2, -2
         Z( 2*K ) = ZERO
         Z( 2*K-1 ) = Z( K )
         Z( 2*K-2 ) = ZERO
         Z( 2*K-3 ) = Z( K-1 )
   30 CONTINUE
C
      I0 = 1
      N0 = N
C
C     Reverse the qd-array, if warranted.
C
      IF( CBIAS*Z( 4*I0-3 ).LT.Z( 4*N0-3 ) ) THEN
         IPN4 = 4*( I0+N0 )
         DO 40 I4 = 4*I0, 2*( I0+N0-1 ), 4
            TEMP = Z( I4-3 )
            Z( I4-3 ) = Z( IPN4-I4-3 )
            Z( IPN4-I4-3 ) = TEMP
            TEMP = Z( I4-1 )
            Z( I4-1 ) = Z( IPN4-I4-5 )
            Z( IPN4-I4-5 ) = TEMP
   40    CONTINUE
      END IF
C
C     Initial split checking via dqd and Li's test.
C
      PP = 0
C
      DO 80 K = 1, 2
C
         D = Z( 4*N0+PP-3 )
         DO 50 I4 = 4*( N0-1 ) + PP, 4*I0 + PP, -4
            IF( Z( I4-1 ).LE.TOL2*D ) THEN
               Z( I4-1 ) = -ZERO
               D = Z( I4-3 )
            ELSE
               D = Z( I4-3 )*( D / ( D+Z( I4-1 ) ) )
            END IF
   50    CONTINUE
C
C        dqd maps Z to ZZ plus Li's test.
C
         EMIN = Z( 4*I0+PP+1 )
         D = Z( 4*I0+PP-3 )
         DO 60 I4 = 4*I0 + PP, 4*( N0-1 ) + PP, 4
            Z( I4-2*PP-2 ) = D + Z( I4-1 )
            IF( Z( I4-1 ).LE.TOL2*D ) THEN
               Z( I4-1 ) = -ZERO
               Z( I4-2*PP-2 ) = D
               Z( I4-2*PP ) = ZERO
               D = Z( I4+1 )
            ELSE IF( SAFMIN*Z( I4+1 ).LT.Z( I4-2*PP-2 ) .AND.
     $               SAFMIN*Z( I4-2*PP-2 ).LT.Z( I4+1 ) ) THEN
               TEMP = Z( I4+1 ) / Z( I4-2*PP-2 )
               Z( I4-2*PP ) = Z( I4-1 )*TEMP
               D = D*TEMP
            ELSE
               Z( I4-2*PP ) = Z( I4+1 )*( Z( I4-1 ) / Z( I4-2*PP-2 ) )
               D = Z( I4+1 )*( D / Z( I4-2*PP-2 ) )
            END IF
            EMIN = MIN( EMIN, Z( I4-2*PP ) )
   60    CONTINUE
         Z( 4*N0-PP-2 ) = D
C
C        Now find qmax.
C
         QMAX = Z( 4*I0-PP-2 )
         DO 70 I4 = 4*I0 - PP + 2, 4*N0 - PP - 2, 4
            QMAX = MAX( QMAX, Z( I4 ) )
   70    CONTINUE
C
C        Prepare for the next iteration on K.
C
         PP = 1 - PP
   80 CONTINUE
C
C     Initialise variables to pass to DLASQ3.
C
      TTYPE = 0
      DMIN1 = ZERO
      DMIN2 = ZERO
      DN    = ZERO
      DN1   = ZERO
      DN2   = ZERO
      G     = ZERO
      TAU   = ZERO
C
      ITER = 2
      NFAIL = 0
      NDIV = 2*( N0-I0 )
C
      DO 160 IWHILA = 1, N + 1
         IF( N0.LT.1 )
     $      GO TO 170
C
C        While array unfinished do
C
C        E(N0) holds the value of SIGMA when submatrix in I0:N0
C        splits from the rest of the array, but is negated.
C
         DESIG = ZERO
         IF( N0.EQ.N ) THEN
            SIGMA = ZERO
         ELSE
            SIGMA = -Z( 4*N0-1 )
         END IF
         IF( SIGMA.LT.ZERO ) THEN
            INFO = 1
            RETURN
         END IF
C
C        Find last unreduced submatrix's top index I0, find QMAX and
C        EMIN. Find Gershgorin-type bound if Q's much greater than E's.
C
         EMAX = ZERO
         IF( N0.GT.I0 ) THEN
            EMIN = ABS( Z( 4*N0-5 ) )
         ELSE
            EMIN = ZERO
         END IF
         QMIN = Z( 4*N0-3 )
         QMAX = QMIN
         DO 90 I4 = 4*N0, 8, -4
            IF( Z( I4-5 ).LE.ZERO )
     $         GO TO 100
            IF( QMIN.GE.FOUR*EMAX ) THEN
               QMIN = MIN( QMIN, Z( I4-3 ) )
               EMAX = MAX( EMAX, Z( I4-5 ) )
            END IF
            QMAX = MAX( QMAX, Z( I4-7 )+Z( I4-5 ) )
            EMIN = MIN( EMIN, Z( I4-5 ) )
   90    CONTINUE
         I4 = 4
C
  100    CONTINUE
         I0 = I4 / 4
         PP = 0
C
         IF( N0-I0.GT.1 ) THEN
            DEE = Z( 4*I0-3 )
            DEEMIN = DEE
            KMIN = I0
            DO 110 I4 = 4*I0+1, 4*N0-3, 4
               DEE = Z( I4 )*( DEE /( DEE+Z( I4-2 ) ) )
               IF( DEE.LE.DEEMIN ) THEN
                  DEEMIN = DEE
                  KMIN = ( I4+3 )/4
               END IF
  110       CONTINUE
            IF( (KMIN-I0)*2.LT.N0-KMIN .AND.
     $         DEEMIN.LE.HALF*Z(4*N0-3) ) THEN
               IPN4 = 4*( I0+N0 )
               PP = 2
               DO 120 I4 = 4*I0, 2*( I0+N0-1 ), 4
                  TEMP = Z( I4-3 )
                  Z( I4-3 ) = Z( IPN4-I4-3 )
                  Z( IPN4-I4-3 ) = TEMP
                  TEMP = Z( I4-2 )
                  Z( I4-2 ) = Z( IPN4-I4-2 )
                  Z( IPN4-I4-2 ) = TEMP
                  TEMP = Z( I4-1 )
                  Z( I4-1 ) = Z( IPN4-I4-5 )
                  Z( IPN4-I4-5 ) = TEMP
                  TEMP = Z( I4 )
                  Z( I4 ) = Z( IPN4-I4-4 )
                  Z( IPN4-I4-4 ) = TEMP
  120          CONTINUE
            END IF
         END IF
C
C        Put -(initial shift) into DMIN.
C
         DMIN = -MAX( ZERO, QMIN-TWO*SQRT( QMIN )*SQRT( EMAX ) )
C
C        Now I0:N0 is unreduced.
C        PP = 0 for ping, PP = 1 for pong.
C        PP = 2 indicates that flipping was applied to the Z array and
C               and that the tests for deflation upon entry in DLASQ3
C               should not be performed.
C
         NBIG = 100*( N0-I0+1 )
         DO 140 IWHILB = 1, NBIG
            IF( I0.GT.N0 )
     $         GO TO 150
C
C           While submatrix unfinished take a good dqds step.
C
            CALL DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL,
     $                   ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1,
     $                   DN2, G, TAU )
C
            PP = 1 - PP
C
C           When EMIN is very small check for splits.
C
            IF( PP.EQ.0 .AND. N0-I0.GE.3 ) THEN
               IF( Z( 4*N0 ).LE.TOL2*QMAX .OR.
     $             Z( 4*N0-1 ).LE.TOL2*SIGMA ) THEN
                  SPLT = I0 - 1
                  QMAX = Z( 4*I0-3 )
                  EMIN = Z( 4*I0-1 )
                  OLDEMN = Z( 4*I0 )
                  DO 130 I4 = 4*I0, 4*( N0-3 ), 4
                     IF( Z( I4 ).LE.TOL2*Z( I4-3 ) .OR.
     $                   Z( I4-1 ).LE.TOL2*SIGMA ) THEN
                        Z( I4-1 ) = -SIGMA
                        SPLT = I4 / 4
                        QMAX = ZERO
                        EMIN = Z( I4+3 )
                        OLDEMN = Z( I4+4 )
                     ELSE
                        QMAX = MAX( QMAX, Z( I4+1 ) )
                        EMIN = MIN( EMIN, Z( I4-1 ) )
                        OLDEMN = MIN( OLDEMN, Z( I4 ) )
                     END IF
  130             CONTINUE
                  Z( 4*N0-1 ) = EMIN
                  Z( 4*N0 ) = OLDEMN
                  I0 = SPLT + 1
               END IF
            END IF
C
  140    CONTINUE
C
         INFO = 2
C
C        Maximum number of iterations exceeded, restore the shift
C        SIGMA and place the new d's and e's in a qd array.
C        This might need to be done for several blocks
C
         I1 = I0
         N1 = N0
 145     CONTINUE
         TEMPQ = Z( 4*I0-3 )
         Z( 4*I0-3 ) = Z( 4*I0-3 ) + SIGMA
         DO K = I0+1, N0
            TEMPE = Z( 4*K-5 )
            Z( 4*K-5 ) = Z( 4*K-5 ) * (TEMPQ / Z( 4*K-7 ))
            TEMPQ = Z( 4*K-3 )
            Z( 4*K-3 ) = Z( 4*K-3 ) + SIGMA + TEMPE - Z( 4*K-5 )
         END DO
C
C        Prepare to do this on the previous block if there is one
C
         IF( I1.GT.1 ) THEN
            N1 = I1-1
            DO WHILE( ( I1.GE.2 ) .AND. ( Z(4*I1-5).GE.ZERO ) )
               I1 = I1 - 1
            END DO
            SIGMA = -Z(4*N1-1)
            GO TO 145
         END IF

         DO K = 1, N
            Z( 2*K-1 ) = Z( 4*K-3 )
C
C        Only the block 1..N0 is unfinished.  The rest of the e's
C        must be essentially zero, although sometimes other data
C        has been stored in them.
C
            IF( K.LT.N0 ) THEN
               Z( 2*K ) = Z( 4*K-1 )
            ELSE
               Z( 2*K ) = 0
            END IF
         END DO
         RETURN
C
C        end IWHILB
C
  150    CONTINUE
C
  160 CONTINUE
C
      INFO = 3
      RETURN
C
C     end IWHILA
C
  170 CONTINUE
C
C     Move q's to the front.
C
      DO 180 K = 2, N
         Z( K ) = Z( 4*K-3 )
  180 CONTINUE
C
C     Sort and compute sum of eigenvalues.
C
      CALL DLASRT( 'D', N, Z, IINFO )
C
      E = ZERO
      DO 190 K = N, 1, -1
         E = E + Z( K )
  190 CONTINUE
C
C     Store trace, sum(eigenvalues) and information on performance.
C
      Z( 2*N+1 ) = TRACE
      Z( 2*N+2 ) = E
      Z( 2*N+3 ) = DBLE( ITER )
      Z( 2*N+4 ) = DBLE( NDIV ) / DBLE( N**2 )
      Z( 2*N+5 ) = HUNDRD*NFAIL / DBLE( ITER )
      RETURN
C
C     End of DLASQ2
C
      END

      SUBROUTINE DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL,
     $                   ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1,
     $                   DN2, G, TAU )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            I0, ITER, N0, NDIV, NFAIL, PP
      DOUBLE PRECISION   DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, G,
     $                   QMAX, SIGMA, TAU
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   CBIAS
      PARAMETER          ( CBIAS = 1.50D0 )
      DOUBLE PRECISION   ZERO, QURTR, HALF, ONE, TWO, HUNDRD
      PARAMETER          ( ZERO = 0.0D0, QURTR = 0.250D0, HALF = 0.5D0,
     $                     ONE = 1.0D0, TWO = 2.0D0, HUNDRD = 100.0D0 )
C     ..
C     .. Local Scalars ..
      INTEGER            IPN4, J4, N0IN, NN, TTYPE
      DOUBLE PRECISION   EPS, S, T, TEMP, TOL, TOL2
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLASQ4, DLASQ5, DLASQ6
C     ..
C     .. External Function ..
      DOUBLE PRECISION   DLAMCH
      LOGICAL            DISNAN
      EXTERNAL           DISNAN, DLAMCH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
C     ..
C     .. Executable Statements ..
C
      N0IN = N0
      EPS = DLAMCH( 'Precision' )
      TOL = EPS*HUNDRD
      TOL2 = TOL**2
C
C     Check for deflation.
C
   10 CONTINUE
C
      IF( N0.LT.I0 )
     $   RETURN
      IF( N0.EQ.I0 )
     $   GO TO 20
      NN = 4*N0 + PP
      IF( N0.EQ.( I0+1 ) )
     $   GO TO 40
C
C     Check whether E(N0-1) is negligible, 1 eigenvalue.
C
      IF( Z( NN-5 ).GT.TOL2*( SIGMA+Z( NN-3 ) ) .AND.
     $    Z( NN-2*PP-4 ).GT.TOL2*Z( NN-7 ) )
     $   GO TO 30
C
   20 CONTINUE
C
      Z( 4*N0-3 ) = Z( 4*N0+PP-3 ) + SIGMA
      N0 = N0 - 1
      GO TO 10
C
C     Check  whether E(N0-2) is negligible, 2 eigenvalues.
C
   30 CONTINUE
C
      IF( Z( NN-9 ).GT.TOL2*SIGMA .AND.
     $    Z( NN-2*PP-8 ).GT.TOL2*Z( NN-11 ) )
     $   GO TO 50
C
   40 CONTINUE
C
      IF( Z( NN-3 ).GT.Z( NN-7 ) ) THEN
         S = Z( NN-3 )
         Z( NN-3 ) = Z( NN-7 )
         Z( NN-7 ) = S
      END IF
      T = HALF*( ( Z( NN-7 )-Z( NN-3 ) )+Z( NN-5 ) )
      IF( Z( NN-5 ).GT.Z( NN-3 )*TOL2.AND.T.NE.ZERO ) THEN
         S = Z( NN-3 )*( Z( NN-5 ) / T )
         IF( S.LE.T ) THEN
            S = Z( NN-3 )*( Z( NN-5 ) /
     $          ( T*( ONE+SQRT( ONE+S / T ) ) ) )
         ELSE
            S = Z( NN-3 )*( Z( NN-5 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
         END IF
         T = Z( NN-7 ) + ( S+Z( NN-5 ) )
         Z( NN-3 ) = Z( NN-3 )*( Z( NN-7 ) / T )
         Z( NN-7 ) = T
      END IF
      Z( 4*N0-7 ) = Z( NN-7 ) + SIGMA
      Z( 4*N0-3 ) = Z( NN-3 ) + SIGMA
      N0 = N0 - 2
      GO TO 10
C
   50 CONTINUE
      IF( PP.EQ.2 )
     $   PP = 0
C
C     Reverse the qd-array, if warranted.
C
      IF( DMIN.LE.ZERO .OR. N0.LT.N0IN ) THEN
         IF( CBIAS*Z( 4*I0+PP-3 ).LT.Z( 4*N0+PP-3 ) ) THEN
            IPN4 = 4*( I0+N0 )
            DO 60 J4 = 4*I0, 2*( I0+N0-1 ), 4
               TEMP = Z( J4-3 )
               Z( J4-3 ) = Z( IPN4-J4-3 )
               Z( IPN4-J4-3 ) = TEMP
               TEMP = Z( J4-2 )
               Z( J4-2 ) = Z( IPN4-J4-2 )
               Z( IPN4-J4-2 ) = TEMP
               TEMP = Z( J4-1 )
               Z( J4-1 ) = Z( IPN4-J4-5 )
               Z( IPN4-J4-5 ) = TEMP
               TEMP = Z( J4 )
               Z( J4 ) = Z( IPN4-J4-4 )
               Z( IPN4-J4-4 ) = TEMP
   60       CONTINUE
            IF( N0-I0.LE.4 ) THEN
               Z( 4*N0+PP-1 ) = Z( 4*I0+PP-1 )
               Z( 4*N0-PP ) = Z( 4*I0-PP )
            END IF
            DMIN2 = MIN( DMIN2, Z( 4*N0+PP-1 ) )
            Z( 4*N0+PP-1 ) = MIN( Z( 4*N0+PP-1 ), Z( 4*I0+PP-1 ),
     $                            Z( 4*I0+PP+3 ) )
            Z( 4*N0-PP ) = MIN( Z( 4*N0-PP ), Z( 4*I0-PP ),
     $                          Z( 4*I0-PP+4 ) )
            QMAX = MAX( QMAX, Z( 4*I0+PP-3 ), Z( 4*I0+PP+1 ) )
            DMIN = -ZERO
         END IF
      END IF
C
C     Choose a shift.
C
      CALL DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1,
     $             DN2, TAU, TTYPE, G )
C
C     Call dqds until DMIN > 0.
C
   70 CONTINUE
C
      CALL DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN,
     $             DN1, DN2, IEEE, EPS )
C
      NDIV = NDIV + ( N0-I0+2 )
      ITER = ITER + 1
C
C     Check status.
C
      IF( DMIN.GE.ZERO .AND. DMIN1.GE.ZERO ) THEN
C
C        Success.
C
         GO TO 90
C
      ELSE IF( DMIN.LT.ZERO .AND. DMIN1.GT.ZERO .AND.
     $         Z( 4*( N0-1 )-PP ).LT.TOL*( SIGMA+DN1 ) .AND.
     $         ABS( DN ).LT.TOL*SIGMA ) THEN
C
C        Convergence hidden by negative DN.
C
         Z( 4*( N0-1 )-PP+2 ) = ZERO
         DMIN = ZERO
         GO TO 90
      ELSE IF( DMIN.LT.ZERO ) THEN
C
C        TAU too big. Select new TAU and try again.
C
         NFAIL = NFAIL + 1
         IF( TTYPE.LT.-22 ) THEN
C
C           Failed twice. Play it safe.
C
            TAU = ZERO
         ELSE IF( DMIN1.GT.ZERO ) THEN
C
C           Late failure. Gives excellent shift.
C
            TAU = ( TAU+DMIN )*( ONE-TWO*EPS )
            TTYPE = TTYPE - 11
         ELSE
C
C           Early failure. Divide by 4.
C
            TAU = QURTR*TAU
            TTYPE = TTYPE - 12
         END IF
         GO TO 70
      ELSE IF( DISNAN( DMIN ) ) THEN
C
C        NaN.
C
         IF( TAU.EQ.ZERO ) THEN
            GO TO 80
         ELSE
            TAU = ZERO
            GO TO 70
         END IF
      ELSE
C
C        Possible underflow. Play it safe.
C
         GO TO 80
      END IF
C
C     Risk of underflow.
C
   80 CONTINUE
      CALL DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DN1, DN2 )
      NDIV = NDIV + ( N0-I0+2 )
      ITER = ITER + 1
      TAU = ZERO
C
   90 CONTINUE
      IF( TAU.LT.SIGMA ) THEN
         DESIG = DESIG + TAU
         T = SIGMA + DESIG
         DESIG = DESIG - ( T-SIGMA )
      ELSE
         T = SIGMA + TAU
         DESIG = SIGMA - ( T-TAU ) + DESIG
      END IF
      SIGMA = T
C
      RETURN
C
C     End of DLASQ3
C
      END

      SUBROUTINE DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN,
     $                   DN1, DN2, TAU, TTYPE, G )
C
C  -- LAPACK computational routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      INTEGER            I0, N0, N0IN, PP, TTYPE
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   CNST1, CNST2, CNST3
      PARAMETER          ( CNST1 = 0.5630D0, CNST2 = 1.010D0,
     $                   CNST3 = 1.050D0 )
      DOUBLE PRECISION   QURTR, THIRD, HALF, ZERO, ONE, TWO, HUNDRD
      PARAMETER          ( QURTR = 0.250D0, THIRD = 0.3330D0,
     $                   HALF = 0.50D0, ZERO = 0.0D0, ONE = 1.0D0,
     $                   TWO = 2.0D0, HUNDRD = 100.0D0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I4, NN, NP
      DOUBLE PRECISION   A2, B1, B2, GAM, GAP1, GAP2, S
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
C     ..
C     .. Executable Statements ..
C
C     A negative DMIN forces the shift to take that absolute value
C     TTYPE records the type of shift.
C
      IF( DMIN.LE.ZERO ) THEN
         TAU = -DMIN
         TTYPE = -1
         RETURN
      END IF
C
      NN = 4*N0 + PP
      IF( N0IN.EQ.N0 ) THEN
C
C        No eigenvalues deflated.
C
         IF( DMIN.EQ.DN .OR. DMIN.EQ.DN1 ) THEN
C
            B1 = SQRT( Z( NN-3 ) )*SQRT( Z( NN-5 ) )
            B2 = SQRT( Z( NN-7 ) )*SQRT( Z( NN-9 ) )
            A2 = Z( NN-7 ) + Z( NN-5 )
C
C           Cases 2 and 3.
C
            IF( DMIN.EQ.DN .AND. DMIN1.EQ.DN1 ) THEN
               GAP2 = DMIN2 - A2 - DMIN2*QURTR
               IF( GAP2.GT.ZERO .AND. GAP2.GT.B2 ) THEN
                  GAP1 = A2 - DN - ( B2 / GAP2 )*B2
               ELSE
                  GAP1 = A2 - DN - ( B1+B2 )
               END IF
               IF( GAP1.GT.ZERO .AND. GAP1.GT.B1 ) THEN
                  S = MAX( DN-( B1 / GAP1 )*B1, HALF*DMIN )
                  TTYPE = -2
               ELSE
                  S = ZERO
                  IF( DN.GT.B1 )
     $               S = DN - B1
                  IF( A2.GT.( B1+B2 ) )
     $               S = MIN( S, A2-( B1+B2 ) )
                  S = MAX( S, THIRD*DMIN )
                  TTYPE = -3
               END IF
            ELSE
C
C              Case 4.
C
               TTYPE = -4
               S = QURTR*DMIN
               IF( DMIN.EQ.DN ) THEN
                  GAM = DN
                  A2 = ZERO
                  IF( Z( NN-5 ) .GT. Z( NN-7 ) )
     $               RETURN
                  B2 = Z( NN-5 ) / Z( NN-7 )
                  NP = NN - 9
               ELSE
                  NP = NN - 2*PP
                  GAM = DN1
                  IF( Z( NP-4 ) .GT. Z( NP-2 ) )
     $               RETURN
                  A2 = Z( NP-4 ) / Z( NP-2 )
                  IF( Z( NN-9 ) .GT. Z( NN-11 ) )
     $               RETURN
                  B2 = Z( NN-9 ) / Z( NN-11 )
                  NP = NN - 13
               END IF
C
C              Approximate contribution to norm squared from I < NN-1.
C
               A2 = A2 + B2
               DO 10 I4 = NP, 4*I0 - 1 + PP, -4
                  IF( B2.EQ.ZERO )
     $               GO TO 20
                  B1 = B2
                  IF( Z( I4 ) .GT. Z( I4-2 ) )
     $               RETURN
                  B2 = B2*( Z( I4 ) / Z( I4-2 ) )
                  A2 = A2 + B2
                  IF( HUNDRD*MAX( B2, B1 ).LT.A2 .OR. CNST1.LT.A2 )
     $               GO TO 20
   10          CONTINUE
   20          CONTINUE
               A2 = CNST3*A2
C
C              Rayleigh quotient residual bound.
C
               IF( A2.LT.CNST1 )
     $            S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )
            END IF
         ELSE IF( DMIN.EQ.DN2 ) THEN
C
C           Case 5.
C
            TTYPE = -5
            S = QURTR*DMIN
C
C           Compute contribution to norm squared from I > NN-2.
C
            NP = NN - 2*PP
            B1 = Z( NP-2 )
            B2 = Z( NP-6 )
            GAM = DN2
            IF( Z( NP-8 ).GT.B2 .OR. Z( NP-4 ).GT.B1 )
     $         RETURN
            A2 = ( Z( NP-8 ) / B2 )*( ONE+Z( NP-4 ) / B1 )
C
C           Approximate contribution to norm squared from I < NN-2.
C
            IF( N0-I0.GT.2 ) THEN
               B2 = Z( NN-13 ) / Z( NN-15 )
               A2 = A2 + B2
               DO 30 I4 = NN - 17, 4*I0 - 1 + PP, -4
                  IF( B2.EQ.ZERO )
     $               GO TO 40
                  B1 = B2
                  IF( Z( I4 ) .GT. Z( I4-2 ) )
     $               RETURN
                  B2 = B2*( Z( I4 ) / Z( I4-2 ) )
                  A2 = A2 + B2
                  IF( HUNDRD*MAX( B2, B1 ).LT.A2 .OR. CNST1.LT.A2 )
     $               GO TO 40
   30          CONTINUE
   40          CONTINUE
               A2 = CNST3*A2
            END IF
C
            IF( A2.LT.CNST1 )
     $         S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )
         ELSE
C
C           Case 6, no information to guide us.
C
            IF( TTYPE.EQ.-6 ) THEN
               G = G + THIRD*( ONE-G )
            ELSE IF( TTYPE.EQ.-18 ) THEN
               G = QURTR*THIRD
            ELSE
               G = QURTR
            END IF
            S = G*DMIN
            TTYPE = -6
         END IF
C
      ELSE IF( N0IN.EQ.( N0+1 ) ) THEN
C
C        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.
C
         IF( DMIN1.EQ.DN1 .AND. DMIN2.EQ.DN2 ) THEN
C
C           Cases 7 and 8.
C
            TTYPE = -7
            S = THIRD*DMIN1
            IF( Z( NN-5 ).GT.Z( NN-7 ) )
     $         RETURN
            B1 = Z( NN-5 ) / Z( NN-7 )
            B2 = B1
            IF( B2.EQ.ZERO )
     $         GO TO 60
            DO 50 I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
               A2 = B1
               IF( Z( I4 ).GT.Z( I4-2 ) )
     $            RETURN
               B1 = B1*( Z( I4 ) / Z( I4-2 ) )
               B2 = B2 + B1
               IF( HUNDRD*MAX( B1, A2 ).LT.B2 )
     $            GO TO 60
   50       CONTINUE
   60       CONTINUE
            B2 = SQRT( CNST3*B2 )
            A2 = DMIN1 / ( ONE+B2**2 )
            GAP2 = HALF*DMIN2 - A2
            IF( GAP2.GT.ZERO .AND. GAP2.GT.B2*A2 ) THEN
               S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
            ELSE
               S = MAX( S, A2*( ONE-CNST2*B2 ) )
               TTYPE = -8
            END IF
         ELSE
C
C           Case 9.
C
            S = QURTR*DMIN1
            IF( DMIN1.EQ.DN1 )
     $         S = HALF*DMIN1
            TTYPE = -9
         END IF
C
      ELSE IF( N0IN.EQ.( N0+2 ) ) THEN
C
C        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.
C
C        Cases 10 and 11.
C
         IF( DMIN2.EQ.DN2 .AND. TWO*Z( NN-5 ).LT.Z( NN-7 ) ) THEN
            TTYPE = -10
            S = THIRD*DMIN2
            IF( Z( NN-5 ).GT.Z( NN-7 ) )
     $         RETURN
            B1 = Z( NN-5 ) / Z( NN-7 )
            B2 = B1
            IF( B2.EQ.ZERO )
     $         GO TO 80
            DO 70 I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
               IF( Z( I4 ).GT.Z( I4-2 ) )
     $            RETURN
               B1 = B1*( Z( I4 ) / Z( I4-2 ) )
               B2 = B2 + B1
               IF( HUNDRD*B1.LT.B2 )
     $            GO TO 80
   70       CONTINUE
   80       CONTINUE
            B2 = SQRT( CNST3*B2 )
            A2 = DMIN2 / ( ONE+B2**2 )
            GAP2 = Z( NN-7 ) + Z( NN-9 ) -
     $             SQRT( Z( NN-11 ) )*SQRT( Z( NN-9 ) ) - A2
            IF( GAP2.GT.ZERO .AND. GAP2.GT.B2*A2 ) THEN
               S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
            ELSE
               S = MAX( S, A2*( ONE-CNST2*B2 ) )
            END IF
         ELSE
            S = QURTR*DMIN2
            TTYPE = -11
         END IF
      ELSE IF( N0IN.GT.( N0+2 ) ) THEN
C
C        Case 12, more than two eigenvalues deflated. No information.
C
         S = ZERO
         TTYPE = -12
      END IF
C
      TAU = S
      RETURN
C
C     End of DLASQ4
C
      END

      SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2,
     $                   DN, DNM1, DNM2, IEEE, EPS )
C
C  -- LAPACK computational routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2017
C
C     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            I0, N0, PP
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU,
     $                   SIGMA, EPS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
C     ..
C
C  =====================================================================
C
C     .. Parameter ..
      DOUBLE PRECISION   ZERO, HALF
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5 )
C     ..
C     .. Local Scalars ..
      INTEGER            J4, J4P2
      DOUBLE PRECISION   D, EMIN, TEMP, DTHRESH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
C
      IF( ( N0-I0-1 ).LE.0 )
     $   RETURN
C
      DTHRESH = EPS*(SIGMA+TAU)
      IF( TAU.LT.DTHRESH*HALF ) TAU = ZERO
      IF( TAU.NE.ZERO ) THEN
      J4 = 4*I0 + PP - 3
      EMIN = Z( J4+4 )
      D = Z( J4 ) - TAU
      DMIN = D
      DMIN1 = -Z( J4 )
C
      IF( IEEE ) THEN
C
C        Code for IEEE arithmetic.
C
         IF( PP.EQ.0 ) THEN
            DO 10 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-2 ) = D + Z( J4-1 )
               TEMP = Z( J4+1 ) / Z( J4-2 )
               D = D*TEMP - TAU
               DMIN = MIN( DMIN, D )
               Z( J4 ) = Z( J4-1 )*TEMP
               EMIN = MIN( Z( J4 ), EMIN )
   10       CONTINUE
         ELSE
            DO 20 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-3 ) = D + Z( J4 )
               TEMP = Z( J4+2 ) / Z( J4-3 )
               D = D*TEMP - TAU
               DMIN = MIN( DMIN, D )
               Z( J4-1 ) = Z( J4 )*TEMP
               EMIN = MIN( Z( J4-1 ), EMIN )
   20       CONTINUE
         END IF
C
C        Unroll last two steps.
C
         DNM2 = D
         DMIN2 = DMIN
         J4 = 4*( N0-2 ) - PP
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM2 + Z( J4P2 )
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
         DMIN = MIN( DMIN, DNM1 )
C
         DMIN1 = DMIN
         J4 = J4 + 4
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM1 + Z( J4P2 )
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
         DMIN = MIN( DMIN, DN )
C
      ELSE
C
C        Code for non IEEE arithmetic.
C
         IF( PP.EQ.0 ) THEN
            DO 30 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-2 ) = D + Z( J4-1 )
               IF( D.LT.ZERO ) THEN
                  RETURN
               ELSE
                  Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
                  D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU
               END IF
               DMIN = MIN( DMIN, D )
               EMIN = MIN( EMIN, Z( J4 ) )
   30       CONTINUE
         ELSE
            DO 40 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-3 ) = D + Z( J4 )
               IF( D.LT.ZERO ) THEN
                  RETURN
               ELSE
                  Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
                  D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU
               END IF
               DMIN = MIN( DMIN, D )
               EMIN = MIN( EMIN, Z( J4-1 ) )
   40       CONTINUE
         END IF
C
C        Unroll last two steps.
C
         DNM2 = D
         DMIN2 = DMIN
         J4 = 4*( N0-2 ) - PP
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM2 + Z( J4P2 )
         IF( DNM2.LT.ZERO ) THEN
            RETURN
         ELSE
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
         END IF
         DMIN = MIN( DMIN, DNM1 )
C
         DMIN1 = DMIN
         J4 = J4 + 4
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM1 + Z( J4P2 )
         IF( DNM1.LT.ZERO ) THEN
            RETURN
         ELSE
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
         END IF
         DMIN = MIN( DMIN, DN )
C
      END IF
      ELSE
C     This is the version that sets d's to zero if they are small enough
         J4 = 4*I0 + PP - 3
         EMIN = Z( J4+4 )
         D = Z( J4 ) - TAU
         DMIN = D
         DMIN1 = -Z( J4 )
         IF( IEEE ) THEN
C
C     Code for IEEE arithmetic.
C
            IF( PP.EQ.0 ) THEN
               DO 50 J4 = 4*I0, 4*( N0-3 ), 4
                  Z( J4-2 ) = D + Z( J4-1 )
                  TEMP = Z( J4+1 ) / Z( J4-2 )
                  D = D*TEMP - TAU
                  IF( D.LT.DTHRESH ) D = ZERO
                  DMIN = MIN( DMIN, D )
                  Z( J4 ) = Z( J4-1 )*TEMP
                  EMIN = MIN( Z( J4 ), EMIN )
 50            CONTINUE
            ELSE
               DO 60 J4 = 4*I0, 4*( N0-3 ), 4
                  Z( J4-3 ) = D + Z( J4 )
                  TEMP = Z( J4+2 ) / Z( J4-3 )
                  D = D*TEMP - TAU
                  IF( D.LT.DTHRESH ) D = ZERO
                  DMIN = MIN( DMIN, D )
                  Z( J4-1 ) = Z( J4 )*TEMP
                  EMIN = MIN( Z( J4-1 ), EMIN )
 60            CONTINUE
            END IF
C
C     Unroll last two steps.
C
            DNM2 = D
            DMIN2 = DMIN
            J4 = 4*( N0-2 ) - PP
            J4P2 = J4 + 2*PP - 1
            Z( J4-2 ) = DNM2 + Z( J4P2 )
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
            DMIN = MIN( DMIN, DNM1 )
C
            DMIN1 = DMIN
            J4 = J4 + 4
            J4P2 = J4 + 2*PP - 1
            Z( J4-2 ) = DNM1 + Z( J4P2 )
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
            DMIN = MIN( DMIN, DN )
C
         ELSE
C
C     Code for non IEEE arithmetic.
C
            IF( PP.EQ.0 ) THEN
               DO 70 J4 = 4*I0, 4*( N0-3 ), 4
                  Z( J4-2 ) = D + Z( J4-1 )
                  IF( D.LT.ZERO ) THEN
                     RETURN
                  ELSE
                     Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
                     D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU
                  END IF
                  IF( D.LT.DTHRESH) D = ZERO
                  DMIN = MIN( DMIN, D )
                  EMIN = MIN( EMIN, Z( J4 ) )
 70            CONTINUE
            ELSE
               DO 80 J4 = 4*I0, 4*( N0-3 ), 4
                  Z( J4-3 ) = D + Z( J4 )
                  IF( D.LT.ZERO ) THEN
                     RETURN
                  ELSE
                     Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
                     D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU
                  END IF
                  IF( D.LT.DTHRESH) D = ZERO
                  DMIN = MIN( DMIN, D )
                  EMIN = MIN( EMIN, Z( J4-1 ) )
 80            CONTINUE
            END IF
C
C     Unroll last two steps.
C
            DNM2 = D
            DMIN2 = DMIN
            J4 = 4*( N0-2 ) - PP
            J4P2 = J4 + 2*PP - 1
            Z( J4-2 ) = DNM2 + Z( J4P2 )
            IF( DNM2.LT.ZERO ) THEN
               RETURN
            ELSE
               Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
               DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
            END IF
            DMIN = MIN( DMIN, DNM1 )
C
            DMIN1 = DMIN
            J4 = J4 + 4
            J4P2 = J4 + 2*PP - 1
            Z( J4-2 ) = DNM1 + Z( J4P2 )
            IF( DNM1.LT.ZERO ) THEN
               RETURN
            ELSE
               Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
               DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
            END IF
            DMIN = MIN( DMIN, DN )
C
         END IF
      END IF
C
      Z( J4+2 ) = DN
      Z( 4*N0-PP ) = EMIN
      RETURN
C
C     End of DLASQ5
C
      END

      SUBROUTINE DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN,
     $                   DNM1, DNM2 )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            I0, N0, PP
      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   Z( * )
C     ..
C
C  =====================================================================
C
C     .. Parameter ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
C     ..
C     .. Local Scalars ..
      INTEGER            J4, J4P2
      DOUBLE PRECISION   D, EMIN, SAFMIN, TEMP
C     ..
C     .. External Function ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
C
      IF( ( N0-I0-1 ).LE.0 )
     $   RETURN
C
      SAFMIN = DLAMCH( 'Safe minimum' )
      J4 = 4*I0 + PP - 3
      EMIN = Z( J4+4 )
      D = Z( J4 )
      DMIN = D
C
      IF( PP.EQ.0 ) THEN
         DO 10 J4 = 4*I0, 4*( N0-3 ), 4
            Z( J4-2 ) = D + Z( J4-1 )
            IF( Z( J4-2 ).EQ.ZERO ) THEN
               Z( J4 ) = ZERO
               D = Z( J4+1 )
               DMIN = D
               EMIN = ZERO
            ELSE IF( SAFMIN*Z( J4+1 ).LT.Z( J4-2 ) .AND.
     $               SAFMIN*Z( J4-2 ).LT.Z( J4+1 ) ) THEN
               TEMP = Z( J4+1 ) / Z( J4-2 )
               Z( J4 ) = Z( J4-1 )*TEMP
               D = D*TEMP
            ELSE
               Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
               D = Z( J4+1 )*( D / Z( J4-2 ) )
            END IF
            DMIN = MIN( DMIN, D )
            EMIN = MIN( EMIN, Z( J4 ) )
   10    CONTINUE
      ELSE
         DO 20 J4 = 4*I0, 4*( N0-3 ), 4
            Z( J4-3 ) = D + Z( J4 )
            IF( Z( J4-3 ).EQ.ZERO ) THEN
               Z( J4-1 ) = ZERO
               D = Z( J4+2 )
               DMIN = D
               EMIN = ZERO
            ELSE IF( SAFMIN*Z( J4+2 ).LT.Z( J4-3 ) .AND.
     $               SAFMIN*Z( J4-3 ).LT.Z( J4+2 ) ) THEN
               TEMP = Z( J4+2 ) / Z( J4-3 )
               Z( J4-1 ) = Z( J4 )*TEMP
               D = D*TEMP
            ELSE
               Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
               D = Z( J4+2 )*( D / Z( J4-3 ) )
            END IF
            DMIN = MIN( DMIN, D )
            EMIN = MIN( EMIN, Z( J4-1 ) )
   20    CONTINUE
      END IF
C
C     Unroll last two steps.
C
      DNM2 = D
      DMIN2 = DMIN
      J4 = 4*( N0-2 ) - PP
      J4P2 = J4 + 2*PP - 1
      Z( J4-2 ) = DNM2 + Z( J4P2 )
      IF( Z( J4-2 ).EQ.ZERO ) THEN
         Z( J4 ) = ZERO
         DNM1 = Z( J4P2+2 )
         DMIN = DNM1
         EMIN = ZERO
      ELSE IF( SAFMIN*Z( J4P2+2 ).LT.Z( J4-2 ) .AND.
     $         SAFMIN*Z( J4-2 ).LT.Z( J4P2+2 ) ) THEN
         TEMP = Z( J4P2+2 ) / Z( J4-2 )
         Z( J4 ) = Z( J4P2 )*TEMP
         DNM1 = DNM2*TEMP
      ELSE
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) )
      END IF
      DMIN = MIN( DMIN, DNM1 )
C
      DMIN1 = DMIN
      J4 = J4 + 4
      J4P2 = J4 + 2*PP - 1
      Z( J4-2 ) = DNM1 + Z( J4P2 )
      IF( Z( J4-2 ).EQ.ZERO ) THEN
         Z( J4 ) = ZERO
         DN = Z( J4P2+2 )
         DMIN = DN
         EMIN = ZERO
      ELSE IF( SAFMIN*Z( J4P2+2 ).LT.Z( J4-2 ) .AND.
     $         SAFMIN*Z( J4-2 ).LT.Z( J4P2+2 ) ) THEN
         TEMP = Z( J4P2+2 ) / Z( J4-2 )
         Z( J4 ) = Z( J4P2 )*TEMP
         DN = DNM1*TEMP
      ELSE
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) )
      END IF
      DMIN = MIN( DMIN, DN )
C
      Z( J4+2 ) = DN
      Z( 4*N0-PP ) = EMIN
      RETURN
C
C     End of DLASQ6
C
      END

      SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
C
C  -- LAPACK driver routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2017
C
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, N, NRHS
C     ..
C     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
C     ..
C
C  =====================================================================
C
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGETRF, ZGETRS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGESV ', -INFO )
         RETURN
      END IF
C
C     Compute the LU factorization of A.
C
      CALL ZGETRF( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN
C
C        Solve the system A*X = B, overwriting B with X.
C
         CALL ZGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
     $                INFO )
      END IF
      RETURN
C
C     End of ZGESV
C
      END

      SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
C     ..
C     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      LOGICAL            NOTRAN
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLASWP, ZTRSM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRS', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
C
      IF( NOTRAN ) THEN
C
C        Solve A * X = B.
C
C        Apply row interchanges to the right hand sides.
C
         CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
C
C        Solve L*X = B, overwriting B with X.
C
         CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
C
C        Solve U*X = B, overwriting B with X.
C
         CALL ZTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
C
C        Solve A**T * X = B  or A**H * X = B.
C
C        Solve U**T *X = B or U**H *X = B, overwriting B with X.
C
         CALL ZTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
C
C        Solve L**T *X = B, or L**H *X = B overwriting B with X.
C
         CALL ZTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,
     $               LDA, B, LDB )
C
C        Apply row interchanges to the solution vectors.
C
         CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
C
      RETURN
C
C     End of ZGETRS
C
      END

      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
C     ..
C     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGETRF2, ZLASWP, ZTRSM
C     ..
C     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRF', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
C
C     Determine the block size for this environment.
C
      NB = ILAENV( 1, 'ZGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
C
C        Use unblocked code.
C
         CALL ZGETRF2( M, N, A, LDA, IPIV, INFO )
      ELSE
C
C        Use blocked code.
C
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
C
C           Factor diagonal and subdiagonal blocks and test for exact
C           singularity.
C
            CALL ZGETRF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
C
C           Adjust INFO and the pivot indices.
C
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
C
C           Apply interchanges to columns 1:J-1.
C
            CALL ZLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
C
            IF( J+JB.LE.N ) THEN
C
C              Apply interchanges to columns J+JB:N.
C
               CALL ZLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
C
C              Compute block row of U.
C
               CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
C
C                 Update trailing submatrix.
C
                  CALL ZGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
C
C     End of ZGETRF
C
      END

      RECURSIVE SUBROUTINE ZGETRF2( M, N, A, LDA, IPIV, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2016
C
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
C     ..
C     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                     ZERO = ( 0.0D+0, 0.0D+0 ) )
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN
      COMPLEX*16         TEMP
      INTEGER            I, IINFO, N1, N2
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            IZAMAX
      EXTERNAL           DLAMCH, IZAMAX
C     ..
C     .. External Subroutines ..
      EXTERNAL           ZGEMM, ZSCAL, ZLASWP, ZTRSM, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRF2', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN

      IF ( M.EQ.1 ) THEN
C
C        Use unblocked code for one row case
C        Just need to handle IPIV and INFO
C
         IPIV( 1 ) = 1
         IF ( A(1,1).EQ.ZERO )
     $      INFO = 1
C
      ELSE IF( N.EQ.1 ) THEN
C
C        Use unblocked code for one column case
C
C
C        Compute machine safe minimum
C
         SFMIN = DLAMCH('S')
C
C        Find pivot and test for singularity
C
         I = IZAMAX( M, A( 1, 1 ), 1 )
         IPIV( 1 ) = I
         IF( A( I, 1 ).NE.ZERO ) THEN
C
C           Apply the interchange
C
            IF( I.NE.1 ) THEN
               TEMP = A( 1, 1 )
               A( 1, 1 ) = A( I, 1 )
               A( I, 1 ) = TEMP
            END IF
C
C           Compute elements 2:M of the column
C
            IF( ABS(A( 1, 1 )) .GE. SFMIN ) THEN
               CALL ZSCAL( M-1, ONE / A( 1, 1 ), A( 2, 1 ), 1 )
            ELSE
               DO 10 I = 1, M-1
                  A( 1+I, 1 ) = A( 1+I, 1 ) / A( 1, 1 )
   10          CONTINUE
            END IF
C
         ELSE
            INFO = 1
         END IF

      ELSE
C
C        Use recursive code
C
         N1 = MIN( M, N ) / 2
         N2 = N-N1
C
C               [ A11 ]
C        Factor [ --- ]
C               [ A21 ]
C
         CALL ZGETRF2( M, N1, A, LDA, IPIV, IINFO )

         IF ( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO
C
C                              [ A12 ]
C        Apply interchanges to [ --- ]
C                              [ A22 ]
C
         CALL ZLASWP( N2, A( 1, N1+1 ), LDA, 1, N1, IPIV, 1 )
C
C        Solve A12
C
         CALL ZTRSM( 'L', 'L', 'N', 'U', N1, N2, ONE, A, LDA,
     $               A( 1, N1+1 ), LDA )
C
C        Update A22
C
         CALL ZGEMM( 'N', 'N', M-N1, N2, N1, -ONE, A( N1+1, 1 ), LDA,
     $               A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA )
C
C        Factor A22
C
         CALL ZGETRF2( M-N1, N2, A( N1+1, N1+1 ), LDA, IPIV( N1+1 ),
     $                 IINFO )
C
C        Adjust INFO and the pivot indices
C
         IF ( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO + N1
         DO 20 I = N1+1, MIN( M, N )
            IPIV( I ) = IPIV( I ) + N1
   20    CONTINUE
C
C        Apply interchanges to A21
C
         CALL ZLASWP( N1, A( 1, 1 ), LDA, N1+1, MIN( M, N), IPIV, 1 )
C
      END IF
      RETURN
C
C     End of ZGETRF2
C
      END

      SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
C
C  -- LAPACK auxiliary routine (version 3.7.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     June 2017
C
C     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
C     ..
C     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
C     ..
C
C =====================================================================
C
C     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      COMPLEX*16         TEMP
C     ..
C     .. Executable Statements ..
C
C     Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
C     K1 through K2.
C
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = K1 + ( K1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
C
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
C
      RETURN
C
C     End of ZLASWP
C
      END

      SUBROUTINE ZLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
C
C  -- LAPACK auxiliary routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   C( * ), S( * )
      COMPLEX*16         A( LDA, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP
      COMPLEX*16         TEMP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLASR ', INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
C
C        Form  P * A
C
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
C
C        Form A * P**T
C
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
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
C
      RETURN
C
C     End of ZLASR
C
      END

      SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
C
C  -- LAPACK computational routine (version 3.7.0) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     December 2016
C
C     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * )
      COMPLEX*16         Z( LDZ, * )
C     ..
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ),
     $                   CONE = ( 1.0D0, 0.0D0 ) )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
C     ..
C     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,
     $                   NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASRT, XERBLA,
     $                   ZLASET, ZLASR, ZSWAP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
C
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSTEQR', -INFO )
         RETURN
      END IF
C
C     Quick return if possible
C
      IF( N.EQ.0 )
     $   RETURN
C
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = CONE
         RETURN
      END IF
C
C     Determine the unit roundoff and over/underflow thresholds.
C
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
C
C     Compute the eigenvalues and eigenvectors of the tridiagonal
C     matrix.
C
      IF( ICOMPZ.EQ.2 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
C
      NMAXIT = N*MAXIT
      JTOT = 0
C
C     Determine where the matrix splits and choose QL or QR iteration
C     for each block, according to whether top or bottom diagonal
C     element is smaller.
C
      L1 = 1
      NM1 = N - 1
C
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
C
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
C
C     Scale submatrix in rows and columns L to LEND
C
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
C
C     Choose between QL and QR iteration
C
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
C
      IF( LEND.GT.L ) THEN
C
C        QL Iteration
C
C        Look for small subdiagonal element.
C
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GO TO 60
   50       CONTINUE
         END IF
C
         M = LEND
C
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
C
C        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
C        to compute its eigensystem.
C
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL ZLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF
C
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
C
C        Form shift.
C
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
C
         S = ONE
         C = ONE
         P = ZERO
C
C        Inner loop
C
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
C
C           If eigenvectors are desired, then save rotations.
C
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
C
   70    CONTINUE
C
C        If eigenvectors are desired, then apply saved rotations.
C
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL ZLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
C
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
C
C        Eigenvalue found.
C
   80    CONTINUE
         D( L ) = P
C
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140
C
      ELSE
C
C        QR Iteration
C
C        Look for small superdiagonal element.
C
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GO TO 110
  100       CONTINUE
         END IF
C
         M = LEND
C
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
C
C        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
C        to compute its eigensystem.
C
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL ZLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF
C
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
C
C        Form shift.
C
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
C
         S = ONE
         C = ONE
         P = ZERO
C
C        Inner loop
C
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
C
C           If eigenvectors are desired, then save rotations.
C
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
C
  120    CONTINUE
C
C        If eigenvectors are desired, then apply saved rotations.
C
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL ZLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
C
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
C
C        Eigenvalue found.
C
  130    CONTINUE
         D( L ) = P
C
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140
C
      END IF
C
C     Undo scaling if necessary
C
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      END IF
C
C     Check for no convergence to an eigenvalue after a total
C     of N*MAXIT iterations.
C
      IF( JTOT.EQ.NMAXIT ) THEN
         DO 150 I = 1, N - 1
            IF( E( I ).NE.ZERO )
     $         INFO = INFO + 1
  150    CONTINUE
         RETURN
      END IF
      GO TO 10
C
C     Order eigenvalues and eigenvectors.
C
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
C
C        Use Quick Sort
C
         CALL DLASRT( 'I', N, D, INFO )
C
      ELSE
C
C        Use Selection Sort to minimize swaps of eigenvectors
C
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
      RETURN
C
C     End of ZSTEQR
C
      END

c-----------------------------------------------------------------------
c---  END
c-----------------------------------------------------------------------