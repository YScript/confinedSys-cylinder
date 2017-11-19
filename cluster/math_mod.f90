MODULE math_mod


CONTAINS
!*************************************************************************************************************
!*************************************************************************************************************


!   IMSL ROUTINE     -EIGCH
!--------------------------------------------------------------------
!   COMPUTER         -CDC/SINGLE
!   LATEST REVISION  -JANUARY 1,1978
!   PURPOSE          -EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
!                       A COMPLEX HERMITIAON MATRIX
!   USAGE          -CALL EIGCH (A,N,IJOB,D,Z,IZ,WK,IER)
!   ARGUMENTS A    -INPUT COMPLEX HERMITIAN MATRIX OF ORDER
!                     N,STORED IN HERMITIAN STORAGE DMODE,
!                     WHOSE EIGENVALUES AND EIGENVECTORS ARE
!                     TO BE COMPUTED. INPUT A IS DESTROYED IF
!                     IJOB IS EQUAL TO 0 OR 1.
!                   NOTE - THE ROUTINE TREATS A AS A REAL VECTOR
!                     OF LENGTH N*(N+1).AN APPROPRIATE
!                     EQUIVALENCE STATEMENT MAY BE REQUIRED.
!                     SEE DOCUMENT EXAMPLE.
!             N    -INPUT ORDER OF THE MATRIX A AND MATRIX Z.
!             IJOB -INPUT OPTION PARAMETER,WHEN
!                     IJOB=0,COMPUTE EIGENVALUES ONLY.
!                     IJOB=1,COMPUTE EIGENVALUES AND EIGEN-
!                       VECTORS.
!                     IJOB=2,COMPUTE EIGENVALUES,EIGENVECTORS
!                       AND PERFORMANCE INDEX.
!                     IJOB=3,COMPUTE PERFORMANCE INDEX ONLY.
!                     IF THE PERFORMANCE INDEX IS COMPUTED, IT IS
!                     RETURNED IN WK(1). THE ROUTINES HAVE
!                     PERFORMED (WELL, SATISFACTORILY,POORLY) IF
!                     WK(1) IS (LESS THAN 1, BETWEEN 1 AND 100,
!                     GREATER THEN 100).
!             D   - OUTPUT VECTOR OF LENGTH N CONTAINING THE
!                     EIGENVALUES OF A.
!             Z   - OUTPUT N BY N COMPLEX MATRIX CONTAINING
!                     THE EIGENVECTORS OF A.
!                     THE EIGENVECTOR IN COLUMN J OF Z CORRES-
!                     PONDS TO THE EIGENVALUE D(J).
!                     IF IJOB =0, Z IS NOT USED.
!                   NOTE - THE ROUTINE TREATS Z AS A REAL VECTOR
!                     OF LENGTH 2*N*N. AN APPROPRIATE EQUIVALENCE
!                     STATEMENT MAY BE REQUIRED. SEE DOCUMENT
!                     EXAMPLE
!             IZ  - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
!                     SPECIFIED IN THE DIMENSION STATEMENT IN THE
!                     CALLING PROGRAM. IZ MUST BE GREATER THAN
!                     OR EQUAL TO N IF IJOB IS NOT EQUAL TO ZERO.
!             WK  - WORK AREA, THE LENGTH OF WK DEPENDS
!                     ON THE VALUE OF IJOB, WHEN
!                     IJOB = 0, THE LENGTH OF WK IS AT LEAST 3N.
!                     IJOB = 1, THE LENGTH OF WK IS AT LEAST 3N.
!                     IJOB = 2, THE LENGTH OF WK IS AT LEAST
!                       N*N+4N.
!                     IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
!             IER - ERROR PARAMETER. (OUTPUT)
!                   TERMINAL ERROR
!                     IER = 128+J, INDICATES THAT EQRT2S
!                       FAILED TO CONVERGE ON EIGENVALUE J.
!                       EIGENVALUES J+1,J+2,...,N HAVE BEEN
!                       COMPUTED CORRECTLY.
!                     THE PERFORMANCE INDEX IS SET TO 1000.0.
!                   WARNING ERROR (WITH FIX)
!                     IER = 66, INDICATES IJOB IS LESS THAN 0 OR
!                       IJOB IS GREATER THAN 3. IJOB IS SET TO 1.
!                     IER = 67, INDICATES IJOB IS NOT EQUAL TO
!                       ZERO, AND IZ IS LESS THAN THE ORDER OF
!                       MATRIX A. IJOB IS SET TO ZERO.
!                     IER = 68, INDICATES THAT MATRIX A IS NOT
!                       HERMITIAN BECAUSE SOME DIAGONAL ELEMENT(S)
!                       ARE NOT REAL. EIGCH SETS THE IMAGINARY
!                       PART OF THESE ELEMENTS TO ZERO AND
!                       PROCEEDS WITH THE COMPUTATIONS.
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - EHBCKH,EHDUSH,EQRT2S,UERTST,UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!

      SUBROUTINE EIGCH  (A,N,IJOB,D,Z,IZ,WK,IER)
!                                  SPECIFICATIONS FOR ARGUMENTS
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER            N,IJOB,IZ,IER
      REAL*8             A(1),D(N),Z(1),WK(1)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            JER,K,I,NE,NTAU,NA,NI,NI2,IM1,J,IIZ,NZ,IIZ1, &
     &                   JZ,JZI,L,M,II,IL,KK,LZ,MZ,LK,KZ
      REAL*8             ANORM,ASUM,PI,SUMZ,SUMR,SUMI,S,TEN,RDELP, &
     &                   ZERO,ONE,THOUS,AN,SIGNA
!     EQUIVALENCE
      DATA               RDELP/1.D-30/
      DATA               ZERO,ONE/0.0D0,1.0D0/,TEN/10.0D0/
      DATA               THOUS/1000.0D0/
!                                  INITIALIZE ERROR PARAMETERS
!                                  FIRST EXECUTABLE STATEMENT
!     DO 988 I=1,90
!     WRITE(8,*)'EIGCH A(89)=',A(89),'    A(90)=',A(90)
!988  CONTINUE
      IER = 0
      JER = 0
      IF (IJOB.GE.0.AND.IJOB.LE.3) GO TO 5
!                                  WARNING ERROR -IJOB IS NOT IN THE
!                                    RANGE
      IER = 66
      IJOB = 1
      GO TO 10
  5   IF (IJOB.EQ.0) GO TO 30
 10   IF (IZ.GE.N) GO TO 15
!                          WARNING ERROR - IZ IS LESS THAN N
!                            EIGENVECTORS CAN NOT BE COMPUTED,
!                            IJOB SET TO ZERO
      IER = 67
      IJOB = 0
 15   K = 2
      DO 25 I=1,N
         IF (A(K).EQ.ZERO) GO TO 20
         A(K) = ZERO
!                          WARNING ERROR - SOME DIAGONAL
!                            ELEMENT(S) NOT REAL
      IER = 68
 20   K = K+I+I+2
 25   CONTINUE
      IF (IJOB.EQ.3) GO TO 95
 30   NE = 1
      NTAU = NE+N
      NA = NTAU+N+N
      NI = (N*(N+1))/2
      NI2 = NI+NI
      IF (IJOB.NE.2) GO TO 40
!                          SAVE INPUT A IF IJOB = 2
      K = NA
      DO 35 I=1,NI2
         WK(K) = A(I)
         K = K+1
 35   CONTINUE
!                          SEPARATE A INTO REAL AND IMAGINARY
!                            PARTS
 40   IF (NI.LT.2) GO TO 55
      IM1 = 1
      DO 50 I=2,NI
         K = IM1+I
         PI = A(K)
         DO 45 J=1,IM1
            A(K) = A(K-1)
            K = K-1
 45      CONTINUE
         A(I) = PI
         IM1 = I
 50   CONTINUE
!                          REDUCE HERMITIAN MATRIX TO A REAL
!                            SYMMETRIC TRIDIAGONAL MATRIX
!     WRITE(6,*)'CALL EHOUSH BEFOR'
 55   CALL EHOUSH (A(1),A(NI+1),N,D,WK(NE),WK(NTAU))
!     WRITE(6,*)'CALL EHOUSH AFTER'
      IIZ = 1
      IF (IJOB.NE.0) IIZ = IZ+IZ
      IF (IIZ.EQ.1) GO TO 70
!                          SET Z TO AN IDENTITY MATRIX
      NZ = (IZ+IZ)*N
      DO 60 I=1,NZ
         Z(I) = ZERO
 60   CONTINUE
      K = 1
      IIZ1 = IIZ+1
      DO 65 I=1,N
         Z(K) = ONE
         K = K+IIZ1
 65   CONTINUE
!                          COMPUTE EIGENVALUES AND EIGENVECTORS
 70   CALL EQRT2S (D,WK(NE),N,Z(1),IIZ,JER)
!     WRITE(6,*)'CALL EQRT2S AFTER'
      IF (IJOB.EQ.0) GO TO 9000
!                          BACKTRANSFORM THE EIGENVECTORS
      CALL EHBCKH (A(1),A(NI+1),WK(NTAU),N,Z(1),Z(IZ+1),IIZ)
!     WRITE(6,*)'CALL EHBCKH AFTER'
!                          CONVERT Z (EIGENVECTORS) TO COMPLEX
!                            FORMAT Z(IZ,N)
      JZ = 0
      DO 85 J=1,N
         JZI = JZ+IZ
         DO 75 I=1,N
            K = JZI+I
            WK(I) = Z(K)
 75      CONTINUE
         K = JZ+N
         L = K+N-1
         M = N
         DO 80 I=1,N
            Z(L) = Z(K)
            Z(L+1) = WK(M)
            K = K-1
            L = L-2
            M = M-1
 80      CONTINUE
         JZ = JZ+IZ+IZ
 85   CONTINUE
!                          Z IS NOW IN COMPLEX FORMAT Z(IZ,N)
      IF (IJOB.NE.2) GO TO 9000
!                          MOVE ORIGINAL MATRIX BACK TO A
      K = NA
      DO 90 I=1,NI2
         A(I) = WK(K)
         K = K+1
 90   CONTINUE
      WK(1) = THOUS
      IF (JER.NE.0) GO TO 9000
!                          COMPUTE 1-NORM OF A
 95   ANORM = ZERO
      II = 1
      DO 105 I=1,N
         ASUM = ZERO
         IL = II
         KK = 2
         DO 100 L=1,N
            ASUM = ASUM+CDABS(DCMPLX(A(IL),A(IL+1)))
            IF (L.GE.I) KK = L+L
            IL = IL+KK
 100  CONTINUE
      ANORM = DMAX1(ANORM,ASUM)
      II = II+I+I
 105  CONTINUE
      IF (ANORM.EQ.ZERO) ANORM = ONE
!                          COMPUTE PERFORMANCE INDEX
      PI = ZERO
      DO 120 I=1,N
         II = 1
         S = ZERO
         SUMZ = ZERO
         LZ = (IZ+IZ)*(I-1)+1
         LZ = IZ*(I-1)*2+1
         MZ = LZ
         DO 115 L=1,N
            LK =II
            KK = 2
            SUMZ = SUMZ+CDABS(DCMPLX(Z(LZ),Z(LZ+1)))
            SUMR = -D(I)*Z(LZ)
            SUMI = -D(I)*Z(LZ+1)
            KZ = MZ
            DO 110 K=1,N
               SIGNA = ONE
               IF (K.GT.L) SIGNA = -ONE
               SUMR = SUMR+A(LK)*Z(KZ)-SIGNA*A(LK+1)*Z(KZ+1)
               SUMI = SUMI+A(LK)*Z(KZ+1)+SIGNA*A(LK+1)*Z(KZ)
               IF (K.GE.L) KK = K+K
               LK = LK+KK
               KZ = KZ+2
 110        CONTINUE
            S = S+CDABS(DCMPLX(SUMR,SUMI))
            LZ = LZ+2
            II = II+L+L
 115     CONTINUE
         IF (SUMZ.EQ.ZERO) GO TO 120
         PI = DMAX1(PI,S/SUMZ)
 120  CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WK(1) = PI
 9000 CONTINUE
!     WRITE(6,*)'IER=',IER,'    JER=',JER
      IF (IER.NE.0) CALL UERTST (IER,'EIGCH')
!     WRITE(6,*)'CALL UERTST AFTER'
      IF (JER.EQ.0) GO TO 9005
      IER = JER
      CALL UERTST (IER,'EIGCH' )
!     WRITE(6,*)'CALL UERTST AFTER'
 9005 RETURN
      END  SUBROUTINE EIGCH
!   IMSL ROUTINE NAME   - EHBCKH
!
!----------------------------------------------------------------------
!
!   COMPUTER            - CDC/SINGLE
!
!   LATEST REVISION     - JANUARY 1,  1978
!
!   PURPOSE             - BACK TRANSFORMATION OF THE EIGENVECTORS OF A
!                           REAL SYMMETRIC TRIDIAGONAL MATRIX OBTAINED
!                           FROM THE HOUSEHOLDER REDUCTION OF A
!                           HERMITIAN MATRIX
!
!   USAGE               - CALL EHBCKH (AR,AI,TAU,N,ZR,ZI,IZ)
!
!   ARGUMENTS    AR     - INPUT VECTOR OF LENGTH N*(N+1)/2 CONTAINING
!                           THE DETAILS OF THE HOUSEHOLDER REDUCTION
!                           AS GENERATED BY IMSL ROUTINE EHOUSH.
!                AI     - INPUT VECTOR OF LENGTH N*(N+1)/2 CONTAINING
!                           THE IMAGINARY COUNTERPARTS TO AR, ABOVE.
!                TAU    - INPUT MATRIX OF DIMENSION 2 BY N CONTAINING
!                           THE REMAINING INFORMATION ABOUT UNITARY
!                           TRANSFORMATION (SEE PARAMETER DESCRIPTION
!                           OF TAU IN IMSL ROUTINE EHOUSH)
!                N      - INPUT SCALAR CONTAINING THE ORDER OF THE
!                           HERMITIAN MATRIX.
!                ZR     - INPUT/OUTPUT MATRIX OF DIMENSION N BY N.
!                           ON INPUT CONTAINS THE REAL COMPONENTS
!                           OF THE EIGENVECTORS TO BE BACK TRANSFORMED
!                           AS GENERATED BY IMSL ROUTINE EQRT2S
!                           (SEE PARAMETER DESCRIPTION OF Z IN EQRT2S).
!                           ON OUTPUT CONTAINS THE REAL COMPONENTS
!                           OF THE TRANSFORMED EIGENVECTORS.
!                ZI     - OUTPUT MATRIX OF DIMENSION N BY N
!                           CONTAINING THE IMAGINARY COMPONENTS
!                           OF THE TRANSFORMED EIGENVECTORS.
!                IZ     - INPUT ROW DIMENSION OF MATRICES ZR AND ZI
!                           EXACTLY AS SPECIFIED IN THE DIMENSION
!                           STATEMENT IN THE CALLING PROGRAM.
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTIOPN OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1978 BY IMSL, INC.ALL RIGHTS RESERVED.
!
!   WARRANTY            -IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                          APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                          EXPRESSED OR IMPLIED,IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE EHBCKH (AR,AI,TAU,N,ZR,ZI,IZ)
!                                  SPECIFICATIONS FOR ARGUMENTS
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER          N,IZ
      REAL*8           AR(1),AI(1),TAU(2,1),ZR(IZ,1),ZI(IZ,1)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER          J,K,NR,L,NRM1,INX1,INX2,K1
      REAL*8           DELTA,ZERO,ALPHA1,ALPHA2
      DATA             ZERO/0.0D0/
!                                  TRANSFORM THE EIGENVECTORS OF THE
!                                    REAL SYMMETRIC TRIDIAGONAL MATRIX
!                                    TO THOSE OF THE HERMITIAN TRIDIA-
!                                    GONAL MATRIX
!                                  FIRST EXECUTABLE STATEMENT
      DO 5 J=1,N
         DO 5 K=1,N
            ZI(J,K)=-ZR(J,K)*TAU(2,J)
            ZR(J,K)=ZR(J,K)*TAU(1,J)
   5  CONTINUE
      IF (N .LE. 2) GOTO 30
!                                 RECOVER THE HOUSEHOLDER MATRICES IN
!                                   REVERSE ORDER
      DO 25 L=3,N
         NR=N-L+2
         NRM1=NR-1
         INX1=(NR*(NRM1))/2+NR
         INX2=INX1-1
         IF (AI(INX1) .EQ. ZERO) GOTO 25
         DELTA=AI(INX1)*DSQRT(AR(INX2)**2+AI(INX2)**2)
      DO 20 J=1,N
         ALPHA1=ZERO
         ALPHA2=ZERO
         DO 10 K=NR,N
            K1=(K*(K-1))/2+NRM1
            ALPHA1=ALPHA1+AR(K1)*ZR(K,J)+AI(K1)*ZI(K,J)
            ALPHA2=ALPHA2-AI(K1)*ZR(K,J)+AR(K1)*ZI(K,J)
  10     CONTINUE
         ALPHA1=ALPHA1/DELTA
         ALPHA2=ALPHA2/DELTA
         DO 15 K=NR,N
            K1=(K*(K-1))/2+NRM1
            ZR(K,J)=ZR(K,J)-AR(K1)*ALPHA1+AI(K1)*ALPHA2
            ZI(K,J)=ZI(K,J)-AR(K1)*ALPHA2-AI(K1)*ALPHA1
  15     CONTINUE
  20  CONTINUE
  25  CONTINUE
  30  RETURN
      END SUBROUTINE EHBCKH
!   IMSL ROUTINE NAME   - EHOUSH
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDC/SINGLE
!
!   LATEST REVISION     - JANUARY 1, 1978
!
!   PURPOSE             - REDUCTION OF A COMPLEX HERMITIAN MATRIX TO
!                           REAL SYMMETRIC TRIDIAGONAL FORM
!
!   USAGE               -CALL EHOUSH (AR,AI,N,D,E,TAU)
!
!   ARGUMENTS    AR     - INPUT/OUTPUT VECTOR OF LENGTH N(N+1)/2.
!                           ON INPUT AR CONTAINS THE REAL COMPONENTS 0F
!                           THE FULL LOWER TRIANGLE OF THE HERMITIAN
!                           MATRIX STORED IN SYMMETRIC STORAGE DMODE.
!                           ON OUTPUT AR CONTAINS THE REAL COMPONENTS
!                           OF THE FULL LOWER TRIANGLE OF THE PRODUCT
!                           OF THE UNITARY HERMITIAN MATRICES STORED IN
!                           SYMMETRIC STORAGE MODE.
!                AI     - INPUT/OUTPUT VECTOR OF LENGTH N(N+1)/2
!                           CONTAINING THE IMAGINARY COUNTERPARTS TO
!                           AR, ABOVE.
!                N      - INPUT SCALAR CONTAINING THE ORDER OF THE
!                           HERMITIAN MATRIX.
!                D      - OUTPUT VECTOR OF LENGTH N CONTAINING THE
!                           DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
!                E      - OUTPUT VECTOR OF LENGTH N CONTAINING THE
!                           SUB-DIAGONAL IN ITS LAST (N-1) ELEMENTS.
!                TAU    - OUTPUT MATRIX OF DIMENSION 2 BY N
!                           CONTAINING THE REMAINING INFORMATION
!                           ABOUT THE UNITARY TRANSFORMATIONS.
!                           THE OUTPUT INFORMATION CONTAINED IN TAU
!                           ALONG WITH THE OUTPUT INFORMATION IN AR
!                           AND AI ARE USED TO BACK-TRANSFORM
!                           EIGENVECTORS OF THE SYMMETRIC TRIDIAGONAL
!                           MATRIX TO EIGENVECTORS OF THE ORIGINAL
!                           INPUT MATRIX (SPECIFIED BY THE INPUT VALUES
!                           OF AR AND AI).
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1978 BY IMSL,INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE EHOUSH (AR,AI,N,D,E,TAU)
!                                  SPECIFICATIONS FOR ARGUMENTS
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER            N
      REAL*8             AR(1),AI(1),D(1),E(1),TAU(2,1)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NM1,NN,I,NR,NRM1,L,INDX,J,JJ,INX1,INX2,JP1,KK, &
     &                   IX,IM1
      REAL*8             RHO,TOLER,ZERO,ONE,T1,T2,TESTBB,VR,ROOT,DELTA, &
     &                   RATIO,RDELP,Q1,Q2,X1,X2,TT1,TT2,BB
      DATA               ZERO/0.0D0/,ONE/1.0D0/
      DATA               RDELP/1.D-30/
!                                  FIRST EXECUTABLE STATEMENT
      NM1=N-1
      TOLER=ZERO
      NN=(N*(N+1))/2
      DO 5 I=1,NN
         T1=DABS(AR(I))
         T2=DABS(AI(I))
         IF (T2.GT.T1) T1=T2
         IF (T1.GT.TOLER) TOLER=T1
   5  CONTINUE
      TESTBB=RDELP*TOLER
      IF (N.LE.2) GO TO 65
!                                 PERFORM N - 2 SIMILARITY
!                                   TRANSFORMATIONS
      DO 60 NR=2,NM1
         NRM1=NR-1
         VR=ZERO
         TAU(1,NR)=ZERO
         TAU(2,NR)=ZERO
         TAU(2,1)=ZERO
         DO 10 L=NR,N
            INDX=(L*(L-1))/2+NRM1
            VR=AR(INDX)**2+AI(INDX)**2+VR
!           WRITE(6,*)'EHOUSH VR=',VR
   10    CONTINUE
         INDX=(NR*NRM1)/2+NRM1
         IF ((TESTBB)**2 .GE. VR) GO TO 60
         ROOT = CDABS(DCMPLX(AR(INDX),AI(INDX)))*DSQRT(VR)
!        WRITE(6,*)'EHOUSH ROOT=',ROOT
         IF (ROOT.NE.ZERO) GO TO 15
         AR(INDX)=DSQRT(VR)
         DELTA=VR
         TAU(1,1)=-AR(INDX)
         GO TO 20
   15    DELTA=VR+ROOT
         RATIO=VR/ROOT
         TAU(1,1)=-RATIO*AR(INDX)
         TAU(2,1)= RATIO*AI(INDX)
         AR(INDX)=(RATIO+ONE)*AR(INDX)
         AI(INDX)=(RATIO+ONE)*AI(INDX)
!                                  THE MATRIX TO BE USED IN THE
!                                    SIMILARITY TRANSFORMATION HAS
!                                    BEEN DETERMINED. THE TRANSFOR-
!                                    MATION FOLLOWS
   20    DO 35 J=NR,N
            JJ=(J*(J-1))/2
            INDX=JJ+NRM1
!        WRITE(6,*)'EHOUSH DELTA=',DELTA
            TAU(1,J)=AR(INDX)/DELTA
            TAU(2,J)=AI(INDX)/DELTA
            D(J)=ZERO
            E(J)=ZERO
            DO 25 L=NR,J
               INX1=(L*(L-1))/2+NRM1
               INX2=JJ+L
               D(J)= D(J)+AR(INX2)*AR(INX1)-AI(INX2)*AI(INX1)
               E(J)= E(J)+AR(INX2)*AI(INX1)+AI(INX2)*AR(INX1)
   25       CONTINUE
            JP1=J+1
            IF (JP1 .GT. N) GO TO 40
            DO 30 L=JP1,N
               KK=(L*(L-1))/2
               INX1=KK+NRM1
               INX2=KK+J
               D(J)=D(J)+AR(INX2)*AR(INX1)+AI(INX2)*AI(INX1)
               E(J)=E(J)+AR(INX2)*AI(INX1)-AI(INX2)*AR(INX1)
   30       CONTINUE
   35    CONTINUE
   40    RHO=ZERO
         DO 45 L=NR,N
            RHO=RHO+D(L)*TAU(1,L)+E(L)*TAU(2,L)
   45    CONTINUE
         IX=(NRM1*(NR-2))/2
         DO 55 I=NR,N
            IX=IX+I-1
            INX2=IX+NRM1
            DO 50 J=NR,I
               INX1=IX+J
               X1=TAU(1,I)*D(J)+TAU(2,I)*E(J)
               X2=TAU(2,I)*D(J)-TAU(1,I)*E(J)
               Q1=D(I)-RHO*AR(INX2)
               Q2=E(I)-RHO*AI(INX2)
               T1=Q1*TAU(1,J)+Q2*TAU(2,J)
               T2=Q2*TAU(1,J)-Q1*TAU(2,J)
               AR(INX1)=AR(INX1)-X1-T1
               AI(INX1)=AI(INX1)-X2-T2
   50       CONTINUE
   55    CONTINUE
         TAU(1,NR)=TAU(1,1)
         TAU(2,NR)=TAU(2,1)
   60 CONTINUE
!                                  THE MATRIX HAS BEEN REDUCED TO TRI-
!                                    DIAGONAL HERMITIAN FORM. THE SUB-
!                                    DIAGONAL HAS BEEN TEMPORARILY
!                                    STORED IN VECTOR TAU. STORE THE
!                                    DIAGONAL OF THE REDUCED MATRIX IN D
   65 INDX=0
      DO 70 I=1,N
         INDX=INDX+I
         D(I)=AR(INDX)
   70 CONTINUE
!                                  PERFORM THE DIAGONAL UNITARY SIMILA-
!                                    RITY TRANSFORMATION
      TAU(1,1)=ONE
      TAU(2,1)=ZERO
      E(1)=ZERO
      IF (N .EQ. 1) GO TO 85
      INDX=(N*NM1)/2+NM1
      TAU(1,N)=AR(INDX)
      TAU(2,N)=-AI(INDX)
!                                  CALCLATE SUBDIAGONAL E OF THE REAL
!                                    SYMMETRIC TRIDIAGONAL MATRIX. CAL-
!                                    CULATE TAU, THE DIAGONAL OF THE
!                                    DIAGONAL UNITARY MATRIX
      INDX=1
      DO 80 I=2,N
         INDX=INDX+I
         IM1=I-1
         BB= DSQRT(TAU(1,I)**2+TAU(2,I)**2)
         E(I)=BB
         AI(INDX)=BB
         IF (TESTBB .LT. BB) GO TO 75
         TAU(1,I)=ONE
         TAU(2,I)=ZERO
         BB=ONE
   75 TT1=TAU(1,I)*TAU(1,IM1)-TAU(2,I)*TAU(2,IM1)
         TT2=TAU(1,I)*TAU(2,IM1)+TAU(2,I)*TAU(1,IM1)
!        WRITE(6,*)'EHOUSH BB=',BB
         TAU(1,I)=TT1/BB
         TAU(2,I)=TT2/BB
   80 CONTINUE
   85 RETURN
      END SUBROUTINE EHOUSH
!   IMSL ROUTINE NAME   - EQRT2S
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDC/SINGLE
!
!   LATEST REVISION     - JANUARY 1, 1978
!
!   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
!                           A SYMMETRIC TRIDIAGONAL MATRIX USING THE
!                           QL METHOD.
!
!   USAGE               - CALL EQRT2S (D,E,N,Z,IZ,IER)
!
!   ARGUMENTS    D      - ON INPUT, THE VECTOR D OF LENGTH N CONTAINS
!                           THE DIAGONAL ELEMENTS OF THE SYMMETRIC
!                           TRIDIAGONAL MATRIX T.
!                           ON OUTPUT, D CONTAINS THE EIGENVALUES OF
!                           T IN ASCENDING ORDER.
!                E      - ON INPUT, THE VECTOR E OF LENGTH N CONTAINS
!                           THE SUB-DIAGONAL ELEMENTS OF T IN POSITION
!                           2,...,N. OUTPUT, E IS DESTROYED.
!                N      - ORDER OF TRIDIAGONAL MATRIX T.(INPUT)
!                Z      - ON INPUT, Z CONTAINS THE IDENTITY MATRIX OF
!                           ORDER N.
!                           ON OUTPUT, Z CONTAINS THE EIGENVECTORS
!                           OF T. THE EIGENVECTOR IN COLUMN J OF Z
!                           CORRESPONDS TO THE EIGENVALUE D(J).
!                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
!                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
!                           CALLING PROGRAM. IF IZ IS LESS THAN N, THE
!                           EIGENVECTORS ARE NOT COMPUTED. IN THIS CASE
!                           Z IS NOT USED.
!                IER    - ERROR PARAMETER
!                         TERMINAL ERROR
!                           IER = 128+J, INDICATES THAT EQRT2S FAILED
!                           TO CONVERGE ON EIGENVALUE J. EIGENVALUES
!                           AND EIGENVECTORS 1,...,J-1 HACE BEEN
!                           COMPUTED CORRECTLY, BUT THE EIGENVALUES
!                           ARE UNORDERED.
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - UERTST,UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   REMARKS      IMSL ROUTINE EQRT2S IS DESIGNED TO ACCEPT OUTPUT
!                  VECTORS D AND E FROM IMSL ROUTINE EHOUSS AS INPUT
!                  D AND E OF EQRT2S. GIVEN A SYMMETRIC TRIDIAGONAL
!                  MATRIX, T, VECTOR D CONTAINS THE DIAGONAL ELEMENTS
!                  OF T AND VECTOR E CONTAINS THE SUBDIAGONAL ELEMENTS
!                  OF T. SEE THE EHOUSS DOCUMENT.
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   EARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE EQRT2S (D,E,N,Z,IZ,IER)
      IMPLICIT REAL*8 (A-H,O-Z)
!
      REAL*8             RDELP,R,B,F,H,S,P,C,G
      REAL*8             D(1),E(1),Z(IZ,1)
!     DIMENSION          D(1),E(1),Z(IZ,1)
      DATA               RDELP/1.D-30/
      DATA               ZERO,ONE/0.0D0,1.0D0/
!                                  MOVE THE LAST N-1 ELEMENTS
!                                  OF E INTO THE FIRST N-1 LOCATIONS
!                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (N .EQ. 1) GO TO 9005
      DO 5 I=2,N
         E(I-1) = E(I)
    5 CONTINUE
      E(N) = ZERO
      B = ZERO
      F = ZERO
      DO  60  L=1,N
         J = 0
         H = RDELP*(DABS(D(L))+DABS(E(L)))
         IF (B.LT.H) B = H
!                                  LOOK FOR SMALL SUB-DIAGONAL ELEMENT
         DO 10  M=L,N
            K=M
            IF (DABS(E(K)) .LE. B) GO TO 15
   10    CONTINUE
   15    M = K
         IF (M.EQ.L) GO TO 55
   20    IF (J .EQ. 30) GO TO 85
         J = J+1
         L1 = L+1
         G = D(L)
         P = (D(L1)-G)/(E(L)+E(L))
         R = DSQRT(P*P+ONE)
         D(L) = E(L)/(P+DSIGN(R,P))
         H = G-D(L)
         DO 25 I = L1,N
            D(I) = D(I)-H
   25    CONTINUE
         F = F+H
!                                   QL TRANSFORMATION
         P = D(M)
         C = ONE
         S = ZERO
         MM1 = M-1
         MM1PL = MM1+L
         IF (L.GT.MM1) GO TO 50
         DO 45 II=L,MM1
            I = MM1PL-II
            G = C*E(I)
            H = C*P
            IF (DABS(P).LT.DABS(E(I))) GO TO 30
            C = E(I)/P
            R = DSQRT(C*C+ONE)
            E(I+1) = S*P*R
            S = C/R
            C = ONE/R
            GO TO 35
   30       C = P/E(I)
            R = DSQRT(C*C+ONE)
            E(I+1) = S*E(I)*R
            S = ONE/R
            C = C*S
   35       P = C*D(I)-S*G
            D(I+1) = H+S*(C*G+S*D(I))
            IF (IZ .LT. N) GO TO 45
!                                  FORM VECTOR
            DO 40 K=1,N
               H = Z(K,I+1)
               Z(K,I+1) = S*Z(K,I)+C*H
               Z(K,I) = C*Z(K,I)-S*H
   40       CONTINUE
   45    CONTINUE
   50    E(L) = S*P
         D(L) = C*P
         IF ( DABS(E(L)) .GT.B) GO TO 20
   55    D(L) = D(L) +F
   60 CONTINUE
!                                  ORDER EIGENVALUES AND EIGENVECTORS
      DO  80  I=1,N
         K = I
         P = D(I)
         IP1 = I+1
         IF (IP1.GT.N) GO TO 70
         DO 65  J=IP1,N
            IF (D(J) .GE. P) GO TO 65
            K = J
            P = D(J)
   65    CONTINUE
   70    IF (K.EQ.I) GO TO 80
         D(K) = D(I)
         D(I) = P
         IF (IZ .LT. N) GO TO 80
         DO 75 J = 1,N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
   75    CONTINUE
   80 CONTINUE
      GO TO 9005
   85 IER =128+L
 9000 CONTINUE
      CALL UERTST(IER,'EQRT2S')
 9005 RETURN
      END SUBROUTINE EQRT2S
!   IMSL ROUTINE NAME   - UERTST
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDC/SINGLE
!
!   LATEST REVISION     - JANUARY 1,1978
!
!   PURPOSE             - PRINT A MESSAGE REFLECTING AN ERROR CONDITION
!
!   USAGE               - CALL UERTST (IER,NAME)
!
!   ARGUMENTS    IER    - ERROR PARAMETER. (INPUT)
!                           IER =I+J WHERE
!                             I = 128 IMPLIES TERMINAL ERROR,
!                             I =  64 IMPLIES WARNING WITH FIX, AND
!                             I =  32 IMPLIES WARNING.
!                             J = ERROR CODE RELEVANT TO CALLING
!                                 ROUTINE.
!                NAME   - A SIX CHARACTER LITERAL STRING GIVING THE
!                           NAME OF THE CALLING ROUTINE. (INPUT)
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   REMARKS      THE ERROR MESSAGE PRODUCED BY UERTST IS WRITTEN
!                ONTO THE STANDARD OUTPUT UNIT. THE OUTPUT UNIT
!                NUMBER CAN BE DETERMINED BY CALLING UGETIO AS
!                FOLLOWS..   CALL UGETIO(1,NIN,NOUT).
!                THE OUTPUT UNIT NUMBER CAN BE CHANGED BY CALLING
!                UGETIO AS FOLLOWS..
!                                NIN = 0
!                                NOUT = NEW OUTPUT UNIT NUMBER
!                                CALL UGETIO(3,NIN,NOUT)
!                SEE THE UGETIO DOCUMENT FOR MORE DETAILS.
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE UERTST (IER,NAME)
      IMPLICIT REAL*8 (A-H,O-Z)
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      CHARACTER*8        NAMSET,NAMEQ,IEQ
      CHARACTER(LEN=*), INTENT(IN) :: NAME
      DATA               NAMSET/'UERSET'/
      DATA               NAMEQ/'      '/
!                                  FIRST EXECUTABLE STATEMENT
      DATA               LEVEL/4/,IEQDF/0/,IEQ/'='/
      IF (IER.GT.999) GO TO 25
      IF (IER.LT.-32) GO TO 55
      IF (IER.LE.128) GO TO 5
      IF (LEVEL.LT.1) GO TO 30
!                                  PRINT TERMINAL MESSAGE
      CALL UGETIO(1,NIN,IOUNIT)
      IF (IEQDF.EQ.1) WRITE(IOUNIT,35) IER,NAMEQ,IEQ,NAME
      IF (IEQDF.EQ.0) WRITE(IOUNIT,36) IER,NAME
      GO TO 30
    5 IF (IER.LE.64) GO TO 10
      IF (LEVEL.LT.2) GO TO 30
!                                  PRINT WARNING WITH FIX MESSAGE
      CALL UGETIO(1,NIN,IOUNIT)
      IF (IEQDF.EQ.1) WRITE(IOUNIT,40) IER,NAMEQ,IEQ,NAME
!W    IF (IEQDF.EQ.0) WRITE(IOUNIT,41) IER,NAME
      GO TO 30
   10 IF (IER.LE.32) GO TO 15
!                                  PRINT WARNING MESSAGE
      IF (LEVEL.LT.3) GO TO 30
      CALL UGETIO(1,NIN,IOUNIT)
      IF (IEQDF.EQ.1) WRITE(IOUNIT,45) IER,NAMEQ,IEQ,NAME
      IF (IEQDF.EQ.0) WRITE(IOUNIT,46) IER,NAME
      GO TO 30
   15 CONTINUE
!                                  CHECK FOR UERSET CALL
      IF (NAME.NE.NAMSET) GO TO 25
      LEVOLD = LEVEL
      LEVEL = IER
      IER = LEVOLD
      IF (LEVEL.LT.0) LEVEL = 4
      IF (LEVEL.GT.4) LEVEL = 4
      GO TO 30
   25 CONTINUE
      IF (LEVEL.LT.4) GO TO 30
!                                  PRINT NON-DEFINED MESSAGE
      CALL UGETIO(1,NIN,IOUNIT)
      IF (IEQDF.EQ.1) WRITE(IOUNIT,50) IER,NAMEQ,IEQ,NAME
      IF (IEQDF.EQ.0) WRITE(IOUNIT,51) IER,NAME
   30 IEQDF = 0
      RETURN
   35 FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3, &
     &       20H) FROM IMSL ROUTINE ,A6,A1,A6)
   36 FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3, &
     &       20H) FROM IMSL ROUTINE ,A6)
   40 FORMAT(36H *** WARNING WITH FIX ERROR  (IER = ,I3, &
     &       20H) FROM IMSL ROUTINE ,A6,A1,A6)
   41 FORMAT(36H *** WARNING WITH FIX ERROR  (IER = ,I3, &
     &       20H) FROM IMSL ROUTINE ,A6)
   45 FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3, &
     &       20H) FROM IMSL ROUTINE ,A6,A1,A6)
   46 FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3, &
     &       20H) FROM IMSL ROUTINE ,A6)
   50 FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5, &
     &       20H) FROM IMSL ROUTINE ,A6,A1,A6)
   51 FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5, &
     &       20H) FROM IMSL ROUTINE ,A6)
!                                  SAVE P FOR P = R CASE
!                                    P IS THE PAGE NAME
!                                    R IS THE ROUTINE NAME
   55 IEQDF = 1
      NAMEQ = NAME
   65 RETURN
      END SUBROUTINE UERTST
!   IMSL ROUTINE NAME   - UGETIO
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDC/SINGLE
!
!   LATEST REVISION     - JANUARY 1, 1978
!
!   PURPOSE             - TO RETRIEVE CURRENT VALUES AND TO SET NEW
!                           VALUES FOR INPUT AND OUTPUT UNIT
!                           IDENTIFIERS.
!
!   USAGE               - CALL UGETIO(IOPT,NIN,NOUT)
!
!   ARGUMENTS    IOPT   - OPTION PARAMETER. (INPUT)
!                           IF IOPT=1, THE CURRENT INPUT AND OUTPUT
!                           UNIT IDENTIFIER VALUES ARE RETURNED IN NIN
!                           AND NOUT, RESPECTIVELY.
!                           IF IOPT=2 (3) THE INTERNAL VALUE OF
!                           NIN (NOUT) IS RESET FOR SUBSEQUENT USE.
!                NIN    - INPUIT UNIT IDENTIFIER.
!                           OUTPUT IF IOPT=1, INPUT IF IOPT=2.
!                NOUT   - OUTPUT UNIT IDENTIFIER.
!                           OUTPUT IF IOPT=1, INPUT IF IOPT=3.
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINE - NONE REQUIRED
!
!   NOTATION           - INFORMATION ON SPECIAL NOTATION AND
!                          CONVENTIONS IS AVAILABLE IN THE MANUAL
!                          INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   REMARKS     EACH IMSL ROUTINE THAT PERFORMS INPUT AND/OR OUTPUT
!               OPERATIONS CALLS UGETIO TO OBTAIN THE CURRENT UNIT
!               IDENTIFIER VALUES. IF UGETIO IS CALLED WITH IOPT=2 OR 3
!               NEW UNIT IDENTIFIER VALUES ARE ESTABLISHED. SUBSEQUENT
!               INPUT/OUTPUT IS PERFORMED ON THE NEW UNITS.
!
!   COPYRIGHT          - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY           - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                          APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                          EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE UGETIO(IOPT,NIN,NOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,NIN,NOUT
!                                  SPRCIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NIND,NOUTD
      DATA               NIND/5/,NOUTD/6/
!                                  FIRST EXECUTABLE STATEMENT
      IF (IOPT.EQ.3) GO TO 10
      IF (IOPT.EQ.2) GO TO 5
      IF (IOPT.NE.1) GO TO 9005
      NIN = NIND
      NOUT = NOUTD
      GO TO 9005
    5 NIND = NIN
      GO TO 9005
   10 NOUTD = NOUT
 9005 RETURN
      END SUBROUTINE UGETIO



END MODULE math_mod
