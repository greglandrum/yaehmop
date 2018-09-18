!*==CBORIS.spg  processed by SPAG 6.72Dc at 06:53 on 17 Sep 2018
      SUBROUTINE CBORIS(N,Nd,A,B,C,D,E,F,Fail)
      IMPLICIT NONE
!*--CBORIS4
!*** Start of declarations inserted by SPAG
      REAL*8 A , B , C , D , E , F
      INTEGER i , ia , ii , j , ja , k , lf , N , Nd
!*** End of declarations inserted by SPAG
      DIMENSION A(Nd,Nd) , B(Nd,Nd) , C(Nd,Nd) , D(Nd) , E(Nd) , F(Nd)
      INTEGER Fail
!
!
!  for YAeHMOP this common block is not needed
!              gL
!      COMMON/PASS/NATOM,NDIM,LB,LE,NVAL,ECUT,EERR,IDEL
!
!
! SUBROUTINE CBORIS SOLVES THE COMPLEX EIGENVALUE
! PROBLEM
!
!      A*X=B*X*LAMBDA
!
! WHERE
!
! A     --- GIVEN MATRIX OF ORDER N
! B     --- GIVEN POSITIVE DEFINITE MATRIX OF ORDER N
! X     --- EIGENVECTORS, COLUMN BY COLUMN,
!           NORMALIZED TO (BX,X)=I
! LAMBDA --- EIGENVALUES IN INCREASING ORDER
!
! MATRICES A AND B SHOULD BE GIVEN IN
! ARRAYS A AND B RESPECTIVELY IN THE FOLLOWING WAY
!
! DIAGONAL ELEMENTS IN THE DIAGONAL
! REAL PARTS OF THE LOWER TRIANGLE IN THE LOWER
! TRIANGLE
! IMAG PARTS OF THE LOWER TRIANGLE IN THE UPPER
! TRIANGLE
!
! REAL PARTS OF EIGENVECTORS ARE STORED IN A,
! IMAG PARTS IN ARRAY C.  ARRAY B IS DESTROYED DURING
! COMPUTATION.  IT HOLDS THE LOWER TRIANGLE OF THE
! CHOLESKI DECOMPOSITION OF MATRIX B (SEE
! SUBROUTINE CCHOL).
!
! ALL MATRICES ARE OF ORDER N, DECLARED IN THE
! CALLING PROGRAM WITH DIMENSION ND WHICH NEED
! NOT TO BE EQUAL TO N.
!
! EIGENVALUES ARE STORED IN ARRAY D.  ARRAYS E
! AND F ARE USED FOR INTERMEDIATE RESULTS.
!
! FAIL GETS THE FOLLOWING VALUES
!
!  0 --- COMPUTATION FINISHED SUCCESFULLY
!  1 --- B NOT POSITIVE DEFINITE
!  2 --- QR ALGORITHM DOES NOT CONVERGE
!
! SUBROUTINES USED
!
!     CCHOL
!     CTRED2
!     CTQL2
!
!
! PROGRAMMED BY E. ZAKRAJSEK
! JUNE 21,1974
!
!
!
!
! DECOMPOSE MATRIX B
!
!      WRITE(4,1918)
!1918  FORMAT(' WELCOME TO CBORIS')
!
!      CALL TIME(0)
!
      Fail = 1
      CALL CCHOL(N,Nd,B,lf)
      IF ( lf==0 ) THEN
!
! MOVE MATRIX A
!
         DO i = 1 , N
            C(i,i) = 0.D0
            IF ( i/=1 ) THEN
               ia = i - 1
               DO j = 1 , ia
                  C(i,j) = A(j,i)
                  C(j,i) = -A(j,i)
                  A(j,i) = A(i,j)
               ENDDO
            ENDIF
         ENDDO
!
!  COMPUTE (L(-1)*A)
!
         DO j = 1 , N
            DO i = 1 , N
               IF ( i/=1 ) THEN
                  ia = i - 1
                  DO k = 1 , ia
                     A(i,j) = A(i,j) - A(k,j)*B(i,k) + C(k,j)*B(k,i)
                     C(i,j) = C(i,j) - A(k,j)*B(k,i) - C(k,j)*B(i,k)
                  ENDDO
               ENDIF
               A(i,j) = A(i,j)/B(i,i)
               C(i,j) = C(i,j)/B(i,i)
            ENDDO
         ENDDO
!
!  COMPUTE  A*L(-H)
!
         DO i = 1 , N
            DO j = 1 , i
               IF ( j/=1 ) THEN
                  ja = j - 1
                  DO k = 1 , ja
                     A(i,j) = A(i,j) - A(i,k)*B(j,k) - C(i,k)*B(k,j)
                     C(i,j) = C(i,j) + A(i,k)*B(k,j) - C(i,k)*B(j,k)
                  ENDDO
               ENDIF
               A(i,j) = A(i,j)/B(j,j)
               C(i,j) = C(i,j)/B(j,j)
            ENDDO
         ENDDO
!
!     PUT MATRIX TOGETHER INTO A
!
         DO i = 1 , N
            IF ( i/=N ) THEN
               ia = i + 1
               DO j = ia , N
                  A(i,j) = C(j,i)
               ENDDO
            ENDIF
         ENDDO
!
!     DIAGONALIZE A
!
         Fail = 2
         CALL CTRED2(N,Nd,A,C,D,E,F)
         CALL CTQL2(N,Nd,D,E,F,A,C,lf)
         IF ( lf==0 ) THEN
!
!     COMPUTE L(-H)*A
!
            DO j = 1 , N
               DO ii = 1 , N
                  i = N - ii + 1
                  IF ( i/=N ) THEN
                     ia = i + 1
                     DO k = ia , N
                        A(i,j) = A(i,j) - A(k,j)*B(k,i) - C(k,j)*B(i,k)
                        C(i,j) = C(i,j) + A(k,j)*B(i,k) - C(k,j)*B(k,i)
                     ENDDO
                  ENDIF
                  A(i,j) = A(i,j)/B(i,i)
                  C(i,j) = C(i,j)/B(i,i)
               ENDDO
            ENDDO
            Fail = 0
         ENDIF
      ENDIF
!  999 WRITE(4,1919)
! 1919 FORMAT(' DONE CBORIS')
!
!      CALL TIME(0)
!
      END
