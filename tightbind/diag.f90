!*==CCHOL.spg  processed by SPAG 6.72Dc at 06:53 on 17 Sep 2018
      SUBROUTINE CCHOL(N,Nd,A,Fail)
      IMPLICIT NONE
!*--CCHOL4
!*** Start of declarations inserted by SPAG
      REAL*8 A
      INTEGER i , ia , j , k , ka , N , Nd
!*** End of declarations inserted by SPAG
      DIMENSION A(Nd,Nd)
      INTEGER Fail
!
!
! SUBROUTINE CCHOL COMPUTES CHOLESKI
! DECOMPOSITION OF GIVEN COMPLEX POSITIVE DEFINITE
! MATRIX A.
!
! INPUT DATA
!
!  N --- ORDER OF MATRIX
!  ND -- DIMENSION OF ARRAY A (IT CAN BE
!        GREATER THAN OR EQUAL TO N)
!  A --- GIVEN MATRIX
!        IT IS SUPPOSED TO BE STORED IN THE
!        FOLLOWING WAY.  DIAGONAL ELEMENTS,
!        BEING REAL, ARE STORED ON THE DIAGONAL,
!        REAL PARTS OF OFFDIAGONAL ELEMENTS
!        ARE STORED IN THE LOWER TRIANGLE OF A,
!        IMAG PARTS OF THE LOWER TRIANGLE ARE
!        STORED IN THE UPPER TRIANGLE OF A.
!
!     EXIT INFORMATION
!
!  A --- COMPLEX ELEMENTS OF MATRIX L, DEFINED BY
!     A=L*L(H)
!        ARE STORED IN THE SAME WAY AS ORIGINAL
!        ELEMENTS OF A, THAT IS REAL PARTS OF THE
!        LOWER TRIANGLE OF L IN THE LOWER TRIANGLE
!        OF A AND THE CORRESPONDING IMAG PARTS IN
!        THE UPPER TRIANGLE OF A.
!  FAIL --- IS SET TO ZERO IF THE DECOMPOSITION WAS
!           SUCCESFUL AND TO NONZERO IF
!           THE MATRIX WAS NOT POSITIVE DEFINITE.
!
!
!    PROGRAMMED BY E. ZAKRAJSEK
!     JUNE 20, 1974
!
!
!
!     SUPPOSE DECOMPOSITION WILL FAIL
!
      Fail = 1
!
      DO i = 1 , N
!
!     TEST FOR POSITIVE DEFINITNESS
!
         IF ( A(i,i)<=0.0D0 ) RETURN
!
!      COMPUTE COLUMN I
!
         A(i,i) = DSQRT(A(i,i))
         IF ( i/=N ) THEN
            ia = i + 1
            DO j = ia , N
               A(j,i) = A(j,i)/A(i,i)
               A(i,j) = A(i,j)/A(i,i)
            ENDDO
!
!     REDUCE REMAINING COLUMNS
!
            DO k = ia , N
               A(k,k) = A(k,k) - A(k,i)*A(k,i) - A(i,k)*A(i,k)
               IF ( k/=N ) THEN
                  ka = k + 1
                  DO j = ka , N
                     A(j,k) = A(j,k) - A(j,i)*A(k,i) - A(i,j)*A(i,k)
                     A(k,j) = A(k,j) - A(i,j)*A(k,i) + A(j,i)*A(i,k)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      Fail = 0
      END
!*==CTRED2.spg  processed by SPAG 6.72Dc at 06:53 on 17 Sep 2018
      SUBROUTINE CTRED2(N,Nd,A,B,D,E,F)
      IMPLICIT NONE
!*--CTRED289
!*** Start of declarations inserted by SPAG
      REAL*8 A , ali , all , alr , B , c , chep , D , E , F , g , r ,   &
           & s , sm , t
      INTEGER i , ia , j , k , kk , l , N , Nd
!*** End of declarations inserted by SPAG
      DIMENSION A(Nd,Nd) , B(Nd,Nd) , D(Nd) , E(Nd) , F(Nd)
!
!
!     SUBROUTINE CTRED2 REDUCES GIVEN COMPLEX
!     HERMITIAN MATRIX TO A TRIDIAGONAL FORM
!
!     PARAMETERS
!
!     N    --- ORDER OF THE MATRIX
!     ND   --- DIMENSION OF ARRAYS A AND B
!     A    --- GIVEN MATRIX, REPLACED BY REAL PART
!              OF THE TRANSFORMATION MATRIX
!     B    --- IMAG PART OF TRANSFORMATION MATRIX
!     D    --- DIAGONAL PART OF THE TRIADIAGONAL MATRIX
!     E    --- REAL PART OF THE CODIAGONAL OF THE
!              TRIDIAGONAL MATRIX
!              (LAST N-1 LOCATIONS)
!     F    --- IMAG PARTS OF THE LOWER CODIAGONAL.
!
!     THE GIVEN MATRIX SHOULD BE STORED IN THE
!     FOLLOWING WAY
!
!          --- DIAGONAL ELEMENTS IN THE DIAGONAL
!          --- REAL PART OF THE LOWER TRIANGLE IN THE
!              LOWER TRIANGLE
!          --- IMAG PARTS OF THE LOWER TRIANGLE
!              IN THE UPPER TRIANGLE
!
!
!     PROGRAMMED BY E. ZAKRAJSEK
!     JUNE 20,1974
!
!
!
      chep = 2.0D0**(-56)
      D(1) = A(1,1)
      IF ( N/=1 ) THEN
!
! MAIN K LOOP
!
         DO k = 2 , N
            l = k - 1
!
!     COMPUTE NORM
!
            all = 0.D0
            DO i = k , N
               all = all + A(i,l)*A(i,l) + A(l,i)*A(l,i)
            ENDDO
            all = DSQRT(all)
!
!     COMPUTE CONSTANTS
!
            c = 1.0D0
            s = 0.D0
            r = DSQRT(A(k,l)*A(k,l)+A(l,k)*A(l,k))
            IF ( DABS(r)<1.D-50 ) r = 0.D0
            IF ( r/=0.0D0 ) THEN
               c = A(k,l)/r
               s = A(l,k)/r
            ENDIF
            alr = all*c
            ali = all*s
            A(l,l) = 0.0D0
!
!     TEST FOR SUPERFLUOUS TRANSFORMATION
!
            sm = all*(all+r)
            IF ( DABS(sm)<1.D-50 ) sm = 0.D0
            IF ( sm/=0.D0 ) THEN
               g = 1.0D0/sm
               A(l,l) = g
!
               A(k,l) = A(k,l) + alr
               A(l,k) = A(l,k) + ali
!
!     NOW COMPUTE U=A*W
!     AND STORE INTO (E,F)
!
               t = 0.0D0
               DO i = k , N
                  c = A(i,i)*A(i,l)
                  s = A(i,i)*A(l,i)
                  IF ( i/=k ) THEN
                     ia = i - 1
                     DO j = k , ia
                        c = c + A(i,j)*A(j,l) - A(j,i)*A(l,j)
                        s = s + A(i,j)*A(l,j) + A(j,i)*A(j,l)
                     ENDDO
                  ENDIF
                  IF ( i/=N ) THEN
                     ia = i + 1
                     DO j = ia , N
                        c = c + A(j,i)*A(j,l) + A(i,j)*A(l,j)
                        s = s + A(j,i)*A(l,j) - A(i,j)*A(j,l)
                     ENDDO
                  ENDIF
                  E(i) = g*c
                  F(i) = g*s
                  t = t + A(i,l)*c + A(l,i)*s
               ENDDO
               t = t*(g*g)
!
!    TRANSFORM  MATRIX
!
               DO i = k , N
                  A(i,i) = A(i,i) - 2.0D0*(A(i,l)*E(i)+A(l,i)*F(i))     &
                         & + t*(A(i,l)*A(i,l)+A(l,i)*A(l,i))
                  IF ( i/=k ) THEN
                     ia = i - 1
                     DO j = k , ia
                        A(i,j) = A(i,j) - A(i,l)*E(j) - A(l,i)*F(j)     &
                               & - A(j,l)*E(i) - A(l,j)*F(i)            &
                               & + t*(A(i,l)*A(j,l)+A(l,i)*A(l,j))
                        A(j,i) = A(j,i) - A(l,i)*E(j) + A(i,l)*F(j)     &
                               & + A(l,j)*E(i) - A(j,l)*F(i)            &
                               & + t*(A(l,i)*A(j,l)-A(i,l)*A(l,j))
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
!
!     STORE DIAGONAL AND CODIAGONAL ELEMENTS
!
            D(k) = A(k,k)
            E(k) = -alr
            F(k) = -ali
         ENDDO
      ENDIF
!
!     NOW ACCUMULATE TRANSFORMATIONS
!
      A(N,N) = 1.D0
      B(N,N) = 0.D0
      IF ( N==1 ) RETURN
      DO kk = 2 , N
         k = N - kk + 2
         l = k - 1
!
!     SKIP TRANSFORMATION IF UNIT
!
         IF ( DABS(A(l,l))<1.D-50 ) A(l,l) = 0.D0
         IF ( A(l,l)/=0.0D0 ) THEN
!
!     COMPUTE PRODUCT
!
            DO j = k , N
               c = 0.D0
               s = 0.D0
               DO i = k , N
                  c = c + A(i,l)*A(i,j) + A(l,i)*B(i,j)
                  s = s + A(i,l)*B(i,j) - A(l,i)*A(i,j)
               ENDDO
               c = c*A(l,l)
               s = s*A(l,l)
               DO i = k , N
                  A(i,j) = A(i,j) - c*A(i,l) + s*A(l,i)
                  B(i,j) = B(i,j) - c*A(l,i) - s*A(i,l)
               ENDDO
            ENDDO
         ENDIF
!
!     MAKE NEW LINE
!
         DO i = k , N
            A(i,l) = 0.D0
            A(l,i) = 0.D0
            B(i,l) = 0.D0
            B(l,i) = 0.D0
         ENDDO
         A(l,l) = 1.D0
         B(l,l) = 0.D0
      ENDDO
      END
!*==CTQL2.spg  processed by SPAG 6.72Dc at 06:53 on 17 Sep 2018
      SUBROUTINE CTQL2(N,Nd,D,E,F,A,B,Fail)
      IMPLICIT NONE
!*--CTQL2272
!*** Start of declarations inserted by SPAG
      REAL*8 A , B , bb , c , chep , D , E , F , ff , g , h , hi , hr , &
           & p , r , s
      INTEGER i , i1 , ia , j , k , l , m , ma , N , Nd
!*** End of declarations inserted by SPAG
      DIMENSION A(Nd,Nd) , B(Nd,Nd) , D(Nd) , E(Nd) , F(Nd)
      INTEGER Fail
!
!
!     SUBROUTINE CTQL2 COMPUTES THE EIGENVALUES AND
!     EIGENVECTORS OF A COMPLEX HERMITIAN TRIDIAGONAL
!     MATRIX
!
!     PARAMETERS
!
!     N    --- ORDER OF MATRIX
!     ND   --- DIMENSION OF A AND B
!     D    --- DIAGONAL GIVEN
!     E    --- REAL PART OF CODIAGONAL GIVEN
!              (LAST N-1 LOCATIONS)
!     F    --- IMAG PART OF THE LOWER CODIAGONAL
!     A    --- REAL PART OF EIGENVECTORS
!     B    --- IMAG PART OF EIGENVECTORS
!     FAIL --- RECEIVES VALUE OF 1 INSTEAD OF ZERO
!              IF SOME EIGENVALUE TAKES MORE THAN 30
!              ITERATIONS.
!
!
!     EIGENVALUES ARE OBTAINED IN INCREASING OF
!     MAGNITUDE IN VECTOR D, EIGENVECTORS ARE STORED
!     BY COLUMNS.  ARRAYS A AND B SHOULD BE PRESET TO
!     SOME UNITARY MATRIX SUCH AS THE IDENTITY MATRIX
!     OR THE MATRIX PRODUCED BY CTRED2.
!
!
!     PROGRAMMED BY E.  ZAKRAJSEK
!     JUNE 21, 1974
!
!
!
!
!     ***************************************
!     *                                     *
!     * NEXT LINE OF PROGRAM DEFINES        *
!     * MACHINE DEPENDENT CONSTANT CHEP     *
!     * DEFINED AS THE SMALLEST REAL        *
!     * NUMBER FOR WHICH                    *
!     *                                     *
!     *        1.0+CHEP .GT. 1.0            *
!     *                                     *
!     ***************************************
!
      chep = 2.0D0**(-56)
!
!     FIRST MAKE REAL CODIAGONAL MOVED DOWN
!     TO FIRST LOCATION
!
      IF ( N/=1 ) THEN
         DO k = 2 , N
            r = DSQRT(E(k)*E(k)+F(k)*F(k))
            IF ( DABS(r)<1.D-50 ) r = 0.D0
            IF ( r/=0.0D0 ) THEN
               c = E(k)/r
               s = F(k)/r
!
!     ACCUMULATE ROTATION
!
               DO i = 1 , N
                  p = A(i,k)*c - B(i,k)*s
                  B(i,k) = A(i,k)*s + B(i,k)*c
                  A(i,k) = p
               ENDDO
!
!     TRANSFORM NEXT E
!
               IF ( k/=N ) THEN
                  l = k + 1
                  p = E(l)*c - F(l)*s
                  F(l) = E(l)*s + F(l)*c
                  E(l) = p
               ENDIF
            ENDIF
            E(k-1) = r
         ENDDO
      ENDIF
      E(N) = 0.D0
!
!     INITIALIZE
!
      bb = 0.D0
      ff = 0.D0
      Fail = 1
!
!     MAIN LOOP
!
      DO l = 1 , N
         j = 0
         h = chep*(DABS(D(l))+DABS(E(l)))
         IF ( bb<h ) bb = h
!
!     LOOK FOR SMALL SUBDIAGONAL ELEMENT
!
         DO m = l , N
            IF ( DABS(E(m))<=bb ) GOTO 50
         ENDDO
 50      IF ( m/=l ) THEN
!
!     NEXT ITERATION
!
 60         IF ( j==30 ) RETURN
            j = j + 1
!
!     FORM SHIFT
!
            p = (D(l+1)-D(l))/(2.0D0*E(l))
            r = DSQRT(1.0D0+p*p)
            h = D(l) - E(l)/(p+DSIGN(r,p))
!
            DO i = l , N
               D(i) = D(i) - h
            ENDDO
            ff = ff + h
!
!     QL TRANSFORMATION
!
            p = D(m)
            c = 1.0D0
            s = 0.0D0
            ma = m - 1
!
            DO ia = l , ma
               i = ma - ia + l
               i1 = i + 1
               g = c*E(i)
               h = c*p
               IF ( DABS(p)<DABS(E(i)) ) THEN
!
                  c = p/E(i)
                  r = DSQRT(c*c+1.D0)
                  E(i1) = s*E(i)*r
                  s = 1.0D0/r
                  c = c/r
               ELSE
!
                  c = E(i)/p
                  r = DSQRT(c*c+1.0D0)
                  E(i1) = s*p*r
                  s = c/r
                  c = 1.0D0/r
               ENDIF
!
               p = c*D(i) - s*g
               D(i1) = h + s*(c*g+s*D(i))
!
!     FORM VECTOR
!
               DO k = 1 , N
                  hr = A(k,i1)
                  hi = B(k,i1)
                  A(k,i1) = s*A(k,i) + c*hr
                  B(k,i1) = s*B(k,i) + c*hi
                  A(k,i) = c*A(k,i) - s*hr
                  B(k,i) = c*B(k,i) - s*hi
               ENDDO
!
            ENDDO
!
            E(l) = s*p
            D(l) = c*p
            IF ( DABS(E(l))>bb ) GOTO 60
         ENDIF
!
!     ROOT FOUND
!
         D(l) = D(l) + ff
      ENDDO
!
!     ORDER EIGENVALUES AND EIGENVECTORS
!
      DO i = 1 , N
         k = i
         DO j = i , N
            IF ( D(j)<D(k) ) k = j
         ENDDO
!
!
         IF ( k/=i ) THEN
            p = D(i)
            D(i) = D(k)
            D(k) = p
!
            DO j = 1 , N
               p = A(j,i)
               A(j,i) = A(j,k)
               A(j,k) = p
!
               p = B(j,i)
               B(j,i) = B(j,k)
               B(j,k) = p
            ENDDO
         ENDIF
      ENDDO
      Fail = 0
      END
