!*==ABFNS.spg  processed by SPAG 6.72Dc at 06:55 on 17 Sep 2018
      SUBROUTINE ABFNS(A,B,Sk1,Sk2,Rr,L1,L2,M,N1,N2,Maxcal)
!
!  ********************************************************************
!  *                                                                  *
!  *     SUBROUTINE ABFNS        CALLED FROM MOV                   *
!  *                                                                  *
!  *                                                                  *
!  *       ABFNS    SUBROUTINE TO CALCULATE THE AB FUNCTIONS          *
!  *                                                                  *
!  *       SUBROUTINES USED:                                          *
!  *                                                                  *
!  *             NONE                                                 *
!  *                                                                  *
!  *                                                                  *
!  *       ORIGIN LOST IN ANTIQUITY                                   *
!
!        modified by greg in modernity (august 1993) so that it doesn't use that
!         damn common block any more.
!
!  *                                                                  *
!  ********************************************************************
!
      IMPLICIT NONE
!*--ABFNS25
!*** Start of declarations inserted by SPAG
      REAL*8 A , B , c , d , h , r , ra , rho1 , rho2 , rho22 , Rr ,    &
           & Sk1 , Sk2 , t , tr
      INTEGER i , il , in , ir , is , ix , j , k , L1 , L2 , M ,        &
            & Maxcal , N1 , N2
!*** End of declarations inserted by SPAG
!
      DIMENSION A(30) , B(30)
!     THIS ONLY WORKS FOR PRINCIPAL QUANTUM # < OR = 7
!
!      COMMON/LOCLAP/SK1,SK2,RR,L1,L2,M,N1,N2,MAXCAL
!
      j = Maxcal + 1
      rho1 = 0.5D0*(Sk1+Sk2)*Rr
      rho2 = 0.5D0*(Sk1-Sk2)*Rr
      IF ( DABS(rho1)>165.D0 ) GOTO 200
      IF ( DABS(rho2)>165.D0 ) GOTO 200
      c = DEXP(-rho1)
      A(1) = c/rho1
      DO i = 2 , j
         A(i) = (DBLE(FLOAT(i-1))*A(i-1)+c)/rho1
      ENDDO
      ix = j
      ir = DABS(2.D0*rho2)
      is = MIN0(ir+1,19)
      IF ( rho2/=0 ) THEN
         d = DEXP(rho2)
         h = 1.D0/d
!
!  IF THE VALUE OF RHO2 IS TOO SMALL THE SINH MUST BE OBTAINED
!  BY SUMMING THE INFINITE SERIES RATHER THAN BY ADDITION OF
!  TWO EXPONENTIALS.
!
         r = d - h
         IF ( DABS(r)<0.1 ) THEN
            ra = rho2
            rho22 = rho2*rho2
            t = rho2
            DO i = 2 , 50 , 2
               t = t*rho22/DBLE(FLOAT(i*i+i))
               ra = ra + t
               IF ( t<1.D-30 ) GOTO 20
            ENDDO
 20         r = ra + ra
         ENDIF
!
!  AS MANY SUCCESSIVE B-FUNCTIONS ARE GENERATED FROM B-0 BY THE
!  RECURRENCE FORMULAE AS ACCURACY WILL PERMIT.
!
         B(1) = r/rho2
         DO i = 2 , ix , is
            IF ( ir/=0 ) THEN
               il = is - 1
               IF ( 1<=il ) THEN
                  DO j = 1 , il
                     k = i + j - 1
                     IF ( (-1)**k<=0 ) THEN
                        B(k) = (r+DBLE(FLOAT(k-1))*B(k-1))/rho2
                     ELSE
                        B(k) = -(d+h-DBLE(FLOAT(k-1))*B(k-1))/rho2
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
            in = i + is - 1
!
!  AFTER THE RECURRENCE FORMULAE HAVE BEEN APPLIED AN APPROPRIATE
!  NUMBER OF TIMES THE NEXT B-FUNCTION IS OBTAINED BY SUMMATION
!  OF THE INFINITE SERIES.
!
            IF ( in>ix ) GOTO 100
            IF ( (-1)**in<=0 ) THEN
               tr = 1.
               B(in) = 2.D0*tr/DBLE(FLOAT(in))
               DO j = 1 , 500
                  tr = tr*rho2**2/DBLE(FLOAT((2*j)*(2*j-1)))
                  IF ( DABS(tr/B(in))<=1.0D-7 ) GOTO 50
                  B(in) = B(in) + 2.D0*tr/DBLE(FLOAT(in+2*j))
               ENDDO
            ELSE
               tr = rho2
               B(in) = -2.D0*tr/DBLE(FLOAT(in+1))
               DO j = 1 , 500
                  tr = tr*rho2**2/DBLE(FLOAT((2*j)*(2*j+1)))
                  IF ( DABS(tr/B(in))<=1.0D-7 ) GOTO 50
                  B(in) = B(in) - 2.D0*tr/DBLE(FLOAT(in+1+2*j))
               ENDDO
            ENDIF
!
!  IF THE ARGUMENT IS ZERO A SEPARATE FORMULA MUST BE USED.
!
 50      ENDDO
      ELSE
         DO i = 1 , ix , 2
            B(i) = 2.D0/DBLE(FLOAT(i))
            B(i+1) = 0.D0
         ENDDO
      ENDIF
 100  RETURN
 200  DO i = 1 , 20
         A(i) = 0.D0
         B(i) = 0.D0
      ENDDO
      GOTO 100
      END
