      SUBROUTINE ABFNS(A,B,SK1,SK2,RR,L1,L2,M,N1,N2,MAXCAL)
C
C  ********************************************************************
C  *                                                                  *
C  *     SUBROUTINE ABFNS        CALLED FROM MOV                   *
C  *                                                                  *
C  *                                                                  *
C  *       ABFNS    SUBROUTINE TO CALCULATE THE AB FUNCTIONS          *
C  *                                                                  *
C  *       SUBROUTINES USED:                                          *
C  *                                                                  *
C  *             NONE                                                 *
C  *                                                                  *
C  *                                                                  *
C  *       ORIGIN LOST IN ANTIQUITY                                   *
C
C        modified by greg in modernity (august 1993) so that it doesn't use that
C         damn common block any more.
C
C  *                                                                  *
C  ********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION A(30),B(30)
C     THIS ONLY WORKS FOR PRINCIPAL QUANTUM # < OR = 7
C
C      COMMON/LOCLAP/SK1,SK2,RR,L1,L2,M,N1,N2,MAXCAL
C
      J=MAXCAL+1
      RHO1=0.5D0*(SK1+SK2)*RR
      RHO2=0.5D0*(SK1-SK2)*RR
      IF(DABS(RHO1).GT.165.D0) GO TO 100
      IF(DABS(RHO2).GT.165.D0) GO TO 100
      C=DEXP(-RHO1)
      A(1)=C/RHO1
      DO 15 I=2,J
 15   A(I)=(DBLE(FLOAT(I-1))*A(I-1)+C)/RHO1
      IX=J
      IR=DABS(2.D0*RHO2)
      IS=MIN0(IR+1,19)
      IF(RHO2) 25,35,25
 25   D=DEXP(RHO2)
      H=1.D0/D
C
C  IF THE VALUE OF RHO2 IS TOO SMALL THE SINH MUST BE OBTAINED
C  BY SUMMING THE INFINITE SERIES RATHER THAN BY ADDITION OF
C  TWO EXPONENTIALS.
C
      R=D-H
      IF(DABS(R)-0.1) 26,28,28
 26   RA=RHO2
      RHO22=RHO2*RHO2
      T=RHO2
      DO 27 I=2,50,2
      T=T*RHO22/DBLE(FLOAT(I*I+I))
      RA=RA+T
      IF(T.LT.1.D-30) GO TO 999
 27   CONTINUE
 999  CONTINUE
      R=RA+RA
C
C  AS MANY SUCCESSIVE B-FUNCTIONS ARE GENERATED FROM B-0 BY THE
C  RECURRENCE FORMULAE AS ACCURACY WILL PERMIT.
C
 28   B(1)=R/RHO2
      DO 51 I=2,IX,IS
      IF(IR.EQ.0) GO TO 40
 32   IL=IS-1
      IF(1.GT.IL) GO TO 9050
      DO 31 J=1,IL
      K=I+J-1
      IF((-1)**K) 29,29,30
 29   B(K)=(R+DBLE(FLOAT(K-1))*B(K-1))/RHO2
      GO TO 31
 30   B(K)=-(D+H-DBLE(FLOAT(K-1))*B(K-1))/RHO2
 31   CONTINUE
 9050 CONTINUE
 40   IN=I+IS-1
C
C  AFTER THE RECURRENCE FORMULAE HAVE BEEN APPLIED AN APPROPRIATE
C  NUMBER OF TIMES THE NEXT B-FUNCTION IS OBTAINED BY SUMMATION
C  OF THE INFINITE SERIES.
C
      IF(IN-IX) 39,39,38
 39   IF((-1)**IN) 44,44,42
 42   TR=RHO2
105   B(IN)=-2.D0*TR/DBLE(FLOAT(IN+1))
      DO 43 J=1,500
      TR=TR*RHO2**2/DBLE(FLOAT((2*J)*(2*J+1)))
      IF(DABS(TR/B(IN))-1.0D-7 ) 51,51,43
 43   B(IN)=B(IN)-2.D0*TR/DBLE(FLOAT(IN+1+2*J))
      GOTO 51
 44   TR=1.
 107  B(IN)=2.D0*TR/DBLE(FLOAT(IN))
      DO 46 J=1,500
      TR=TR*RHO2**2/DBLE(FLOAT((2*J)*(2*J-1)))
      IF(DABS(TR/B(IN))-1.0D-7 ) 51,51,46
 46   B(IN)=B(IN)+2.D0*TR/DBLE(FLOAT(IN+2*J))
 51   CONTINUE
C
C  IF THE ARGUMENT IS ZERO A SEPARATE FORMULA MUST BE USED.
C
      GO TO 38
 35   DO 36 I=1,IX,2
      B(I)=2.D0/DBLE(FLOAT(I))
 36   B(I+1)=0.D0
 38   RETURN
 100  DO 101 I=1,20
      A(I)=0.D0
 101  B(I)=0.D0
      GO TO 38
      END
