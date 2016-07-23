      SUBROUTINE LOVLAP (STRAD,A,B,SK1,SK2,R,L1,L2,M1,N1,N2,MAX)
C
C  ********************************************************************
C  *                                                                  *
C  *     SUBROUTINE LOVLAP       CALLED FROM MOV                   *
C  *                                                                  *
C  *                                                                  *
C  *       LOVLAP   SUBROUTINE TO CALCULATE THE OVERLAP               *
C  *                COMPONENT INDEPENDENT OF THE ANGLE BETWEEN        *
C  *                THE ATOMS                                         *
C  *                                                                  *
C  *       SUBROUTINES USED:                                          *
C  *                                                                  *
C  *             NONE                                                 *
C  *                                                                  *
C  *                                                                  *
C  *       ORIGIN LOST IN ANTIQUITY                                   *
C
C    modified by greg in modernity (august 1993) so that it doesn't use that
C         damn common block any more.
C
C  *                                                                  *
C  ********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION A(30),B(30)
      DIMENSION FACT(25)
      DIMENSION BINCOE(8,8)
C
C      COMMON/LOCLAP/SK1,SK2,R,L1,L2,M1,N1,N2,MAX
C
      DATA BINCOE/8*1.D0,   0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,
     1 2*0.D0,1.D0,3.D0,6.D0,10.D0,15.D0,21.D0,   3*0.D0,1.D0,4.D0,
     2 10.D0,20.D0,35.D0,   4*0.D0,1.D0,5.D0,15.D0,35.D0,   5*0.D0,
     3 1.D0,6.D0,21.D0,   6*0.D0,1.D0,7.D0,   7*0.D0,1.D0/
C

C      write (*,*) 'Lovlap: ',sk1,sk2,r,l1,l2,m1,n1,n2,max

      MAXXA=1
      MAXXB=1
      MAXXC=1
      FACT(1)=1.D0
      DO 5 I=2,25
 5    FACT(I)=FACT(I-1)*DBLE(FLOAT(I-1))
      M2=M1
      STRAD=0.D0
      RHOA=R*SK1
      RHOB=R*SK2
      RHOAP=RHOA**N1
      RHOAP=RHOAP*RHOAP
      RHOAP=RHOAP*RHOA
      RHOAB=RHOB**N2
      RHOAB=RHOAB*RHOAB
      RHOAB=RHOAB*RHOB
      RHOPO=RHOAP*RHOAB
      TERMA=0.5D0**(L1+L2+1)*DSQRT(DBLE(FLOAT((L1+L1+1)*(L2+L2+1)))*
     1 FACT(L1-M1+1)*FACT(L2-M1+1)/(FACT(N1+N1+1)*FACT(N2+N2+1)*
     2 FACT(L1+M1+1)*FACT(L2+M1+1))*RHOPO)
      JEND=1+((L1-M1)/2)
      KEND=1+((L2-M2)/2)
      IEB=M1+1
      DO 50 J=1,JEND
      JU=J-1
      IAB=N1-L1+JU+JU+1
      ICB=L1-M1-JU-JU+1
      CON1=FACT(L1+L1-JU-JU+1)/(FACT(L1-M1-JU-JU+1)*FACT(JU+1)*
     1 FACT(L1-JU+1))
      DO 50 K=1,KEND
      KU=K-1
      CON12=CON1*FACT(L2+L2-KU-KU+1)/(FACT(L2-M2-KU-KU+1)*FACT(KU+1)*
     1 FACT(L2-KU+1))
      IEV=JU+KU+L2
      IF(2*(IEV/2).NE.IEV) CON12=-CON12
      IBB=N2-L2+KU+KU+1
      IDB=L2-M2-KU-KU+1
      VALUE=0.D0
      DO 90 I6=1,IEB
      DO 90 I5=1,IEB
      VALUE1=BINCOE(IEB,I6)*BINCOE(IEB,I5)
      IEV=I5+I6
      IF(2*(IEV/2).NE.IEV) VALUE1=-VALUE1
      DO 90 I4=1,IDB
      VALUE1=-VALUE1
      VALUE2=BINCOE(IDB,I4)*VALUE1
      DO 90 I3=1,ICB
      VALUE3=BINCOE(ICB,I3)*VALUE2
      DO 90 I2=1,IBB
      VALUE3=-VALUE3
      VALUE4=BINCOE(IBB,I2)*VALUE3
      DO 90 I1=1,IAB
      TERM=VALUE4*BINCOE(IAB,I1)
      IR=I1+I2+IEB+IEB-I6-I6-I3+IDB-I4+ICB-1
      IP=IAB-I1+IBB-I2+IEB+IEB-I5-I5+ICB-I3+IDB-I4+1
90    VALUE=VALUE+A(IP)*B(IR)*TERM
50    STRAD=STRAD+VALUE*CON12
      STRAD=STRAD*TERMA

      RETURN
      END
