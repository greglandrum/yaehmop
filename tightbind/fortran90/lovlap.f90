!*==LOVLAP.spg  processed by SPAG 6.72Dc at 06:56 on 17 Sep 2018
      SUBROUTINE LOVLAP(Strad,A,B,Sk1,Sk2,R,L1,L2,M1,N1,N2,Max)
!
!  ********************************************************************
!  *                                                                  *
!  *     SUBROUTINE LOVLAP       CALLED FROM MOV                   *
!  *                                                                  *
!  *                                                                  *
!  *       LOVLAP   SUBROUTINE TO CALCULATE THE OVERLAP               *
!  *                COMPONENT INDEPENDENT OF THE ANGLE BETWEEN        *
!  *                THE ATOMS                                         *
!  *                                                                  *
!  *       SUBROUTINES USED:                                          *
!  *                                                                  *
!  *             NONE                                                 *
!  *                                                                  *
!  *                                                                  *
!  *       ORIGIN LOST IN ANTIQUITY                                   *
!
!    modified by greg in modernity (august 1993) so that it doesn't use that
!         damn common block any more.
!
!  *                                                                  *
!  ********************************************************************
!
      IMPLICIT NONE
!*--LOVLAP27
!*** Start of declarations inserted by SPAG
      REAL*8 A , B , bincoe , con1 , con12 , fact , R , rhoa , rhoab ,  &
           & rhoap , rhob , rhopo , Sk1 , Sk2 , Strad , term , terma ,  &
           & value , value1 , value2
      REAL*8 value3 , value4
      INTEGER i , i1 , i2 , i3 , i4 , i5 , i6 , iab , ibb , icb , idb , &
            & ieb , iev , ip , ir , j , jend , ju , k , kend
      INTEGER ku , L1 , L2 , M1 , m2 , Max , maxxa , maxxb , maxxc ,    &
            & N1 , N2
!*** End of declarations inserted by SPAG
!
      DIMENSION A(30) , B(30)
      DIMENSION fact(25)
      DIMENSION bincoe(8,8)
!
!      COMMON/LOCLAP/SK1,SK2,R,L1,L2,M1,N1,N2,MAX
!
      DATA bincoe/8*1.D0 , 0.D0 , 1.D0 , 2.D0 , 3.D0 , 4.D0 , 5.D0 ,    &
         & 6.D0 , 7.D0 , 2*0.D0 , 1.D0 , 3.D0 , 6.D0 , 10.D0 , 15.D0 ,  &
         & 21.D0 , 3*0.D0 , 1.D0 , 4.D0 , 10.D0 , 20.D0 , 35.D0 ,       &
         & 4*0.D0 , 1.D0 , 5.D0 , 15.D0 , 35.D0 , 5*0.D0 , 1.D0 , 6.D0 ,&
         & 21.D0 , 6*0.D0 , 1.D0 , 7.D0 , 7*0.D0 , 1.D0/
!
 
!      write (*,*) 'Lovlap: ',sk1,sk2,r,l1,l2,m1,n1,n2,max
 
      maxxa = 1
      maxxb = 1
      maxxc = 1
      fact(1) = 1.D0
      DO i = 2 , 25
         fact(i) = fact(i-1)*DBLE(FLOAT(i-1))
      ENDDO
      m2 = M1
      Strad = 0.D0
      rhoa = R*Sk1
      rhob = R*Sk2
      rhoap = rhoa**N1
      rhoap = rhoap*rhoap
      rhoap = rhoap*rhoa
      rhoab = rhob**N2
      rhoab = rhoab*rhoab
      rhoab = rhoab*rhob
      rhopo = rhoap*rhoab
      terma = 0.5D0**(L1+L2+1)                                          &
            & *DSQRT(DBLE(FLOAT((L1+L1+1)*(L2+L2+1)))*fact(L1-M1+1)     &
            & *fact(L2-M1+1)/(fact(N1+N1+1)*fact(N2+N2+1)*fact(L1+M1+1) &
            & *fact(L2+M1+1))*rhopo)
      jend = 1 + ((L1-M1)/2)
      kend = 1 + ((L2-m2)/2)
      ieb = M1 + 1
      DO j = 1 , jend
         ju = j - 1
         iab = N1 - L1 + ju + ju + 1
         icb = L1 - M1 - ju - ju + 1
         con1 = fact(L1+L1-ju-ju+1)                                     &
              & /(fact(L1-M1-ju-ju+1)*fact(ju+1)*fact(L1-ju+1))
         DO k = 1 , kend
            ku = k - 1
            con12 = con1*fact(L2+L2-ku-ku+1)                            &
                  & /(fact(L2-m2-ku-ku+1)*fact(ku+1)*fact(L2-ku+1))
            iev = ju + ku + L2
            IF ( 2*(iev/2)/=iev ) con12 = -con12
            ibb = N2 - L2 + ku + ku + 1
            idb = L2 - m2 - ku - ku + 1
            value = 0.D0
            DO i6 = 1 , ieb
               DO i5 = 1 , ieb
                  value1 = bincoe(ieb,i6)*bincoe(ieb,i5)
                  iev = i5 + i6
                  IF ( 2*(iev/2)/=iev ) value1 = -value1
                  DO i4 = 1 , idb
                     value1 = -value1
                     value2 = bincoe(idb,i4)*value1
                     DO i3 = 1 , icb
                        value3 = bincoe(icb,i3)*value2
                        DO i2 = 1 , ibb
                           value3 = -value3
                           value4 = bincoe(ibb,i2)*value3
                           DO i1 = 1 , iab
                              term = value4*bincoe(iab,i1)
                              ir = i1 + i2 + ieb + ieb - i6 - i6 - i3 + &
                                 & idb - i4 + icb - 1
                              ip = iab - i1 + ibb - i2 + ieb + ieb -    &
                                 & i5 - i5 + icb - i3 + idb - i4 + 1
                              value = value + A(ip)*B(ir)*term
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            Strad = Strad + value*con12
         ENDDO
      ENDDO
      Strad = Strad*terma
 
      END
