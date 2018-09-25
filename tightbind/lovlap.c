/*
Produced by gmFortran V30.59(10/26/17) on 9/17/18 at 9:45:37
*/
#define LPROTOTYPE
#include "fortran.h"
/*
*/
void lovlap(double *strad,double *a,double *b,double *sk1,double *sk2,double *r,
    int *l1,int *l2,int *m1,int *n1,int *n2,int *max)
{
/*
  ********************************************************************
  *                                                                  *
  *     SUBROUTINE LOVLAP       CALLED FROM MOV                   *
  *                                                                  *
  *                                                                  *
  *       LOVLAP   SUBROUTINE TO CALCULATE THE OVERLAP               *
  *                COMPONENT INDEPENDENT OF THE ANGLE BETWEEN        *
  *                THE ATOMS                                         *
  *                                                                  *
  *       SUBROUTINES USED:                                          *
  *                                                                  *
  *             NONE                                                 *
  *                                                                  *
  *                                                                  *
  *       ORIGIN LOST IN ANTIQUITY                                   *

    modified by greg in modernity (august 1993) so that it doesn't use that
         damn common block any more.

  *                                                                  *
  ********************************************************************

*/
/*
      COMMON/LOCLAP/SK1,SK2,R,L1,L2,M1,N1,N2,MAX

*/
static double bincoe[8][8] = {
    1.e0,1.e0,1.e0,1.e0,1.e0,1.e0,1.e0,1.e0,0.e0,1.e0,2.e0,3.e0,4.e0,5.e0,6.e0,
    7.e0,0.e0,0.e0,1.e0,3.e0,6.e0,10.e0,15.e0,21.e0,0.e0,0.e0,0.e0,1.e0,4.e0,
    10.e0,20.e0,35.e0,0.e0,0.e0,0.e0,0.e0,1.e0,5.e0,15.e0,35.e0,0.e0,0.e0,0.e0,
    0.e0,0.e0,1.e0,6.e0,21.e0,0.e0,0.e0,0.e0,0.e0,0.e0,0.e0,1.e0,7.e0,0.e0,0.e0,
    0.e0,0.e0,0.e0,0.e0,0.e0,1.e0
};
static double fact[25];
static int maxxa,maxxb,maxxc,i,m2;
static double rhoa,rhob,rhoap,rhoab,rhopo,terma;
static int jend,kend,ieb,j,ju,iab,icb;
static double con1;
static int k,ku;
static double con12;
static int iev,ibb,idb;
static double value;
static int i6,i5;
static double value1;
static int i4;
static double value2;
static int i3;
static double value3;
static int i2;
static double value4;
static int i1;
static double term;
static int ir,ip;
/*
      write (*,*) 'Lovlap: ',sk1,sk2,r,l1,l2,m1,n1,n2,max
*/
    maxxa = 1;
    maxxb = 1;
    maxxc = 1;
    fact[0] = 1.e0;
    for(i=1; i<25; i++) {
        fact[i] = fact[i-1]*(double)(float)i;
    }
    m2 = *m1;
    *strad = 0.e0;
    rhoa = *r**sk1;
    rhob = *r**sk2;
    rhoap = pow(rhoa,(double)*n1);
    rhoap = rhoap*rhoap;
    rhoap = rhoap*rhoa;
    rhoab = pow(rhob,(double)*n2);
    rhoab = rhoab*rhoab;
    rhoab = rhoab*rhob;
    rhopo = rhoap*rhoab;
    terma = pow(0.5e0,(double)(*l1+*l2+1))*sqrt((double)(float)((*l1+*l1+1)*(*
      l2+*l2+1))*fact[*l1-*m1]*fact[*l2-*m1]/(fact[*n1+*n1]*fact[*n2+*n2]*fact[*
      l1+*m1]*fact[*l2+*m1])*rhopo);
    jend = 1+(*l1-*m1)/2;
    kend = 1+(*l2-m2)/2;
    ieb = *m1+1;
    for(j=0; j<jend; j++) {
        ju = j;
        iab = *n1-*l1+ju+ju+1;
        icb = *l1-*m1-ju-ju+1;
        con1 = fact[*l1+*l1-ju-ju]/(fact[*l1-*m1-ju-ju]*fact[ju]*fact[*l1-ju]);
        for(k=0; k<kend; k++) {
            ku = k;
            con12 = con1*fact[*l2+*l2-ku-ku]/(fact[*l2-m2-ku-ku]*fact[ku]*fact[*
              l2-ku]);
            iev = ju+ku+*l2;
            if(2*(iev/2) != iev) con12 = -con12;
            ibb = *n2-*l2+ku+ku+1;
            idb = *l2-m2-ku-ku+1;
            value = 0.e0;
            for(i6=1; i6<=ieb; i6++) {
                for(i5=1; i5<=ieb; i5++) {
                    value1 = bincoe[i6-1][ieb-1]*bincoe[i5-1][ieb-1];
                    iev = i5+i6;
                    if(2*(iev/2) != iev) value1 = -value1;
                    for(i4=1; i4<=idb; i4++) {
                        value1 = -value1;
                        value2 = bincoe[i4-1][idb-1]*value1;
                        for(i3=1; i3<=icb; i3++) {
                            value3 = bincoe[i3-1][icb-1]*value2;
                            for(i2=1; i2<=ibb; i2++) {
                                value3 = -value3;
                                value4 = bincoe[i2-1][ibb-1]*value3;
                                for(i1=1; i1<=iab; i1++) {
                                    term = value4*bincoe[i1-1][iab-1];
                                    ir = i1+i2+ieb+ieb-i6-i6-i3+idb-i4+icb-1;
                                    ip = iab-i1+ibb-i2+ieb+ieb-i5-i5+icb-i3+idb-
                                      i4+1;
                                    value = value+a[ip-1]*b[ir-1]*term;
                                }
                            }
                        }
                    }
                }
            }
            *strad = *strad+value*con12;
        }
    }
    *strad = *strad*terma;
    return;
}
