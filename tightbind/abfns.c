/*
Produced by gmFortran V30.59(10/26/17) on 9/17/18 at 9:45:37
*/
#define LPROTOTYPE
#include "fortran.h"
/*
*/
void abfns(double *a,double *b,double *sk1,double *sk2,double *rr,int *l1,
    int *l2,int *m,int *n1,int *n2,int *maxcal)
{
/*
  ********************************************************************
  *                                                                  *
  *     SUBROUTINE ABFNS        CALLED FROM MOV                   *
  *                                                                  *
  *                                                                  *
  *       ABFNS    SUBROUTINE TO CALCULATE THE AB FUNCTIONS          *
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
static int D1,D2;
static int j;
static double rho1,rho2,c;
static int i,ix,ir,is;
static double d,h,r,ra,rho22,t;
static int il,k,in;
static double tr;
/*
     THIS ONLY WORKS FOR PRINCIPAL QUANTUM # < OR = 7

      COMMON/LOCLAP/SK1,SK2,RR,L1,L2,M,N1,N2,MAXCAL

*/
    j = *maxcal+1;
    rho1 = 0.5e0*(*sk1+*sk2)**rr;
    rho2 = 0.5e0*(*sk1-*sk2)**rr;
    if(fabs(rho1) > 165.e0) goto S100;
    if(fabs(rho2) > 165.e0) goto S100;
    c = exp(-rho1);
    a[0] = c/rho1;
    for(i=2; i<=j; i++) {
        a[i-1] = ((double)(float)(i-1)*a[i-2]+c)/rho1;
    }
    ix = j;
    ir = fabs(2.e0*rho2);
    is = fifmin0(ir+1,19);
    if(rho2 == 0) goto S35;
    d = exp(rho2);
    h = 1.e0/d;
/*
  IF THE VALUE OF RHO2 IS TOO SMALL THE SINH MUST BE OBTAINED
  BY SUMMING THE INFINITE SERIES RATHER THAN BY ADDITION OF
  TWO EXPONENTIALS.

*/
    r = d-h;
    if(fabs(r) >= 0.1) goto S28;
    ra = rho2;
    rho22 = rho2*rho2;
    t = rho2;
    for(i=2; i<=50; i+=2) {
        t = t*rho22/(double)(float)(i*i+i);
        ra = ra+t;
        if(t < 1.e-30) goto S999;
    }
S999:
    r = ra+ra;
/*
  AS MANY SUCCESSIVE B-FUNCTIONS ARE GENERATED FROM B-0 BY THE
  RECURRENCE FORMULAE AS ACCURACY WILL PERMIT.

*/
S28:
    b[0] = r/rho2;
    for(i=2,D1=is,D2=(ix-i+D1)/D1; D2>0; D2--,i+=D1) {
        if(ir == 0) goto S40;
        il = is-1;
        if(1 > il) goto S9050;
        for(j=1; j<=il; j++) {
            k = i+j-1;
            if(pow(-1,k) > 0) goto S30; // fifipow replaced by pow
            b[k-1] = (r+(double)(float)(k-1)*b[k-2])/rho2;
            goto S31;
S30:
            b[k-1] = -((d+h-(double)(float)(k-1)*b[k-2])/rho2);
S31:;
        }
S9050:
S40:
        in = i+is-1;
/*
  AFTER THE RECURRENCE FORMULAE HAVE BEEN APPLIED AN APPROPRIATE
  NUMBER OF TIMES THE NEXT B-FUNCTION IS OBTAINED BY SUMMATION
  OF THE INFINITE SERIES.

*/
        if(in-ix > 0) goto S38;
        if(pow(-1,in) <= 0) goto S44; // fifipow replaced by pow cause in is an integer
        tr = rho2;
        b[in-1] = -(2.e0*tr/(double)(float)(in+1));
        for(j=1; j<=500; j++) {
            tr = tr*pow(rho2,2.0)/(double)(float)(2*j*(2*j+1));
            if(fabs(tr/b[in-1]) <= 1.0e-7) goto S51;
            b[in-1] = b[in-1]-2.e0*tr/(double)(float)(in+1+2*j);
        }
        goto S51;
S44:
        tr = 1.;
        b[in-1] = 2.e0*tr/(double)(float)in;
        for(j=1; j<=500; j++) {
            tr = tr*pow(rho2,2.0)/(double)(float)(2*j*(2*j-1));
            if(fabs(tr/b[in-1]) <= 1.0e-7) goto S51;
            b[in-1] = b[in-1]+2.e0*tr/(double)(float)(in+2*j);
        }
S51:;
    }
/*
  IF THE ARGUMENT IS ZERO A SEPARATE FORMULA MUST BE USED.

*/
    goto S38;
S35:
    for(i=1; i<=ix; i+=2) {
        b[i-1] = 2.e0/(double)(float)i;
        b[i] = 0.e0;
    }
S38:
    return;
S100:
    for(i=1; i<=20; i++) {
        a[i-1] = 0.e0;
        b[i-1] = 0.e0;
    }
    goto S38;
}
