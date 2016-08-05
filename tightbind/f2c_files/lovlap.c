/* lovlap.f -- translated by f2c (version 19950602).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

#include "f2c.h"
#include <math.h>


#ifdef abs
#undef abs
#endif
#define abs(x)  fabs((x))


/* Table of constant values */

static doublereal c_b3 = .5;

/* Subroutine */ int lovlap_(strad, a, b, sk1, sk2, r, l1, l2, m1, n1, n2,
        max__)
doublereal *strad, *a, *b, *sk1, *sk2, *r;
integer *l1, *l2, *m1, *n1, *n2, *max__;
{
    /* Initialized data */

    static doublereal bincoe[64]        /* was [8][8] */ = { 1.,1.,1.,1.,1.,
            1.,1.,1.,0.,1.,2.,3.,4.,5.,6.,7.,0.,0.,1.,3.,6.,10.,15.,21.,0.,0.,
            0.,1.,4.,10.,20.,35.,0.,0.,0.,0.,1.,5.,15.,35.,0.,0.,0.,0.,0.,1.,
            6.,21.,0.,0.,0.,0.,0.,0.,1.,7.,0.,0.,0.,0.,0.,0.,0.,1. };

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;

    /* Builtin functions */
    double pow_di(), sqrt();

    /* Local variables */
    static doublereal fact[25];
    static integer jend, kend;
    static doublereal con12, rhoa, rhob, term;
    static integer i, j, k;
    static doublereal rhoab, terma, rhoap, value;
    static integer maxxa, maxxb, maxxc;
    static doublereal rhopo;
    static integer i1, i2, i3, m2, i6, i5, i4;
    static doublereal value1, value2, value3, value4;
    static integer ip, ir, ju, ku, iab, ibb, icb, idb, ieb, iev;
    static doublereal con1;


/*  ********************************************************************
*/
/*  *                                                                  *
*/
/*  *     SUBROUTINE LOVLAP       CALLED FROM MOV                   * */
/*  *                                                                  *
*/
/*  *                                                                  *
*/
/*  *       LOVLAP   SUBROUTINE TO CALCULATE THE OVERLAP               *
*/
/*  *                COMPONENT INDEPENDENT OF THE ANGLE BETWEEN        *
*/
/*  *                THE ATOMS                                         *
*/
/*  *                                                                  *
*/
/*  *       SUBROUTINES USED:                                          *
*/
/*  *                                                                  *
*/
/*  *             NONE                                                 *
*/
/*  *                                                                  *
*/
/*  *                                                                  *
*/
/*  *       ORIGIN LOST IN ANTIQUITY                                   *
*/

/*   modified by greg in modernity (august 1993) so that it doesn't use th
at*/
/*         damn common block any more. */

/*  *                                                                  *
*/
/*  ********************************************************************
*/



/*      COMMON/LOCLAP/SK1,SK2,R,L1,L2,M1,N1,N2,MAX */

    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */

/*      write (*,*) 'Lovlap: ',sk1,sk2,r,l1,l2,m1,n1,n2,max */
    maxxa = 1;
    maxxb = 1;
    maxxc = 1;
    fact[0] = 1.;
    for (i = 2; i <= 25; ++i) {
/* L5: */
        fact[i - 1] = fact[i - 2] * (doublereal) ((real) (i - 1));
    }
    m2 = *m1;
    *strad = 0.;
    rhoa = *r * *sk1;
    rhob = *r * *sk2;
    rhoap = pow_di(&rhoa, n1);
    rhoap *= rhoap;
    rhoap *= rhoa;
    rhoab = pow_di(&rhob, n2);
    rhoab *= rhoab;
    rhoab *= rhob;
    rhopo = rhoap * rhoab;
    i__1 = *l1 + *l2 + 1;
    terma = pow_di(&c_b3, &i__1) * sqrt((doublereal) ((real) ((*l1 + *l1 + 1)
            * (*l2 + *l2 + 1))) * fact[*l1 - *m1] * fact[*l2 - *m1] / (fact[*
            n1 + *n1] * fact[*n2 + *n2] * fact[*l1 + *m1] * fact[*l2 + *m1]) *
             rhopo);
    jend = (*l1 - *m1) / 2 + 1;
    kend = (*l2 - m2) / 2 + 1;
    ieb = *m1 + 1;
    i__1 = jend;
    for (j = 1; j <= i__1; ++j) {
        ju = j - 1;
        iab = *n1 - *l1 + ju + ju + 1;
        icb = *l1 - *m1 - ju - ju + 1;
        con1 = fact[*l1 + *l1 - ju - ju] / (fact[*l1 - *m1 - ju - ju] * fact[
                ju] * fact[*l1 - ju]);
        i__2 = kend;
        for (k = 1; k <= i__2; ++k) {
            ku = k - 1;
            con12 = con1 * fact[*l2 + *l2 - ku - ku] / (fact[*l2 - m2 - ku -
                    ku] * fact[ku] * fact[*l2 - ku]);
            iev = ju + ku + *l2;
            if (iev / 2 << 1 != iev) {
                con12 = -con12;
            }
            ibb = *n2 - *l2 + ku + ku + 1;
            idb = *l2 - m2 - ku - ku + 1;
            value = 0.;
            i__3 = ieb;
            for (i6 = 1; i6 <= i__3; ++i6) {
                i__4 = ieb;
                for (i5 = 1; i5 <= i__4; ++i5) {
                    value1 = bincoe[ieb + (i6 << 3) - 9] * bincoe[ieb + (i5 <<
                             3) - 9];
                    iev = i5 + i6;
                    if (iev / 2 << 1 != iev) {
                        value1 = -value1;
                    }
                    i__5 = idb;
                    for (i4 = 1; i4 <= i__5; ++i4) {
                        value1 = -value1;
                        value2 = bincoe[idb + (i4 << 3) - 9] * value1;
                        i__6 = icb;
                        for (i3 = 1; i3 <= i__6; ++i3) {
                            value3 = bincoe[icb + (i3 << 3) - 9] * value2;
                            i__7 = ibb;
                            for (i2 = 1; i2 <= i__7; ++i2) {
                                value3 = -value3;
                                value4 = bincoe[ibb + (i2 << 3) - 9] * value3;
                                i__8 = iab;
                                for (i1 = 1; i1 <= i__8; ++i1) {
                                    term = value4 * bincoe[iab + (i1 << 3) -
                                            9];
                                    ir = i1 + i2 + ieb + ieb - i6 - i6 - i3 +
                                            idb - i4 + icb - 1;
                                    ip = iab - i1 + ibb - i2 + ieb + ieb - i5
                                            - i5 + icb - i3 + idb - i4 + 1;
/* L90: */
                                    value += a[ip] * b[ir] * term;
                                }
                            }
                        }
                    }
                }
            }
/* L50: */
            *strad += value * con12;
        }
    }
    *strad *= terma;
    return 0;
} /* lovlap_ */

