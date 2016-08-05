/* abfns.f -- translated by f2c (version 19950602).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

#include "f2c.h"
#include <math.h>

#ifdef abs
#undef abs
#endif
#define abs(x)  fabs((x))

//#include "bind.h"

/* Table of constant values */

static integer c_n1 = -1;

/* Subroutine */ int abfns_(a, b, sk1, sk2, rr, l1, l2, m, n1, n2, maxcal)
doublereal *a, *b, *sk1, *sk2, *rr;
integer *l1, *l2, *m, *n1, *n2, *maxcal;
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double exp();
    integer pow_ii();

    /* Local variables */
    static doublereal rho22, c, d, h;
    static integer i, j, k;
    static doublereal r, t, ra;
    static integer il, in, ir, is, ix;
    static doublereal tr, rho1, rho2;


/*  ********************************************************************
*/
/*  *                                                                  *
*/
/*  *     SUBROUTINE ABFNS        CALLED FROM MOV                   * */
/*  *                                                                  *
*/
/*  *                                                                  *
*/
/*  *       ABFNS    SUBROUTINE TO CALCULATE THE AB FUNCTIONS          *
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

/*       modified by greg in modernity (august 1993) so that it doesn't us
e that*/
/*         damn common block any more. */

/*  *                                                                  *
*/
/*  ********************************************************************
*/


/*     THIS ONLY WORKS FOR PRINCIPAL QUANTUM # < OR = 7 */

/*      COMMON/LOCLAP/SK1,SK2,RR,L1,L2,M,N1,N2,MAXCAL */

    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    j = *maxcal + 1;
    rho1 = (*sk1 + *sk2) * .5 * *rr;
    rho2 = (*sk1 - *sk2) * .5 * *rr;
    if (abs(rho1) > 165.) {
        goto L100;
    }
    if (abs(rho2) > 165.) {
        goto L100;
    }
    c = exp(-rho1);
    a[1] = c / rho1;
    i__1 = j;
    for (i = 2; i <= i__1; ++i) {
/* L15: */
        a[i] = ((doublereal) ((real) (i - 1)) * a[i - 1] + c) / rho1;
    }
    ix = j;
    ir = (d__1 = rho2 * 2., (integer) abs(d__1));
/* Computing MIN */
    i__1 = ir + 1;
    is = min(i__1,19);
    if (rho2 != 0.) {
        goto L25;
    } else {
        goto L35;
    }
L25:
    d = exp(rho2);
    h = 1. / d;

/*  IF THE VALUE OF RHO2 IS TOO SMALL THE SINH MUST BE OBTAINED */
/*  BY SUMMING THE INFINITE SERIES RATHER THAN BY ADDITION OF */
/*  TWO EXPONENTIALS. */

    r = d - h;
    if (abs(r) - (float).1 >= 0.) {
        goto L28;
    } else {
        goto L26;
    }
L26:
    ra = rho2;
    rho22 = rho2 * rho2;
    t = rho2;
    for (i = 2; i <= 50; i += 2) {
        t = t * rho22 / (doublereal) ((real) (i * i + i));
        ra += t;
        if (t < 1e-30) {
            goto L999;
        }
/* L27: */
    }
L999:
    r = ra + ra;

/*  AS MANY SUCCESSIVE B-FUNCTIONS ARE GENERATED FROM B-0 BY THE */
/*  RECURRENCE FORMULAE AS ACCURACY WILL PERMIT. */

L28:
    b[1] = r / rho2;
    i__1 = ix;
    i__2 = is;
    for (i = 2; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2) {
        if (ir == 0) {
            goto L40;
        }
/* L32: */
        il = is - 1;
        if (1 > il) {
            goto L9050;
        }
        i__3 = il;
        for (j = 1; j <= i__3; ++j) {
            k = i + j - 1;
            if (pow_ii(&c_n1, &k) <= 0) {
                goto L29;
            } else {
                goto L30;
            }
L29:
            b[k] = (r + (doublereal) ((real) (k - 1)) * b[k - 1]) / rho2;
            goto L31;
L30:
            b[k] = -(d + h - (doublereal) ((real) (k - 1)) * b[k - 1]) / rho2;
L31:
            ;
        }
L9050:
L40:
        in = i + is - 1;

/*  AFTER THE RECURRENCE FORMULAE HAVE BEEN APPLIED AN APPROPRIATE */
/*  NUMBER OF TIMES THE NEXT B-FUNCTION IS OBTAINED BY SUMMATION */
/*  OF THE INFINITE SERIES. */

        if (in - ix <= 0) {
            goto L39;
        } else {
            goto L38;
        }
L39:
        if (pow_ii(&c_n1, &in) <= 0) {
            goto L44;
        } else {
            goto L42;
        }
L42:
        tr = rho2;
/* L105: */
        b[in] = tr * -2. / (doublereal) ((real) (in + 1));
        for (j = 1; j <= 500; ++j) {
/* Computing 2nd power */
            d__1 = rho2;
            tr = tr * (d__1 * d__1) / (doublereal) ((real) ((j << 1) * ((j <<
                    1) + 1)));
            if ((d__1 = tr / b[in], abs(d__1)) - 1e-7 <= 0.) {
                goto L51;
            } else {
                goto L43;
            }
L43:
            b[in] -= tr * 2. / (doublereal) ((real) (in + 1 + (j << 1)));
        }
        goto L51;
L44:
        tr = (float)1.;
/* L107: */
        b[in] = tr * 2. / (doublereal) ((real) in);
        for (j = 1; j <= 500; ++j) {
/* Computing 2nd power */
            d__1 = rho2;
            tr = tr * (d__1 * d__1) / (doublereal) ((real) ((j << 1) * ((j <<
                    1) - 1)));
            if ((d__1 = tr / b[in], abs(d__1)) - 1e-7 <= 0.) {
                goto L51;
            } else {
                goto L46;
            }
L46:
            b[in] += tr * 2. / (doublereal) ((real) (in + (j << 1)));
        }
L51:
        ;
    }

/*  IF THE ARGUMENT IS ZERO A SEPARATE FORMULA MUST BE USED. */

    goto L38;
L35:
    i__2 = ix;
    for (i = 1; i <= i__2; i += 2) {
        b[i] = 2. / (doublereal) ((real) i);
/* L36: */
        b[i + 1] = 0.;
    }
L38:
    return 0;
L100:
    for (i = 1; i <= 20; ++i) {
        a[i] = 0.;
/* L101: */
        b[i] = 0.;
    }
    goto L38;
} /* abfns_ */

