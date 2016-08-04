/* cboris.f -- translated by f2c (version 19950602).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

#include "f2c.h"
#include <math.h>
//#include "bind.h"

extern int cchol_(integer *n, integer *nd, doublereal *a, integer *fail);
extern int ctred2_(integer *n, integer *nd, doublereal *a, doublereal *b,
                                        doublereal *d,doublereal *e,doublereal *f);

extern int ctql2_(integer *n, integer *nd, doublereal *d, doublereal *e,
                                        doublereal *f,doublereal *a,doublereal *b,integer *fail);

int cboris_(integer *,integer *,doublereal *,doublereal *,doublereal *,doublereal *,
                        doublereal *,doublereal *, integer *);

/* Subroutine */ int cboris_(n, nd, a, b, c, d, e, f, fail)
integer *n, *nd;
doublereal *a, *b, *c, *d, *e, *f;
integer *fail;
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
            i__3;

    /* Local variables */
    integer i, j, k;
    integer ia, ja, lf, ii;



/*  for YAeHMOP this common block is not needed */
/*              gL */
/*      COMMON/PASS/NATOM,NDIM,LB,LE,NVAL,ECUT,EERR,IDEL */


/* SUBROUTINE CBORIS SOLVES THE COMPLEX EIGENVALUE */
/* PROBLEM */

/*      A*X=B*X*LAMBDA */

/* WHERE */

/* A     --- GIVEN MATRIX OF ORDER N */
/* B     --- GIVEN POSITIVE DEFINITE MATRIX OF ORDER N */
/* X     --- EIGENVECTORS, COLUMN BY COLUMN, */
/*           NORMALIZED TO (BX,X)=I */
/* LAMBDA --- EIGENVALUES IN INCREASING ORDER */

/* MATRICES A AND B SHOULD BE GIVEN IN */
/* ARRAYS A AND B RESPECTIVELY IN THE FOLLOWING WAY */

/* DIAGONAL ELEMENTS IN THE DIAGONAL */
/* REAL PARTS OF THE LOWER TRIANGLE IN THE LOWER */
/* TRIANGLE */
/* IMAG PARTS OF THE LOWER TRIANGLE IN THE UPPER */
/* TRIANGLE */

/* REAL PARTS OF EIGENVECTORS ARE STORED IN A, */
/* IMAG PARTS IN ARRAY C.  ARRAY B IS DESTROYED DURING */
/* COMPUTATION.  IT HOLDS THE LOWER TRIANGLE OF THE */
/* CHOLESKI DECOMPOSITION OF MATRIX B (SEE */
/* SUBROUTINE CCHOL). */

/* ALL MATRICES ARE OF ORDER N, DECLARED IN THE */
/* CALLING PROGRAM WITH DIMENSION ND WHICH NEED */
/* NOT TO BE EQUAL TO N. */

/* EIGENVALUES ARE STORED IN ARRAY D.  ARRAYS E */
/* AND F ARE USED FOR INTERMEDIATE RESULTS. */

/* FAIL GETS THE FOLLOWING VALUES */

/*  0 --- COMPUTATION FINISHED SUCCESFULLY */
/*  1 --- B NOT POSITIVE DEFINITE */
/*  2 --- QR ALGORITHM DOES NOT CONVERGE */

/* SUBROUTINES USED */

/*     CCHOL */
/*     CTRED2 */
/*     CTQL2 */


/* PROGRAMMED BY E. ZAKRAJSEK */
/* JUNE 21,1974 */




/* DECOMPOSE MATRIX B */

/*      WRITE(4,1918) */
/* 1918  FORMAT(' WELCOME TO CBORIS') */

/*      CALL TIME(0) */

    /* Parameter adjustments */
    --f;
    --e;
    --d;
    c_dim1 = *nd;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    b_dim1 = *nd;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *nd;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *fail = 1;
    cchol_(n, nd, &b[b_offset], &lf);
    if (lf != 0) {
        goto L999;
    }

/* MOVE MATRIX A */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
        c[i + i * c_dim1] = 0.;
        if (i == 1) {
            goto L11;
        }
        ia = i - 1;
        i__2 = ia;
        for (j = 1; j <= i__2; ++j) {
            c[i + j * c_dim1] = a[j + i * a_dim1];
            c[j + i * c_dim1] = -a[j + i * a_dim1];
/* L10: */
            a[j + i * a_dim1] = a[i + j * a_dim1];
        }
L11:
        ;
    }

/*  COMPUTE (L(-1)*A) */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
        i__2 = *n;
        for (i = 1; i <= i__2; ++i) {
            if (i == 1) {
                goto L21;
            }
            ia = i - 1;
            i__3 = ia;
            for (k = 1; k <= i__3; ++k) {
                a[i + j * a_dim1] = a[i + j * a_dim1] - a[k + j * a_dim1] * b[
                        i + k * b_dim1] + c[k + j * c_dim1] * b[k + i *
                        b_dim1];
/* L20: */
                c[i + j * c_dim1] = c[i + j * c_dim1] - a[k + j * a_dim1] * b[
                        k + i * b_dim1] - c[k + j * c_dim1] * b[i + k *
                        b_dim1];
            }
L21:
            a[i + j * a_dim1] /= b[i + i * b_dim1];
/* L22: */
            c[i + j * c_dim1] /= b[i + i * b_dim1];
        }
    }

/*  COMPUTE  A*L(-H) */

    i__2 = *n;
    for (i = 1; i <= i__2; ++i) {
        i__1 = i;
        for (j = 1; j <= i__1; ++j) {
            if (j == 1) {
                goto L31;
            }
            ja = j - 1;
            i__3 = ja;
            for (k = 1; k <= i__3; ++k) {
                a[i + j * a_dim1] = a[i + j * a_dim1] - a[i + k * a_dim1] * b[
                        j + k * b_dim1] - c[i + k * c_dim1] * b[k + j *
                        b_dim1];
/* L30: */
                c[i + j * c_dim1] = c[i + j * c_dim1] + a[i + k * a_dim1] * b[
                        k + j * b_dim1] - c[i + k * c_dim1] * b[j + k *
                        b_dim1];
            }
L31:
            a[i + j * a_dim1] /= b[j + j * b_dim1];
/* L32: */
            c[i + j * c_dim1] /= b[j + j * b_dim1];
        }
    }

/*     PUT MATRIX TOGETHER INTO A */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
        if (i == *n) {
            goto L41;
        }
        ia = i + 1;
        i__2 = *n;
        for (j = ia; j <= i__2; ++j) {
/* L40: */
            a[i + j * a_dim1] = c[j + i * c_dim1];
        }
L41:
        ;
    }

/*     DIAGONALIZE A */

    *fail = 2;
    ctred2_(n, nd, &a[a_offset], &c[c_offset], &d[1], &e[1], &f[1]);
    ctql2_(n, nd, &d[1], &e[1], &f[1], &a[a_offset], &c[c_offset], &lf);
    if (lf != 0) {
        goto L999;
    }

/*     COMPUTE L(-H)*A */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
        i__2 = *n;
        for (ii = 1; ii <= i__2; ++ii) {
            i = *n - ii + 1;
            if (i == *n) {
                goto L51;
            }
            ia = i + 1;
            i__3 = *n;
            for (k = ia; k <= i__3; ++k) {
                a[i + j * a_dim1] = a[i + j * a_dim1] - a[k + j * a_dim1] * b[
                        k + i * b_dim1] - c[k + j * c_dim1] * b[i + k *
                        b_dim1];
/* L50: */
                c[i + j * c_dim1] = c[i + j * c_dim1] + a[k + j * a_dim1] * b[
                        i + k * b_dim1] - c[k + j * c_dim1] * b[k + i *
                        b_dim1];
            }
L51:
            a[i + j * a_dim1] /= b[i + i * b_dim1];
            c[i + j * c_dim1] /= b[i + i * b_dim1];
/* L52: */
        }
    }
    *fail = 0;
/*  999 WRITE(4,1919) */
/* 1919 FORMAT(' DONE CBORIS') */

/*      CALL TIME(0) */

L999:
    return 0;
} /* cboris_ */

