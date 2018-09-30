/*
Produced by gmFortran V30.59(10/26/17) on 9/18/18 at 9:02:32
*/
#define LPROTOTYPE
#include "fortran.h"
void cchol(int *n,int *nd,double *a,int *fail)
{
static int i,ia,j,k,ka;
/*

 SUBROUTINE CCHOL COMPUTES CHOLESKI
 DECOMPOSITION OF GIVEN COMPLEX POSITIVE DEFINITE
 MATRIX A.

 INPUT DATA

  N --- ORDER OF MATRIX
  ND -- DIMENSION OF ARRAY A (IT CAN BE
        GREATER THAN OR EQUAL TO N)
  A --- GIVEN MATRIX
        IT IS SUPPOSED TO BE STORED IN THE
        FOLLOWING WAY.  DIAGONAL ELEMENTS,
        BEING REAL, ARE STORED ON THE DIAGONAL,
        REAL PARTS OF OFFDIAGONAL ELEMENTS
        ARE STORED IN THE LOWER TRIANGLE OF A,
        IMAG PARTS OF THE LOWER TRIANGLE ARE
        STORED IN THE UPPER TRIANGLE OF A.

     EXIT INFORMATION

  A --- COMPLEX ELEMENTS OF MATRIX L, DEFINED BY
     A=L*L(H)
        ARE STORED IN THE SAME WAY AS ORIGINAL
        ELEMENTS OF A, THAT IS REAL PARTS OF THE
        LOWER TRIANGLE OF L IN THE LOWER TRIANGLE
        OF A AND THE CORRESPONDING IMAG PARTS IN
        THE UPPER TRIANGLE OF A.
  FAIL --- IS SET TO ZERO IF THE DECOMPOSITION WAS
           SUCCESFUL AND TO NONZERO IF
           THE MATRIX WAS NOT POSITIVE DEFINITE.


    PROGRAMMED BY E. ZAKRAJSEK
     JUNE 20, 1974



     SUPPOSE DECOMPOSITION WILL FAIL

*/
    *fail = 1;
/*
*/
    for(i=1; i<=*n; i++) {
/*
     TEST FOR POSITIVE DEFINITNESS

*/
        if(*(a+i-1+(i-1)**nd) <= 0.0e0) return;
/*
      COMPUTE COLUMN I

*/
        *(a+i-1+(i-1)**nd) = sqrt(*(a+i-1+(i-1)**nd));
        if(i == *n) goto S13;
        ia = i+1;
        for(j=ia-1; j<*n; j++) {
            *(a+j+(i-1)**nd) = *(a+j+(i-1)**nd)/ *(a+i-1+(i-1)**nd);
            *(a+i-1+j**nd) = *(a+i-1+j**nd)/ *(a+i-1+(i-1)**nd);
        }
/*
     REDUCE REMAINING COLUMNS

*/
        for(k=ia; k<=*n; k++) {
            *(a+k-1+(k-1)**nd) = *(a+k-1+(k-1)**nd)-*(a+k-1+(i-1)**nd)**(a+k-1+
              (i-1)**nd)-*(a+i-1+(k-1)**nd)**(a+i-1+(k-1)**nd);
            if(k == *n) goto S12;
            ka = k+1;
            for(j=ka-1; j<*n; j++) {
                *(a+j+(k-1)**nd) = *(a+j+(k-1)**nd)-*(a+j+(i-1)**nd)**(a+k-1+(i
                  -1)**nd)-*(a+i-1+j**nd)**(a+i-1+(k-1)**nd);
                *(a+k-1+j**nd) = *(a+k-1+j**nd)-*(a+i-1+j**nd)**(a+k-1+(i-1)**
                  nd)+*(a+j+(i-1)**nd)**(a+i-1+(k-1)**nd);
            }
S12:;
        }
S13:;
    }
    *fail = 0;
    return;
}
void ctred2(int *n,int *nd,double *a,double *b,double *d,double *e,double *f)
{
static double chep;
static int k,l;
static double all;
static int i;
static double c,s,r,alr,ali,sm,g,t;
static int ia,j,kk;
/*

     SUBROUTINE CTRED2 REDUCES GIVEN COMPLEX
     HERMITIAN MATRIX TO A TRIDIAGONAL FORM

     PARAMETERS

     N    --- ORDER OF THE MATRIX
     ND   --- DIMENSION OF ARRAYS A AND B
     A    --- GIVEN MATRIX, REPLACED BY REAL PART
              OF THE TRANSFORMATION MATRIX
     B    --- IMAG PART OF TRANSFORMATION MATRIX
     D    --- DIAGONAL PART OF THE TRIADIAGONAL MATRIX
     E    --- REAL PART OF THE CODIAGONAL OF THE
              TRIDIAGONAL MATRIX
              (LAST N-1 LOCATIONS)
     F    --- IMAG PARTS OF THE LOWER CODIAGONAL.

     THE GIVEN MATRIX SHOULD BE STORED IN THE
     FOLLOWING WAY

          --- DIAGONAL ELEMENTS IN THE DIAGONAL
          --- REAL PART OF THE LOWER TRIANGLE IN THE
              LOWER TRIANGLE
          --- IMAG PARTS OF THE LOWER TRIANGLE
              IN THE UPPER TRIANGLE


     PROGRAMMED BY E. ZAKRAJSEK
     JUNE 20,1974



*/
    chep = pow(2.0e0,(double)-56);
    d[0] = *(a+0+0**nd);
    if(*n == 1) goto S31;
/*
 MAIN K LOOP

*/
    for(k=2; k<=*n; k++) {
        l = k-1;
/*
     COMPUTE NORM

*/
        all = 0.e0;
        for(i=k-1; i<*n; i++) {
            all = all+*(a+i+(l-1)**nd)**(a+i+(l-1)**nd)+*(a+l-1+i**nd)**(a+l-1+
              i**nd);
        }
        all = sqrt(all);
/*
     COMPUTE CONSTANTS

*/
        c = 1.0e0;
        s = 0.e0;
        r = sqrt(*(a+k-1+(l-1)**nd)**(a+k-1+(l-1)**nd)+*(a+l-1+(k-1)**nd)**(a+l
          -1+(k-1)**nd));
        if(fabs(r) < 1.e-50) r = 0.e0;
        if(r == 0.0e0) goto S11;
        c = *(a+k-1+(l-1)**nd)/r;
        s = *(a+l-1+(k-1)**nd)/r;
S11:
        alr = all*c;
        ali = all*s;
        *(a+l-1+(l-1)**nd) = 0.0e0;
/*
     TEST FOR SUPERFLUOUS TRANSFORMATION

*/
        sm = all*(all+r);
        if(fabs(sm) < 1.e-50) sm = 0.e0;
        if(sm == 0.e0) goto S20;
        g = 1.0e0/sm;
        *(a+l-1+(l-1)**nd) = g;
/*
*/
        *(a+k-1+(l-1)**nd) = *(a+k-1+(l-1)**nd)+alr;
        *(a+l-1+(k-1)**nd) = *(a+l-1+(k-1)**nd)+ali;
/*
     NOW COMPUTE U=A*W
     AND STORE INTO (E,F)

*/
        t = 0.0e0;
        for(i=k; i<=*n; i++) {
            c = *(a+i-1+(i-1)**nd)**(a+i-1+(l-1)**nd);
            s = *(a+i-1+(i-1)**nd)**(a+l-1+(i-1)**nd);
            if(i == k) goto S13;
            ia = i-1;
            for(j=k-1; j<ia; j++) {
                c = c+*(a+i-1+j**nd)**(a+j+(l-1)**nd)-*(a+j+(i-1)**nd)**(a+l-1+
                  j**nd);
                s = s+*(a+i-1+j**nd)**(a+l-1+j**nd)+*(a+j+(i-1)**nd)**(a+j+(l-1)
                  **nd);
            }
S13:
            if(i == *n) goto S15;
            ia = i+1;
            for(j=ia-1; j<*n; j++) {
                c = c+*(a+j+(i-1)**nd)**(a+j+(l-1)**nd)+*(a+i-1+j**nd)**(a+l-1+
                  j**nd);
                s = s+*(a+j+(i-1)**nd)**(a+l-1+j**nd)-*(a+i-1+j**nd)**(a+j+(l-1)
                  **nd);
            }
S15:
            e[i-1] = g*c;
            f[i-1] = g*s;
            t = t+*(a+i-1+(l-1)**nd)*c+*(a+l-1+(i-1)**nd)*s;
        }
        t = t*(g*g);
/*
    TRANSFORM  MATRIX

*/
        for(i=k; i<=*n; i++) {
            *(a+i-1+(i-1)**nd) = *(a+i-1+(i-1)**nd)-2.0e0*(*(a+i-1+(l-1)**nd)*e
              [i-1]+*(a+l-1+(i-1)**nd)*f[i-1])+t*(*(a+i-1+(l-1)**nd)**(a+i-1+(l
              -1)**nd)+*(a+l-1+(i-1)**nd)**(a+l-1+(i-1)**nd));
            if(i == k) goto S18;
            ia = i-1;
            for(j=k-1; j<ia; j++) {
                *(a+i-1+j**nd) = *(a+i-1+j**nd)-*(a+i-1+(l-1)**nd)*e[j]-*(a+l-1+
                  (i-1)**nd)*f[j]-*(a+j+(l-1)**nd)*e[i-1]-*(a+l-1+j**nd)*f[i-1]+
                  t*(*(a+i-1+(l-1)**nd)**(a+j+(l-1)**nd)+*(a+l-1+(i-1)**nd)**(a+
                  l-1+j**nd));
                *(a+j+(i-1)**nd) = *(a+j+(i-1)**nd)-*(a+l-1+(i-1)**nd)*e[j]+*(a+
                  i-1+(l-1)**nd)*f[j]+*(a+l-1+j**nd)*e[i-1]-*(a+j+(l-1)**nd)*f
                  [i-1]+t*(*(a+l-1+(i-1)**nd)**(a+j+(l-1)**nd)-*(a+i-1+(l-1)**
                  nd)**(a+l-1+j**nd));
            }
S18:;
        }
/*
     STORE DIAGONAL AND CODIAGONAL ELEMENTS

*/
S20:
        d[k-1] = *(a+k-1+(k-1)**nd);
        e[k-1] = -alr;
        f[k-1] = -ali;
    }
/*
     NOW ACCUMULATE TRANSFORMATIONS

*/
S31:
    *(a+*n-1+(*n-1)**nd) = 1.e0;
    *(b+*n-1+(*n-1)**nd) = 0.e0;
    if(*n == 1) return;
    for(kk=2; kk<=*n; kk++) {
        k = *n-kk+2;
        l = k-1;
/*
     SKIP TRANSFORMATION IF UNIT

*/
        if(fabs(*(a+l-1+(l-1)**nd)) < 1.e-50) *(a+l-1+(l-1)**nd) = 0.e0;
        if(*(a+l-1+(l-1)**nd) == 0.0e0) goto S36;
/*
     COMPUTE PRODUCT

*/
        for(j=k-1; j<*n; j++) {
            c = 0.e0;
            s = 0.e0;
            for(i=k; i<=*n; i++) {
                c = c+*(a+i-1+(l-1)**nd)**(a+i-1+j**nd)+*(a+l-1+(i-1)**nd)**(b+
                  i-1+j**nd);
                s = s+*(a+i-1+(l-1)**nd)**(b+i-1+j**nd)-*(a+l-1+(i-1)**nd)**(a+
                  i-1+j**nd);
            }
            c = c**(a+l-1+(l-1)**nd);
            s = s**(a+l-1+(l-1)**nd);
            for(i=k; i<=*n; i++) {
                *(a+i-1+j**nd) = *(a+i-1+j**nd)-c**(a+i-1+(l-1)**nd)+s**(a+l-1+
                  (i-1)**nd);
                *(b+i-1+j**nd) = *(b+i-1+j**nd)-c**(a+l-1+(i-1)**nd)-s**(a+i-1+
                  (l-1)**nd);
            }
        }
/*
     MAKE NEW LINE

*/
S36:
        for(i=k; i<=*n; i++) {
            *(a+i-1+(l-1)**nd) = 0.e0;
            *(a+l-1+(i-1)**nd) = 0.e0;
            *(b+i-1+(l-1)**nd) = 0.e0;
            *(b+l-1+(i-1)**nd) = 0.e0;
        }
        *(a+l-1+(l-1)**nd) = 1.e0;
        *(b+l-1+(l-1)**nd) = 0.e0;
    }
    return;
}
void ctql2(int *n,int *nd,double *d,double *e,double *f,double *a,double *b,
    int *fail)
{
static double chep;
static int k;
static double r,c,s;
static int i;
static double p;
static int l;
static double bb,ff;
static int j;
static double h;
static int m,ma,ia,i1;
static double g,hr,hi;
/*

     SUBROUTINE CTQL2 COMPUTES THE EIGENVALUES AND
     EIGENVECTORS OF A COMPLEX HERMITIAN TRIDIAGONAL
     MATRIX

     PARAMETERS

     N    --- ORDER OF MATRIX
     ND   --- DIMENSION OF A AND B
     D    --- DIAGONAL GIVEN
     E    --- REAL PART OF CODIAGONAL GIVEN
              (LAST N-1 LOCATIONS)
     F    --- IMAG PART OF THE LOWER CODIAGONAL
     A    --- REAL PART OF EIGENVECTORS
     B    --- IMAG PART OF EIGENVECTORS
     FAIL --- RECEIVES VALUE OF 1 INSTEAD OF ZERO
              IF SOME EIGENVALUE TAKES MORE THAN 30
              ITERATIONS.


     EIGENVALUES ARE OBTAINED IN INCREASING OF
     MAGNITUDE IN VECTOR D, EIGENVECTORS ARE STORED
     BY COLUMNS.  ARRAYS A AND B SHOULD BE PRESET TO
     SOME UNITARY MATRIX SUCH AS THE IDENTITY MATRIX
     OR THE MATRIX PRODUCED BY CTRED2.


     PROGRAMMED BY E.  ZAKRAJSEK
     JUNE 21, 1974




     ***************************************
     *                                     *
     * NEXT LINE OF PROGRAM DEFINES        *
     * MACHINE DEPENDENT CONSTANT CHEP     *
     * DEFINED AS THE SMALLEST REAL        *
     * NUMBER FOR WHICH                    *
     *                                     *
     *        1.0+CHEP .GT. 1.0            *
     *                                     *
     ***************************************

*/
    chep = pow(2.0e0,(double)-56);
/*
     FIRST MAKE REAL CODIAGONAL MOVED DOWN
     TO FIRST LOCATION

*/
    if(*n == 1) goto S12;
    for(k=2; k<=*n; k++) {
        r = sqrt(e[k-1]*e[k-1]+f[k-1]*f[k-1]);
        if(fabs(r) < 1.e-50) r = 0.e0;
        if(r == 0.0e0) goto S11;
        c = e[k-1]/r;
        s = f[k-1]/r;
/*
     ACCUMULATE ROTATION

*/
        for(i=1; i<=*n; i++) {
            p = *(a+i-1+(k-1)**nd)*c-*(b+i-1+(k-1)**nd)*s;
            *(b+i-1+(k-1)**nd) = *(a+i-1+(k-1)**nd)*s+*(b+i-1+(k-1)**nd)*c;
            *(a+i-1+(k-1)**nd) = p;
        }
/*
     TRANSFORM NEXT E

*/
        if(k == *n) goto S11;
        l = k+1;
        p = e[l-1]*c-f[l-1]*s;
        f[l-1] = e[l-1]*s+f[l-1]*c;
        e[l-1] = p;
S11:
        e[k-2] = r;
    }
S12:
    e[*n-1] = 0.e0;
/*
     INITIALIZE

*/
    bb = 0.e0;
    ff = 0.e0;
    *fail = 1;
/*
     MAIN LOOP

*/
    for(l=1; l<=*n; l++) {
        j = 0;
        h = chep*(fabs(d[l-1])+fabs(e[l-1]));
        if(bb < h) bb = h;
/*
     LOOK FOR SMALL SUBDIAGONAL ELEMENT

*/
        for(m=l; m<=*n; m++) {
            if(fabs(e[m-1]) <= bb) goto S21;
        }
S21:
        if(m == l) goto S31;
/*
     NEXT ITERATION

*/
S24:
        if(j == 30) return;
        j = j+1;
/*
     FORM SHIFT

*/
        p = (d[l]-d[l-1])/(2.0e0*e[l-1]);
        r = sqrt(1.0e0+p*p);
        h = d[l-1]-e[l-1]/(p+fdsign(r,p));
/*
*/
        for(i=l; i<=*n; i++) {
            d[i-1] = d[i-1]-h;
        }
        ff = ff+h;
/*
     QL TRANSFORMATION

*/
        p = d[m-1];
        c = 1.0e0;
        s = 0.0e0;
        ma = m-1;
/*
*/
        for(ia=l; ia<=ma; ia++) {
            i = ma-ia+l;
            i1 = i+1;
            g = c*e[i-1];
            h = c*p;
            if(fabs(p) < fabs(e[i-1])) goto S26;
/*
*/
            c = e[i-1]/p;
            r = sqrt(c*c+1.0e0);
            e[i1-1] = s*p*r;
            s = c/r;
            c = 1.0e0/r;
            goto S27;
/*
*/
S26:
            c = p/e[i-1];
            r = sqrt(c*c+1.e0);
            e[i1-1] = s*e[i-1]*r;
            s = 1.0e0/r;
            c = c/r;
/*
*/
S27:
            p = c*d[i-1]-s*g;
            d[i1-1] = h+s*(c*g+s*d[i-1]);
/*
     FORM VECTOR

*/
            for(k=1; k<=*n; k++) {
                hr = *(a+k-1+(i1-1)**nd);
                hi = *(b+k-1+(i1-1)**nd);
                *(a+k-1+(i1-1)**nd) = s**(a+k-1+(i-1)**nd)+c*hr;
                *(b+k-1+(i1-1)**nd) = s**(b+k-1+(i-1)**nd)+c*hi;
                *(a+k-1+(i-1)**nd) = c**(a+k-1+(i-1)**nd)-s*hr;
                *(b+k-1+(i-1)**nd) = c**(b+k-1+(i-1)**nd)-s*hi;
            }
/*
*/
        }
/*
*/
        e[l-1] = s*p;
        d[l-1] = c*p;
        if(fabs(e[l-1]) > bb) goto S24;
/*
     ROOT FOUND

*/
S31:
        d[l-1] = d[l-1]+ff;
    }
/*
     ORDER EIGENVALUES AND EIGENVECTORS

*/
    for(i=1; i<=*n; i++) {
        k = i;
        for(j=i; j<=*n; j++) {
            if(d[j-1] < d[k-1]) k = j;
        }
/*

*/
        if(k == i) goto S42;
        p = d[i-1];
        d[i-1] = d[k-1];
        d[k-1] = p;
/*
*/
        for(j=1; j<=*n; j++) {
            p = *(a+j-1+(i-1)**nd);
            *(a+j-1+(i-1)**nd) = *(a+j-1+(k-1)**nd);
            *(a+j-1+(k-1)**nd) = p;
/*
*/
            p = *(b+j-1+(i-1)**nd);
            *(b+j-1+(i-1)**nd) = *(b+j-1+(k-1)**nd);
            *(b+j-1+(k-1)**nd) = p;
        }
S42:;
    }
    *fail = 0;
    return;
}
