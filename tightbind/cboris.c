/*
Produced by gmFortran V30.59(10/26/17) on 9/18/18 at 9:02:32
*/
#define LPROTOTYPE
#include "fortran.h"
void cboris(int *n,int *nd,double *a,double *b,double *c,double *d,double *e,
    double *f,int *fail)
{
extern void cchol(),ctred2(),ctql2();
static int lf,i,ia,j,k,ja,ii;
/*

  for YAeHMOP this common block is not needed
              gL
      COMMON/PASS/NATOM,NDIM,LB,LE,NVAL,ECUT,EERR,IDEL


 SUBROUTINE CBORIS SOLVES THE COMPLEX EIGENVALUE
 PROBLEM

      A*X=B*X*LAMBDA

 WHERE

 A     --- GIVEN MATRIX OF ORDER N
 B     --- GIVEN POSITIVE DEFINITE MATRIX OF ORDER N
 X     --- EIGENVECTORS, COLUMN BY COLUMN,
           NORMALIZED TO (BX,X)=I
 LAMBDA --- EIGENVALUES IN INCREASING ORDER

 MATRICES A AND B SHOULD BE GIVEN IN
 ARRAYS A AND B RESPECTIVELY IN THE FOLLOWING WAY

 DIAGONAL ELEMENTS IN THE DIAGONAL
 REAL PARTS OF THE LOWER TRIANGLE IN THE LOWER
 TRIANGLE
 IMAG PARTS OF THE LOWER TRIANGLE IN THE UPPER
 TRIANGLE

 REAL PARTS OF EIGENVECTORS ARE STORED IN A,
 IMAG PARTS IN ARRAY C.  ARRAY B IS DESTROYED DURING
 COMPUTATION.  IT HOLDS THE LOWER TRIANGLE OF THE
 CHOLESKI DECOMPOSITION OF MATRIX B (SEE
 SUBROUTINE CCHOL).

 ALL MATRICES ARE OF ORDER N, DECLARED IN THE
 CALLING PROGRAM WITH DIMENSION ND WHICH NEED
 NOT TO BE EQUAL TO N.

 EIGENVALUES ARE STORED IN ARRAY D.  ARRAYS E
 AND F ARE USED FOR INTERMEDIATE RESULTS.

 FAIL GETS THE FOLLOWING VALUES

  0 --- COMPUTATION FINISHED SUCCESFULLY
  1 --- B NOT POSITIVE DEFINITE
  2 --- QR ALGORITHM DOES NOT CONVERGE

 SUBROUTINES USED

     CCHOL
     CTRED2
     CTQL2


 PROGRAMMED BY E. ZAKRAJSEK
 JUNE 21,1974




 DECOMPOSE MATRIX B

      WRITE(4,1918)
1918  FORMAT(' WELCOME TO CBORIS')

      CALL TIME(0)

*/
    *fail = 1;
    cchol(n,nd,b,&lf);
    if(lf != 0) goto S999;
/*
 MOVE MATRIX A

*/
    for(i=1; i<=*n; i++) {
        *(c+i-1+(i-1)**nd) = 0.e0;
        if(i == 1) goto S11;
        ia = i-1;
        for(j=0; j<ia; j++) {
            *(c+i-1+j**nd) = *(a+j+(i-1)**nd);
            *(c+j+(i-1)**nd) = -*(a+j+(i-1)**nd);
            *(a+j+(i-1)**nd) = *(a+i-1+j**nd);
        }
S11:;
    }
/*
  COMPUTE (L(-1)*A)

*/
    for(j=0; j<*n; j++) {
        for(i=1; i<=*n; i++) {
            if(i == 1) goto S21;
            ia = i-1;
            for(k=0; k<ia; k++) {
                *(a+i-1+j**nd) = *(a+i-1+j**nd)-*(a+k+j**nd)**(b+i-1+k**nd)+*(c+
                  k+j**nd)**(b+k+(i-1)**nd);
                *(c+i-1+j**nd) = *(c+i-1+j**nd)-*(a+k+j**nd)**(b+k+(i-1)**nd)-*
                  (c+k+j**nd)**(b+i-1+k**nd);
            }
S21:
            *(a+i-1+j**nd) = *(a+i-1+j**nd)/ *(b+i-1+(i-1)**nd);
            *(c+i-1+j**nd) = *(c+i-1+j**nd)/ *(b+i-1+(i-1)**nd);
        }
    }
/*
  COMPUTE  A*L(-H)

*/
    for(i=1; i<=*n; i++) {
        for(j=1; j<=i; j++) {
            if(j == 1) goto S31;
            ja = j-1;
            for(k=0; k<ja; k++) {
                *(a+i-1+(j-1)**nd) = *(a+i-1+(j-1)**nd)-*(a+i-1+k**nd)**(b+j-1+
                  k**nd)-*(c+i-1+k**nd)**(b+k+(j-1)**nd);
                *(c+i-1+(j-1)**nd) = *(c+i-1+(j-1)**nd)+*(a+i-1+k**nd)**(b+k+(j
                  -1)**nd)-*(c+i-1+k**nd)**(b+j-1+k**nd);
            }
S31:
            *(a+i-1+(j-1)**nd) = *(a+i-1+(j-1)**nd)/ *(b+j-1+(j-1)**nd);
            *(c+i-1+(j-1)**nd) = *(c+i-1+(j-1)**nd)/ *(b+j-1+(j-1)**nd);
        }
    }
/*
     PUT MATRIX TOGETHER INTO A

*/
    for(i=1; i<=*n; i++) {
        if(i == *n) goto S41;
        ia = i+1;
        for(j=ia; j<=*n; j++) {
            *(a+i-1+(j-1)**nd) = *(c+j-1+(i-1)**nd);
        }
S41:;
    }
/*
     DIAGONALIZE A

*/
    *fail = 2;
    ctred2(n,nd,a,c,d,e,f);
    ctql2(n,nd,d,e,f,a,c,&lf);
    if(lf != 0) goto S999;
/*
     COMPUTE L(-H)*A

*/
    for(j=1; j<=*n; j++) {
        for(ii=1; ii<=*n; ii++) {
            i = *n-ii+1;
            if(i == *n) goto S51;
            ia = i+1;
            for(k=ia-1; k<*n; k++) {
                *(a+i-1+(j-1)**nd) = *(a+i-1+(j-1)**nd)-*(a+k+(j-1)**nd)**(b+k+
                  (i-1)**nd)-*(c+k+(j-1)**nd)**(b+i-1+k**nd);
                *(c+i-1+(j-1)**nd) = *(c+i-1+(j-1)**nd)+*(a+k+(j-1)**nd)**(b+i
                  -1+k**nd)-*(c+k+(j-1)**nd)**(b+k+(i-1)**nd);
            }
S51:
            *(a+i-1+(j-1)**nd) = *(a+i-1+(j-1)**nd)/ *(b+i-1+(i-1)**nd);
            *(c+i-1+(j-1)**nd) = *(c+i-1+(j-1)**nd)/ *(b+i-1+(i-1)**nd);
        }
    }
    *fail = 0;
/*
  999 WRITE(4,1919)
 1919 FORMAT(' DONE CBORIS')

      CALL TIME(0)

*/
S999:
    return;
}
