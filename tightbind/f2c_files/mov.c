/*******************************************************
*      Copyright (C) 1995 Greg Landrum
*
*  This file is part of yaehmop.
*
*   This is free software.
*
*  Permission is granted to modify, or otherwise fold, spindle, and mutilate this
*    code provided all copyright notices are left intact.
*
*  This code may be distributed to your heart's content, in whatever form,
*    provided no fee is charged for the distribution, all copyright notices are
*    left intact, and the source is distributed (without fee) along with any
*    binaries to anyone who requests it.
*
*  There are, of course, no warranties at all on this program.
*
********************************************************************/



/****************************************************************************
*
*     this file contains stuff for evaluating overlaps
*
*  created:  greg landrum  August 1993
*
*****************************************************************************/
#include "bind.h"

#ifdef COMMON_BLK_HACK
/*********
  this is a global variable which is used to fill a FORTRAN common block
  used by the f77 subroutines that this procedure calls.

  This may or may not be a total hack....

  suffice it to say that this probably isn't the most stable thing to
   do and you should most definitely test it with your compiler before
   trusting it!
*********/

struct {
  real sk1,sk2,r;
  int l1,l2,m,n1,n2,max;
} loclap;
#endif


#if 0
#define lovlap lovlap_
#define abfns abfns_
#endif

/****************************************************************************
*
*                   Procedure mov
*
* Arguments: sigma,pi,delta,phi: pointers to reals
*                 which1,which2: int
*                          dist: real
*        q_num1,q_num2,l1,l2,nn: int
*                         atoms: atom_type
*
*
* Returns: none
*
* Action:  does whatever MOV did in the original program
*
*   comments will follow the clue when I get one
*
****************************************************************************/
void mov(sigma,pi,delta,phi,which1,which2,dist,q_num1,q_num2,
         l1,l2,nn,atoms)
  real *sigma,*pi,*delta,*phi,dist;
  int which1,which2,q_num1,q_num2,l1,l2,nn;
  atom_type *atoms;
{
  real oldsk1=-1,oldsk2=-1,oldr=-1;

  int i,j,num_zeta1,num_zeta2;
  real coeff1,coeff2;

#ifndef COMMON_BLK_HACK
  real sk1,sk2,r;
  int m,n1,n2,max;
#endif

  real A_fn_values[30], B_fn_values[30];
  real ang_ind_overlap[4];


/*
fprintf(stderr,"MOV: (i,j,dist) %d %d %lf\n",which1,which2,dist);
*/
  /* initialize the components of the overlap to zero */
  *sigma=*pi=*delta=*phi=0.0;
  ang_ind_overlap[0]=ang_ind_overlap[1]=ang_ind_overlap[2]=
    ang_ind_overlap[3]=0.0;

  /* figure out whether or not we are using double zeta f'ns */
  if( l1 == 2 && atoms[which1].coeff2 != 0 ) num_zeta1 = 2;
  else num_zeta1 = 1;
  if( l2 == 2 && atoms[which2].coeff2 != 0 ) num_zeta2 = 2;
  else num_zeta2 = 1;

  max = q_num1+q_num2;

  /* I unrolled the loop in the original routine */

  /* get the coefficients */
  switch(l1){
  case 0:
    sk1 = atoms[which1].exp_s;
    break;
  case 1:
    sk1 = atoms[which1].exp_p;
    break;
  case 2:
    sk1 = atoms[which1].exp_d;
    break;
  case 3:
    sk1 = atoms[which1].exp_f;
    break;
  }

  switch(l2){
  case 0:
    sk2 = atoms[which2].exp_s;
    break;
  case 1:
    sk2 = atoms[which2].exp_p;
    break;
  case 2:
    sk2 = atoms[which2].exp_d;
    break;
  case 3:
    sk2 = atoms[which2].exp_f;
    break;
  }

  if(!(atoms[which1].coeff1) || l1 <= 1) coeff1 = 1;
  else coeff1 = atoms[which1].coeff1;
  if(!(atoms[which2].coeff1) || l2 <= 1) coeff2 = 1;
  else coeff2 = atoms[which2].coeff1;


  /* call the routine to evaluate the A & B functions */
  abfns(A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,
        &q_num2,&max);

  /* I don't know what nn is yet, but I'll figure it out eventually */
  for(i=0;i<nn;i++){
    m = i;
    lovlap(&(ang_ind_overlap[i]),A_fn_values,B_fn_values,&sk1,&sk2,
           &dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
  }


/*
fprintf(stderr,"COEFFS: %lf %lf \n",coeff1,coeff2);
*/
  /* add in the contributions we have found thus far */
  *sigma += coeff1*coeff2*ang_ind_overlap[0];
  *pi += coeff1*coeff2*ang_ind_overlap[1];
  *delta += coeff1*coeff2*ang_ind_overlap[2];
  *phi += coeff1*coeff2*ang_ind_overlap[3];

/*
fprintf(stderr,"sigma,pi: %6.4lf %6.4lf\n",*sigma,*pi);
*/
  if( num_zeta2 == 2){

    sk2 = atoms[which2].exp_d2;
    coeff2 = atoms[which2].coeff2;
    abfns(A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);

    for(i=0;i<nn;i++){
      m = i;
      lovlap(&(ang_ind_overlap[i]),A_fn_values,B_fn_values,&sk1,&sk2,
             &dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
    }

    /* add in the contributions we have found thus far */
/*
fprintf(stderr,"COEFFS: %lf %lf \n",coeff1,coeff2);
*/
    *sigma += coeff1*coeff2*ang_ind_overlap[0];
    *pi += coeff1*coeff2*ang_ind_overlap[1];
    *delta += coeff1*coeff2*ang_ind_overlap[2];
    *phi += coeff1*coeff2*ang_ind_overlap[3];
  }

  if( num_zeta1 == 2 ){
    sk1 = atoms[which1].exp_d2;
    coeff1 = atoms[which1].coeff2;
    abfns(A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);

    for(i=0;i<nn;i++){
      m = i;
      lovlap(&(ang_ind_overlap[i]),A_fn_values,B_fn_values,&sk1,&sk2,
             &dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
    }

    /* add in the contributions we have found thus far */
/*
fprintf(stderr,"COEFFS: %lf %lf \n",coeff1,coeff2);
*/
    *sigma += coeff1*coeff2*ang_ind_overlap[0];
    *pi += coeff1*coeff2*ang_ind_overlap[1];
    *delta += coeff1*coeff2*ang_ind_overlap[2];
    *phi += coeff1*coeff2*ang_ind_overlap[3];


    if( num_zeta2 == 2 ){
      /* we've already done the second zeta for orbital 2, so now do the first */
      sk2 = atoms[which2].exp_d;
      coeff2 = atoms[which2].coeff1;
      abfns(A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);

      for(i=0;i<nn;i++){
        m = i;
        lovlap(&(ang_ind_overlap[i]),A_fn_values,B_fn_values,&sk1,&sk2,
               &dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
      }

      /* add in the contributions we have found thus far */
/*
fprintf(stderr,"COEFFS: %lf %lf \n",coeff1,coeff2);
*/
      *sigma += coeff1*coeff2*ang_ind_overlap[0];
      *pi += coeff1*coeff2*ang_ind_overlap[1];
      *delta += coeff1*coeff2*ang_ind_overlap[2];
      *phi += coeff1*coeff2*ang_ind_overlap[3];

    }
  }
}
