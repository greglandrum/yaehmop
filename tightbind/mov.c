/*******************************************************

Copyright (C) 1995 Greg Landrum
All rights reserved

This file is part of yaehmop.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

********************************************************************/

/********************************************************************************
*
*     this file contains stuff for evaluating overlaps
*
*  created:  greg landrum  August 1993
*
********************************************************************************/

/***
  Edit History:

  March '98: WG
  - f orbital overlaps added

***/
#include "bind.h"

/********************************************************************************
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
********************************************************************************/
void mov(sigma,pi,delta,phi,which1,which2,dist,q_num1,q_num2,
         l1,l2,atoms)
  real *sigma,*pi,*delta,*phi,dist;
  int which1,which2,q_num1,q_num2,l1,l2;
  atom_type *atoms;
{
  int i,j,num_zeta1,num_zeta2;
  real coeff_1,coeff_2,sk1,sk2,r;

  real A_fn_values[30], B_fn_values[30];
  real ang_ind_overlap[4];

  int loopvar,max, m=0, nn;  /* definitions for abfns.f and lovlap.f */
  max = q_num1 + q_num2;

  if(l1>l2){
    nn=l2;
  }
  else{
    nn=l1;
  }

  /* initialize the components of the overlap to zero */
  *sigma=*pi=*delta=*phi=0.0;
  ang_ind_overlap[0]=ang_ind_overlap[1]=ang_ind_overlap[2]=ang_ind_overlap[3]=0.0;

  /* figure out whether or not we are using double zeta f'ns */

  if( (l1 == 2 && atoms[which1].coeff_d2 != 0) || (l1 == 3 && atoms[which1].coeff_f2 != 0) ) num_zeta1 = 2;
  else num_zeta1 = 1;
  if( (l2 == 2 && atoms[which2].coeff_d2 != 0) || (l2 == 3 && atoms[which2].coeff_f2 != 0) ) num_zeta2 = 2;
  else num_zeta2 = 1;

  /* get the exponents */

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

  /* deal with zeta1 - zeta1 overlap */

  if(!(atoms[which1].coeff_d1 || atoms[which1].coeff_f1) || l1 <= 1){
    coeff_1 = 1;
  }
  else{
    if(l1 == 2) coeff_1 = atoms[which1].coeff_d1;
    if(l1 == 3) coeff_1 = atoms[which1].coeff_f1;
  }

  if(!(atoms[which2].coeff_d1 || atoms[which2].coeff_f1) || l2 <= 1){
    coeff_2 = 1;
  }
  else{
    if(l2 == 2) coeff_2 = atoms[which2].coeff_d1;
    if(l2 == 3) coeff_2 = atoms[which2].coeff_f1;
  }

  /* call the routine to evaluate the A & B functions */

  abfns(A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);

  /* test print of A and B functions */

/*    fprintf(stderr,"zeta1[%f]-zeta1[%f] call ...\n",sk1,sk2);
    loopvar=0;
    while(loopvar <= max)
      {
        fprintf(stderr,"A(%d)= %.16f \t B(%d)= %f\n",loopvar,A_fn_values[loopvar],loopvar,B_fn_values[loopvar]);
        loopvar++;
      }
*/

  /* call the routine to evaluate sigma,pi,... overlaps in the local ref frame */

  for(i=0;i<=nn;i++)
    {
      m=i;
      lovlap(&(ang_ind_overlap[i]),A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
    }

  /* add in the contributions we have found thus far */

  *sigma += coeff_1*coeff_2*ang_ind_overlap[0];
  *pi += coeff_1*coeff_2*ang_ind_overlap[1];
  *delta += coeff_1*coeff_2*ang_ind_overlap[2];
  *phi += coeff_1*coeff_2*ang_ind_overlap[3];

  /* now do zeta1 - zeta2 overlap if applicable */

  if( num_zeta2 == 2){
    if( l2 == 2){
      sk2 = atoms[which2].exp_d2;
      coeff_2 = atoms[which2].coeff_d2;
    }
    if( l2 == 3){
      sk2 = atoms[which2].exp_f2;
      coeff_2 = atoms[which2].coeff_f2;
    }

    abfns(A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);

    for(i=0;i<=nn;i++){
      m=i;
      lovlap(&(ang_ind_overlap[i]),A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
    }

    *sigma += coeff_1*coeff_2*ang_ind_overlap[0];
    *pi += coeff_1*coeff_2*ang_ind_overlap[1];
    *delta += coeff_1*coeff_2*ang_ind_overlap[2];
    *phi += coeff_1*coeff_2*ang_ind_overlap[3];
  }

  /* now do zeta2 - zeta2 */

  if( num_zeta1 == 2 ){
    if( l1 == 2 ){
      sk1 = atoms[which1].exp_d2;
      coeff_1 = atoms[which1].coeff_d2;
    }
    if( l1 == 3 ){
      sk1 = atoms[which1].exp_f2;
      coeff_1 = atoms[which1].coeff_f2;
    }

    abfns(A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);

    for(i=0;i<=nn;i++){
      m=i;
      lovlap(&(ang_ind_overlap[i]),A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
    }

    *sigma += coeff_1*coeff_2*ang_ind_overlap[0];
    *pi += coeff_1*coeff_2*ang_ind_overlap[1];
    *delta += coeff_1*coeff_2*ang_ind_overlap[2];
    *phi += coeff_1*coeff_2*ang_ind_overlap[3];

    /* finally do zeta2 - zeta1 */

    if( num_zeta2 == 2 ){
      if( l2 == 2){
        sk2 = atoms[which2].exp_d;
        coeff_2 = atoms[which2].coeff_d1;
      }
      if( l2 == 3){
        sk2 = atoms[which2].exp_f;
        coeff_2 = atoms[which2].coeff_f1;
      }

      abfns(A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);

      for(i=0;i<=nn;i++){
        m=i;
        lovlap(&(ang_ind_overlap[i]),A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
      }
      *sigma += coeff_1*coeff_2*ang_ind_overlap[0];
      *pi += coeff_1*coeff_2*ang_ind_overlap[1];
      *delta += coeff_1*coeff_2*ang_ind_overlap[2];
      *phi += coeff_1*coeff_2*ang_ind_overlap[3];
    }
  }
}
