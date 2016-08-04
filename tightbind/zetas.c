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

/****************************************************************************
*
*   This is the stuff for doing self consistent zeta modifications
*
*  created:  greg landrum  January 1994
*
*****************************************************************************/
#include "bind.h"




/****************************************************************************
 *
 *                   Procedure update_zetas
 *
 * Arguments:  cell: pointer to cell type
 *         net_chgs: pointer to real
 *         zeta_tol: a real
 *        converged: pointer to int
 *            reset: a char
 *
 *
 * Returns: none
 *
 * Action:  This uses the charges in 'net_chgs to update the atomic
 *       zeta in the atom list stored in 'cell.
 *
 *   if the total change in zeta values is less than zeta_tol then 'converged
 *     is set to 1.
 *
 *   if 'reset is nonzero then the last_charge array is zeroed and no other action.
 *    is taken.
 *
 ****************************************************************************/
void update_zetas(cell_type *cell,real *net_chgs,real zeta_tol,int *converged,char reset)
{
  static real *last_chgs=0;
  static int num_calls=0;
  static int max_calls=100;
  atom_type *atom;
  int i,num_atoms;
  real total_delta;
  real delta_q;
  real delta;

  if( reset ) num_calls = 0;
  num_calls++;

  num_atoms = cell->num_atoms;

  /* get space for the array to store the last charge values */
  if( !last_chgs ){
    last_chgs = (real *)calloc(num_atoms,sizeof(real));
    if(!last_chgs)fatal("Can't allocate last_charge array in update_zetas.");
  }

  /* zero out the last charge array if that's needed */
  if( reset ){
    bzero((char *)last_chgs,num_atoms*sizeof(real));
    return;
  }


  total_delta = 0.0;

  fprintf(output_file,"\n;      ->>>>>>>> Updating zeta values based on net charges <<<<<<<-\n");
  for(i=0;i<num_atoms;i++){
    atom = &(cell->atoms[i]);

    delta_q = net_chgs[i] - last_chgs[i];
    if(atom->ns != 1) delta = ZETA_SCALE * delta_q;
    else delta = ZETA_SCALE_1S * delta_q;
    last_chgs[i] = net_chgs[i];

    fprintf(output_file,"Atom: %d  %s     Delta Charge = %lg    Delta zeta = %lg\n",i,atom->symb,delta_q,delta);

    total_delta += fabs(delta);

    if( atom->ns != 0){
      fprintf(output_file,"\t %dS: %lg -> ",atom->ns,atom->exp_s);
      atom->exp_s += .1*delta/atom->ns;
      fprintf(output_file,"%lg\n ",atom->exp_s);
    }
    if( atom->np != 0){
      fprintf(output_file,"\t %dP: %lg -> ",atom->np,atom->exp_p);
      atom->exp_p += .1*delta/atom->np;
      fprintf(output_file,"%lg\n ",atom->exp_p);
    }
    /* this needs to be changed.... */
    if( atom->nd != 0){
/*      atom->exp_d += delta;*/
    }
  }

  fprintf(output_file,"Total zeta change this cycle: %lg  tolerance: %lg\n",total_delta,
          zeta_tol);

fprintf(stderr,"Total zeta change this cycle: %lg  tolerance: %lg\n",total_delta,
          zeta_tol);
  fprintf(output_file,"\n\n ;-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n\n");

  if( total_delta < zeta_tol || num_calls == max_calls) {
    *converged  = 1;
  }

}
