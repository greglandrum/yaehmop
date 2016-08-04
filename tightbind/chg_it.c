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
*   This is the stuff for doing charge iteration
*
*  created:  greg landrum  October 1995
*
*****************************************************************************/
#include "bind.h"




/****************************************************************************
 *
 *                   Procedure update_chg_it_parms
 *
 * Arguments: details: pointer to detail_type
 *               cell: pointer to cell type
 *          AO_occups: pointer to real
 *          converged: pointer to int
 *           num_orbs: int
 *      orbital_lookup_table: pointer to int
 *
 *
 * Returns: none
 *
 * Action:  This uses the AO occupations in 'AO_occups to update
 *   atomic Hii's using the charge iteration formula.
 *
 ****************************************************************************/
void update_chg_it_parms(details,cell,AO_occups,converged,num_orbs,
                         orbital_lookup_table)
  detail_type *details;
  cell_type *cell;
  real *AO_occups;
  int *converged;
  int num_orbs;
  int *orbital_lookup_table;
{
  static real *AO_store=0;
  static int num_calls=0;
  atom_type *atom;
  chg_it_parm_type *parms;
  int i,j,num_atoms;
  real new_Hss, new_Hpp, new_Hdd;
  real old_Hss, old_Hpp, old_Hdd;
  real tot_s_occup, tot_p_occup, tot_d_occup;
  real old_s_occup, old_p_occup, old_d_occup;
  real damped_s_occup, damped_p_occup, damped_d_occup;
  real tot_chg,old_chg;
  real denom,adjust,lambda;
  int orb_tab;
  int begin_atom,end_atom;

  /* get storage space if this is the first call */
  if( !AO_store ){
    AO_store = (real *)calloc(num_orbs,sizeof(real));
    if( !AO_store ) fatal("Can't get AO_store memory.");
  }
  parms = &(details->chg_it_parms);

  num_calls++;

  num_atoms = cell->num_atoms;
  /* loop over atoms */
  *converged = 1;

  fprintf(output_file,";Charge Iteration Step... New parms:\n");
  for(i=0;i<num_atoms;i++){
    /* do we need to update this atom? */
    if( cell->atoms[i].chg_it_vary ){
      atom = &(cell->atoms[i]);
      /* find the orbitals */
      find_atoms_orbs(num_orbs,cell->num_atoms,i,orbital_lookup_table,&begin_atom,
                      &end_atom);
      if( begin_atom >= 0 ){
        orb_tab = begin_atom;
        /********

          figure out the occupations

        *********/
        tot_s_occup = 0;
        old_s_occup = 0;
        if(atom->ns){
          tot_s_occup = AO_occups[orb_tab];
          old_s_occup = AO_store[orb_tab];
          orb_tab++;
          old_Hss = atom->coul_s;
        }
        tot_p_occup = 0;
        old_p_occup = 0;
        if(atom->np) {
          for(j=0;j<3;j++){
            tot_p_occup += AO_occups[orb_tab];
            old_p_occup += AO_store[orb_tab];
            orb_tab++;
          }
          old_Hpp = atom->coul_p;
        }
        tot_d_occup = 0;
        old_d_occup = 0;
        if(atom->nd) {
          for(j=0;j<5;j++){
            tot_d_occup += AO_occups[orb_tab];
            old_d_occup += AO_store[orb_tab];
            orb_tab++;
          }
          old_Hdd = atom->coul_d;
        }
      }


      /* find denom.  this is used for convergence checking */
      denom = 0;
      if( atom->ns && fabs(old_s_occup - tot_s_occup) > denom)
        denom = fabs(old_s_occup - tot_s_occup);
      if( atom->np && fabs(old_p_occup - tot_p_occup) > denom)
        denom = fabs(old_p_occup - tot_p_occup);
      if( atom->nd && fabs(old_d_occup - tot_d_occup) > denom)
        denom = fabs(old_d_occup - tot_d_occup);

      /* figure out the lambda */
      if( parms->variable_step ){
        if( num_calls == 1 ) adjust = denom;
        else{
          adjust = parms->damp1;
        }

        lambda = adjust / denom;
      }else {
        lambda = parms->lambda;
      }

fprintf(stderr,"It: %d, denom: %lf tol: %lf\n",num_calls,denom,parms->tolerance);
      if( denom > parms->tolerance ){
        *converged = 0;

        /* make sure that the lambda isn't too big, so we don't explode */
        if( lambda > parms->lampri ) lambda = parms->lampri;

#if 0
        /* damp the AO occupations */
        if( num_calls != 1 ){
          damped_s_occup = old_s_occup + lambda * (tot_s_occup - old_s_occup);
          damped_p_occup = old_p_occup + lambda * (tot_p_occup - old_p_occup);
          damped_d_occup = old_d_occup + lambda * (tot_d_occup - old_d_occup);
        } else{
          damped_s_occup = tot_s_occup;
          damped_p_occup = tot_p_occup;
          damped_d_occup = tot_d_occup;
        }
#endif


        /* figure out the (damped) net charge */
        old_chg = atom->num_valence - (old_s_occup + old_p_occup +
                                       old_d_occup);
        tot_chg = atom->num_valence - (tot_s_occup + tot_p_occup +
                                       tot_d_occup);
#if 0
        if( num_calls != 1 ){
          tot_chg = old_chg + lambda*(tot_chg - old_chg);
        } else{
          tot_chg = lambda*tot_chg;
        }
#endif

        /* now update the Hii's */
        new_Hss = tot_chg*tot_chg*atom->s_A + tot_chg*atom->s_B
          + atom->s_C;
        new_Hpp = tot_chg*tot_chg*atom->p_A + tot_chg*atom->p_B
          + atom->p_C;
        new_Hdd = tot_chg*tot_chg*atom->d_A + tot_chg*atom->d_B
          + atom->d_C;

        atom->coul_s = atom->coul_s + lambda*(new_Hss - atom->coul_s);
        atom->coul_p = atom->coul_p + lambda*(new_Hpp - atom->coul_p);
        atom->coul_d = atom->coul_d + lambda*(new_Hdd - atom->coul_d);


        /* write out the parameters that we're varying */
        fprintf(output_file,"%d %s ",i,atom->symb);
        if( atom->ns ) fprintf(output_file,"S:% -6.4lf ",new_Hss);
        if( atom->np ) fprintf(output_file,"P:% -6.4lf ",new_Hpp);
        if( atom->nd ) fprintf(output_file,"D:% -6.4lf ",new_Hdd);
        fprintf(output_file,"\n");
      }
    }
  }

  /* copy over the AO occupations */
  bcopy((char *)AO_occups,(char *)AO_store,num_orbs*sizeof(real));


  if( num_calls == parms->max_it ) *converged = 1;


  if( *converged )
    fprintf(stderr,"Charge iteration converged after %d (of %d max) iterations\n",
                           num_calls,parms->max_it);

}
