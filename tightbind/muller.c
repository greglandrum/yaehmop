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
 *     this file contains everything needed to do an iterative eHT
 *      calculation using Edgar Muller's technique
 *
 *
 *  created:  greg landrum  July 1996
 *
 *****************************************************************************/
/***
  Edit History:

  March '98: WG
  - zeta coefficients for d's renamed (coeff1,coeff2 -> coeff_d1,coeff_d2)
***/

#include "bind.h"




/****************************************************************************
 *
 *                   Procedure calc_muller_init_parms
 *
 * Arguments: atom: pointer to atom type
 *
 *
 * Returns: none
 *
 * Action:  This uses the specified initial occupations of 'atom
 *   to generate an initial parameter set.
 *
 ****************************************************************************/
void calc_muller_init_parms(atom_type *atom)
{
  real *E_vals,*Z_vals;

  if( atom->nd ){
    E_vals = atom->muller_d_E;
    Z_vals = atom->muller_d_Z;
    atom->coul_d = E_vals[0] + E_vals[1] * atom->init_s_occup +
      E_vals[2] * atom->init_p_occup + E_vals[3] * atom->init_d_occup +
        E_vals[4] * (atom->init_s_occup+atom->init_p_occup)*atom->init_d_occup +
          E_vals[5]*atom->init_d_occup*atom->init_d_occup +
            E_vals[6]*(atom->init_s_occup+atom->init_p_occup)*
              (atom->init_s_occup+atom->init_p_occup);
    atom->exp_d = Z_vals[0] + Z_vals[1]*atom->init_s_occup +
      Z_vals[2]*atom->init_p_occup + Z_vals[3]*atom->init_d_occup;
    atom->exp_d2 = 0;
    atom->coeff_d1 = 1;
    atom->coeff_d2 = 0;

    if( atom->np ){
      E_vals = atom->muller_p_E;
      Z_vals = atom->muller_p_Z;
      atom->coul_p = E_vals[0] + E_vals[1] * atom->init_s_occup +
        E_vals[2] * atom->init_p_occup + E_vals[3] * atom->init_d_occup +
          E_vals[4] * (atom->init_s_occup+atom->init_p_occup)*atom->init_d_occup +
            E_vals[5]*atom->init_d_occup*atom->init_d_occup +
              E_vals[6]*(atom->init_s_occup+atom->init_p_occup)*
                (atom->init_s_occup+atom->init_p_occup);
      atom->exp_p = Z_vals[0] + Z_vals[1]*atom->init_s_occup +
        Z_vals[2]*atom->init_p_occup + Z_vals[3]*atom->init_d_occup;
    }
    if( atom->ns ){
      E_vals = atom->muller_s_E;
      Z_vals = atom->muller_s_Z;
      atom->coul_s = E_vals[0] + E_vals[1] * atom->init_s_occup +
        E_vals[2] * atom->init_p_occup + E_vals[3] * atom->init_d_occup +
          E_vals[4] * (atom->init_s_occup+atom->init_p_occup)*atom->init_d_occup +
            E_vals[5]*atom->init_d_occup*atom->init_d_occup +
              E_vals[6]*(atom->init_s_occup+atom->init_p_occup)*
                (atom->init_s_occup+atom->init_p_occup);
      atom->exp_s = Z_vals[0] + Z_vals[1]*atom->init_s_occup +
        Z_vals[2]*atom->init_p_occup + Z_vals[3]*atom->init_d_occup;
    }
  } else if(atom->np){
    E_vals = atom->muller_p_E;
    Z_vals = atom->muller_p_Z;
    atom->coul_p = E_vals[0] + E_vals[1] * atom->init_s_occup +
      E_vals[2] * atom->init_p_occup +
        E_vals[3]*(atom->init_s_occup+atom->init_p_occup)*
          (atom->init_s_occup+atom->init_p_occup);
    atom->exp_p = Z_vals[0] + Z_vals[1]*atom->init_s_occup +
      Z_vals[2]*atom->init_p_occup;

    if( atom->ns ){
      E_vals = atom->muller_s_E;
      Z_vals = atom->muller_s_Z;
      atom->coul_s = E_vals[0] + E_vals[1] * atom->init_s_occup +
        E_vals[2] * atom->init_p_occup +
          E_vals[3]*(atom->init_s_occup+atom->init_p_occup)*
            (atom->init_s_occup+atom->init_p_occup);
      atom->exp_s = Z_vals[0] + Z_vals[1]*atom->init_s_occup +
        Z_vals[2]*atom->init_p_occup;
    }
  } else if(atom->ns){
    E_vals = atom->muller_s_E;
    Z_vals = atom->muller_s_Z;
    atom->coul_s = E_vals[0] + E_vals[1] * atom->init_s_occup +
      E_vals[2] * atom->init_s_occup * atom->init_s_occup;
    atom->exp_s = Z_vals[0] + Z_vals[1]*atom->init_s_occup +
      Z_vals[2]*atom->init_s_occup*atom->init_s_occup;
  }
}

/****************************************************************************
 *
 *                   Procedure calc_muller_parms
 *
 * Arguments: atom: pointer to atom type
 *  s_occup,p_occup,d_occup: real
 *  s_E,s_Z,p_E,p_Z,d_E,d_Z: pointers to real
 *
 *
 * Returns: none
 *
 * Action:  This uses the specified occupations of 'atom
 *   to generate a new parameter set.
 *
 ****************************************************************************/
void calc_muller_parms(atom_type *atom, real s_occup,real p_occup,real d_occup,
                       real *s_E,real *s_Z,real *p_E,real *p_Z,real *d_E,
                       real *d_Z)
{
  real *E_vals,*Z_vals;

  if( atom->nd ){
    E_vals = atom->muller_d_E;
    Z_vals = atom->muller_d_Z;
    *d_E = E_vals[0] +
      E_vals[1] * s_occup +
        E_vals[2] * p_occup +
          E_vals[3] * d_occup +
            E_vals[4] * (s_occup+p_occup)*d_occup +
              E_vals[5]*d_occup*d_occup +
                E_vals[6]*(s_occup+p_occup)*(s_occup+p_occup);
    *d_Z = Z_vals[0] +
      Z_vals[1]*s_occup +
        Z_vals[2]*p_occup +
          Z_vals[3]*d_occup;

    if( atom->np ){
      E_vals = atom->muller_p_E;
      Z_vals = atom->muller_p_Z;
      *p_E = E_vals[0] +
        E_vals[1] * s_occup +
          E_vals[2] * p_occup +
            E_vals[3] * d_occup +
              E_vals[4] * (s_occup+p_occup)*d_occup +
                E_vals[5]*d_occup*d_occup +
                  E_vals[6]*(s_occup+p_occup)*(s_occup+p_occup);
      *p_Z = Z_vals[0] +
        Z_vals[1]*s_occup +
          Z_vals[2]*p_occup +
            Z_vals[3]*d_occup;
    }
    if( atom->ns ){
      E_vals = atom->muller_s_E;
      Z_vals = atom->muller_s_Z;
      *s_E = E_vals[0] +
        E_vals[1] * s_occup +
          E_vals[2] * p_occup +
            E_vals[3] * d_occup +
              E_vals[4] * (s_occup+p_occup)*d_occup +
                E_vals[5]*d_occup*d_occup +
                  E_vals[6]*(s_occup+p_occup)*(s_occup+p_occup);
      *s_Z = Z_vals[0] +
        Z_vals[1]*s_occup +
          Z_vals[2]*p_occup +
            Z_vals[3]*d_occup;
    }
  } else if(atom->np){
    E_vals = atom->muller_p_E;
    Z_vals = atom->muller_p_Z;
    *p_E = E_vals[0] +
      E_vals[1] * s_occup +
        E_vals[2] * p_occup +
          E_vals[3]*(s_occup+p_occup)*(s_occup+p_occup);
    *p_Z = Z_vals[0] +
      Z_vals[1]*s_occup +
        Z_vals[2]*p_occup;

    if( atom->ns ){
      E_vals = atom->muller_s_E;
      Z_vals = atom->muller_s_Z;
      *s_E = E_vals[0] +
        E_vals[1] * s_occup +
          E_vals[2] * p_occup +
            E_vals[3]*(s_occup+p_occup)*(s_occup+p_occup);
      *s_Z = Z_vals[0] +
        Z_vals[1]*s_occup +
          Z_vals[2]*p_occup;
    }
  } else if(atom->ns){
    E_vals = atom->muller_s_E;
    Z_vals = atom->muller_s_Z;
    *s_E = E_vals[0] +
      E_vals[1] * s_occup +
        E_vals[2] * s_occup * s_occup;
    *s_Z = Z_vals[0] +
      Z_vals[1]*s_occup +
        Z_vals[2]*s_occup*s_occup;
  }
}



/****************************************************************************
 *
 *                   Procedure update_muller_it_parms
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
 *   atomic Hii's and zetas using Muller's iteration technique
 *
 ****************************************************************************/
void update_muller_it_parms(details,cell,AO_occups,converged,num_orbs,
                      orbital_lookup_table)
  detail_type *details;
  cell_type *cell;
  real *AO_occups;
  int *converged;
  int num_orbs;
  int *orbital_lookup_table;
{
  static int num_its = 0;
  static char *atoms_done=0;
  static int num_atoms;

  real *E_vals,*Z_vals;
  real new_s_Hii,new_s_zeta,new_p_Hii,new_p_zeta,new_d_Hii,new_d_zeta;
  real temp_s_Hii,temp_s_zeta,temp_p_Hii,temp_p_zeta,temp_d_Hii,temp_d_zeta;
  real dH,dZ,max_dH=0.0,max_dZ=0.0;

  atom_type *atom;
  equiv_atom_type *equiv_list;
  int begin_atom,end_atom;
  int orb_tab;
  real num_s,num_p,num_d;
  real occups[9];

  int i,j,k;

  /* first get space for the atoms_done array, if we need it */
  if( !atoms_done || cell->num_atoms > num_atoms ){
    num_atoms = cell->num_atoms;
    atoms_done = (char *)malloc(num_atoms*sizeof(char));
    if( !atoms_done ) fatal("can't get memory for atoms_done array");
  }
  /* zero out the atoms_done array */
  bzero(atoms_done,num_atoms*sizeof(char));

  /*******

    first loop through the equivalent atoms,

    after we finish, tag each atom in
     the atoms_done array so that we don't end up doing an
     atom more than once.

  ********/
  equiv_list = cell->equiv_atoms;
  while(equiv_list){
    atom = &cell->atoms[equiv_list->equiv_atoms[0]];

    if( atom->chg_it_vary ){
      /* zero out the occups array, then do the average */
      bzero(occups,9*sizeof(real));
      for(i=0;i<equiv_list->num_equiv;i++){
        atom = &cell->atoms[equiv_list->equiv_atoms[i]];
        atoms_done[equiv_list->equiv_atoms[i]] = 1;
        find_atoms_orbs(num_orbs,cell->num_atoms,equiv_list->equiv_atoms[i],
                        orbital_lookup_table,&begin_atom,&end_atom);
        if( begin_atom >= 0){
          for(j=0;j<end_atom-begin_atom;j++){
            occups[j] +=
              AO_occups[j+begin_atom]/(real)equiv_list->num_equiv;
          }
        }
      }

      /****

        okay, we've got the average occupations, do the update

      ****/

      /* start by determing the total occupation of each type of orbital */
      atom = &cell->atoms[equiv_list->equiv_atoms[0]];
      num_s = 0.0;
      num_p = 0.0;
      num_d = 0.0;
      orb_tab = 0;
      if( atom->ns ){
        num_s = occups[orb_tab];
        orb_tab++;
      }
      if(atom->np){
        num_p = occups[orb_tab]+occups[orb_tab+1]+occups[orb_tab+2];
        orb_tab += 3;
      }
      if(atom->nd){
        num_d = occups[orb_tab]+occups[orb_tab+1]+occups[orb_tab+2]+
          occups[orb_tab+3]+occups[orb_tab+4];
      }

      /********

        check for negative AO occupations (these arise from
        counterintuitive orbital mixing and are not real)

      *********/
#ifdef ZERO_BOGUS_OCCUPS
      if( num_s < 0.0 ){
        fprintf(stderr,"Negative s occupation (%6.4lf) zeroed\n",num_s);
        num_s = 0.0;
      }
      if( num_p < 0.0 ){
        fprintf(stderr,"Negative p occupation (%6.4lf) zeroed\n",num_p);
        num_p = 0.0;
      }
      if( num_d < 0.0 ){
        fprintf(stderr,"Negative d occupation (%6.4lf) zeroed\n",num_d);
        num_d = 0.0;
      }
#else
      if( num_s < 0.0 ){
        fprintf(stderr,"Negative s occupation (%6.4lf) found\n",num_s);
      }
      if( num_p < 0.0 ){
        fprintf(stderr,"Negative p occupation (%6.4lf) found\n",num_p);
      }
      if( num_d < 0.0 ){
        fprintf(stderr,"Negative d occupation (%6.4lf) found\n",num_d);
      }
#endif


      /**********

        now use those total occupations to update the parameters

      **********/
      calc_muller_parms(atom,num_s,num_p,num_d,&temp_s_Hii,&temp_s_zeta,
                        &temp_p_Hii,&temp_p_zeta,&temp_d_Hii,&temp_d_zeta);
      if( atom->ns ){
        new_s_Hii = (1.0 - details->muller_mix) * atom->coul_s +
          details->muller_mix * temp_s_Hii;
        new_s_zeta = (1.0 - details->muller_mix) * atom->exp_s +
          details->muller_mix * temp_s_zeta;
        dH = fabs(atom->coul_s-new_s_Hii);
        dZ = fabs(atom->exp_s-new_s_zeta);
        if( dH > max_dH ) max_dH = dH;
        if( dZ > max_dZ ) max_dZ = dZ;
      }
      if( atom->np ){
        new_p_Hii = (1.0 - details->muller_mix) * atom->coul_p +
          details->muller_mix * temp_p_Hii;
        new_p_zeta = (1.0 - details->muller_mix) * atom->exp_p +
          details->muller_mix * temp_p_zeta;
        dH = fabs(atom->coul_p-new_p_Hii);
        dZ = fabs(atom->exp_p-new_p_zeta);
        if( dH > max_dH ) max_dH = dH;
        if( dZ > max_dZ ) max_dZ = dZ;
      }
      if(atom->nd){
        new_d_Hii = (1.0 - details->muller_mix) * atom->coul_d +
          details->muller_mix * temp_d_Hii;
        new_d_zeta = (1.0 - details->muller_mix) * atom->exp_d +
          details->muller_mix * temp_d_zeta;
        dH = fabs(atom->coul_d-new_d_Hii);
        dZ = fabs(atom->exp_d-new_d_zeta);
        if( dH > max_dH ) max_dH = dH;
        if( dZ > max_dZ ) max_dZ = dZ;
      }

      /* we've got the new parameters, update the equivalent atoms */
      for(i=0;i<equiv_list->num_equiv;i++){
        atom = &cell->atoms[equiv_list->equiv_atoms[i]];
        if( atom->nd ){
          atom->exp_d = new_d_zeta;
          atom->coul_d = new_d_Hii;
        }
        if( atom->np ){
          atom->exp_p = new_p_zeta;
          atom->coul_p = new_p_Hii;
        }
        if( atom->ns ){
          atom->exp_s = new_s_zeta;
          atom->coul_s = new_s_Hii;
        }
      }
    }
    equiv_list = equiv_list->next;
  }

  /******

    okay, that takes care of the equivalent atoms list, do any
    leftover atoms

    This is done exactly the same way as the update for the
     equivalent atoms.

  *******/

  for(i=0;i<cell->num_atoms;i++){
    atom = &cell->atoms[i];
    find_atoms_orbs(num_orbs,cell->num_atoms,i,
                    orbital_lookup_table,&begin_atom,&end_atom);

    if( !atoms_done[i] && atom->chg_it_vary && begin_atom >= 0){
      num_s = 0.0;
      num_p = 0.0;
      num_d = 0.0;
      orb_tab = begin_atom;
      if( atom->ns ){
        num_s = AO_occups[orb_tab];
        orb_tab++;
      }
      if(atom->np){
        num_p = AO_occups[orb_tab]+AO_occups[orb_tab+1]+AO_occups[orb_tab+2];
        orb_tab += 3;
      }
      if(atom->nd){
        num_d = AO_occups[orb_tab]+AO_occups[orb_tab+1]+AO_occups[orb_tab+2]+
          AO_occups[orb_tab+3]+AO_occups[orb_tab+4];
      }
#ifdef ZERO_BOGUS_OCCUPS
      if( num_s < 0.0 ){
        fprintf(stderr,"Negative s occupation (%6.4lf) zeroed\n",num_s);
        num_s = 0.0;
      }
      if( num_p < 0.0 ){
        fprintf(stderr,"Negative p occupation (%6.4lf) zeroed\n",num_p);
        num_p = 0.0;
      }
      if( num_d < 0.0 ){
        fprintf(stderr,"Negative d occupation (%6.4lf) zeroed\n",num_d);
        num_d = 0.0;
      }
#else
      if( num_s < 0.0 ){
        fprintf(stderr,"Negative s occupation (%6.4lf) found\n",num_s);
      }
      if( num_p < 0.0 ){
        fprintf(stderr,"Negative p occupation (%6.4lf) found\n",num_p);
      }
      if( num_d < 0.0 ){
        fprintf(stderr,"Negative d occupation (%6.4lf) found\n",num_d);
      }
#endif

      calc_muller_parms(atom,num_s,num_p,num_d,&temp_s_Hii,&temp_s_zeta,
                        &temp_p_Hii,&temp_p_zeta,&temp_d_Hii,&temp_d_zeta);
      if( atom->ns ){
        new_s_Hii = (1.0 - details->muller_mix) * atom->coul_s +
          details->muller_mix * temp_s_Hii;
        new_s_zeta = (1.0 - details->muller_mix) * atom->exp_s +
          details->muller_mix * temp_s_zeta;
        dH = fabs(atom->coul_s-new_s_Hii);
        dZ = fabs(atom->exp_s-new_s_zeta);
        if( dH > max_dH ) max_dH = dH;
        if( dZ > max_dZ ) max_dZ = dZ;
        atom->exp_s = new_s_zeta;
        atom->coul_s = new_s_Hii;
      }
      if( atom->np ){
        new_p_Hii = (1.0 - details->muller_mix) * atom->coul_p +
          details->muller_mix * temp_p_Hii;
        new_p_zeta = (1.0 - details->muller_mix) * atom->exp_p +
          details->muller_mix * temp_p_zeta;
        dH = fabs(atom->coul_p-new_p_Hii);
        dZ = fabs(atom->exp_p-new_p_zeta);
        if( dH > max_dH ) max_dH = dH;
        if( dZ > max_dZ ) max_dZ = dZ;
        atom->exp_p = new_p_zeta;
        atom->coul_p = new_p_Hii;
      }
      if(atom->nd){
        new_d_Hii = (1.0 - details->muller_mix) * atom->coul_d +
          details->muller_mix * temp_d_Hii;
        new_d_zeta = (1.0 - details->muller_mix) * atom->exp_d +
          details->muller_mix * temp_d_zeta;
        dH = fabs(atom->coul_d-new_d_Hii);
        dZ = fabs(atom->exp_d-new_d_zeta);
        if( dH > max_dH ) max_dH = dH;
        if( dZ > max_dZ ) max_dZ = dZ;
        atom->exp_d = new_d_zeta;
        atom->coul_d = new_d_Hii;
      }
    }
  }

  num_its++;
fprintf(stderr,"Muller it %d: dH = %lf dZ = %lf\n",num_its,max_dH,max_dZ);
write_atom_parms(details,cell->atoms,cell->num_atoms,1);

  /* now check convergence */
  if( max_dH <= details->muller_E_tol && max_dZ <= details->muller_Z_tol ){
fprintf(stderr,"Muller iteration converged after %d steps\n",num_its);
    *converged = 1;
  }else{
    *converged = 0;
  }

}


