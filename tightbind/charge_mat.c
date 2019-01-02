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
*   This is the charge matrix information
*
*  created:  greg landrum  January 1994
*
*****************************************************************************/
#include "bind.h"


/****************************************************************************
 *
 *                   Procedure eval_charge_matrix
 *
 * Arguments:  cell: pointer to cell type
 *         eigenset: eigenset_type
 *          overlap: hermetian_matrix_type
 *         num_orbs: int
 * orbital_lookup_table: pointer to int.
 *         chg_matrix: pointer to real
 *             accum: pointer to real
 *
 *
 * Returns: none
 *
 * Action:  Evaluates the charge_matrix at the current k point
 *
 *  'chg_matrix is an array which is used to build the charge matrix... it should
 *    be num_orbs*num_orbs in size.
 *  'accum is an array which should be at least (number of atoms) long... it is
 *    just used to make this more efficient.
 *
 *  NOTE:  In the interests of efficiency, this function make assumptions
 *    about the form in which Hermetian matrices are stored.  If the
 *    representation ever changes, then this will need to be rewritten...
 *    ahhhhh, the suffering one must do for execution speed. :-)
 *
 ****************************************************************************/
void eval_charge_matrix(cell,eigenset,overlap,num_orbs,
                          orbital_lookup_table,chg_matrix,accum)
  cell_type *cell;
  eigenset_type eigenset;
  hermetian_matrix_type overlap;
  int num_orbs,*orbital_lookup_table;
  real *chg_matrix,*accum;
{
  int num_atoms;
  int i,j,k,l;
  int itab,jtab,ktab;
  int start_orb,end_orb;
  int electrons_done,num_electrons;
  int top_of_degeneracy;
  real weight;
  real net_chg;
  real AO_chg,AO_chgI;
  real Sjk_R,Sjk_I,Cik_R,Cik_I;

  num_atoms = cell->num_atoms;
  num_electrons = cell->num_electrons;

  /* now, loop over crystal orbitals, then atomic orbitals */
  electrons_done = 0;

  for(i=0;i<num_orbs;i++){
    itab = i*num_orbs;
    for(j=0;j<num_orbs;j++){
      jtab = j*num_orbs;

      /* loop over the other AO's in _THIS_ MO */
      AO_chg = AO_chgI = 0;
      AO_chg = HERMETIAN_R(overlap,j,j) * EIGENVECT_R(eigenset,i,j);
      AO_chgI = HERMETIAN_I(overlap,j,j) * EIGENVECT_I(eigenset,i,j);

      for( k=j+1; k<num_orbs; k++){
        ktab = k*num_orbs;
        /****

          to save a ton of pointer math, set some temporary
          variables here.

        *****/
        Sjk_R = overlap.mat[jtab+k];
        Sjk_I = overlap.mat[ktab+j];
        Cik_R = eigenset.vectR[itab+k];
        Cik_I = eigenset.vectI[itab+k];

        AO_chg +=  Sjk_R * Cik_R + Sjk_I * Cik_I;
        AO_chgI += Sjk_R * Cik_I - Sjk_I * Cik_R;
      }

      for( k=0; k<j; k++){
        ktab = k*num_orbs;
        Sjk_R = overlap.mat[ktab+j];
        Sjk_I = overlap.mat[jtab+k];
        Cik_R = eigenset.vectR[itab+k];
        Cik_I = eigenset.vectI[itab+k];

        AO_chg +=  Sjk_R * Cik_R - Sjk_I * Cik_I;
        AO_chgI += Sjk_R * Cik_I + Sjk_I * Cik_R;
      }

      chg_matrix[itab+j] = EIGENVECT_R(eigenset,i,j)*AO_chg +
        EIGENVECT_I(eigenset,i,j)*AO_chgI;

    }

  }
}

/****************************************************************************
 *
 *                   Procedure reduced_charge_mat
 *
 * Arguments: num_atoms: int
 *         num_orbs: int
 * orbital_lookup_table: pointer to int.
 *         Chg_matrix: pointer to real
 *        RChg_matrix: pointer to real
 *
 *
 * Returns: none
 *
 * Action:  Generates the reduced charge matrix by summing
 *    up the individual atomic contributions in Chg_matrix.
 *
 *   the results are stored in 'RChg_matrix which should be at least
 *     num_atoms * num_orbs
 *
 ****************************************************************************/
void reduced_charge_mat(num_atoms,num_orbs,orbital_lookup_table,Chg_matrix,RChg_matrix)
  int num_atoms,num_orbs,*orbital_lookup_table;
  real *Chg_matrix;
  real *RChg_matrix;
{
  int i,j,k,l;
  int ktab;
  int i_orb_tab, j_orb_tab;
  int i_orb_end,j_orb_end;
  int num_elements;
  real temp;

  num_elements = 0;

  /**********
    loop over the atoms in the unit cell
  **********/
  for(i=0;i<num_atoms;i++){
    /* find this atom's orbitals */
    find_atoms_orbs(num_orbs,num_atoms,i,orbital_lookup_table,&i_orb_tab,&i_orb_end);
    if( i_orb_tab >= 0 ){
      /* first do the cross terms involving this atom */
      for(j=0;j<num_orbs;j++){
        temp = 0;
        for(k=i_orb_tab; k<i_orb_end; k++){
          temp += Chg_matrix[j*num_orbs+k];
        }
        /*********
          store the element....
          **********/
        RChg_matrix[i*num_orbs+j] = temp;
      }
    }
  }
}