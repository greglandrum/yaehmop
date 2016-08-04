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
*   This is the stuff for doing modified Mulliken analysis
*
*  created:  greg landrum  April 1996
*
*****************************************************************************/
#include "bind.h"





/****************************************************************************
 *
 *                   Procedure modified_mulliken
 *
 * Arguments: details: pointer to detail_type
 *             cell: pointer to cell type
 *         eigenset: eigenset_type
 *          overlap: hermetian_matrix_type
 *         num_orbs: int
 *      occupations: pointer to type real.
 * orbital_lookup_table: pointer to int.
 *         OP_matrix: pointer to real
 *     mod_OP_matrix: pointer to real
 *         net_chgs: pointer to real
 *             accum: pointer to real
 *
 *
 * Returns: none
 *
 * Action:  This is the driver function for doing modified Mulliken analysis.
 *
 *  'OP_matrix is an array which contains the overlap population matrix...
 *    it should be num_orbs*num_orbs in size.
 *  'mod_OP_matrix is an array which is used to build the modified
 *    overlap population matrix...
 *    it should be num_orbs*num_orbs in size.
 *
 *  'accum is an array which should be at least num_orbs long... it is
 *    just used to make this more efficient.
 *  'net_chgs should be at least cell->num_atoms long, it is used to store and
 *      return the net atomic charges
 *
 *  'occupations is an array which should be at least num_orbs in length. It should
 *     be filled with the occupation numbers of the various orbitals.
 *
 ****************************************************************************/
void modified_mulliken(cell,eigenset,overlap,num_orbs,
                       occupations,orbital_lookup_table,OP_matrix,
                       mod_OP_matrix,net_chgs,accum)
  cell_type *cell;
  eigenset_type eigenset;
  hermetian_matrix_type overlap;
  int num_orbs,*orbital_lookup_table;
  real *occupations,*OP_matrix,*mod_OP_matrix,*accum,*net_chgs;
{
  int num_atoms;
  int i,j,k;
  int itab,jtab,ktab;
  int begin_of_atom,end_of_atom;
  real OP_accum,net_chg,tot_chg;
  real OP_accumI;

  num_atoms = cell->num_atoms;

  /* zero out arrays which will be used */
  bzero((char *)accum,num_orbs*sizeof(real));


  /********

    this is using the following formula:



    i != j
                          OP[i][j]
    mod_OP[i][j] =  -------------------
                    OP[i][i] + OP[j][j]

    i == j

    mod_OP[i][i] = OP[i][i]

  *********/

  /* i indexes atomic orbitals */
  for(i=0;i<num_orbs;i++){
    itab = i*num_orbs;

    /* j indexes atomic orbitals */
    for(j=0;j<num_orbs;j++){
      jtab = j*num_orbs;

      if( i != j ){
        mod_OP_matrix[itab+j] = OP_matrix[itab+j]/
          (OP_matrix[itab+i] + OP_matrix[jtab+j]);
      }
      else{
        mod_OP_matrix[itab+i] =  OP_matrix[itab+i];
      }
    }
  }

  /*******

    generate the net charges

  ********/
  for(i=0;i<cell->num_atoms;i++){
    find_atoms_orbs(num_orbs,cell->num_atoms,i,orbital_lookup_table,
                    &begin_of_atom,&end_of_atom);
    if( begin_of_atom >= 0 ){
      net_chg = 0.0;

      for(j = begin_of_atom;j<end_of_atom;j++){
        jtab = j*num_orbs;
        for(k=0;k<num_orbs;k++){
          if( k==j ) net_chg -= mod_OP_matrix[jtab+k];
          else net_chg -= mod_OP_matrix[jtab+j]*mod_OP_matrix[jtab+k];
        }
      }
      /* subtract off the number of valence electrons */
      net_chg += cell->atoms[i].num_valence;

      net_chgs[i] = net_chg;
    }
  }
}
