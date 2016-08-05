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
*     this file contains stuff for building the R space hamiltonian
*
*  created:  greg landrum  September 1993
*
*****************************************************************************/
#include "bind.h"

/***
  Edit History:

  March '98: WG
  - f orbitals added

***/
/****************************************************************************
 *
 *                   Procedure R_space_Hamiltonian
 *
 * Arguments:  cell: pointer to cell type
 *          details: pointer to detail type
 *          overlap: hermetian_matrix_type
 *            hamil: hermetian_matrix_type
 *         num_orbs: int
 * orbital_lookup_table: pointer to int.
 *
 * Returns: none
 *
 * Action:  Builds the R space hamiltonian... At the moment this is done by just
 *     filling the diagonal elements.  i.e. the hamiltonian is just an array
 *     of diagonal elements... the rest of the work can be done later when the
 *     K space hamiltonian (the important one) is built.
 *
 ****************************************************************************/
void R_space_Hamiltonian(cell,details,overlap,hamil,num_orbs,orbital_lookup_table)
  cell_type *cell;
  detail_type *details;
  hermetian_matrix_type overlap,hamil;
  int num_orbs;
  int *orbital_lookup_table;
{
  int i,j;
  int orb_tab1;
  atom_type *atom_ptr1;

  /* put in the diagonal elements */
  for(i=0;i<cell->num_atoms;i++){
    atom_ptr1 = &(cell->atoms[i]);
    orb_tab1 = orbital_lookup_table[i];

    /* trap dummy atoms */
    if( orb_tab1 >= 0 ){
      if( atom_ptr1->ns != 0 ) hamil.mat[orb_tab1+BEGIN_S]
        = atom_ptr1->coul_s;
      if( atom_ptr1->np != 0 )
        for(j=BEGIN_P;j<=END_P;j++){
          hamil.mat[orb_tab1+j] = atom_ptr1->coul_p;
        }
      if( atom_ptr1->nd != 0 )
        for(j=BEGIN_D;j<=END_D;j++){
          hamil.mat[orb_tab1+j] = atom_ptr1->coul_d;
        }
      if(atom_ptr1->nf != 0)
        for(j=BEGIN_F;j<END_F;j++){
          hamil.mat[orb_tab1+j] = atom_ptr1->coul_f;
        }
    }
  }
}

/****************************************************************************
 *
 *                   Procedure full_R_space_Hamiltonian
 *
 * Arguments:  cell: pointer to cell type
 *          details: pointer to detail type
 *          overlap: hermetian_matrix_type
 *            hamil: hermetian_matrix_type
 *         num_orbs: int
 * orbital_lookup_table: pointer to int.
 *      mult_by_overlap: char
 *
 * Returns: none
 *
 * Action:  Builds the full R space hamiltonian...
 *
 *    This includes all the elements.
 *
 ****************************************************************************/
void full_R_space_Hamiltonian(cell_type *cell,detail_type *details,
                              hermetian_matrix_type overlap,hermetian_matrix_type hamil,
                              int num_orbs,
                              int *orbital_lookup_table, char mult_by_overlap)
{
  int i,j,k;
  int orb_tab1,orb_tab2;
  atom_type *atom_ptr1,*atom_ptr2;
  real temp,temp2;
  static real *diagonal_elements=0;
#if 0
real *rham;

rham = (real *)calloc(num_orbs*num_orbs,sizeof(real));
#endif

  /* get space for the diagonal elements if it hasn't been done already */
  if( !diagonal_elements ){
    diagonal_elements = (real *)calloc(num_orbs,sizeof(real));
    if( !diagonal_elements )fatal("Can't get memory to build hamiltonian.");
  }

  /******
    put in the diagonal elements. These are just the coulomb
    integrals, which are given in the input file.
  *******/
  for(i=0;i<cell->num_atoms;i++){
    atom_ptr1 = &(cell->atoms[i]);
    orb_tab1 = orbital_lookup_table[i];

    /* trap dummy atoms */
    if( orb_tab1 >= 0 ){
      if( atom_ptr1->ns != 0 ){
        diagonal_elements[orb_tab1+BEGIN_S] =
          hamil.mat[(orb_tab1+BEGIN_S)*num_orbs+orb_tab1+BEGIN_S] =
            atom_ptr1->coul_s;
      }
      if( atom_ptr1->np != 0 ){
        for(j=BEGIN_P;j<=END_P;j++){
          diagonal_elements[orb_tab1+j] =
            hamil.mat[(orb_tab1+j)*num_orbs+orb_tab1+j] =
              atom_ptr1->coul_p;
        }
      }
      if( atom_ptr1->nd != 0 ){
        for(j=BEGIN_D;j<=END_D;j++){
          diagonal_elements[orb_tab1+j] =
            hamil.mat[(orb_tab1+j)*num_orbs+orb_tab1+j] =
              atom_ptr1->coul_d;
        }
      }

      if( atom_ptr1->nf != 0 ){
        for(j=BEGIN_F;j<=END_F;j++){
          diagonal_elements[orb_tab1+j] =
            hamil.mat[(orb_tab1+j)*num_orbs+orb_tab1+j] =
              atom_ptr1->coul_f;
        }
      }
    }
  }
  /* now do the off diagonals */
  for(i=1;i<num_orbs;i++){
    for(j=0;j<i;j++){
      /******
        The form of these is stolen straight from the new3 code
      *******/
      if( details->weighted_Hij ){
        /* use the weighted Hij formula */
        temp = diagonal_elements[i]+diagonal_elements[j];
        temp2 = (diagonal_elements[i]-diagonal_elements[j])/temp;
        temp2 = temp2*temp2;
        temp2 = .5*(details->the_const + temp2 + temp2*temp2*(1-details->the_const));
        hamil.mat[j*num_orbs+i] = temp2*temp;
      }
      else{
        /* don't use the weighted Hij formula */
        hamil.mat[j*num_orbs+i] = .5*details->the_const*(diagonal_elements[i] +
                                                         diagonal_elements[j]);
      }
      hamil.mat[i*num_orbs+j] = hamil.mat[j*num_orbs+i];
      if(details->Execution_Mode == MOLECULAR || mult_by_overlap){
        hamil.mat[j*num_orbs+i] *= overlap.mat[j*num_orbs+i];
        hamil.mat[i*num_orbs+j] *= overlap.mat[i*num_orbs+j];
      }

    }
  }

#if 0
  fprintf(output_file,"\n\n----------------------------- OVERLAP:\n");
printmat(overlap.mat,num_orbs,num_orbs,output_file,1e-6,details->line_width);

  fprintf(output_file,"\n\n----------------------------------Hamiltonian:\n");
printmat(rham,num_orbs,num_orbs,output_file,1e-6,details->line_width);
#endif
}


