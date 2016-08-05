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
*   This is the stuff for doing Mulliken analysis
*
*  created:  greg landrum  January 1994
*
*****************************************************************************/
#include "bind.h"

/****************************************************************************
 *
 *                   Procedure calc_occupations
 *
 * Arguments:        details: pointer to detail_type
 *             num_electrons: real
 *                  num_orbs: int
 *               occupations: pointer to type real
 *                  eigenset: eigenset_type
 *
 * Returns: none
 *
 * Action:  calculates the number of electrons in each orbital... Takes
 *    degeneracies into account.
 *
 *       If details is nonzero and details->num_orbital_occups is specified then
 *         the information in
 *         details->orbital_occups will be used to CHANGE the occupations of the
 *         specified orbitals.
 *
 ****************************************************************************/
void calc_occupations(details,num_electrons,num_orbs,occupations,eigenset)
  detail_type *details;
  real num_electrons;
  int num_orbs;
  real *occupations;
  eigenset_type eigenset;
{
  int i,begin_degen,end_degen,num_degen_levels,last_occup;
  real num_degen_electrons;
  real electrons_left,electrons_per_level;

  electrons_left = num_electrons;

  /* just loop until we run out of electrons */
  for(i=0;i<num_orbs && electrons_left > 0.0;i++){
    if( electrons_left >= 2.0){
      occupations[i] = 2.0;
      electrons_left -= 2.0;
    }
    else if( electrons_left > 0.0 ){
      occupations[i] = electrons_left;
      electrons_left = 0.0;
    }
  }

  /* now check for degeneracies */
  last_occup = i-1;
  end_degen = i;
  while( fabs(EIGENVAL(eigenset,last_occup) - EIGENVAL(eigenset,end_degen)) < DEGEN_TOL
        && end_degen < num_orbs ){
    end_degen++;
  }
  begin_degen = last_occup - 1;
  while( fabs(EIGENVAL(eigenset,last_occup) - EIGENVAL(eigenset,begin_degen)) < DEGEN_TOL
        && begin_degen > -1){
    begin_degen--;
  }
  /* we went one step too far, so increment begin_degen */
  begin_degen++;

  num_degen_levels = end_degen-begin_degen;
  if( num_degen_levels > 1){
    fprintf(output_file,"; >>>>> The HOMO was found to be %d-fold degenerate.\n",num_degen_levels);

    /********
      now adjust the occupation numbers for the degenerate orbitals
        do this by finding the number of electrons in degenerate levels and
        distributing them equally amongst those levels.
    ******/
    num_degen_electrons = 0;
    for(i=begin_degen;i<end_degen;i++){
      num_degen_electrons += occupations[i];
    }
    electrons_per_level = num_degen_electrons / (real)num_degen_levels;
    for(i=begin_degen;i<end_degen;i++){
      occupations[i] = electrons_per_level;
    }
  }

  /* now check to see if we need adjust any occupancies */
  if( details && details->num_orbital_occups > 0 ){
    for(i=0;i<details->num_orbital_occups;i++){
      occupations[details->orbital_occups[i].orb] = details->orbital_occups[i].occup;
    }
  }
}


/****************************************************************************
 *
 *                   Procedure reduced_mulliken
 *
 * Arguments: num_atoms: int
 *         num_orbs: int
 * orbital_lookup_table: pointer to int.
 *         OP_matrix: pointer to real
 *        ROP_matrix: pointer to real
 *
 *
 * Returns: none
 *
 * Action:  Generates the reduced overlap population matrix by summing
 *    up the individual atomic contributions in OP_matrix.
 *
 *   the results are stored in 'ROP_matrix which should be at least
 *     num_atoms * num_atoms
 *
 ****************************************************************************/
void reduced_mulliken(num_atoms,num_orbs,orbital_lookup_table,OP_matrix,ROP_matrix)
  int num_atoms,num_orbs,*orbital_lookup_table;
  real *OP_matrix;
  real *ROP_matrix;
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
      for(j=0;j<i;j++){
        find_atoms_orbs(num_orbs,num_atoms,j,orbital_lookup_table,&j_orb_tab,
                        &j_orb_end);
        if( j_orb_tab >= 0 ){
          temp = 0;
          for(k=i_orb_tab; k<i_orb_end; k++){
            ktab = k*num_orbs;
            for(l=j_orb_tab; l<j_orb_end; l++){
              temp += OP_matrix[ktab+l];
            }
          }

          /*********
            store the element....
            treat the matrix as a symmetrical matrix (see notes.outl for the format)
            **********/
          ROP_matrix[num_elements++] = temp;
        }
      }

      /* now add up the diagonal elements which arise from this atom */
      temp = 0;
      for(j=i_orb_tab; j<i_orb_end; j++){
        temp += OP_matrix[j*num_orbs+j];
      }

      ROP_matrix[num_elements++] = temp;
    }
  }
}


/****************************************************************************
 *
 *                   Procedure FMO_reduced_mulliken
 *
 * Arguments: details: pointer to detail_type
 *   num_atoms,num_orbs: integers
 *         OP_matrix: pointer to real
 *        ROP_matrix: pointer to real
 *
 *
 * Returns: none
 *
 * Action:  Generates the reduced overlap population matrix by summing
 *    up the individual fragment contributions in OP_matrix.
 *
 *   the results are stored in 'ROP_matrix which should be at least
 *     details->num_FMO_frags * details->num_FMO_frags
 *
 ****************************************************************************/
void FMO_reduced_mulliken(details,num_atoms,num_orbs,OP_matrix,ROP_matrix)
  detail_type *details;
  int num_atoms,num_orbs;
  real *OP_matrix;
  real *ROP_matrix;
{
  FMO_frag_type *FMO_frag1,*FMO_frag2;
  int frag1,frag2;
  int orbs_so_far1,orbs_so_far2;
  int i,j,k,l;
  int ktab;
  int i_orb_tab, j_orb_tab;
  int i_orb_end,j_orb_end;
  int num_elements;
  real temp;

  num_elements = 0;

  /**********

    loop over all of the fragments

  **********/
  orbs_so_far1 = 0;
  for(frag1=0;frag1<details->num_FMO_frags;frag1++){
    FMO_frag1 = &(details->FMO_frags[frag1]);

    /* do the cross terms involving this fragment */
    orbs_so_far2 = 0;
    for(frag2=0;frag2<frag1;frag2++){
      FMO_frag2 = &(details->FMO_frags[frag2]);

      temp = 0;

#if 0
      /* loop over the orbitals on each fragment */
      for(atom1 = 0; atom1<FMO_frag1->num_atoms; atom1++){
        find_atoms_orbs(num_orbs,num_atoms,FMO_frag1->atoms_in_frag[atom1],
                        FMO_frag1->orbital_lookup_table,
                        &begin_atom1,&end_atom1);

        if( begin_atom1 >= 0 ){
          for(atom2 = 0; atom2<FMO_frag2->num_atoms; atom2++){
            find_atoms_orbs(num_orbs,num_atoms,FMO_frag2->atoms_in_frag[atom2],
                            FMO_frag2->orbital_lookup_table,
                            &begin_atom2,&end_atom2);

            if( begin_atom2 >= 0 ){

              for(k=begin_atom1; k<end_atom1; k++){
                ktab = k*num_orbs;
                for(l=begin_atom2; l<end_atom2; l++){
                  temp += OP_matrix[ktab+l];
                }
              }
            }
          }
        }
      }
#endif

      /* loop over the orbitals on each fragment */
      for(k=0; k<FMO_frag1->num_orbs; k++){
        ktab = (k+orbs_so_far1)*num_orbs;
        for(l=0; l<FMO_frag2->num_orbs; l++){
          temp += OP_matrix[ktab+orbs_so_far2+l];
        }
      }

      /*********
        store the element....
        treat the matrix as a symmetric matrix
        (see notes.outl for the format)
      **********/
      ROP_matrix[num_elements++] = temp;
      orbs_so_far2 += FMO_frag2->num_orbs;
    }
#if 0
    /* now add up the diagonal elements which arise from this fragment */
    temp = 0;
    for(atom1 = 0; atom1<FMO_frag1->num_atoms; atom1++){
      find_atoms_orbs(num_orbs,num_atoms,FMO_frag1->atoms_in_frag[atom1],
                      orbital_lookup_table,&begin_atom1,&end_atom1);
      if( begin_atom1 >= 0 ){
        for(j=begin_atom1; j<end_atom1; j++){
          temp += OP_matrix[j*num_orbs+j];
        }
      }
    }
#endif
    /* now add up the diagonal elements which arise from this fragment */
    temp = 0;
    for(j=0; j<FMO_frag1->num_orbs; j++,orbs_so_far1++){
      temp += OP_matrix[orbs_so_far1*num_orbs+orbs_so_far1];
    }

    ROP_matrix[num_elements++] = temp;
  }
}




/****************************************************************************
 *
 *                   Procedure eval_mulliken
 *
 * Arguments: details: pointer to detail_type
               cell: pointer to cell type
 *         eigenset: eigenset_type
 *          overlap: hermetian_matrix_type
 *         num_orbs: int
 *      occupations: pointer to type real.
 * orbital_lookup_table: pointer to int.
 *         OP_matrix: pointer to real
 *             accum: pointer to real
 *         net_chgs: pointer to real
 *
 *
 * Returns: none
 *
 * Action:  This is the driver function for doing Mulliken analysis.
 *
 *  'OP_matrix is an array which is used to build the overlap population matrix...
 *    it should be num_orbs*num_orbs in size.
 *  'accum is an array which should be at least num_orbs long... it is
 *    just used to make this more efficient.
 *  'net_chgs should be at least cell->num_atoms long, it is used to store and
 *      return the net atomic charges
 *
 *  'occupations is an array which should be at least num_orbs in length. It should
 *     be filled with the occupation numbers of the various orbitals.
 *
 ****************************************************************************/
void eval_mulliken(cell,eigenset,overlap,num_orbs,
                   occupations,orbital_lookup_table,OP_matrix,net_chgs,accum)
  cell_type *cell;
  eigenset_type eigenset;
  hermetian_matrix_type overlap;
  int num_orbs,*orbital_lookup_table;
  real *occupations,*OP_matrix,*accum,*net_chgs;
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
  bzero((char *)OP_matrix,num_orbs*num_orbs*sizeof(real));

  /********

    this is using the following formula:

                   MO's
                   ---\
                   \
    OP[i][j] = 2 * /     occupations[k] * c[k][i] * c[k][j] * S[i][j]
                   ---/
                     k


   for off diagonal elements.... diagonal elements are obtained by using the above
     and dividing by 2
  *********/

  /* i indexes atomic orbitals */
  for(i=0;i<num_orbs;i++){
    itab = i*num_orbs;

    /* j indexes atomic orbitals */
    for(j=0;j<num_orbs;j++){
      jtab = j*num_orbs;

      /* k indexes the molecular orbitals */
      OP_accum = 0.0;
      OP_accumI = 0.0;
      for(k=0;k<num_orbs;k++){
        ktab = k;

        OP_accum += occupations[k]*EIGENVECT_R(eigenset,k,i)*EIGENVECT_R(eigenset,k,j);
        OP_accum += occupations[k]*EIGENVECT_I(eigenset,k,i)*EIGENVECT_I(eigenset,k,j);
        OP_accumI -= occupations[k]*EIGENVECT_R(eigenset,k,i)*EIGENVECT_I(eigenset,k,j);
        OP_accumI += occupations[k]*EIGENVECT_I(eigenset,k,i)*EIGENVECT_R(eigenset,k,j);
      }

      if( i != j ){
        OP_matrix[itab+j] = 2.0 * OP_accum * HERMETIAN_R(overlap,i,j);

      }
      else{
        OP_matrix[itab+j] =  OP_accum * HERMETIAN_R(overlap,i,j);
      }

      /* this is used later to find net charges */
      accum[i] +=   OP_accum * HERMETIAN_R(overlap,i,j);


    }
  }

  /*******
    generate the net charges by summing up the diagonal elements of the overlap
    population matrix for each atom.
  ********/
  for(i=0;i<cell->num_atoms;i++){
    find_atoms_orbs(num_orbs,cell->num_atoms,i,orbital_lookup_table,
                    &begin_of_atom,&end_of_atom);
    if( begin_of_atom >= 0 ){
      net_chg = 0.0;

      for(j = begin_of_atom;j<end_of_atom;j++){
        net_chg -= accum[j];
      }

      /* subtract off the number of valence electrons */
      net_chg += cell->atoms[i].num_valence;

      net_chgs[i] = net_chg;
    }
  }
}

