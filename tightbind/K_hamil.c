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
*     this file contains stuff for building the K space hamiltonian
*
*  created:  greg landrum  September 1993
*
*****************************************************************************/
#include "bind.h"


/****************************************************************************
 *
 *                   Procedure build_k_hamil_FAT
 *
 * Arguments:  cell: pointer to cell type
 *    hamilR,hamilK: hermetian_matrix_type
 *         overlapK: hermetian_matrix_type
 *         num_orbs: int
 * orbital_lookup_table: pointer to int.
 *
 * Returns: none
 *
 * Action:  This just forms the k space hamiltonian for this point:
 *
 *     Hij(K) = Hij(R)*Sij(K)
 *     Hii(K) = Hii(R)(7/3 + 1.75*Sii(K))  <- I don't know what the hell this
 *                                            is all about.
 *
 *
 ****************************************************************************/
void build_k_hamil_FAT(cell,hamilR,hamilK,overlapK,num_orbs)
  cell_type *cell;
  hermetian_matrix_type hamilR,hamilK;
  hermetian_matrix_type overlapK;
  int num_orbs;
{
  int i,j;
  int itab,jtab;

  for(i=0;i<num_orbs;i++){
    itab = i*num_orbs;


    /* do the diagonal elements */
    hamilK.mat[itab+i] = hamilR.mat[itab+i]*(1-THE_CONST+
                                     THE_CONST*overlapK.mat[itab+i]);

    /* now the off diagonals */
    for(j=0;j<i;j++){
      jtab = j*num_orbs;

      hamilK.mat[jtab+i] = hamilR.mat[itab+j]*overlapK.mat[jtab+i];
      hamilK.mat[itab+j] = hamilR.mat[jtab+i]*overlapK.mat[itab+j];
    }
  }


#ifdef PRINTMAT
fprintf(output_file,"\n\n^^^^^^H(K)^^^^^^\n");
printmat(hamilK.mat,num_orbs,num_orbs,output_file,1e-6,details->line_width);
#endif
}

/****************************************************************************
 *
 *                   Procedure build_k_hamil_THIN
 *
 * Arguments:  cell: pointer to cell type
 *    hamilR,hamilK: pointer to real type
 *         overlapK: pointer to real type
 *         num_orbs: int
 * orbital_lookup_table: pointer to int.
 *
 * Returns: none
 *
 * Action:  This just forms the k space hamiltonian for this point:
 *
 *  This is the same as the FAT version (above) but the off diagonal
 *   H(R) elements must be generated....
 *
 *  unless I'm very confused, the time spent here is gonna be TOTALLY
 *    insignificant in terms of the diagonalization and building S(K).
 *
 ****************************************************************************/
void build_k_hamil_THIN(cell,hamilR,hamilK,overlapK,num_orbs)
  cell_type *cell;
  hermetian_matrix_type hamilR,hamilK;
  hermetian_matrix_type overlapK;
  int num_orbs;
{
  int i,j;
  int itab,jtab;
  real temp,temp2;

  hamilK.mat[0] = hamilR.mat[0]*(1-THE_CONST+THE_CONST*overlapK.mat[0]);
  for(i=1;i<num_orbs;i++){
    itab = i*num_orbs;
    hamilK.mat[itab+i] = hamilR.mat[i]*(1-THE_CONST+THE_CONST*overlapK.mat[itab+i]);
    for(j=0;j<i;j++){
      jtab = j*num_orbs;
      /******
        The form of these is stolen straight from the new3 code
      *******/
      temp = hamilR.mat[i]+hamilR.mat[j];
      temp2 = (hamilR.mat[i]-hamilR.mat[j])/temp;
      temp2 = temp2*temp2;
      temp2 = .5*(THE_CONST + temp2 + temp2*temp2*(1-THE_CONST));
      temp2 *= temp;
      hamilK.mat[jtab+i] = temp2*overlapK.mat[jtab+i];
      hamilK.mat[itab+j] = temp2*overlapK.mat[itab+j];
    }
  }

#ifdef PRINTMAT
  fprintf(output_file,"\n\n^^^^^^H(K)^^^^^^\n");
printmat(hamilK.mat,num_orbs,num_orbs,output_file,1e-6,details->line_width);
#endif
}

