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
*   This is the stuff required to find the inertia tensor and its principle
*    axes.
*
*  NOTE: in order to avoid having to write my own matrix diagonalization code
*   the routines in this file are dependant on the meschach matrix library
*   which is freely available from netlib.  If you do not understand what this
*   means or how to get things from netlib feel free to ask me.
*
*  created:  greg landrum  February 1994
*
*****************************************************************************/
#include "bind.h"
#include "symmetry.h"
#ifdef REAL_SYMM_ANALYSIS
#include <meschach/matrix.h>
#endif

/****************************************************************************
*
*                   Procedure find_princ_axes
*
* Arguments:  atoms:  pointer to atom_type
*              locs:  pointer to point_type
*        princ_axes:  a 3x3 matrix of reals
*           moments:  a 3-vector of reals
*         num_atoms:  integer
*
* Returns: none
*
* Action:  builds the inertia tensor and diagonalizes it to find the principle
*    axes, which are returned in 'princ_axes.
*
*  NOTE: in order for this to give physically meaningful results, it is important
*   to translate the atoms to center of mass coordinates first.
*
*****************************************************************************/
void find_princ_axes(atoms,locs,princ_axes,moments,num_atoms)
  atom_type *atoms;
  point_type *locs;
  real moments[3],princ_axes[3][3];
  int num_atoms;
{
  int i,j;
  atom_type *atom;
  real mass;

#ifdef REAL_SYMM_ANALYSIS
  MAT *MOI_tensor,*evects;
  VEC *evals;

  /********

    get space for the matrices, using meschach's routines

  ********/
  MOI_tensor = m_get(3,3);
  evects = m_get(3,3);
  evals = v_get(3);

  /********

    now build the moment of inertia tensor....
    we only have to build half, 'cause it's symmetrical.

  ********/
  for(i=0;i<num_atoms;i++){
    atom = &atoms[i];

    mass = atom->at_number;

    /* we don't want dummy atoms contributing */
    if( mass >= 0 ){

      MOI_tensor->me[0][0] += mass*(atom->loc.z*atom->loc.z +
                                    atom->loc.y*atom->loc.y);
      MOI_tensor->me[0][1] -= mass*(atom->loc.x*atom->loc.y);
      MOI_tensor->me[0][2] -= mass*(atom->loc.x*atom->loc.z);

      MOI_tensor->me[1][1] += mass*(atom->loc.z*atom->loc.z +
                                    atom->loc.x*atom->loc.x);
      MOI_tensor->me[1][2] -= mass*(atom->loc.y*atom->loc.z);

      MOI_tensor->me[2][2] += mass*(atom->loc.y*atom->loc.y +
                                    atom->loc.x*atom->loc.x);
    }
  }

  /* copy across the diagonal */
  MOI_tensor->me[1][0] = MOI_tensor->me[0][1];
  MOI_tensor->me[2][0] = MOI_tensor->me[0][2];
  MOI_tensor->me[2][1] = MOI_tensor->me[1][2];


  /*******

    diagonalize it and copy the elements into the princ_axes matrix
    and the moments array

  ********/
  symmeig(MOI_tensor,evects,evals);

  for(i=0;i<3;i++){
    moments[i] = evals->ve[i];
    for(j=0;j<3;j++){
      princ_axes[i][j] = evects->me[i][j];
    }
  }
#else
  bzero(princ_axes,9*sizeof(real));
  princ_axes[0][0] = 1.0;
  princ_axes[1][1] = 1.0;
  princ_axes[2][2] = 1.0;
#endif

  /* that's it, return now... */
}
