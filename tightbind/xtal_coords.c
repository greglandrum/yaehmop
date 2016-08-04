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
*     this file has the routines for dealing with crystallographic
*   coordinates
*
*   This code is based on some FORTRAN code provided by Edgar Muller.
*
*  created:  greg landrum  June 1995
*
*****************************************************************************/
#include "bind.h"


/****************************************************************************
*
*                   Procedure eval_xtal_coord_locs
*
* Arguments: cell: pointer to cell_type
*       printing: char
*
* Returns: none
*
* Action: converts the atoms stored in 'cell from crystallographic
*   coordinates to cartesian coordinates.
*
*   'printing is used to control whether or not all the coordinates
*    are printed out.
*
*   NOTE:  If this function is called more than once, it is
*     absolutely essential that the atomic locations in 'cell->atoms
*     be updated from 'cell->raw_atoms between calls.  Otherwise the
*     cartesian positions will be converted, which will totally
*     screw things up.  It is theoretically possible for the
*     function to catch this, but it involves way too much work
*     for me.  Just be careful if you use this yourself.
*
*****************************************************************************/
void eval_xtal_coord_locs(cell_type *cell,char printing)
{
  int i;
  real cos_alpha,cos_beta,cos_gamma,sin_gamma;
  real weird_term,cell_volume;
  real tform_mat[3][3];

  /* write the crystallographic information if we need to */
  if( printing ){
    fprintf(output_file,"# CRYSTAL SPECIFICATION:\n");
    fprintf(output_file,"#Cell Constants: \n");
    for(i=0;i<cell->dim;i++)
      fprintf(output_file,"%6.4lf ",cell->xtal_defn.axis_lengths[i]);
    fprintf(output_file,"\n");

    fprintf(output_file,"# alpha: % 6.4lf\n",cell->xtal_defn.angles[0]);
    fprintf(output_file,"# beta: % 6.4lf\n",cell->xtal_defn.angles[1]);
    fprintf(output_file,"# gamma: % 6.4lf\n",cell->xtal_defn.angles[2]);
  }


  /* set up the transformation matrix */
  cos_alpha = cos(cell->xtal_defn.angles[0]*PI/180.0);
  cos_beta = cos(cell->xtal_defn.angles[1]*PI/180.0);
  cos_gamma = cos(cell->xtal_defn.angles[2]*PI/180.0);
  sin_gamma = sin(cell->xtal_defn.angles[2]*PI/180.0);
  weird_term = sqrt(1-cos_alpha*cos_alpha-cos_beta*cos_beta-
                    cos_gamma*cos_gamma + 2.0*cos_alpha*cos_beta*cos_gamma);
  cell_volume = weird_term*cell->xtal_defn.axis_lengths[0]*
    cell->xtal_defn.axis_lengths[1]*cell->xtal_defn.axis_lengths[2];

  if( printing ){
    fprintf(output_file,"#Cell Volume: %6.4lf cubic Angstroms\n",cell_volume);
  }


  /* first do the diagonal terms */
  tform_mat[0][0] = cell->xtal_defn.axis_lengths[0];
  tform_mat[1][1] = cell->xtal_defn.axis_lengths[1]*sin_gamma;
  tform_mat[2][2] = cell->xtal_defn.axis_lengths[2]*weird_term/sin_gamma;

  /* now the off diagonals */
  tform_mat[0][1] = cell->xtal_defn.axis_lengths[1]*cos_gamma;
  tform_mat[0][2] = cell->xtal_defn.axis_lengths[2]*cos_beta;
  tform_mat[1][2] = cell->xtal_defn.axis_lengths[2]*(cos_alpha-cos_beta*cos_gamma) /
    sin_gamma;
  tform_mat[1][0] = 0.0;
  tform_mat[2][0] = 0.0;
  tform_mat[2][1] = 0.0;

  /* now transform all the atomic locations */
  transform_atoms(cell->atoms,tform_mat,cell->num_raw_atoms);

  /* transform the lattice vectors (if necessary) */
  if( cell->dim > 0 )
    transform_atoms(&cell->atoms[cell->num_atoms],tform_mat,cell->dim);

  /* print out the new locations */
  if( printing ){

    fprintf(output_file,"# Positions of atoms from crystal coordinates\n");
    for(i=0;i<cell->num_raw_atoms;i++){
      fprintf(output_file,"%4d %4s %8.4lf %8.4lf %8.4lf\n",
              i+1,cell->atoms[i].symb,cell->atoms[i].loc.x,
              cell->atoms[i].loc.y,cell->atoms[i].loc.z);
    }
  }
}


