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
*     this file has the routines for dealing with Z matrices
*
*  created:  greg landrum  April 1994
*
*****************************************************************************/
#include "bind.h"

#define SMALL_X 1e-05


/* helper function to find a particular numbered atom in an array of atoms */
int find_atom(atoms,num_atoms,which)
  atom_type *atoms;
  int num_atoms,which;
{
  static char err_string[120];
  int i;

  for(i=0;i<num_atoms;i++){
    if( atoms[i].which_atom == which ) return(i);
  }

  /* we didn't find it... something's wrong, error out */
  sprintf(err_string,"Can't find atom %d (of %d) in find_atom.\n",
          which,num_atoms);
  FATAL_BUG(err_string);
}




/****************************************************************************
*
*                   Procedure eval_Zmat_locs
*
* Arguments:   atoms: pointer to atom_type
*          num_atoms: int
*           cell_dim: int
*           printing: char
*
* Returns: none
*
* Action: This converts the Z matrix coordinates in the atom list
*   into cartesian coordinates.
*
*   The procedure for doing this is adopted from the source code for GAMESS
*     (in mbldr.f), which was snarfed from one of the early version of gaussian.
*
*
*   The convention used is:
*     the first atom goes at the origin.
*     the second atom is along the positive Z direction
*     the third atom is in the XZ plane
*
*  NOTE: this code only deals with bond-length,angle,dihedral angle specifications.
*   it is not set up to deal with the second bond angle type Z matrix.
*
*  NOTE/APOLOGY:  This code is neither straightforward nor easy to read.  That's
*   because I haven't really bothered to figure it out to any great degree.
*   It should work.
*
*
*  Further explanation:  because of the fact that the user may have
*   specified the atom numbers differently from where they appeared
*   in the output file (to, for example, put all their dummy atoms
*   at the end of the array of atoms), it is necessary to use the
*   which_atom field of each atom structure in order to determine in
*   what order the atoms were really specified.
*
*****************************************************************************/
void eval_Zmat_locs(atom_type *atoms,int num_atoms, int cell_dim,char printing)
{
  int i,j;
  int which_atom,which_atom2;
  atom_type *atom;
  point_type *loc;
  Z_mat_type *Zmat_loc;
  point_type V_ref1_ref2,V_ref2_ref3,N_V_ref1_ref2,N_V_ref2_ref3;
  point_type V1_X_V2,N_V1_X_V2,V_last;
  point_type V_ref3_newpos;
  real bond_length,alpha,beta;
  real dot_p,denom;


  /*******
    The first atom is at the origin.
  *******/
  which_atom = find_atom(atoms,num_atoms+cell_dim,0);
  atom = &(atoms[which_atom]);
  loc = &(atom->loc);
  loc->x = loc->y = loc->z = 0.0;

  if( num_atoms+cell_dim > 1 ){
    /********
      the second atom is on the z axis
      ********/
    which_atom = find_atom(atoms,num_atoms+cell_dim,1);
    atom = &(atoms[which_atom]);
    loc = &(atom->loc);
    loc->x = loc->y = 0.0;
    loc->z = atom->Zmat_loc.bond_length;

  }
  if( num_atoms+cell_dim > 2){
    /********
      the third atom is in the xz plane
      ********/
    which_atom = find_atom(atoms,num_atoms+cell_dim,2);
    atom = &(atoms[which_atom]);
    loc = &(atom->loc);
    bond_length = atom->Zmat_loc.bond_length;
    loc->y = 0.0;

    /* convert the bond angle to degrees */
    alpha = PI*atom->Zmat_loc.alpha/180.0;

    /* now figure out the location of the atom */
    loc->x = bond_length * sin(alpha);

    which_atom = find_atom(atoms,num_atoms+cell_dim,atom->Zmat_loc.ref1);
    if( which_atom == 0 ){
      loc->z = bond_length*cos(alpha);
    }
    else if( which_atom == 1){
      which_atom = find_atom(atoms,num_atoms+cell_dim,1);
      loc->z = atoms[which_atom].loc.z - bond_length*cos(alpha);
    }
    else{
      fatal("Invalid angle reference in Z matrix for atom 3.");
    }
  }


  /* check for atoms lying in a line */
  for(i=3; i<num_atoms+cell_dim && loc->x <= SMALL_X; i++){
    which_atom = find_atom(atoms,num_atoms+cell_dim,i);
    atom = &(atoms[which_atom]);
    loc = &(atom->loc);
    Zmat_loc = &(atom->Zmat_loc);
    bond_length = atom->Zmat_loc.bond_length;
    alpha = PI*atom->Zmat_loc.alpha/180.0;
    beta = PI*atom->Zmat_loc.beta/180.0;

    loc->x = bond_length*sin(alpha);
    loc->y = 0.0;
    loc->z = atoms[which_atom].loc.z -
      bond_length*cos(alpha)*(atoms[atom->Zmat_loc.ref1].loc.z >
                              atoms[atom->Zmat_loc.ref2].loc.z ?
                              1 : -1);
  }

  /* deal with any atoms which are left over */
  for( j=i; j<num_atoms+cell_dim; j++ ){
    which_atom = find_atom(atoms,num_atoms+cell_dim,j);
    atom = &(atoms[which_atom]);
    loc = &(atom->loc);
    Zmat_loc = &(atom->Zmat_loc);
    bond_length = Zmat_loc->bond_length;
    alpha = PI*Zmat_loc->alpha/180.0;
    beta = PI*Zmat_loc->beta/180.0;

    /* determine the vector between ref1 and ref2 */
    vector_diff(&(atoms[Zmat_loc->ref1].loc),&(atoms[Zmat_loc->ref2].loc),
              &V_ref1_ref2);
    /* normalize that vector */
    normalize_vector(&V_ref1_ref2,&N_V_ref1_ref2);

    /* determine the vector between ref2 and ref3 */
    vector_diff(&(atoms[Zmat_loc->ref2].loc),&(atoms[Zmat_loc->ref3].loc),
              &V_ref2_ref3);
    /* normalize that vector */
    normalize_vector(&V_ref2_ref3,&N_V_ref2_ref3);

    /* figure out the cross product between the normalized vectors */
    cross_prod(&N_V_ref2_ref3,&N_V_ref1_ref2,&V1_X_V2);

    /* "normalize" it */
    dot_p = dot_prod(&N_V_ref1_ref2,&N_V_ref2_ref3);
    denom = 1.0-dot_p*dot_p;
    denom = sqrt(denom);

    N_V1_X_V2.x = V1_X_V2.x / denom;
    N_V1_X_V2.y = V1_X_V2.y / denom;
    N_V1_X_V2.z = V1_X_V2.z / denom;

    /* take the cross prod of that and assign it a dumb name */
    cross_prod(&(N_V1_X_V2),&(N_V_ref1_ref2),&(V_last));

    V_ref3_newpos.x = bond_length*(-N_V_ref1_ref2.x*cos(alpha)
                                   +V_last.x*sin(alpha)*cos(beta)
                                   +N_V1_X_V2.x*sin(alpha)*sin(beta));
    V_ref3_newpos.y = bond_length*(-N_V_ref1_ref2.y*cos(alpha)
                                   +V_last.y*sin(alpha)*cos(beta)
                                   +N_V1_X_V2.y*sin(alpha)*sin(beta));
    V_ref3_newpos.z = bond_length*(-N_V_ref1_ref2.z*cos(alpha)
                                   +V_last.z*sin(alpha)*cos(beta)
                                   +N_V1_X_V2.z*sin(alpha)*sin(beta));

    /* now figure out the actual coordinates of the new atom */
    atom->loc.x = atoms[Zmat_loc->ref1].loc.x + V_ref3_newpos.x;
    atom->loc.y = atoms[Zmat_loc->ref1].loc.y + V_ref3_newpos.y;
    atom->loc.z = atoms[Zmat_loc->ref1].loc.z + V_ref3_newpos.z;
  }



  if( printing ){
    /* now print out the positions of the atoms */
    fprintf(output_file,"# Positions of atoms from Z-matrix\n");
    for(i=0;i<num_atoms;i++){
      fprintf(output_file,"% 4d % 4s % 8.4lf % 8.4lf % 8.4lf\n",
              i+1,atoms[i].symb,atoms[i].loc.x,
              atoms[i].loc.y,atoms[i].loc.z);
    }
  }
}


