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
*     this is the stuff for doing transformations of atoms/vectors/etc.
*
*  created:  greg landrum  February 1994
*
*****************************************************************************/
#include "bind.h"
#include "symmetry.h"



/****************************************************************************
*
*                   Procedure vector_diff
*
* Arguments: vect1,vect2: pointers to point_type
*                   diff: pointer to point_type
*
* Returns: none
*
* Action:  This evaluates the difference between the two vectors 'vect1 and 'vect2.
*   the resulting difference is returned in 'diff.
*
*****************************************************************************/
void vector_diff(vect1,vect2,diff)
  point_type *vect1,*vect2,*diff;
{
  diff->x = vect1->x-vect2->x;
  diff->y = vect1->y-vect2->y;
  diff->z = vect1->z-vect2->z;
}

/****************************************************************************
*
*                   Procedure normalize_vector
*
* Arguments: vect, norm_vect: pointers to point_type
*
* Returns: none
*
* Action:  This normalizes vector 'vect and returns the result in 'norm_vect
*
*****************************************************************************/
void normalize_vector(vect,norm_vect)
  point_type *vect,*norm_vect;
{
  real norm_fact;

  norm_fact = vect->x*vect->x + vect->y*vect->y + vect->z*vect->z;
  norm_fact = 1/sqrt(norm_fact);

  norm_vect->x = vect->x * norm_fact;
  norm_vect->y = vect->y * norm_fact;
  norm_vect->z = vect->z * norm_fact;
}


/****************************************************************************
*
*                   Function dot_prod
*
* Arguments: vect1,vect2: pointers to point_type
*
* Returns: real
*
* Action:  Returns the dot product of 'vect1 and 'vect2
*
*****************************************************************************/
real dot_prod(vect1,vect2)
  point_type *vect1,*vect2;
{
  return(vect1->x*vect2->x + vect1->y*vect2->y + vect1->z*vect2->z);
}

/****************************************************************************
*
*                  Procedure scale_vector
*
* Arguments: vect: pointer to point_type
*           scalar: real
*
* Returns: none
*
* Action: multiplies each element of 'vect by 'scalar
*
*****************************************************************************/
void scale_vector(point_type *vect,real scalar)
{
  vect->x *= scalar;
  vect->y *= scalar;
  vect->z *= scalar;
}

/****************************************************************************
*
*                   Procedure cross_prod
*
* Arguments: vect1,vect2: pointers to point_type
*                 result: pointer to point_type
*
* Returns: none
*
* Action:  evaluates the cross product of 'vect1 and 'vect2, returns the resulting
*   vector in 'result
*
*****************************************************************************/
void cross_prod(vect1,vect2,result)
  point_type *vect1,*vect2,*result;
{

  result->x = vect1->y*vect2->z - vect1->z*vect2->y;
  result->y = vect1->z*vect2->x - vect1->x*vect2->z;
  result->z = vect1->x*vect2->y - vect1->y*vect2->x;
}

/****************************************************************************
*
*                   Procedure mult_matrices
*
* Arguments: mat1,mat2,result: pointers to reals
*               dim: integer
*
* Returns: none
*
* Action:  sets 'result = ('mat1)('mat2)
*
*****************************************************************************/
void mult_matrices(real *mat1,real *mat2,real *result,int dim)
{
  int i,j,k,itab;

  /* start by zeroing the result */
  bzero(result,dim*dim*sizeof(real));

  for(i=0;i<dim;i++){
    itab = i*dim;
    for(j=0;j<dim;j++){
      for(k=0;k<dim;k++){
        result[itab+j] += mat1[itab+k]*mat2[k*dim+j];
      }
    }
  }
}


/****************************************************************************
*
*                   Procedure translate_atoms
*
* Arguments: atom_locs: pointer to point_type
*                  pos: point_type
*            num_atoms: int
*
* Returns: none
*
* Action:  Translates all the atomic positions in 'atom_locs by 'pos
*
*****************************************************************************/
void translate_atoms(atom_locs,pos,num_atoms)
  point_type *atom_locs,pos;
  int num_atoms;
{
  int i;

  /* This is pretty simple */
  for(i=0;i<num_atoms;i++){
    atom_locs[i].x += pos.x;
    atom_locs[i].y += pos.y;
    atom_locs[i].z += pos.z;
  }
}



/****************************************************************************
*
*                   Procedure transform_atomic_locs
*
* Arguments: atom_locs: pointer to point_type
*                t_mat: a square T_MAT_DIM matrix of reals
*            num_atoms: int
*
* Returns: none
*
* Action:  Transforms all the atom positions using the matrix passed in.
*
*
*  NOTE: this is written to be more or less transparent for use with homogeneous
*    coordinates (T_MAT_DIM > 3), but you have to put in code to
*    fill in the values for the homogeneous variables.
*
*****************************************************************************/
void transform_atomic_locs(atom_locs,t_mat,num_atoms)
  point_type *atom_locs;
  real t_mat[T_MAT_DIM][T_MAT_DIM];
  int num_atoms;
{
  int atom,i,j;
  static real loc[T_MAT_DIM],new_loc[T_MAT_DIM];


  for(atom=0;atom<num_atoms;atom++){
    loc[0] = atom_locs[atom].x;
    loc[1] = atom_locs[atom].y;
    loc[2] = atom_locs[atom].z;
    /* insert code for homogeneous coordinates here! */

    for(i=0;i<T_MAT_DIM;i++){
      new_loc[i] = 0.0;
      for(j=0; j<T_MAT_DIM; j++){
        new_loc[i] += loc[j]*t_mat[i][j];
      }
    }

    /* copy those coordinates in */
    atom_locs[atom].x = new_loc[0];
    atom_locs[atom].y = new_loc[1];
    atom_locs[atom].z = new_loc[2];
    /* insert code for homogeneous coordinates here! */
  }
  /* that's all there is... */
}


/****************************************************************************
*
*                   Procedure transform_one_point
*
* Arguments:  the_point: pointer to point_type
*                 t_mat: a square T_MAT_DIM matrix of reals
*
* Returns: none
*
* Action:  Transforms 'the_point using the matrix passed in.
*
*
*  NOTE: this is written to be more or less transparent for use with homogeneous
*    coordinates (T_MAT_DIM > 3), but you have to put in code to
*    fill in the values for the homogeneous variables.
*
*****************************************************************************/
void transform_one_point(point_type *the_point,real t_mat[T_MAT_DIM][T_MAT_DIM])
{
  static real loc[T_MAT_DIM],new_loc[T_MAT_DIM];
  int i,j;

  loc[0] = the_point->x;
  loc[1] = the_point->y;
  loc[2] = the_point->z;

  /* insert code for homogeneous coordinates here! */
  for(i=0;i<T_MAT_DIM;i++){
    new_loc[i] = 0.0;
    for(j=0; j<T_MAT_DIM; j++){
      new_loc[i] += loc[j]*t_mat[i][j];
    }
  }

  /* copy those coordinates in */
  the_point->x = new_loc[0];
  the_point->y = new_loc[1];
  the_point->z = new_loc[2];
  /* insert code for homogeneous coordinates here! */

  /* that's all there is... */
}


/****************************************************************************
*
*                   Procedure transform_3x3_transpose
*
* Arguments: atom_locs: pointer to point_type
*                t_mat: a 3x3 matrix of reals
*            num_atoms: int
*
* Returns: none
*
* Action:  Transforms the coordinates of the atoms passed in using the
*   transpose of the transformation matrix 't_mat.
*
*
*  NOTE: this routine exists because of the form the meschach libraries
*    use for the eigenvectors of a matrix.
*
*****************************************************************************/
void transform_3x3_transpose(atom_locs,t_mat,num_atoms)
  point_type *atom_locs;
  real t_mat[3][3];
  int num_atoms;
{
  int atom,i,j;
  static real loc[3],new_loc[3];


  for(atom=0;atom<num_atoms;atom++){
    loc[0] = atom_locs[atom].x;
    loc[1] = atom_locs[atom].y;
    loc[2] = atom_locs[atom].z;
    /* insert code for homogeneous coordinates here! */

    for(i=0;i<3;i++){
      new_loc[i] = 0.0;
      for(j=0; j<3; j++){
        new_loc[i] += loc[j]*t_mat[j][i];
      }
    }

    /* copy those coordinates in */
    atom_locs[atom].x = new_loc[0];
    atom_locs[atom].y = new_loc[1];
    atom_locs[atom].z = new_loc[2];
    /* insert code for homogeneous coordinates here! */
  }
  /* that's all there is... */
}


/****************************************************************************
*
*                   Procedure transform_p_orbs
*
* Arguments:    coeffs: pointer to real
*                t_mat: pointer to pointer to real
*
* Returns: none
*
* Action:  Multiplies the vector of 'coeffs by the transformation matrix
*  't_mat
*
*  NOTE: this is written to be more or less transparent for use with homogeneous
*    coordinates (T_MAT_DIM > 3), but you have to put in code to
*    fill in the values for the homogeneous variables.
*
*****************************************************************************/
void transform_p_orbs(coeffs,t_mat)
  real *coeffs,t_mat[T_MAT_DIM][T_MAT_DIM];
{
  static real result[T_MAT_DIM];
  int i,j;

  for( i=0; i<T_MAT_DIM; i++ ){
    result[i] = 0.0;
    for( j=0; j<T_MAT_DIM; j++ ){
      result[i] += t_mat[i][j]*coeffs[j];
    }
  }
  /* copy the result into the coeffs array */
  bcopy(result,coeffs,3*sizeof(real));

  /* that's it! */
}

/****************************************************************************
*
*                   Procedure transform_d_orbs
*
* Arguments:    coeffs: pointer to real
*              d_t_mat: pointer to pointer to real
*
* Returns: none
*
* Action:  Multiplies the vector of 'coeffs by the transformation matrix
*  't_mat
*
*  NOTE: this is written to be more or less transparent for use with homogeneous
*    coordinates (D_T_MAT_DIM > 3), but you have to put in code to
*    fill in the values for the homogeneous variables.
*
*****************************************************************************/
void transform_d_orbs(coeffs,d_t_mat)
  real *coeffs,d_t_mat[D_T_MAT_DIM][D_T_MAT_DIM];
{
  static real result[D_T_MAT_DIM];
  int i,j;

  bzero(result,D_T_MAT_DIM*sizeof(real));
  for( i=0; i<D_T_MAT_DIM; i++ ){
    for( j=0; j<D_T_MAT_DIM; j++ ){
      result[i] += d_t_mat[i][j]*coeffs[j];
    }
  }
  /* copy the result into the coeffs array */
  bcopy(result,coeffs,5*sizeof(real));

  /* that's it! */
}


/****************************************************************************
*
*                   Procedure transform_orbitals
*
* Arguments:      atom: pointer to atom_type
*               coeffs: pointer to real
*               sym_op: pointer to sym_op_type
*
* Returns: none
*
* Action:   This transforms the atomic orbitals in 'coeffs according to
*   the matrices in 'sym_op.
*
*****************************************************************************/
void transform_orbitals(atom,coeffs,sym_op)
  atom_type *atom;
  real *coeffs;
  sym_op_type *sym_op;
{
  real *coeff_tab;

  coeff_tab = coeffs;

  /* we don't need to do anything to s orbitals */
  if( atom->ns > 0 ) coeff_tab++;

  /* check to see if there are p orbitals to transform */
  if( atom->np > 0 ){
    /* there are, deal with them */
    transform_p_orbs(coeff_tab,sym_op->t_mat);
    coeff_tab += 3;
  }
  if( atom->nd > 0 ){
    transform_d_orbs(coeff_tab,sym_op->d_t_mat);
    coeff_tab += 5;
  }

  /* that's all she wrote */
}

/****************************************************************************
*
*                   Procedure full_transform
*
* Arguments: atoms: pointer to atom_type
*                  COM: point_type
*                t_mat: a 3x3 matrix of reals
*            num_atoms: int
*
* Returns: none
*
* Action:  Transforms the coordinates of the atoms passed in using the
*   transpose of the transformation matrix 't_mat, and the center of mass
*   'COM
*
*
*  NOTE: this routine exists because of the form the meschach libraries
*    use for the eigenvectors of a matrix.
*
*****************************************************************************/
void full_transform(atoms,COM,t_mat,num_atoms)
  atom_type *atoms;
  point_type COM;
  real t_mat[3][3];
  int num_atoms;
{
  int atom,i,j;
  static real loc[3],new_loc[3];


  for(atom=0;atom<num_atoms;atom++){
    loc[0] = atoms[atom].loc.x + COM.x;
    loc[1] = atoms[atom].loc.y + COM.y;
    loc[2] = atoms[atom].loc.z + COM.z;
    /* insert code for homogeneous coordinates here! */

    for(i=0;i<3;i++){
      new_loc[i] = 0.0;
      for(j=0; j<3; j++){
        new_loc[i] += loc[j]*t_mat[j][i];
      }
    }

    /* copy those coordinates in */
    atoms[atom].loc.x = new_loc[0];
    atoms[atom].loc.y = new_loc[1];
    atoms[atom].loc.z = new_loc[2];
    /* insert code for homogeneous coordinates here! */
  }
  /* that's all there is... */
}



/****************************************************************************
*
*                   Procedure transform_atoms
*
* Arguments: atom_locs: pointer to point_type
*                t_mat: a square T_MAT_DIM matrix of reals
*            num_atoms: int
*
* Returns: none
*
* Action:  Transforms all the atom positions using the matrix passed in.
*
*
*       This is the procedure called to transform the atoms from
*         crystallographic->cartestan coordinates.
*
*  NOTE: this is written to be more or less transparent for use with homogeneous
*    coordinates (T_MAT_DIM > 3), but you have to put in code to
*    fill in the values for the homogeneous variables.
*
*****************************************************************************/
void transform_atoms(atoms,t_mat,num_atoms)
  atom_type *atoms;
  real t_mat[T_MAT_DIM][T_MAT_DIM];
  int num_atoms;
{
  int atom,i,j;
  static real loc[T_MAT_DIM],new_loc[T_MAT_DIM];


  for(atom=0;atom<num_atoms;atom++){
    loc[0] = atoms[atom].loc.x;
    loc[1] = atoms[atom].loc.y;
    loc[2] = atoms[atom].loc.z;
    /* insert code for homogeneous coordinates here! */

    for(i=0;i<T_MAT_DIM;i++){
      new_loc[i] = 0.0;
      for(j=0; j<T_MAT_DIM; j++){
        new_loc[i] += loc[j]*t_mat[i][j];
      }
    }

    /* copy those coordinates in */
    atoms[atom].loc.x = new_loc[0];
    atoms[atom].loc.y = new_loc[1];
    atoms[atom].loc.z = new_loc[2];
    /* insert code for homogeneous coordinates here! */
  }
  /* that's all there is... */
}

