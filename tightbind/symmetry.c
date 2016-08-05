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
*     this file contains everything needed to deal with finding symmetry
*      and applying symmetry operations to molecules.
*
*  created:  greg landrum  February 1994
*
*****************************************************************************/
#include "bind.h"
#include "symmetry.h"

/***
  Recent edit history:

  26.09.98 gL:
   store characters of MOs in details structure

***/


/****************************************************************************
*
*                   Procedure name_sym_element
*
* Arguments: elem: pointer to sym_op_type
*        the_file: pointer to FILE
*       num_atoms: int
*      show_equiv: int
*
* Returns: none
*
* Action:  prints out a description of 'elem on 'the_file
*
*  if show_equiv == SHOW_EQUIV, then the atoms which are interconverted by
*    this symmetry op will be printed out as well.
*
*****************************************************************************/
void name_sym_element(elem,the_file,num_atoms,show_equiv)
  sym_op_type *elem;
  FILE *the_file;
  int num_atoms;
  int show_equiv;
{
  int i;

  switch(elem->type){
  case Inversion:
    fprintf(the_file,"Inversion Point\n");
    break;
  case Rotation:
    fprintf(the_file,"%d-Fold Rotation about (%5.3lf,%5.3lf,%5.3lf)\n",
            (int)(TWOPI/elem->angle), elem->axis.x,elem->axis.y,elem->axis.z);
    break;
  case Improper_Rotation:
    fprintf(the_file,"%d-Fold Improper Rotation about (%5.3lf,%5.3lf,%5.3lf)\n",
            (int)(TWOPI/elem->angle), elem->axis.x,elem->axis.y,elem->axis.z);
    break;
  case Mirror:
    fprintf(the_file,"Mirror Plane with normal: (%5.3lf,%5.3lf,%5.3lf)\n",
            elem->axis.x,elem->axis.y,elem->axis.z);
    break;
  case Identity:
    fprintf(the_file,"Identity.\n");
    break;
  }
  if( show_equiv == SHOW_EQUIV ){
    fprintf(the_file,"   Equivalent Atoms: ");
    for(i=0;i<num_atoms;i++){
      /* don't print for dummy atoms */
      if( elem->equiv_atoms[i] >= 0 ){
        fprintf(the_file,"%d -> %d, ",i+1,elem->equiv_atoms[i]+1);
      }
    }
    fprintf(the_file,"\n");
  }
}



/****************************************************************************
*
*                   Procedure compare_molecules
*
* Arguments:  atoms: pointer to atom_type
*       locs1,locs2: pointers to point_type
*         num_atoms: int
*       equiv_atoms: pointer to int
*           present: pointer to char
*          symm_tol: real
*
* Returns: none
*
* Action:  checks to see if the two molecules with atom type described in
*   'atoms  and atomic positions in 'locs1 and 'locs2 are the same.
*
*   if they are then 'present is  set to be nonzero
*
*  This also constructs the list of equivalent atoms: 'equiv_atoms
*
*****************************************************************************/
void compare_molecules(atoms,locs1,locs2,num_atoms,equiv_atoms,present,symm_tol)
  atom_type *atoms;
  point_type *locs1,*locs2;
  int num_atoms;
  int *equiv_atoms;
  char *present;
  real symm_tol;
{
  int i,j;
  char found;

  /* initialize present to 0 so that we can bomb out at any time */
  *present = 0;

  /*********

    loop over each of the atoms in locs1, checking to see if it is present in locs2

  *********/
  for(i=0;i<num_atoms;i++){
    found = 0;
    /* don't do dummy atoms */
    if( atoms[i].at_number >= 0 ){
      for(j=0;j<num_atoms && !found;j++){
        if( atoms[j].at_number >= 0 ){
          /* check to see if they are the same atom type */
          if( atoms[i].at_number == atoms[j].at_number ){
            /* now check to see if they are in the same location */
            if( fabs(locs1[i].x - locs2[j].x) < symm_tol &&
               fabs(locs1[i].y - locs2[j].y) < symm_tol &&
               fabs(locs1[i].z - locs2[j].z) < symm_tol   ){
              /* they are */
              found = 1;
              equiv_atoms[j] = i;
            }
          }
        }
      }
    } else{
      /* it's a dummy atom, pretend we found it and put a -1 in the equiv_atoms array */
      found = 1;
      equiv_atoms[i] = -1;
    }

    /* if we didn't find this atom, we might as well go ahead and return */
    if( !found ) return;
  }
  /* we found every atom, so set present to 1 and return */
  *present = 1;
}


/****************************************************************************
*
*                   Procedure make_new_sym_op
*
* Arguments: None
*
* Returns: pointer to sym_op_type
*
* Action:  Gets the memory for a new symmetry op and initializes the op
*   to be the identity matrix
*
*****************************************************************************/
sym_op_type *make_new_sym_op()
{
  int i;
  sym_op_type *the_op;

  /* get the memory */
  the_op = (sym_op_type *)calloc(1,sizeof(sym_op_type));
  if( !the_op )fatal("Can't allocate space for a symmetry operation in make_new_sym_op.");

  /* set the diagonal elements to 1 */
  for(i=0;i<T_MAT_DIM;i++){
    the_op->t_mat[i][i] = 1.0;
  }

  for(i=0;i<D_T_MAT_DIM;i++){
    the_op->d_t_mat[i][i] = 1.0;
  }

  /* return the element */
  return( the_op );
}


/****************************************************************************
*
*                   Procedure construct_rotn_mats
*
* Arguments:  axis: enum possible_axis
*            angle: real
*           axis_v: pointer to point_type
*            t_mat: real[T_MAT_DIM][T_MAT_DIM]
*          d_t_mat: real[T_MAT_DIM][T_MAT_DIM]
* Returns: none
*
* Action: sets up the transformation matrices for a rotation.
*
*****************************************************************************/
void construct_rotn_mats(enum possible_axis axis,real angle,
                         point_type *axis_v,
                         real t_mat[T_MAT_DIM][T_MAT_DIM],
                         real d_t_mat[D_T_MAT_DIM][D_T_MAT_DIM])
{
  switch(axis){
  case X_Ax:
    /* cartesian transformations */
    t_mat[1][1] = cos(angle);
    t_mat[2][2] = cos(angle);
    t_mat[1][2] = -sin(angle);
    t_mat[2][1] = sin(angle);

    /* d orbital transformations */
    d_t_mat[0][0] = (1.0/4.0)*cos(2.0*angle) + 3.0/4.0;
    d_t_mat[2][2] = cos(angle);
    d_t_mat[3][3] = cos(angle);
    d_t_mat[4][4] = cos(2.0*angle);
    d_t_mat[1][1] = (3.0/4.0)*cos(2.0*angle) + 1.0/4.0;

    d_t_mat[0][4] = (1.0/2.0)*sin(2.0*angle);
    d_t_mat[0][1] = (sqrt(3.0)/4.0)*cos(2.0*angle) - sqrt(3.0)/4.0;

    d_t_mat[2][3] = -sin(angle);

    d_t_mat[3][2] = sin(angle);

    d_t_mat[4][3] = -(1.0/2.0)*sin(2.0*angle);
    d_t_mat[4][1] = -(sqrt(3.0)/2.0)*sin(2.0*angle);

    d_t_mat[1][0] = (sqrt(3.0)/4.0)*cos(2.0*angle) - sqrt(3.0)/4.0;
    d_t_mat[1][4] = (sqrt(3.0)/2.0)*sin(2.0*angle);

    axis_v->x = 1.0;
    break;
  case Y_Ax:
    /* cartesian transformations */
    t_mat[0][0] = cos(angle);
    t_mat[2][2] = cos(angle);
    t_mat[0][2] = sin(angle);
    t_mat[2][0] = -sin(angle);

    /* d orbital transformations */
    d_t_mat[0][0] = (1.0/4.0)*cos(2.0*angle) + 3.0/4.0;
    d_t_mat[2][2] = cos(angle);
    d_t_mat[3][3] = cos(2.0*angle);
    d_t_mat[4][4] = cos(angle);
    d_t_mat[1][1] = (3.0/4.0)*cos(2.0*angle) + 1.0/4.0;

    d_t_mat[0][3] = -(1.0/2.0)*sin(2.0*angle);
    d_t_mat[0][1] = -(sqrt(3.0)/4.0)*cos(2.0*angle) + sqrt(3.0)/4.0;

    d_t_mat[2][4] = sin(angle);

    d_t_mat[3][2] = (1.0/2.0)*sin(2.0*angle);
    d_t_mat[3][1] = -(sqrt(3.0)/2.0)*sin(2.0*angle);

    d_t_mat[4][2] = -sin(angle);

    d_t_mat[1][0] = -(sqrt(3.0)/4.0)*cos(2.0*angle) + sqrt(3.0)/4.0;
    d_t_mat[1][3] = (sqrt(3.0)/2.0)*sin(2.0*angle);

    axis_v->y = 1.0;
    break;
  case Z_Ax:
    /* cartesian transformations */
    t_mat[0][0] = cos(angle);
    t_mat[1][1] = cos(angle);
    t_mat[0][1] = -sin(angle);
    t_mat[1][0] = sin(angle);

    /* d orbital transformations */
    d_t_mat[0][0] = cos(2.0*angle);
    d_t_mat[2][2] = cos(2.0*angle);
    d_t_mat[3][3] = cos(angle);
    d_t_mat[4][4] = cos(angle);

    d_t_mat[0][2] = -sin(2.0*angle);
    d_t_mat[2][0] = sin(2.0*angle);

    d_t_mat[3][4] = -sin(angle);
    d_t_mat[4][3] = sin(angle);

    axis_v->z = 1.0;
    break;
  }
}

/****************************************************************************
*
*                   Procedure gen_sym_ops
*
* Arguments:  the_ops: pointer to pointer to sym_op_type
*             num_ops: pointer to int
*
* Returns: none
*
* Action:  Gets the memory and sets up the matrices for the various symmetry
*           operations.
*
*  NOTE: this is intended to be easy to read... not efficient, this function
*   is only called once so who cares?
*
*****************************************************************************/
void gen_sym_ops(the_ops,num_ops)
  sym_op_type **the_ops;
  int *num_ops;
{
  sym_op_type *op_ptr;
  enum possible_sym_op op;
  enum possible_axis axis;
  real angle;
  int i;

  /**********
    get space for the first operation
    *********/
  op_ptr = make_new_sym_op();

  /* set it to the head of the list */
  *the_ops = op_ptr;

  *num_ops = 1;

  /* loop over operations */
  for( op = 0; op < Identity; op++ ){

    if( op != Inversion ){

      /* loop over axes (if needed) */
      for( axis = 0; axis < No_Axis; axis++ ){

        switch(op){
        case Rotation:
          /******
            generate MAX_ORDER-1 rotations about this axis
            ******/
          for( i=2; i<=MAX_ORDER; i++){
            angle = TWOPI / (real)i;
            op_ptr->angle = angle;
            op_ptr->order = i;
            op_ptr->type = Rotation;

            construct_rotn_mats(axis,angle,&(op_ptr->axis),op_ptr->t_mat,
                                op_ptr->d_t_mat);
            /******
              we're done with that element, get space for and put it at the tail of
              the linked list
            ******/
            op_ptr->next = make_new_sym_op();
            op_ptr = op_ptr->next;
            *num_ops++;
          }
          break;
        case Improper_Rotation:
          /*****
            these are standard rotation matrices with a sign flip for
            the reflection
          *****/
          for( i=3; i<=MAX_ORDER; i++){
            angle = TWOPI / (real)i;
            switch(axis){
            case X_Ax:
              op_ptr->t_mat[1][1] = cos(angle);
              op_ptr->t_mat[2][2] = cos(angle);
              op_ptr->t_mat[1][2] = -sin(angle);
              op_ptr->t_mat[2][1] = sin(angle);
              op_ptr->t_mat[0][0] = -1.0;
              op_ptr->axis.x = 1.0;

              /* d orbital transformations */
              op_ptr->d_t_mat[0][0] = (1.0/4.0)*cos(2.0*angle) + 3.0/4.0;
              op_ptr->d_t_mat[2][2] = -cos(angle);
              op_ptr->d_t_mat[3][3] = -cos(angle);
              op_ptr->d_t_mat[4][4] = cos(2.0*angle);
              op_ptr->d_t_mat[1][1] = (3.0/4.0)*cos(2.0*angle) + 1.0/4.0;

              op_ptr->d_t_mat[0][4] = (1.0/2.0)*sin(2.0*angle);
              op_ptr->d_t_mat[0][1] = (sqrt(3.0)/4.0)*cos(2.0*angle) - sqrt(3.0)/4.0;

              op_ptr->d_t_mat[2][3] = sin(angle);

              op_ptr->d_t_mat[3][2] = -sin(angle);

              op_ptr->d_t_mat[4][3] = -(1.0/2.0)*sin(2.0*angle);
              op_ptr->d_t_mat[4][1] = -(sqrt(3.0)/2.0)*sin(2.0*angle);

              op_ptr->d_t_mat[1][0] = (sqrt(3.0)/4.0)*cos(2.0*angle) - sqrt(3.0)/4.0;
              op_ptr->d_t_mat[1][4] = (sqrt(3.0)/2.0)*sin(2.0*angle);

              break;
            case Y_Ax:
              op_ptr->t_mat[0][0] = cos(angle);
              op_ptr->t_mat[2][2] = cos(angle);
              op_ptr->t_mat[0][2] = sin(angle);
              op_ptr->t_mat[2][0] = -sin(angle);
              op_ptr->t_mat[1][1] = -1.0;

              /* d orbital transformations */
              op_ptr->d_t_mat[0][0] = (1.0/4.0)*cos(2.0*angle) + 3.0/4.0;
              op_ptr->d_t_mat[2][2] = -cos(angle);
              op_ptr->d_t_mat[3][3] = cos(2.0*angle);
              op_ptr->d_t_mat[4][4] = -cos(angle);
              op_ptr->d_t_mat[1][1] = (3.0/4.0)*cos(2.0*angle) + 1.0/4.0;

              op_ptr->d_t_mat[0][3] = -(1.0/2.0)*sin(2.0*angle);
              op_ptr->d_t_mat[0][1] = -(sqrt(3.0)/4.0)*cos(2.0*angle) + sqrt(3.0)/4.0;

              op_ptr->d_t_mat[2][4] = -sin(angle);

              op_ptr->d_t_mat[3][2] = (1.0/2.0)*sin(2.0*angle);
              op_ptr->d_t_mat[3][1] = -(sqrt(3.0)/2.0)*sin(2.0*angle);

              op_ptr->d_t_mat[4][2] = sin(angle);

              op_ptr->d_t_mat[1][0] = -(sqrt(3.0)/4.0)*cos(2.0*angle) + sqrt(3.0)/4.0;
              op_ptr->d_t_mat[1][3] = (sqrt(3.0)/2.0)*sin(2.0*angle);

              op_ptr->axis.y = 1.0;
              break;
            case Z_Ax:
              op_ptr->t_mat[0][0] = cos(angle);
              op_ptr->t_mat[1][1] = cos(angle);
              op_ptr->t_mat[0][1] = -sin(angle);
              op_ptr->t_mat[1][0] = sin(angle);
              op_ptr->t_mat[2][2] = -1.0;

              /* d orbital transformations */
              op_ptr->d_t_mat[0][0] = cos(2.0*angle);
              op_ptr->d_t_mat[2][2] = cos(2.0*angle);
              op_ptr->d_t_mat[3][3] = -cos(angle);
              op_ptr->d_t_mat[4][4] = -cos(angle);

              op_ptr->d_t_mat[0][2] = -sin(2.0*angle);
              op_ptr->d_t_mat[2][0] = sin(2.0*angle);

              op_ptr->d_t_mat[3][4] = sin(angle);
              op_ptr->d_t_mat[4][3] = -sin(angle);

              op_ptr->axis.z = 1.0;
              break;
            }
            op_ptr->angle = angle;
            op_ptr->order = i;
            op_ptr->type = Improper_Rotation;
            /******
              we're done with that element, get space for and put it at the tail of
              the linked list
              ******/
            op_ptr->next = make_new_sym_op();
            op_ptr = op_ptr->next;
            *num_ops++;
          }
          break;
        case Mirror:

          op_ptr->type = Mirror;
          switch(axis){
          case X_Ax:
            op_ptr->t_mat[0][0] = -1.0;

            op_ptr->d_t_mat[2][2] = -1.0;
            op_ptr->d_t_mat[3][3] = -1.0;

            op_ptr->axis.x = 1.0;
            break;
          case Y_Ax:
            op_ptr->t_mat[1][1] = -1.0;

            op_ptr->d_t_mat[2][2] = -1.0;
            op_ptr->d_t_mat[4][4] = -1.0;

            op_ptr->axis.y = 1.0;
            break;
          case Z_Ax:
            op_ptr->t_mat[2][2] = -1.0;

            op_ptr->d_t_mat[3][3] = -1.0;
            op_ptr->d_t_mat[4][4] = -1.0;

            op_ptr->axis.z = 1.0;
            break;
          }
          op_ptr->next = make_new_sym_op();
          op_ptr = op_ptr->next;
          *num_ops++;
          break;
        }
      }
    }
    else{
      /* inversion */
      op_ptr->type = Inversion;
      op_ptr->t_mat[0][0] = -1.0;
      op_ptr->t_mat[1][1] = -1.0;
      op_ptr->t_mat[2][2] = -1.0;

      /********

        Since the d orbs are symmetric with respect to inversion,
        they don't need anything special here, we can just keep the
        identity matrix.

      *********/
      op_ptr->next = make_new_sym_op();
      op_ptr = op_ptr->next;
      *num_ops++;
    }
  }

  /* put in the identity operation */
  op_ptr->type = Identity;
}


/****************************************************************************
*
*                   Procedure find_off_axis_sym_ops
*
* Arguments:  the_ops: pointer to sym_op_type
*             num_ops: pointer to int
*
* Returns: none
*
* Action:  Gets the memory, sets up the matrices for, and test for the
*    presence of symmetry operations which do not lie on the cartesian axes.
*    This looks for rotation axes and mirror planes which are perpendicular
*    to rotation axes of order higher than 2.
*
*****************************************************************************/
void find_off_axis_sym_ops(detail_type *details,cell_type *cell,
                           point_type *locs,point_type *new_locs,
                           point_type cell_dim[3],
                           sym_op_type *the_ops,int *num_ops)
{
  sym_op_type *op_ptr,*this_op,*curr_op,*last_op,*orig_last;
  enum possible_axis temp_ax;
  char killed_it;
  real temp_mat[T_MAT_DIM][T_MAT_DIM];
  real temp_mat2[T_MAT_DIM][T_MAT_DIM];
  real t_mat[T_MAT_DIM][T_MAT_DIM];
  real neg_t_mat[T_MAT_DIM][T_MAT_DIM];
  real temp_d_mat[D_T_MAT_DIM][D_T_MAT_DIM];
  real temp_d_mat2[D_T_MAT_DIM][D_T_MAT_DIM];
  real d_t_mat[D_T_MAT_DIM][D_T_MAT_DIM];
  real neg_d_t_mat[D_T_MAT_DIM][D_T_MAT_DIM];
  point_type tformed_cell[3];
  point_type temp_pt;
  real base_angle;
  real angle;
  real mapping[9];

  char try_x,try_y,try_z,order_is_odd;
  char present_for_basis,present_for_lattice;
  int i,k,num_rots;

  /* move to the last operation in the current list */
  if(!the_ops)return;
  last_op = the_ops;
  while(last_op->next) last_op = last_op->next;
  orig_last = last_op;

  /********

    loop through the current operations, when we find a rotation
    axis of order higher than two, try some perpendicular operations.

  ********/
  curr_op = the_ops;

  /* get space for a new operation  */
  op_ptr = make_new_sym_op();
  last_op->next = op_ptr;

  while(curr_op != orig_last ){
    if( !curr_op->redundant && curr_op->type == Rotation &&
       curr_op->order > 2 ){

      /* don't even ask me why this is needed */
      this_op = curr_op;

      /****

        determine which perpendicular axes we should start from.
        if this is an even order axis, then we only need to start
        from one, otherwise we should try starting from both

      ****/
      order_is_odd = curr_op->order % 2;
      try_x = try_y = try_z = 1;
      if( curr_op->axis.x != 0.0 ) {
        try_x = 0;
        if( !order_is_odd ) try_z = 0;
      }
      if( curr_op->axis.y != 0.0 ){
        try_y = 0;
        if( !order_is_odd ) try_z = 0;
      }
      if( curr_op->axis.z != 0.0 ){
        try_z = 0;
        if( !order_is_odd ) try_y = 0;
      }
      /* error checking */
      if( (!order_is_odd && try_x + try_y + try_z != 1) ||
         (order_is_odd && try_x + try_y + try_z != 2)){
        error("sum of try values is bogus.  this is bad round-off error.");
        return;
      }

      /******

        okay, we are going to need to generate rotations
        corresponding to twice the order of this operation

      ******/
      for(num_rots = 1;num_rots<curr_op->order*2;num_rots++){
        base_angle = (real)num_rots*curr_op->angle/2.0;
        bzero(t_mat,T_MAT_DIM*T_MAT_DIM*sizeof(real));
        for(k=0;k<T_MAT_DIM;k++)t_mat[k][k]=1;
        bzero(d_t_mat,D_T_MAT_DIM*D_T_MAT_DIM*sizeof(real));
        for(k=0;k<D_T_MAT_DIM;k++)d_t_mat[k][k]=1;
        bzero(neg_t_mat,T_MAT_DIM*T_MAT_DIM*sizeof(real));
        for(k=0;k<T_MAT_DIM;k++)neg_t_mat[k][k]=1;
        bzero(neg_d_t_mat,D_T_MAT_DIM*D_T_MAT_DIM*sizeof(real));
        for(k=0;k<D_T_MAT_DIM;k++)neg_d_t_mat[k][k]=1;
        if(curr_op->axis.x != 0)temp_ax = X_Ax;
        else if(curr_op->axis.y != 0)temp_ax = Y_Ax;
        else if(curr_op->axis.z != 0)temp_ax = Z_Ax;

        construct_rotn_mats(temp_ax,base_angle,&temp_pt,
                           t_mat,d_t_mat);
        construct_rotn_mats(temp_ax,-base_angle,&temp_pt,
                           neg_t_mat,neg_d_t_mat);

        /******
          generate MAX_OFF_AXIS_ORDER-1 rotations about this axis
          ******/
        for( i=2; i<=MAX_OFF_AXIS_ORDER; i++){
          angle = TWOPI / (real)i;
          if( try_x ){
            if(!op_ptr) FATAL_BUG("Whoops... op_ptr is zero");
            bzero(temp_mat,T_MAT_DIM*T_MAT_DIM*sizeof(real));
            for(k=0;k<T_MAT_DIM;k++)temp_mat[k][k]=1;
            bzero(temp_d_mat,D_T_MAT_DIM*D_T_MAT_DIM*sizeof(real));
            for(k=0;k<D_T_MAT_DIM;k++)temp_d_mat[k][k]=1;
            temp_ax = X_Ax;
            construct_rotn_mats(temp_ax,angle,&(op_ptr->axis),temp_mat,
                                temp_d_mat);
            /* do the matrix multiplication to form the composite operations */
            mult_matrices((real *)temp_mat,(real *)t_mat,
                          (real *)temp_mat2,T_MAT_DIM);
            /* transfom the axis to construct the new axis */
            op_ptr->axis.x = 1.0;
            transform_one_point(&(op_ptr->axis),temp_mat2);
            mult_matrices((real *)neg_t_mat,(real *)temp_mat2,
                          (real *)op_ptr->t_mat,T_MAT_DIM);
            mult_matrices((real *)temp_d_mat,(real *)d_t_mat,
                          (real *)temp_d_mat2,D_T_MAT_DIM);
            mult_matrices((real *)neg_d_t_mat,(real *)temp_d_mat2,
                          (real *)op_ptr->d_t_mat,D_T_MAT_DIM);

            op_ptr->angle = angle;
            op_ptr->order = i;
            op_ptr->type = Rotation;

            /******
              we're done with that element, get space for and put it at the tail of
              the linked list
              ******/
            op_ptr->next = make_new_sym_op();
            last_op = op_ptr;
            op_ptr = op_ptr->next;
          }
          if( try_y ){
            bzero(temp_mat,T_MAT_DIM*T_MAT_DIM*sizeof(real));
            for(k=0;k<T_MAT_DIM;k++)temp_mat[k][k]=1;
            bzero(temp_d_mat,D_T_MAT_DIM*D_T_MAT_DIM*sizeof(real));
            for(k=0;k<D_T_MAT_DIM;k++)temp_d_mat[k][k]=1;
            temp_ax = Y_Ax;
            construct_rotn_mats(temp_ax,angle,&(op_ptr->axis),temp_mat,
                                temp_d_mat);
            mult_matrices((real *)temp_mat,(real *)t_mat,
                          (real *)temp_mat2,T_MAT_DIM);
            op_ptr->axis.y = 1.0;
            transform_one_point(&(op_ptr->axis),temp_mat2);
            mult_matrices((real *)neg_t_mat,(real *)temp_mat2,
                          (real *)op_ptr->t_mat,T_MAT_DIM);
            mult_matrices((real *)temp_d_mat,(real *)d_t_mat,
                          (real *)temp_d_mat2,D_T_MAT_DIM);
            mult_matrices((real *)neg_d_t_mat,(real *)temp_d_mat2,
                          (real *)op_ptr->d_t_mat,D_T_MAT_DIM);
            op_ptr->angle = angle;
            op_ptr->order = i;
            op_ptr->type = Rotation;
            op_ptr->next = make_new_sym_op();
            last_op = op_ptr;
            op_ptr = op_ptr->next;
          }
          if( try_z ){
            bzero(temp_mat,T_MAT_DIM*T_MAT_DIM*sizeof(real));
            for(k=0;k<T_MAT_DIM;k++)temp_mat[k][k]=1;
            bzero(temp_d_mat,D_T_MAT_DIM*D_T_MAT_DIM*sizeof(real));
            for(k=0;k<D_T_MAT_DIM;k++)temp_d_mat[k][k]=1;
            temp_ax = Z_Ax;
            construct_rotn_mats(temp_ax,angle,&(op_ptr->axis),temp_mat,
                                temp_d_mat);
            mult_matrices((real *)temp_mat,(real *)t_mat,
                          (real *)temp_mat2,T_MAT_DIM);
            op_ptr->axis.z = 1.0;
            transform_one_point(&(op_ptr->axis),temp_mat2);
            mult_matrices((real *)neg_t_mat,(real *)temp_mat2,
                          (real *)op_ptr->t_mat,T_MAT_DIM);
            mult_matrices((real *)temp_d_mat,(real *)d_t_mat,
                          (real *)temp_d_mat2,D_T_MAT_DIM);
            mult_matrices((real *)neg_d_t_mat,(real *)temp_d_mat2,
                          (real *)op_ptr->d_t_mat,D_T_MAT_DIM);
            op_ptr->angle = angle;
            op_ptr->order = i;
            op_ptr->type = Rotation;
            op_ptr->next = make_new_sym_op();
            last_op = op_ptr;
            op_ptr = op_ptr->next;
          }
        }
        /* add the mirrors now */
        if( try_x ){
          bzero(temp_mat,T_MAT_DIM*T_MAT_DIM*sizeof(real));
          for(k=0;k<T_MAT_DIM;k++)temp_mat[k][k]=1;
          bzero(temp_d_mat,D_T_MAT_DIM*D_T_MAT_DIM*sizeof(real));
          for(k=0;k<D_T_MAT_DIM;k++)temp_d_mat[k][k]=1;
          temp_mat[0][0] = -1.0;
          temp_d_mat[2][2] = -1.0;
          temp_d_mat[3][3] = -1.0;
          mult_matrices((real *)temp_mat,(real *)t_mat,
                        (real *)temp_mat2,T_MAT_DIM);
          op_ptr->axis.x = 1.0;
            transform_one_point(&(op_ptr->axis),temp_mat2);
          mult_matrices((real *)neg_t_mat,(real *)temp_mat2,
                        (real *)op_ptr->t_mat,T_MAT_DIM);
          mult_matrices((real *)temp_d_mat,(real *)d_t_mat,
                        (real *)temp_d_mat2,D_T_MAT_DIM);
          mult_matrices((real *)neg_d_t_mat,(real *)temp_d_mat2,
                        (real *)op_ptr->d_t_mat,D_T_MAT_DIM);
          op_ptr->type = Mirror;
          op_ptr->next = make_new_sym_op();
          last_op = op_ptr;
          op_ptr = op_ptr->next;
        }
        if( try_y ){
          bzero(temp_mat,T_MAT_DIM*T_MAT_DIM*sizeof(real));
          for(k=0;k<T_MAT_DIM;k++)temp_mat[k][k]=1;
          bzero(temp_d_mat,D_T_MAT_DIM*D_T_MAT_DIM*sizeof(real));
          for(k=0;k<D_T_MAT_DIM;k++)temp_d_mat[k][k]=1;
          temp_mat[1][1] = -1.0;
          temp_d_mat[2][2] = -1.0;
          temp_d_mat[4][4] = -1.0;
          mult_matrices((real *)temp_mat,(real *)t_mat,
                        (real *)temp_mat2,T_MAT_DIM);
          op_ptr->axis.y = 1.0;
          transform_one_point(&(op_ptr->axis),temp_mat2);
          mult_matrices((real *)neg_t_mat,(real *)temp_mat2,
                        (real *)op_ptr->t_mat,T_MAT_DIM);
          mult_matrices((real *)temp_d_mat,(real *)d_t_mat,
                        (real *)temp_d_mat2,D_T_MAT_DIM);
          mult_matrices((real *)neg_d_t_mat,(real *)temp_d_mat2,
                        (real *)op_ptr->d_t_mat,D_T_MAT_DIM);
          op_ptr->type = Mirror;
          op_ptr->next = make_new_sym_op();
          last_op = op_ptr;
          op_ptr = op_ptr->next;
        }
        if( try_z ){
          bzero(temp_mat,T_MAT_DIM*T_MAT_DIM*sizeof(real));
          for(k=0;k<T_MAT_DIM;k++)temp_mat[k][k]=1;
          bzero(temp_d_mat,D_T_MAT_DIM*D_T_MAT_DIM*sizeof(real));
          for(k=0;k<D_T_MAT_DIM;k++)temp_d_mat[k][k]=1;
          temp_mat[2][2] = -1.0;
          temp_d_mat[3][3] = -1.0;
          temp_d_mat[4][4] = -1.0;
          mult_matrices((real *)temp_mat,(real *)t_mat,
                        (real *)temp_mat2,T_MAT_DIM);
          op_ptr->axis.z = 1.0;
          transform_one_point(&(op_ptr->axis),temp_mat2);
          mult_matrices((real *)neg_t_mat,(real *)temp_mat2,
                        (real *)op_ptr->t_mat,T_MAT_DIM);
          mult_matrices((real *)temp_d_mat,(real *)d_t_mat,
                        (real *)temp_d_mat2,D_T_MAT_DIM);
          mult_matrices((real *)neg_d_t_mat,(real *)temp_d_mat2,
                        (real *)op_ptr->d_t_mat,D_T_MAT_DIM);
          op_ptr->type = Mirror;
          op_ptr->next = make_new_sym_op();
          last_op = op_ptr;
          op_ptr = op_ptr->next;
        }
      }
    }
    curr_op = curr_op->next;
  }
  /* there's an extra element allocated at the end, ditch that now */
  last_op->next = 0;
  free(op_ptr);

  /********

    okay, we just generated a whole slew of potential operations,
    now check to see if they are present (joy!)

  ********/
  this_op = orig_last->next;
  last_op = orig_last;
  while(this_op){
    bcopy((char *)locs,(char *)new_locs,
          cell->num_atoms*sizeof(point_type));
    transform_atomic_locs(new_locs,this_op->t_mat,cell->num_atoms);
    if(cell->dim >= 1){
      bcopy((char *)cell_dim,(char *)tformed_cell,3*sizeof(point_type));
      transform_atomic_locs(tformed_cell,this_op->t_mat,3);
    }
    this_op->equiv_atoms = (int *)calloc(cell->num_atoms,sizeof(int));
    if( !(this_op->equiv_atoms) ) fatal("Can't get memory for equiv_atom array.");

    if( cell->dim == 0 ){
      compare_molecules(cell->atoms,locs,new_locs,
                        cell->num_atoms,this_op->equiv_atoms,
                        &(present_for_basis),details->symm_tol);
      present_for_lattice = 1;
    }else{
      compare_crystal_basis(cell,locs,new_locs,cell_dim,
                            cell->num_atoms,this_op->equiv_atoms,
                            &(present_for_basis),details->symm_tol);
#ifdef GAG_ME_WITH_SYMMETRY
      if( present_for_basis ){
        name_sym_element(this_op,stdout,cell->num_atoms,SHOW_EQUIV);
        printf("Is present for the basis.\n");
      }
#endif
      compare_crystal_lattice(cell,cell_dim,tformed_cell,&(present_for_lattice),
                              details->symm_tol,mapping);

#ifdef GAG_ME_WITH_SYMMETRY
      if( present_for_lattice ){
        name_sym_element(this_op,stdout,cell->num_atoms,NO_EQUIV);
        printf("Is present for the lattice.\n");
        printf("The mapping vector is:\n");
        printf("a -> (%d %d %d)\n",(int)mapping[0],(int)mapping[1],(int)mapping[2]);
        printf("b -> (%d %d %d)\n",(int)mapping[3],(int)mapping[4],(int)mapping[5]);
        printf("c -> (%d %d %d)\n",(int)mapping[6],(int)mapping[7],(int)mapping[8]);
      }
#endif
    }
    if( !present_for_basis || !present_for_lattice ){
      /* it's not, so  remove it from the list */
      last_op->next = this_op->next;
      /* free up all the memory used by this element */
      if(this_op->equiv_atoms) free(this_op->equiv_atoms);
      free(this_op);
      this_op = last_op->next;
    }else{
      last_op = this_op;
      this_op = this_op->next;
    }
  }

  /******

    okay, the brute force generation procedure we just used
    is pretty much guaranteed to have generated some duplicate
    operators, so now let's go back through the list and
    remove the duplicates.

  ******/
  last_op = orig_last;
  this_op  = last_op->next;
  while( this_op ){
    op_ptr = the_ops;
    killed_it = 0;
    while(op_ptr && !killed_it ){
      if( op_ptr != this_op ){
        if( op_ptr->type == this_op->type ){
          temp_pt.x = -op_ptr->axis.x;
          temp_pt.y = -op_ptr->axis.y;
          temp_pt.z = -op_ptr->axis.z;
          if( POINTS_ARE_THE_SAME(&(op_ptr->axis),&(this_op->axis),
                                  details->symm_tol) ||
              POINTS_ARE_THE_SAME(&(temp_pt),&(this_op->axis),
                                  details->symm_tol) ){
            if( op_ptr->type != Rotation || op_ptr->order == this_op->order ){
              last_op->next = this_op->next;
              if(this_op->equiv_atoms)free(this_op->equiv_atoms);
              free(this_op);
              if(last_op) this_op = last_op->next;
              killed_it = 1;
            }
          }
        }
      }
      op_ptr = op_ptr->next;
    }
    /* if we removed the element, then we don't need to advance the pointers */
    if( !killed_it ){
            last_op = this_op;
            this_op = this_op->next;
    }
  }

  /* now count the operators */
  *num_ops = 0;
  this_op = the_ops;
  while(this_op){
    (*num_ops)++;
    this_op = this_op->next;
  }
  /* we're done... whew! */
}


/****************************************************************************
*
*                   Procedure find_sym_ops
*
* Arguments:  cell:  pointer to unit_cell type
*
* Returns: none
*
* Action:  finds the symmetry elements possessed by the atoms described in
*   'cell.
*
*   a linked list of the operations which are present is stored in the global
*    variable sym_ops_present
*
*****************************************************************************/
void find_sym_ops(details,cell)
  detail_type *details;
  cell_type *cell;
{
  int i,j,itab,jtab;
  static point_type *COM_locs,*new_locs;
  static int num_ops=0;
  point_type cell_dim[3],tformed_cell[3];
  atom_type *atom;
  sym_op_type *last_op,*this_op;
  int num_atoms;
  char present_for_basis,present_for_lattice;
  real moments[3];
  real mapping[9];
  long int total_mass;
  int num_ops_present=0;


  /* find the dimensions of the unit cell */
  for(i=0;i<cell->dim;i++){
    itab = cell->tvects[i].begin;
    jtab = cell->tvects[i].end;
    cell_dim[i].x = cell->atoms[jtab].loc.x-cell->atoms[itab].loc.x;
    cell_dim[i].y = cell->atoms[jtab].loc.y-cell->atoms[itab].loc.y;
    cell_dim[i].z = cell->atoms[jtab].loc.z-cell->atoms[itab].loc.z;
  }

  fprintf(output_file,
          "\n\n#---------------------- SYMMETRY ANALYSIS ----------------------\n");


  /********
    there's always at least one operation (Identity), so it's safe to
    use num_ops=0 as the condition for initializing everything
  *********/
  num_atoms = cell->num_atoms;
  if( !num_ops ){
    gen_sym_ops(&sym_ops_present,&num_ops);

    /* get space for the arrays used to store locations */
    COM_locs = (point_type *)calloc(num_atoms,sizeof(point_type));
    new_locs = (point_type *)calloc(num_atoms,sizeof(point_type));
    if(!new_locs) fatal("Can't allocate atomic location storage in find_sym_ops.");
  }

  /**********

    find the centre of mass if this is a molecular problem

  ***********/
  if( cell->dim == 0 ){
    cell->COM.x = cell->COM.y = cell->COM.z = 0.0;
    total_mass=0;
    for(i=0;i<num_atoms;i++){
      atom = &(cell->atoms[i]);
      if(atom->at_number >= 0){
        cell->COM.x += atom->loc.x * atom->at_number;
        cell->COM.y += atom->loc.y * atom->at_number;
        cell->COM.z += atom->loc.z * atom->at_number;
        total_mass += atom->at_number;
      }

      /* copy this atoms coordinates into the COM array */
      bcopy((char *)(&atom->loc),&(COM_locs[i]),sizeof(point_type));
    }
    cell->COM.x /= -(real)total_mass;
    cell->COM.y /= -(real)total_mass;
    cell->COM.z /= -(real)total_mass;

    /********

      translate the atoms to the COM

    *********/
    translate_atoms(COM_locs,cell->COM,num_atoms);

    if( details->find_princ_axes ){
      /*********

        find the principle axes of the inertia tensor, then use
        these to transform the atoms into the frame in which this tensor
        is diagonal.

      **********/
      find_princ_axes(cell->atoms,COM_locs,cell->princ_axes,moments,num_atoms);
      transform_3x3_transpose(COM_locs,cell->princ_axes,num_atoms);

    /*********

      display the principle axes and moments of inertia

      *********/
      fprintf(output_file,"\n;  Principle Axes and Moments of Inertia: \n");
      for( i=0;i<3; i++){
        fprintf(output_file,"\tAxis:  %lf %lf %lf   Moment: %lf\n",
                cell->princ_axes[0][i],cell->princ_axes[1][i],
                cell->princ_axes[2][i],moments[i]);
      }

      fprintf(output_file,"\n; Positions of Atoms in Principle Axis Frame:\n");
      for(i=0;i<num_atoms;i++){
        fprintf(output_file,"%d \t%lf %lf %lf \n",i+1,
                COM_locs[i].x,COM_locs[i].y,COM_locs[i].z);
        atom = &(cell->atoms[i]);
        atom->loc.x = COM_locs[i].x;
        atom->loc.y = COM_locs[i].y;
        atom->loc.z = COM_locs[i].z;
      }
    }else{


      bzero(cell->princ_axes,9*sizeof(real));
      cell->princ_axes[0][0] = 1.0;
      cell->princ_axes[1][1] = 1.0;
      cell->princ_axes[2][2] = 1.0;
    }
  }else{
    for(i=0;i<num_atoms;i++){
      atom = &(cell->atoms[i]);
      COM_locs[i].x=atom->loc.x;
      COM_locs[i].y=atom->loc.y;
      COM_locs[i].z=atom->loc.z;
    }
  }

  /********

    now check for the presence of symmetry elements by looping through the
    list of symmetry operations... remove any element which is not present
    from the list.

  *********/
  last_op = 0;
  this_op = sym_ops_present;
  num_ops = 0;

  while(this_op){

    /* for each operation, get a fresh set of coordinates */
    bcopy((char *)COM_locs,(char *)new_locs,num_atoms*sizeof(point_type));

    /* applying the symmetry operation is just a transformation */
    transform_atomic_locs(new_locs,this_op->t_mat,num_atoms);

    /* if it is a crystal, we also need to transform the lattice */
    if(cell->dim >= 1){
      bcopy((char *)cell_dim,(char *)tformed_cell,3*sizeof(point_type));
      transform_atomic_locs(tformed_cell,this_op->t_mat,3);
    }

    /* get space for the equivalent atoms array */
    this_op->equiv_atoms = (int *)calloc(num_atoms,sizeof(int));
    if( !(this_op->equiv_atoms) ) fatal("Can't get memory for equiv_atom array.");

    /* check to see if the transformed molecule is equivalent */
    if( cell->dim == 0 ){
      compare_molecules(cell->atoms,COM_locs,new_locs,num_atoms,this_op->equiv_atoms,
                        &(present_for_basis),details->symm_tol);
      present_for_lattice = 1;
    }else{
      compare_crystal_basis(cell,COM_locs,new_locs,cell_dim,
                            num_atoms,this_op->equiv_atoms,
                            &(present_for_basis),details->symm_tol);
#ifdef GAG_ME_WITH_SYMMETRY
      if( present_for_basis ){
        name_sym_element(this_op,stdout,num_atoms,SHOW_EQUIV);
        printf("Is present for the basis.\n");
      }
#endif
      compare_crystal_lattice(cell,cell_dim,tformed_cell,&(present_for_lattice),
                              details->symm_tol,mapping);

#ifdef GAG_ME_WITH_SYMMETRY
      if( present_for_lattice ){
        name_sym_element(this_op,stdout,num_atoms,NO_EQUIV);
        printf("Is present for the lattice.\n");
        printf("The mapping vector is:\n");
        printf("a -> (%d %d %d)\n",(int)mapping[0],(int)mapping[1],(int)mapping[2]);
        printf("b -> (%d %d %d)\n",(int)mapping[3],(int)mapping[4],(int)mapping[5]);
        printf("c -> (%d %d %d)\n",(int)mapping[6],(int)mapping[7],(int)mapping[8]);
      }
#endif

    }
    /* was this element present? */
    if( !present_for_basis || !present_for_lattice ){
      /* it's not, check to see if it was present before */
      if( this_op->present ){
        /* it was, display some error messages */
        error("A symmetry element has disappeared.");
        fprintf(output_file,
                "; !!!!!!! A symmetry element has disappeared.  This may cause strange results.\n");
        fprintf(output_file,";         The missing element is:  ");
        name_sym_element(this_op,output_file,num_atoms,NO_EQUIV);
      }

      /* now remove it from the list */
      if( last_op ){
        last_op->next = this_op->next;
        /* free up all the memory used by this element */
        if(this_op->equiv_atoms) free(this_op->equiv_atoms);
        free(this_op);

        this_op = last_op->next;
      }
      else{
        /****
          things are slightly different if the first element of the list isn't
           present
        *****/
        sym_ops_present = this_op->next;
  if(this_op->equiv_atoms)free(this_op->equiv_atoms);
        free(this_op);
        this_op = sym_ops_present;
      }
    }
    else{
      /* yeah, it's present, just advance the pointers */
      last_op = this_op;
      this_op = this_op->next;
      num_ops_present++;
    }
  }
  /* find the minimal rotations */
  last_op = sym_ops_present;
  while(last_op){
    if( last_op->type == Rotation ){
      this_op = sym_ops_present;
      while(this_op){
        if( this_op != last_op && this_op->type == Rotation &&
            this_op->axis.x == last_op->axis.x &&
            this_op->axis.y == last_op->axis.y &&
            this_op->axis.z == last_op->axis.z){
         if( this_op->angle / last_op->angle -
             floor(this_op->angle / last_op->angle) < 0.0001 ){
           this_op->redundant = 1;
         }
        }
        this_op = this_op->next;
      }
    }
    last_op = last_op->next;
  }

  find_off_axis_sym_ops(details,cell,COM_locs,new_locs,cell_dim,
                        sym_ops_present,&num_ops_present);

  /* now print out the symmetry operations that are present */
  fprintf(output_file,"\n%d Elements are present. They are:\n",num_ops_present);
  this_op = sym_ops_present;
  i = 1;
  while(this_op){
    fprintf(output_file,"% 3d ",i);
    i++;
    name_sym_element(this_op,output_file,num_atoms,SHOW_EQUIV);
    this_op = this_op->next;
  }
  details->num_sym_ops = num_ops_present;

  fprintf(output_file,"\nHere are the  Non-Redundant Elements which are present:\n");
  this_op = sym_ops_present;
  while(this_op){
    if( !(this_op->redundant) ){
      name_sym_element(this_op,output_file,num_atoms,SHOW_EQUIV);
    }
    this_op = this_op->next;
  }
  fprintf(output_file,
          "\n;-----------------------------------------------------------------\n");

  /* free up the memory we allocated here */
  if(COM_locs) free(COM_locs);
  if(new_locs) free(new_locs);
}


/****************************************************************************
*
*                   Procedure find_walsh_sym_ops
*
* Arguments:  cell:  pointer to unit_cell type
*          details: pointer to detail_type
*
* Returns: none
*
* Action:  This is used to determine which symmetry elements are conserved
*   along a Walsh distortion.
*  Due to the fact that the center of mass or principle axes may change
*  along the distortion, the principle axis frame is only determined ONCE,
*  then each geometry along the distortion is transformed into this basis.
*
*  Any symmetry element which disappears along the distortion is removed
*  from the list, and any element which appears will never be found.
*
*****************************************************************************/
void find_walsh_sym_ops(cell,details)
  cell_type *cell;
  detail_type *details;
{
  int i,j;
  point_type *COM_locs,*new_locs;
  int num_ops=0;
  atom_type *atom;
  sym_op_type *last_op,*this_op;
  int num_atoms;
  char present;
  long int total_mass;
  real moments[3];
  int num_ops_present=0;
  int num_steps;

  /* set up the first geometry */
  walsh_update(cell,details,0,0);

  num_atoms = cell->num_raw_atoms;
  num_steps = details->walsh_details.num_steps;

  /* generate the initial list of symmetry elements */
  gen_sym_ops(&sym_ops_present,&num_ops);

  /* get space for the arrays used to store locations */
  COM_locs = (point_type *)calloc(num_atoms,sizeof(point_type));
  new_locs = (point_type *)calloc(num_atoms,sizeof(point_type));
  if(!new_locs) fatal("Can't allocate atomic location storage in find_sym_ops.");

  /**********

    find the centre of mass

  ***********/
  cell->COM.x = cell->COM.y = cell->COM.z = 0.0;
  total_mass=0;
  for(i=0;i<num_atoms;i++){
    atom = &(cell->atoms[i]);
    if(atom->at_number >= 0){
      cell->COM.x += atom->loc.x * atom->at_number;
      cell->COM.y += atom->loc.y * atom->at_number;
      cell->COM.z += atom->loc.z * atom->at_number;
      total_mass += atom->at_number;
    }

    /* copy this atom's coordinates into the COM array */
    bcopy((char *)(&atom->loc),&(COM_locs[i]),sizeof(point_type));
  }
  cell->COM.x /= -(real)total_mass;
  cell->COM.y /= -(real)total_mass;
  cell->COM.z /= -(real)total_mass;

  /********

    translate the atoms to the COM

  *********/
  translate_atoms(COM_locs,cell->COM,num_atoms);

  if( details->find_princ_axes ){
    /*********

      find the principle axes of the inertia tensor, then use
      these to transform the atoms into the frame in which this tensor
      is diagonal.

      **********/
    find_princ_axes(cell->atoms,COM_locs,cell->princ_axes,moments,num_atoms);
    transform_3x3_transpose(COM_locs,cell->princ_axes,num_atoms);
  }else{
    bzero(cell->princ_axes,9*sizeof(real));
    cell->princ_axes[0][0] = 1.0;
    cell->princ_axes[1][1] = 1.0;
    cell->princ_axes[2][2] = 1.0;
  }

  /********

    now check for the presence of symmetry elements by looping through the
    list of symmetry operations... remove any element which is not present
    from the list.

  *********/
  last_op = 0;
  this_op = sym_ops_present;
  num_ops = 0;

  while(this_op){

    /* for each operation, get a fresh set of coordinates */
    bcopy((char *)COM_locs,(char *)new_locs,num_atoms*sizeof(point_type));

    /* applying the symmetry operation is just a transformation */
    transform_atomic_locs(new_locs,this_op->t_mat,num_atoms);

    /* get space for the equivalent atoms array */
    this_op->equiv_atoms = (int *)calloc(num_atoms,sizeof(int));
    if( !(this_op->equiv_atoms) ) fatal("Can't get memory for equiv_atom array.");

    /* check to see if the transformed molecule is equivalent */
    compare_molecules(cell->atoms,COM_locs,new_locs,num_atoms,this_op->equiv_atoms,
                      &(present),details->symm_tol);

    /* was this element present? */
    if( !present ){
      /* now remove it from the list */
      if( last_op ){
        last_op->next = this_op->next;
        /* free up all the memory used by this element */
        free(this_op->equiv_atoms);
        free(this_op);

        this_op = last_op->next;
      }
      else{
        /****
          things are slightly different if the first element of the list isn't
           present
        *****/
        sym_ops_present = this_op->next;
        free(this_op);
        free(this_op->equiv_atoms);
        this_op = sym_ops_present;
      }
    }
    else{
      /* yeah, it's present, just advance the pointers */
      last_op = this_op;
      this_op = this_op->next;
      num_ops_present++;
    }
  }

  /* loop through the rest of the walsh steps and find the elements present */
  for( i=1;i<num_steps;i++){
    walsh_update(cell,details,i,0);

    for(j=0;j<num_atoms;j++){
      atom = &(cell->atoms[j]);

      /* copy this atom's coordinates into the COM array */
      bcopy((char *)(&atom->loc),&(COM_locs[j]),sizeof(point_type));
    }
    translate_atoms(COM_locs,cell->COM,num_atoms);
    transform_3x3_transpose(COM_locs,cell->princ_axes,num_atoms);

    last_op = 0;
    this_op = sym_ops_present;
    num_ops = 0;
    num_ops_present = 0;
    while(this_op){

      /* for each operation, get a fresh set of coordinates */
      bcopy((char *)COM_locs,(char *)new_locs,num_atoms*sizeof(point_type));

      /* applying the symmetry operation is just a transformation */
      transform_atomic_locs(new_locs,this_op->t_mat,num_atoms);

      /* get space for the equivalent atoms array */
      this_op->equiv_atoms = (int *)calloc(num_atoms,sizeof(int));
      if( !(this_op->equiv_atoms) ) fatal("Can't get memory for equiv_atom array.");

      /* check to see if the transformed molecule is equivalent */
      compare_molecules(cell->atoms,COM_locs,new_locs,num_atoms,this_op->equiv_atoms,
                        &(present),details->symm_tol);

      /* was this element present? */
      if( !present ){
        /* now remove it from the list */
        if( last_op ){
          last_op->next = this_op->next;
          /* free up all the memory used by this element */
          free(this_op->equiv_atoms);
          free(this_op);

          this_op = last_op->next;
        }
        else{
          /****
            things are slightly different if the first element of the list isn't
            present
            *****/
          sym_ops_present = this_op->next;
          free(this_op);
          free(this_op->equiv_atoms);
          this_op = sym_ops_present;
        }
      }
      else{
        /* yeah, it's present, just advance the pointers */
        last_op = this_op;
        this_op = this_op->next;
        num_ops_present++;
      }
    }
  }

  fprintf(output_file,"# Walsh symmetry analysis\n");
  fprintf(output_file,"%d Elements are present throughout the Walsh distortion.\n",
          num_ops_present);
  details->num_sym_ops = num_ops_present;
  fprintf(output_file,"; They are:\n");
  this_op = sym_ops_present;
  while(this_op){
    name_sym_element(this_op,output_file,num_atoms,NO_EQUIV);
    this_op = this_op->next;
  }
  fprintf(output_file,
          "\n;-----------------------------------------------------------------\n");
  /* free up the memory we allocated here */
  if(COM_locs) free(COM_locs);
  if(new_locs) free(new_locs);
}


/****************************************************************************
*
*                   Procedure find_MO_symmetries
*
* Arguments: num_orbs: int
*             details: pointer to detail_type
*                cell: pointer to cell_type
*            eigenset: eigenset_type
*             overlap: hermetian_matrix_type
*orbital_lookup_table: pointer to int
*
* Returns: none
*
* Action:  Finds the characters of the wavefunctions contained in 'eigenset
*  with respect to the previously determined list of symmetry operations:
*  sym_ops_present (a global variable)
*
*  This is done by normalizing each MO, transforming each atom's contribution
*    using the appropriate matrix (standard 3x3 transformation matrices for
*    p orbitals, special 5x5 matrices for d orbitals, no matrix for s orbitals),
*    then multiplying the transformed coefficients by the corresponding coefficients
*    on the interconverted atoms.
*
*   NOTE: this doesn't deal with f orbitals, it needs to be changed to deal
*    with them.
*
*****************************************************************************/
void find_MO_symmetries(num_orbs,details,cell,eigenset,overlap,orbital_lookup_table)
  int num_orbs;
  detail_type *details;
  cell_type *cell;
  eigenset_type eigenset;
  hermetian_matrix_type overlap;
  int *orbital_lookup_table;
{
  static real *AO_coeffs=0;
  static real *norm_fact=0;
  int i,j,k,ops_so_far;
  int num_atoms;
  int atom1,atom2;
  int begin_atom1,end_atom1,begin_atom2,end_atom2;

  real MO_character;
  sym_op_type *sym_op;

  details->characters = (real *)calloc(num_orbs*details->num_sym_ops,sizeof(real));
  if( !details->characters ) fatal("can't allocate details->characters");

  /* check to see if we need to get space for the AO_coeffs array */
  if( !AO_coeffs ){
    /* this needs to be larger for f orbitals */
    AO_coeffs = (real *)calloc(END_D + 1,sizeof(real));
    if( !AO_coeffs ) fatal("Can't get space for AO_coeffs in find_MO_symmetries.");
    norm_fact = (real *)calloc(num_orbs,sizeof(real));
    if( !norm_fact ) fatal("Can't get memory for norm_fact in find_MO_symmetries.");
  }

  num_atoms = cell->num_atoms;

  /********

    Dump some info to the output file

  *********/
  fprintf(output_file,"\n# Characters of wavefunctions with respect to symmetry \
operations present.\n");
        fprintf(output_file,";         ");
        i=1;
    sym_op = sym_ops_present;
    while(sym_op){
          fprintf(output_file,"% 6d ",i);
          i++;
          sym_op=sym_op->next;
          }

        fprintf(output_file,"\n");

  /* determine the normalization constants */
  for(i=0;i<num_orbs;i++){
    norm_fact[i] = 0.0;
    for(j=0;j<num_orbs;j++){
      norm_fact[i] += EIGENVECT_R(eigenset,i,j)*EIGENVECT_R(eigenset,i,j);
    }
    norm_fact[i] = 1.0/norm_fact[i];
  }


  /* first loop over MO's */
  for(i=0;i<num_orbs;i++){
    fprintf(output_file,"% 4d:---> ",i+1);


    /******

      Walk through the list of symmetry elements

    *******/
    ops_so_far = 0;
    sym_op = sym_ops_present;
    while(sym_op){
      /* loop over the list of atoms which are transformed into each other */
      MO_character = 0.0;
      for(atom1=0;atom1<num_atoms;atom1++){
        /* don't worry about dummy atoms */
        if( cell->atoms[atom1].at_number >= 0 ){
          atom2 = sym_op->equiv_atoms[atom1];
          if( cell->atoms[atom2].at_number >= 0 ){

            /* find the orbitals for the 2 atoms */
            find_atoms_orbs(num_orbs,num_atoms,atom1,orbital_lookup_table,
                            &begin_atom1,&end_atom1);
            find_atoms_orbs(num_orbs,num_atoms,atom2,orbital_lookup_table,
                            &begin_atom2,&end_atom2);

            /***** THIS DOESN'T DEAL WITH IMAGINARY COEFFICIENTS *****/

            /**********
              copy the AO's from atom 1 into a vector so that they can be
              transformed.
            ***********/
            if( begin_atom1 != -1 && begin_atom2 != -1 ){
              for(j=begin_atom1;j<end_atom1;j++){
                AO_coeffs[j-begin_atom1] = EIGENVECT_R(eigenset,i,j);
              }

              /* now transform the orbitals */
              transform_orbitals(&(cell->atoms[atom1]),AO_coeffs,sym_op);

              /* loop through the AO's and multiply coefficients */
              for(j=0,k=begin_atom2;j<(end_atom1-begin_atom1);j++,k++){
                MO_character += AO_coeffs[j]*EIGENVECT_R(eigenset,i,k);
              }
            }
          }
        }
      }
      /* normalize */
      MO_character *= norm_fact[i];

      /* report the result */
      fprintf(output_file,"% -4.3lf ",MO_character);

      details->characters[ops_so_far*num_orbs+i] = MO_character;
      ops_so_far++;

      /* go on to the next symmetry element */
      sym_op = sym_op->next;
    }
    fprintf(output_file,"\n");
  }
}
