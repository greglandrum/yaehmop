/*******************************************************
*      Copyright (C) 1995 Greg Landrum
*
*  This file is part of yaehmop.
*
*   This is free software.
* 
*  Permission is granted to modify, or otherwise fold, spindle, and mutilate this
*    code provided all copyright notices are left intact.
*
*  This code may be distributed to your heart's content, in whatever form,
*    provided no fee is charged for the distribution, all copyright notices are
*    left intact, and the source is distributed (without fee) along with any
*    binaries to anyone who requests it.
*
*  There are, of course, no warranties at all on this program.
*
********************************************************************/


/****************************************************************************
*
*     this file contains the defines and data structures for dealing with
*      symmetry.
*
*     this is kept separate from other files, since it may be generally useful
*
*     it requires that some kind of point_type is defined in a previously
*       #included file.
*
*  created:  greg landrum  February 1994
*
*****************************************************************************/


#ifndef SYM_OPS_DEFINED
#define SYM_OPS_DEFINED


/************
  the dimensionality of the transformation matrices (this should be 3 so
  long as no translations are being looked for)
*************/
#define T_MAT_DIM 3


/************
  the dimensionality of the d orbital transformation matrices (this should be
  5 so long as no translations are being looked for)
*************/
#define D_T_MAT_DIM 5

/* maximum order of a rotation axis */
#define MAX_ORDER 8

/* maximum order of an off axis rotation axis */
#define MAX_OFF_AXIS_ORDER 2

/* tolerance for declaring positions the same */
#define SYMM_TOL 1e-3


/******
  used to indicate whether or not equivalent atoms should be printed when a
  symmetry element is displayed.
*******/
#define NO_EQUIV 0
#define SHOW_EQUIV 1

/* determine if two points are the same to within a tolerance */
#define POINTS_ARE_THE_SAME(a,b,tol) \
  (fabs((a)->x-(b)->x)<tol && fabs((a)->y-(b)->y)<tol && fabs((a)->z-(b)->z)<tol?1:0)


/**********
  list names of all the symmetry operations that will be searched for here

  to facilitate counting the number of operations, make sure that Identity
   is last in this list.
***********/
enum possible_sym_op { Inversion, Rotation, Improper_Rotation,
			  Mirror, Identity };

/************
  the axes that will be checked (note that this will mean different things
    for different sym_ops).

  again, to facilitate the counting of operations, please make sure that
    No_Axis is the last element in this list
************/
enum possible_axis { X_Ax, Y_Ax, Z_Ax, No_Axis };


/********
  the structure used to store the individual symmetry operations and their
  matrix representations.

  the symmetry ops are gonna be stored as a linked list for simplicity's and
  efficiency's sake.
*********/
typedef struct sym_op_type_def{
  enum possible_sym_op type;
  int order;
  char present,redundant;
  real t_mat[T_MAT_DIM][T_MAT_DIM];
  real d_t_mat[D_T_MAT_DIM][D_T_MAT_DIM];
  point_type axis;
  real angle;

  /******
    this is an array, num_atoms long, that is used to keep track of which
    atoms are interconverted by this symmetry operation

    so, if equiv_atoms[1] = 3, then atom 1 is transformed to the coordinates
     of atom 3.
  *******/
  int *equiv_atoms;
  
  struct sym_op_type_def *next;
} sym_op_type;

extern sym_op_type *sym_ops_present;


#endif
