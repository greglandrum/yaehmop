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
*               stack.c
*
* These are the stack manipulation routines
*
*****************************************************************************/
#include "viewkel.h"

/***
  Recent Edit History:
   18.05.98 gL:
    Added transform_norm to do normalized transformations of points

   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
***/

/****************************************************************************
*
*                   Procedure loadmatrix
*
* Arguments: mat: a pointer to matrix_type
* Returns: none
*
* Action:  REPLACES the current matrix with mat.
*
*****************************************************************************/
void loadmatrix( matrix_type *mat )
{
  if( stack ){
    mat->next = stack->next;

    /* deallocate the old current matrix */
    D_FREE( stack );
  }
  stack = mat;
/*fprintf(stderr,"load\n");*/
}

/****************************************************************************
*
*                    Procedure pushmatrix
*
* Arguments: mat: pointer to matrix_type
* Returns: none
*
* Action:  puts a copy of the current matrix on the top of the stack
*   if 'mat is nonzero, then the matrix which is passed in will
*    be pushed onto the stack instead.
*
*****************************************************************************/
void pushmatrix(matrix_type *mat)
{
  matrix_type *newmat;

  newmat = (matrix_type *)D_CALLOC(1,sizeof(matrix_type));
  if( !newmat )fatal( "Memory Allocation" );

  /* fill in the new matrix */
  if( mat == 0 ){
    bcopy(stack->matrix,newmat->matrix,DIM*DIM*sizeof(float));
  } else{
    bcopy(mat->matrix,newmat->matrix,DIM*DIM*sizeof(float));
  }

  /* put it on top */
  newmat->next = stack;
  stack = newmat;
/*fprintf(stderr,"push\n");*/
}

/****************************************************************************
*
*                   Procedure popmatrix
*
* Arguments: none
* Returns: none
*
* Action:  discards the top matrix from the stack.
*
*****************************************************************************/
void popmatrix(void)
{
  matrix_type *temp;

  if( stack == ident ){
    printf("Why are you trying to free ident?\n");
    stack = 0;
    return;
  }

  if( stack ){
    temp = stack->next;
    /* free the old matrix */
    D_FREE( stack );
    stack = temp;
  }
/*fprintf(stderr,"pop\n");*/
}

/****************************************************************************
*
*                   Procedure multmatrix
*
* Arguments: mat: a pointer to matrix_type
* Returns: none
*
* Action:  multiplies the current matrix by mat and replaces it with the
*          result.
*
*****************************************************************************/
void multmatrix( matrix_type *mat )
{
  int i,j;
  float storage[DIM][DIM];

  /* do the multiplication */
  for(i=0; i<DIM; i++ )
    for( j=0; j<DIM; j++ )
      storage[i][j]=stack->matrix[i][0]*mat->matrix[0][j]+
                    stack->matrix[i][1]*mat->matrix[1][j]+
                    stack->matrix[i][2]*mat->matrix[2][j]+
                    stack->matrix[i][3]*mat->matrix[3][j];

  /* put the result at the top of the stack */
#if 0
  for( i=0; i<DIM; i++ )
    for( j=0; j<DIM; j++ )
      stack->matrix[i][j]=storage[i][j];
#endif
  bcopy(storage,stack->matrix,DIM*DIM*sizeof(float));
/*fprintf(stderr,"mult\n");*/
}
/****************************************************************************
*
*                   Procedure premultmatrix
*
* Arguments: mat: a pointer to matrix_type
* Returns: none
*
* Action:  premultiplies the current matrix by mat and replaces it with the
*          result.
*
*****************************************************************************/
void premultmatrix( matrix_type *mat )
{
  int i,j;
  float storage[DIM][DIM];

  /* do the multiplication */
  for(i=0; i<DIM; i++ )
    for( j=0; j<DIM; j++ )
      storage[i][j]=mat->matrix[i][0]*stack->matrix[0][j]+
                    mat->matrix[i][1]*stack->matrix[1][j]+
                    mat->matrix[i][2]*stack->matrix[2][j]+
                    mat->matrix[i][3]*stack->matrix[3][j];

  /* put the result at the top of the stack */
#if 0
  for( i=0; i<DIM; i++ )
    for( j=0; j<DIM; j++ )
      stack->matrix[i][j]=storage[i][j];
#endif
  bcopy(storage,stack->matrix,DIM*DIM*sizeof(float));
/*fprintf(stderr,"premult\n");*/
}

/****************************************************************************
*
*                   Procedure transform
*
* Arguments: point: a pointer to type point
* Returns: none
*
* Action:  multiplies the current matrix by point and replaces point with the
*          result. (the assumption is made that the last element of point is
*          1.0 ).
*
*****************************************************************************/
void transform( point_type *point )
{
  int i;
  float storage[DIM];

  /* first calculate what the new column vector is */
  for( i=0;i<DIM;i++ )
    storage[i] = stack->matrix[i][0]*point->x+
                   stack->matrix[i][1]*point->y+
                   stack->matrix[i][2]*point->z+
                   stack->matrix[i][3];

  /* now copy the result into point */
  point->x=storage[0]/storage[3];
  point->y=storage[1]/storage[3];
  point->z=storage[2];
}

/****************************************************************************
*
*                   Procedure transform_norm
*
* Arguments: point: a pointer to type point
* Returns: none
*
* Action:  multiplies the current matrix by point and replaces point with the
*          *normalized* result.
*          (the assumption is made that the last element of point is
*          1.0 ).
*
*****************************************************************************/
void transform_norm( point_type *point )
{
  int i;
  float storage[DIM];

  /* first calculate what the new column vector is */
  for( i=0;i<DIM;i++ )
    storage[i] = stack->matrix[i][0]*point->x+
                   stack->matrix[i][1]*point->y+
                   stack->matrix[i][2]*point->z+
                   stack->matrix[i][3];

  /* now copy the raw result into point */
  point->x=storage[0];
  point->y=storage[1];
  point->z=storage[2];
  /* normalize */
  V3Normalize(point);
  /* and do the homogeneous transformation */
  point->x /= storage[3];
  point->y /= storage[3];
}

/****************************************************************************
*
*                   Procedure getmatrix
*
* Arguments: mat: a pointer to matrix_type
* Returns: none
*
* Action:  copies the current matrix into mat
*
*****************************************************************************/
void getmatrix(  matrix_type *mat )
{
#if 0
  for( i=0; i<DIM; i++ )
    for( j=0; j<DIM; j++ )
      mat->matrix[i][j]=stack->matrix[i][j];
#endif
  bcopy(stack->matrix,mat->matrix,DIM*DIM*sizeof(float));
/*fprintf(stderr,"get\n");*/
}

/****************************************************************************
*
*                   Procedure xrot
*
* Arguments: theta: a float
* Returns: none
*
* Action:  rotates the current matrix by the euler angle passed in
*
*****************************************************************************/
void xrot( float theta )
{
  int i;
  float temp;
  float cosvar, sinvar;

  /* these variables save calls to cos and sin */
  cosvar=(float)cos((double)theta);
  sinvar=(float)sin((double)theta);

  for(i=0;i<4;i++){
    temp = stack->matrix[i][1];
    stack->matrix[i][1] = temp*cosvar - stack->matrix[i][2]*sinvar;
    stack->matrix[i][2] = temp*sinvar + stack->matrix[i][2]*cosvar;
  }
}

/****************************************************************************
*
*                   Procedure yrot
*
* Arguments: theta: a float
* Returns: none
*
* Action:  rotates the current matrix by the euler angle passed in
*
*****************************************************************************/
void yrot( float theta )
{
  int i;
  float cosvar, sinvar, temp;

  /* these variables save calls to cos and sin */
  cosvar=(float)cos((double)theta);
  sinvar=(float)sin((double)theta);

  for(i=0;i<4;i++){
    temp = stack->matrix[i][0];
    stack->matrix[i][0] = temp*cosvar + stack->matrix[i][2]*sinvar;
    stack->matrix[i][2] = stack->matrix[i][2]*cosvar - temp*sinvar;
  }
}

/****************************************************************************
*
*                   Procedure zrot
*
* Arguments: theta: a float
* Returns: none
*
* Action:  rotates the current matrix by the euler angle passed in
*
*****************************************************************************/
void zrot( float theta )
{
  int i;
  float cosvar, sinvar, temp;

  /* these variables save calls to cos and sin */
  cosvar=(float)cos((double)theta);
  sinvar=(float)sin((double)theta);

  for(i=0;i<4;i++){
    temp = stack->matrix[i][0];
    stack->matrix[i][0] = temp*cosvar - stack->matrix[i][1]*sinvar;
    stack->matrix[i][1] = stack->matrix[i][1]*cosvar + temp*sinvar;
  }
}

/****************************************************************************
*
*                   Procedure translate
*
* Arguments: xtrans, ytrans, ztrans: floats
* Returns: none
*
* Action:  translates the current matrix xtrans, ytrans, and ztrans units
*
*****************************************************************************/
void translate( float xtrans,float ytrans,float ztrans )
{
  int i;

  for(i=0;i<4;i++){
    stack->matrix[i][3] += (stack->matrix[i][0]*xtrans +
                            stack->matrix[i][1]*ytrans +
                            stack->matrix[i][2]*ztrans);
  }
}

/****************************************************************************
*
*                   Procedure scale
*
* Arguments: xscale, yscale, zscale: floats
* Returns: none
*
* Action:  scales the current matrix by xscale, yscale and zscale
*
*****************************************************************************/
void scale( float xscale, float yscale, float zscale )
{
  int i;
  for( i=0; i<4; i++){
    stack->matrix[i][0] *= xscale;
    stack->matrix[i][1] *= yscale;
    stack->matrix[i][2] *= zscale;
  }
}

/****************************************************************************
*
*                   Procedure shear
*
* Arguments: sh1,sh2: floats
*              which: int
* Returns: none
*
* Action:  shears the current matrix
*
*****************************************************************************/
void shear(float sh1,float sh2,int which)
{
  static matrix_type *shearmat=0;

  if( !shearmat ){
    shearmat = (matrix_type *)D_CALLOC(1,sizeof(matrix_type));
    if( !shearmat ) fatal( "Memory allocation" );
  }
  else bzero(shearmat,sizeof(matrix_type));


  /* set up the diagonal elements of the shear matrix */
  shearmat->matrix[0][0] = 1;
  shearmat->matrix[1][1] = 1;
  shearmat->matrix[2][2] = 1;
  shearmat->matrix[3][3] = 1;

  /* now do the elements which are specific to which type shearing is being
     done */
  if( which == X_AX ){
    shearmat->matrix[1][0]=sh1;
    shearmat->matrix[2][0]=sh2;
  }
  else if( which == Y_AX ){
    shearmat->matrix[0][1]=sh1;
    shearmat->matrix[2][1]=sh2;
  }
  else if( which == Z_AX ){
    shearmat->matrix[0][2]=sh1;
    shearmat->matrix[1][2]=sh2;
  }
  else{
    error( "Bogus call to shear!" );
    return;
  }
  multmatrix( shearmat );

}


/****************************************************************************
*
*                   Procedure ortho2
*
* Arguments: left, right, bottom, top
* Returns: none
*
* Action:  does the matrix manipulations to scale the viewport to accomodate
*  world coordinates in the region left, right, bottom, top
*
*****************************************************************************/
void ortho2( float left,float right,float bottom,float top )
{
  float xtemp, ytemp;

  /* do the scaling */
  xtemp = g_xmax / (right - left);
  ytemp = g_ymax / (top - bottom);
  scale( xtemp, -ytemp, 0.0 );

 /* do the translation */
  translate( -left, top, 0.0 );
}


