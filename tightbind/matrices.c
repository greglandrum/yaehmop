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
*     this file contains stuff for the overlap matrices in K space
*
*  created:  greg landrum  September 1995
*
*****************************************************************************/
#include "bind.h"

/****************************************************************************
*
*                   Procedure sparsify_hermetian_matrix
*
* Arguments:  value: real
*           mat: hermetian_matrix_type
*          num_orbs: int
*
* Returns: none
*
* Action: zeroes all the elements of 'mat which are smaller than 'value
*
****************************************************************************/
void sparsify_hermetian_matrix(value,mat,num_orbs)
  real value;
  hermetian_matrix_type mat;
  int num_orbs;
{
  int i,j,itab,jtab;
  int num_zeroed;
  real mag,val_squared;
  real real_part,imag_part;
  
  fprintf(stderr,"sparsify_hermetian_matrix: %lf\n",value);
  val_squared = value*value;

  num_zeroed = 0;
  for(i=0;i<num_orbs;i++){
    itab = i*num_orbs;
    for(j=i+1;j<num_orbs;j++){
      jtab = j*num_orbs;
      real_part = mat.mat[itab+j];
      imag_part = mat.mat[jtab+i];
      mag = real_part*real_part + imag_part*imag_part;
      if( mag <= val_squared ){
	mat.mat[itab+j] = 0.0;
	mat.mat[jtab+i] = 0.0;
	num_zeroed++;
      }
    }
    real_part = mat.mat[itab+i];
    if( real_part*real_part <= val_squared ){
      mat.mat[itab+i] = 0.0;
      num_zeroed++;
    }
    
  }
  fprintf(stderr,"\t %d of %d zeroed (%4.2lf%% nonzero)\n",num_zeroed,num_orbs*num_orbs,
	  100.0*((real)(num_orbs*num_orbs)-(real)2*num_zeroed)/(real)(num_orbs*num_orbs));
}
  

/****************************************************************************
*
*                   Procedure sparsify_matrix
*
* Arguments:  value: real
*           mat_R,mat_I: pointers to int
*          num_orbs: int
*
* Returns: none
*
* Action: zeroes all the elements of 'mat (where 'mat = 'mat_R + 'mat_I) 
*         which are smaller than 'value
*
****************************************************************************/
void sparsify_matrix(value,mat_R,mat_I,num_orbs)
  real value;
  real *mat_R,*mat_I;
  int num_orbs;
{
  int i,j,itab,jtab;
  int num_zeroed;
  real mag,val_squared;
  real real_part,imag_part;
  
  fprintf(stderr,"sparsify_matrix: %lf\n",value);
  val_squared = value*value;

  num_zeroed = 0;
  for(i=0;i<num_orbs;i++){
    itab = i*num_orbs;
    for(j=0;j<num_orbs;j++){
      jtab = j*num_orbs;
      real_part = mat_R[itab+j];
      if( mat_I ){
	imag_part = mat_I[jtab+i];
      }else{
	imag_part = 0.0;
      }
      mag = real_part*real_part + imag_part*imag_part;
      if( mag <= val_squared ){
	mat_R[itab+j] = 0.0;
	if( mat_I ) mat_I[itab+j] = 0.0;
	num_zeroed++;
      }
    }
  }
  fprintf(stderr,"\t %d of %d zeroed (%4.2lf%% nonzero)\n",num_zeroed,num_orbs*num_orbs,
	  100.0*((real)(num_orbs*num_orbs)-(real)num_zeroed)/(real)(num_orbs*num_orbs));
}
  
