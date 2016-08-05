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

