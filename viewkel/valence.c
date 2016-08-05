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

/******************************

  valence.c

  This file contains the implementation of the
  bond length->bond valence scheme based upon

  This file implements the bond valence method using
   the parameters and equations given in:
   M. O'Keeffe and N. E. Brese, JACS 113, 3226-3229 (1991)

  Created by gL November 1997

*******************************/
#include "viewkel.h"
#include "valence.h"

#ifdef INCLUDE_BOND_VALENCE
float calc_R0(float c1,float r1,float c2,float r2)
{
  float result;

  result = sqrt(c1)-sqrt(c2);
  result *= result;
  result *= r1*r2;
  result /= c1*r1 + c2*r2;
  result += r1+r2;

  return(result);
}

/****************************************************************************
 *
 *                   Procedure bond_length_to_bond_valence
 *
 * Arguments:  atom1,atom2: pointers to atom_type
 *                 length: float
 *        R0_val,valence: pointers to float
 *
 * Returns: none
 *
 * Action: Calculates the bond valence between atoms 'atom1 and 'atom2
 *   using the bond length 'length.
 *
 *  If *R0_val is nonzero, it will be used to calculate the valence, which
 *   is returned in 'valence, otherwise the ci and ri parameters from the
 *   atoms are used.
 *
 *  On Return: 'RO_val contains the R0 value (calculated or passed in)
 *   and 'valence contains the calculated value of the valence.
 *
 ****************************************************************************/
void bond_length_to_bond_valence(atom_type *atom1,atom_type *atom2,
                                 float length,float *R0_val,
                                 float *valence)
{

  /* do we need to calculate R0? */
  if(*R0_val == 0.0){

    if( atom1->ci == 0.0 || atom2->ci == 0.0 ) *R0_val = 0.0;
    else{
      *R0_val = calc_R0(atom1->ci,atom1->ri,atom2->ci,atom2->ri);
    }
  }


  /******

    if R0 is zero at this point, there is a problem.
    return a valence of 0

  *******/
  if(*R0_val == 0.0 ){
    *valence = 0.0;
    return;
  }
  *valence = (float)exp((double)(*R0_val-length)/(double)BVAL);

}
#endif

